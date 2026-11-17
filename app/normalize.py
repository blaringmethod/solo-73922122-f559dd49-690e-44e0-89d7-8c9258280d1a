"""FASTA reference parsing and reference-aware VCF allele normalization.

Normalization follows the usual VCF rule: trim common suffix/prefix bases
while keeping at least one base per allele, then shift an indel left while
the last base of the two alleles still agrees.  The result is the left-most
(minimal) representation; all equivalent representations of the same change
in a repeat therefore collapse to one canonical key.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Iterator

import pysam

SEQUENCE_RE = re.compile(r"[ACGT]+")


class FastaError(ValueError):
    """The uploaded reference FASTA violates the acceptance rules."""


class NormalizationError(ValueError):
    """A VCF record does not match the supplied reference sequence."""


@dataclass(frozen=True)
class Reference:
    chrom: str
    seq: str


def read_fasta(path: str) -> Reference:
    """Read exactly one FASTA record whose sequence is uppercase A/C/G/T.

    The record name is the first whitespace-delimited token of the header
    (text after it is allowed as a description).
    """
    name = None
    parts: list[str] = []
    try:
        with open(path, "r", encoding="utf-8", errors="strict", newline=None) as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        raise FastaError(
                            f"line {lineno}: multiple FASTA records; exactly one "
                            "chromosome is required"
                        )
                    header = line[1:].strip()
                    if not header:
                        raise FastaError(f"line {lineno}: empty FASTA header")
                    name = header.split(None, 1)[0]
                    continue
                if name is None:
                    raise FastaError(
                        f"line {lineno}: sequence data before the FASTA header"
                    )
                parts.append(line.strip())
    except UnicodeDecodeError as exc:
        raise FastaError(f"upload is not valid UTF-8 text: {exc}")

    if name is None:
        raise FastaError("missing FASTA header line starting with '>'")
    seq = "".join(parts)
    if not seq:
        raise FastaError("FASTA contains no sequence")
    if not SEQUENCE_RE.fullmatch(seq):
        bad = next(ch for ch in seq if ch not in "ACGT")
        raise FastaError(
            f"reference sequence must contain only A/C/G/T; found {bad!r}"
        )
    return Reference(chrom=name, seq=seq)


def normalize_alleles(pos: int, ref: str, alt: str, seq: str) -> tuple[int, str, str]:
    """Return the left-aligned, base-trimmed (POS, REF, ALT) for one record.

    `seq` is the 1-based reference and `ref` is assumed to match it at `pos`.
    Alleles always stay non-empty, as required by VCF.
    """
    if ref == alt:  # a no-op record has no alternative representation
        return pos, ref, alt
    while True:
        changed = False
        # Redundant common suffix (keep one base in each allele).
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        # Redundant common prefix (keep one base in each allele).
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1
            changed = True
        # Left-align indels: while the alleles end in the same base, prepend
        # the reference base to the left and drop the now-redundant last base.
        while pos > 1 and len(ref) != len(alt) and ref[-1] == alt[-1]:
            prev = seq[pos - 2]
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos -= 1
            changed = True
        if not changed:
            return pos, ref, alt


def iter_records(gz_path: str, chrom: str) -> Iterator[tuple[int, str, str, str]]:
    """Stream the records of `chrom` through the Tabix index.

    Yields (POS, REF, ALT, raw_line) tuples in file (sorted) order.  Reading
    through the index means normalization never performs a full-file scan.
    """
    try:
        tb = pysam.TabixFile(gz_path)
    except Exception as exc:
        raise RuntimeError(f"cannot open published VCF: {exc}")
    try:
        try:
            lines = tb.fetch(chrom)
        except Exception as exc:
            raise RuntimeError(f"indexed fetch failed: {exc}")
        try:
            for lineno, line in enumerate(lines, 1):
                fields = line.split("\t")
                if len(fields) < 5:
                    raise RuntimeError(f"malformed published record {lineno}")
                try:
                    pos = int(fields[1])
                except ValueError:
                    raise RuntimeError(f"malformed POS in published record {lineno}")
                yield pos, fields[3], fields[4], line
        except RuntimeError:
            raise
        except Exception as exc:
            raise RuntimeError(f"indexed fetch failed: {exc}")
    finally:
        tb.close()


def _check_ref(idx: int, pos: int, ref: str, seq: str) -> None:
    end = pos + len(ref) - 1
    if pos < 1 or end > len(seq):
        raise NormalizationError(
            f"record {idx}: REF at POS {pos} spans {pos}-{end}, beyond "
            f"reference length {len(seq)}"
        )
    if seq[pos - 1:end] != ref:
        raise NormalizationError(
            f"record {idx}: REF '{ref}' at POS {pos} does not match the "
            f"reference sequence '{seq[pos - 1:end]}'"
        )


def normalize_records(
    records: Iterable[tuple[int, str, str, str]], reference: Reference
) -> list[dict]:
    """Verify every REF against the reference, then normalize every record.

    Every record is checked before the result list is returned, so a mismatch
    or out-of-range REF fails the whole request rather than yielding partial
    output.  Records are kept one-for-one, including duplicates.
    """
    seq = reference.seq
    out = []
    for idx, (pos, ref, alt, raw) in enumerate(records, 1):
        _check_ref(idx, pos, ref, seq)
        # The upload path only admits ACGT biallelic alleles; guard against a
        # corrupted or externally-modified published file.
        if not SEQUENCE_RE.fullmatch(alt):
            raise RuntimeError(
                f"published record {idx} has an unsupported ALT allele {alt!r}"
            )
        npos, nref, nalt = normalize_alleles(pos, ref, alt, seq)
        out.append(
            {
                "original": raw,
                "original_pos": pos,
                "original_ref": ref,
                "original_alt": alt,
                "pos": npos,
                "ref": nref,
                "alt": nalt,
                "key": f"{reference.chrom}:{npos}:{nref}:{nalt}",
            }
        )
    return out


def normalize_dataset(gz_path: str, reference: Reference) -> list[dict]:
    """Indexed read of the pinned dataset + reference-checked normalization."""
    return normalize_records(iter_records(gz_path, reference.chrom), reference)
