"""Property tests for the normalization algorithm itself.

For a small repeat-rich reference, enumerate every possible biallelic
substitution/indel representation with alleles up to length 4, group the
representations by the haplotype they produce, and require that all
representations of the same change normalize to one canonical allele.
"""
import itertools
from collections import defaultdict

import pytest

from app.normalize import normalize_alleles


def apply(seq: str, pos: int, ref: str, alt: str) -> str:
    assert seq[pos - 1:pos - 1 + len(ref)] == ref
    return seq[:pos - 1] + alt + seq[pos - 1 + len(ref):]


def all_representations(seq: str, max_len: int = 4):
    """Yield (pos, ref, alt, resulting_haplotype) for every simple edit."""
    n = len(seq)
    for pos in range(1, n + 1):
        for rlen in range(1, min(max_len, n - pos + 1) + 1):
            ref = seq[pos - 1:pos - 1 + rlen]
            for alen in range(1, max_len + 1):
                for alt in map("".join, itertools.product("ACGT", repeat=alen)):
                    yield pos, ref, alt, apply(seq, pos, ref, alt)


REFERENCE = "CA" + "A" * 5 + "GT" + "ACGT" * 3 + "T" * 4 + "CAG"


def test_all_equivalent_representations_share_one_canonical_key():
    groups = defaultdict(set)
    for pos, ref, alt, hap in all_representations(REFERENCE):
        if hap == REFERENCE:
            continue  # no-op records are not variants; normalizer leaves them be
        groups[hap].add(normalize_alleles(pos, ref, alt, REFERENCE))

    assert groups, "test generated no variants"
    for hap, normalized in groups.items():
        assert len(normalized) == 1, (
            f"equivalent representations diverged for haplotype {hap}: "
            f"{sorted(normalized)}"
        )


def test_normalized_allele_reproduces_the_same_haplotype():
    for pos, ref, alt, hap in all_representations(REFERENCE):
        if hap == REFERENCE:
            continue
        npos, nref, nalt = normalize_alleles(pos, ref, alt, REFERENCE)
        assert REFERENCE[npos - 1:npos - 1 + len(nref)] == nref
        assert apply(REFERENCE, npos, nref, nalt) == hap


def test_normalized_alleles_are_minimal_and_left_aligned():
    for pos, ref, alt, hap in all_representations(REFERENCE):
        if hap == REFERENCE:
            continue
        npos, nref, nalt = normalize_alleles(pos, ref, alt, REFERENCE)
        # no redundant common suffix / prefix
        assert not (len(nref) > 1 and len(nalt) > 1 and nref[-1] == nalt[-1])
        assert not (len(nref) > 1 and len(nalt) > 1 and nref[0] == nalt[0])
        # no further left shift is possible
        if len(nref) != len(nalt) and npos > 1:
            assert nref[-1] != nalt[-1], (npos, nref, nalt)


@pytest.mark.parametrize("seq", ["A", "AC", "AAAA", "ATATAT", "GCGCGCGC"])
def test_edge_sequences_do_not_crash_or_leave_bounds(seq):
    for pos, ref, alt, hap in all_representations(seq, max_len=2):
        npos, nref, nalt = normalize_alleles(pos, ref, alt, seq)
        assert npos >= 1
        assert nref and nalt
        assert seq[npos - 1:npos - 1 + len(nref)] == nref
        assert apply(seq, npos, nref, nalt) == hap
