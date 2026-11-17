"""Tests for /normalize: reference checks, left-alignment, repeat keys."""
import os
import threading

import pytest
from conftest import make_config, vcf

from app import normalize as norm
from app.main import create_app

# Deliberately repeat-rich reference.  Positions are 1-based:
#   1       C
#   2-11    AAAAAAAAAA            (homopolymer)
#   12      G
#   13-32   ACGTACGTACGTACGTACGT  (tandem repeat)
#   33-40   TTTTTTTT
#   41-47   GATTACA
#   48-53   CCCCCC
#   54-99   A padding (fill)
#   100-103 ACGT
#   105     G
#   110     C
#   200-204 TAAAA
# The fill and overrides keep the sample VCF's REF alleles valid.
_PREFIX = "C" + "A" * 10 + "G" + "ACGT" * 5 + "T" * 8 + "GATTACA" + "C" * 6


def reference_seq(length=300):
    seq = ["A"] * length
    for i, ch in enumerate(_PREFIX, 0):
        seq[i] = ch
    for i, ch in enumerate("ACGT", 99):        # 100-103
        seq[i] = ch
    seq[104] = "G"                         # 105
    seq[109] = "C"                         # 110
    for i, ch in enumerate("TAAAA", 199):      # 200-204
        seq[i] = ch
    return "".join(seq)


REFSEQ = reference_seq()


def fasta(seq=REFSEQ, name="chr1", description="", records=1):
    head = ">" + name + ((" " + description) if description else "")
    body = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
    lines = [head, body]
    for k in range(1, records):
        lines += [f">chr{k + 1}", "ACGT"]
    return "\n".join(lines) + "\n"


def upload_fasta(client, body):
    return client.post(
        "/normalize",
        files={"file": ("reference.fasta", body if isinstance(body, bytes) else body.encode())},
    )


def upload_vcf(client, body, base_version):
    return client.post(
        "/vcf",
        params={"base_version": base_version},
        files={"file": ("variants.vcf", body if isinstance(body, bytes) else body.encode())},
    )


# --------------------------------------------------------------------- happy path

def test_normalize_returns_original_normalized_and_version(client):
    body = vcf([(100, "ACGT", "A"), (105, "G", "T"), (110, "C", "CT"), (200, "TAAAA", "T")])
    assert upload_vcf(client, body, 0).status_code == 201

    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["version"] == 1
    assert out["chrom"] == "chr1"
    assert out["count"] == 4

    originals = [rec for rec in out["records"]]
    assert [rec["original"] for rec in originals] == [
        "chr1\t100\t.\tACGT\tA\t.\t.\t.",
        "chr1\t105\t.\tG\tT\t.\t.\t.",
        "chr1\t110\t.\tC\tCT\t.\t.\t.",
        "chr1\t200\t.\tTAAAA\tT\t.\t.\t.",
    ]
    # designed to already be left-aligned / trimmed
    assert [(rec["pos"], rec["ref"], rec["alt"]) for rec in originals] == [
        (100, "ACGT", "A"),
        (105, "G", "T"),
        (110, "C", "CT"),
        (200, "TAAAA", "T"),
    ]
    for rec in originals:
        assert rec["original_pos"] == int(rec["original"].split("\t")[1])
        assert rec["key"] == f"chr1:{rec['pos']}:{rec['ref']}:{rec['alt']}"


def test_normalize_accepts_multiline_fasta_and_description(client):
    assert upload_vcf(client, vcf([(105, "G", "T")]), 0).status_code == 201
    r = upload_fasta(client, fasta(description="chr1 reference build test"))
    assert r.status_code == 200, r.text
    assert r.json()["records"][0]["pos"] == 105


def test_normalize_preserves_duplicate_records(client):
    body = vcf([(110, "C", "CT"), (110, "C", "CT")])
    assert upload_vcf(client, body, 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["count"] == 2
    assert [rec["original"] for rec in out["records"]] == [
        "chr1\t110\t.\tC\tCT\t.\t.\t.",
        "chr1\t110\t.\tC\tCT\t.\t.\t.",
    ]
    assert len({rec["key"] for rec in out["records"]}) == 1


# ------------------------------------------------------------------ repeat handling

DELETION_REPS = [(1, "CA", "C"), (2, "AA", "A"), (11, "AG", "G")]
INSERTION_REPS = [(1, "C", "CA"), (11, "A", "AA")]


def test_equivalent_deletions_in_homopolymer_get_one_key(client):
    assert upload_vcf(client, vcf(DELETION_REPS), 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["count"] == len(DELETION_REPS), "records must not be merged"
    normalized = [(rec["pos"], rec["ref"], rec["alt"]) for rec in out["records"]]
    assert normalized == [(1, "CA", "C")] * len(DELETION_REPS)
    assert len({rec["key"] for rec in out["records"]}) == 1


def test_equivalent_insertions_in_homopolymer_get_one_key(client):
    assert upload_vcf(client, vcf(INSERTION_REPS), 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["count"] == len(INSERTION_REPS)
    normalized = [(rec["pos"], rec["ref"], rec["alt"]) for rec in out["records"]]
    assert normalized == [(1, "C", "CA")] * len(INSERTION_REPS)
    assert len({rec["key"] for rec in out["records"]}) == 1


def test_equivalent_deletions_in_tandem_repeat_get_one_key(client):
    # Delete one "ACGT" unit, anchored at different units of the 13-32 repeat.
    reps = [(16, "TACGT", "T"), (20, "TACGT", "T"), (24, "TACGT", "T"), (28, "TACGT", "T")]
    assert upload_vcf(client, vcf(reps), 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    out = r.json()
    assert out["count"] == len(reps)
    normalized = {(rec["pos"], rec["ref"], rec["alt"]) for rec in out["records"]}
    assert normalized == {(12, "GACGT", "G")}


def test_right_aligned_insertion_shifts_left(client):
    assert upload_vcf(client, vcf([(11, "A", "AA")]), 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 200, r.text
    rec = r.json()["records"][0]
    assert (rec["pos"], rec["ref"], rec["alt"]) == (1, "C", "CA")


# ------------------------------------------------------------------ failure modes

def test_ref_mismatch_rejects_whole_request(client):
    # 110 is 'C' in the reference, so the second record's REF is wrong.
    body = vcf([(105, "G", "T"), (110, "A", "G")])
    assert upload_vcf(client, body, 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 400
    assert "reference mismatch" in r.json()["detail"]
    assert "does not match" in r.json()["detail"]
    assert "records" not in r.json(), "no partial result may be returned"


def test_ref_span_beyond_reference_end_rejected(client):
    seq = REFSEQ
    pos = len(seq) - 1
    assert upload_vcf(client, vcf([(pos, "AAA", "A")]), 0).status_code == 201
    r = upload_fasta(client, fasta())
    assert r.status_code == 400
    assert "beyond reference length" in r.json()["detail"]


def test_reference_chromosome_must_match_published(client):
    assert upload_vcf(client, vcf([(105, "G", "T")]), 0).status_code == 201
    r = upload_fasta(client, fasta(name="chr2"))
    assert r.status_code == 400
    assert "chr2" in r.json()["detail"] and "chr1" in r.json()["detail"]


def test_normalize_without_published_data_is_404(client):
    r = upload_fasta(client, fasta())
    assert r.status_code == 404
    assert "no VCF" in r.json()["detail"]


BAD_FASTAS = {
    "empty": "",
    "no_header": "ACGT\n",
    "multiple_records": fasta(records=2),
    "lowercase": fasta(seq="acgt"),
    "contains_n": fasta(seq="ACGTN"),
    "iupac": fasta(seq="ACGRT"),
    "whitespace_in_sequence": fasta(seq="AC GT"),
    "header_only": ">chr1\n",
    "empty_header": ">\nACGT\n",
    "sequence_before_header": "ACGT\n>chr1\nACGT\n",
}


@pytest.mark.parametrize("name", sorted(BAD_FASTAS))
def test_invalid_fasta_rejected(client, name):
    assert upload_vcf(client, vcf([(105, "G", "T")]), 0).status_code == 201
    r = upload_fasta(client, BAD_FASTAS[name])
    assert r.status_code == 400, f"{name}: {r.status_code} {r.text}"
    assert "invalid FASTA" in r.json()["detail"]
    # the published data is untouched
    assert client.get("/version").json()["current_version"] == 1


def test_oversize_fasta_rejected(client, app):
    limit = app.state.cfg.max_upload_bytes
    head = ">chr1\n"
    body = head + "A" * (limit + 1)
    assert upload_vcf(client, vcf([(105, "G", "T")]), 0).status_code == 201
    r = upload_fasta(client, body)
    assert r.status_code == 413
    assert client.get("/version").json()["current_version"] == 1


def test_binary_fasta_rejected(client):
    import gzip

    assert upload_vcf(client, vcf([(105, "G", "T")]), 0).status_code == 201
    r = upload_fasta(client, gzip.compress(fasta().encode()))
    assert r.status_code == 400
    assert "invalid FASTA" in r.json()["detail"]


# ------------------------------------------------------------------ pinning / concurrency

def test_normalize_pins_version_and_does_not_block_other_traffic(app, monkeypatch):
    """A slow normalize keeps its version; uploads/queries proceed meanwhile."""
    from fastapi.testclient import TestClient

    started = threading.Event()
    release = threading.Event()
    real = norm.normalize_dataset

    def slow_normalize(gz_path, reference):
        started.set()
        assert release.wait(10), "test did not release the slow normalizer"
        return real(gz_path, reference)

    monkeypatch.setattr(norm, "normalize_dataset", slow_normalize)

    # Context managers are required so lifespan (and thus the instance
    # heartbeat that keeps the in-flight pin alive) actually runs.
    with TestClient(app) as client1, TestClient(app) as client2:
        v1 = vcf([(105, "G", "T")])
        assert upload_vcf(client1, v1, 0).status_code == 201
        v1_dir = os.path.join(app.state.cfg.data_dir, app.state.store.version_info()["path"])

        results = {}

        def run_normalize():
            results["response"] = upload_fasta(client1, fasta())

        t = threading.Thread(target=run_normalize)
        t.start()
        assert started.wait(10), "normalize never started"

        try:
            # While normalization is stuck: interval query and a full
            # replacement (the correction) must both complete.
            q = client2.get("/variants", params={"chrom": "chr1", "start": 105, "end": 105})
            assert q.status_code == 200 and q.json()["version"] == 1
            v2 = vcf([(110, "C", "A")])
            assert upload_vcf(client2, v2, 1).status_code == 201
            assert client2.get("/version").json()["current_version"] == 2

            # The version being normalized is superseded but pinned: the
            # reaper must leave its data/index on disk for the in-flight work.
            app.state.reap_once()
            assert os.path.isdir(v1_dir)

            release.set()
            t.join(20)
            assert not t.is_alive()

            resp = results["response"]
            assert resp.status_code == 200, resp.text
            # Normalization fixed the version it started with, not the new one.
            assert resp.json()["version"] == 1
            assert [r["original"] for r in resp.json()["records"]] == [
                "chr1\t105\t.\tG\tT\t.\t.\t."
            ]
            # Current manifest is unchanged by normalization.
            assert client2.get("/version").json()["current_version"] == 2
        finally:
            release.set()
            t.join(20)
        # Once the pin is gone, the superseded version is reclaimed.
        app.state.reap_once()
        assert not os.path.exists(v1_dir)


def test_two_instances_concurrent_correction_while_normalizing(tmp_path):
    """Concurrent corrections across two instances: exactly one wins."""
    from fastapi.testclient import TestClient

    app1 = create_app(make_config(tmp_path, instance_id="inst-1"))
    app2 = create_app(make_config(tmp_path, instance_id="inst-2"))
    c1, c2 = TestClient(app1), TestClient(app2)
    try:
        assert upload_vcf(c1, vcf([(105, "G", "T")]), 0).status_code == 201
        # normalizing on instance 1 while both instances race to correct v1
        baseline = upload_fasta(c1, fasta())
        assert baseline.status_code == 200 and baseline.json()["version"] == 1

        barrier = threading.Barrier(2)
        codes = []
        lock = threading.Lock()

        def correct(client, body):
            barrier.wait()
            r = upload_vcf(client, body, 1)
            with lock:
                codes.append(r.status_code)

        threads = [
            threading.Thread(target=correct, args=(c1, vcf([(110, "C", "A")]))),
            threading.Thread(target=correct, args=(c2, vcf([(200, "TAAAA", "T")]))),
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert sorted(codes) == [201, 409]
        assert c1.get("/version").json()["current_version"] == 2

        # a subsequent normalize sees the winning manifest and its version
        after = upload_fasta(c2, fasta())
        assert after.status_code == 200
        assert after.json()["version"] == 2
    finally:
        c1.close()
        c2.close()


# ------------------------------------------------------------------ docs / schema

def test_openapi_and_docs_expose_normalize(client):
    spec = client.get("/openapi.json").json()
    op = spec["paths"]["/normalize"]["post"]
    assert "multipart/form-data" in op["requestBody"]["content"]
    assert "200" in op["responses"] and "400" in op["responses"]
    schema = spec["components"]["schemas"]["NormalizeResponse"]
    assert {"version", "chrom", "count", "records"} <= set(schema["properties"])
    docs = client.get("/docs")
    assert docs.status_code == 200
    assert "swagger-ui" in docs.text
