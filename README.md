# Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement and safe concurrent access from two app instances. It additionally
normalizes the published records against an uploaded reference FASTA (indel
left-alignment and redundant-base trimming) without modifying the published
manifest.

## Run

```bash
docker compose up --build
```

Two instances come up sharing one data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published data; uploads may go to either one.

## API

### `POST /vcf?base_version=<n>` (multipart file field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T` (uppercase; no `N`, no
  symbolic alleles, no `*`, no multi-allelic `A,G`).

Records are sorted by position, BGZF-compressed and Tabix-indexed. The upload
carries the version the client last saw (`base_version`, `0` when the store is
empty). If two replacements race — even across the two instances — exactly one
wins; the loser gets `409` and must re-read `/version` and retry.
Returns `201` with `{"version", "chrom", "n_records"}`.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the
`version` that was actually read; that version's data and index are pinned for
the duration of the query, so version, data and index always match.

### `GET /version`

Current published version, chromosome and record count.

### `POST /normalize` (multipart file field `file`)

Normalizes the **currently published** VCF records against an uploaded
reference FASTA. The request does not change the published manifest or
version — it is a read of the current version plus a computation.

The FASTA must satisfy:

- at most 10 MiB (same `VCF_MAX_UPLOAD_BYTES` limit as `/vcf`),
- exactly one record (single chromosome),
- sequence containing only uppercase `A/C/G/T`,
- its record name (first whitespace-delimited token of the header) equal to the
  published chromosome.

Every published record's `REF` is checked against the reference sequence at
`POS`. If any record's `REF` does not match, or its span
`[POS, POS+len(REF)-1]` runs past the end of the reference, the request fails
as a whole with `400`; no partial results are returned. Records are then
left-aligned and stripped of redundant common bases using the reference, so
equivalent representations of the same indel in a repeat region produce the
same `key`. Records are **not** merged or deduplicated — every original record
appears exactly once in the output, in published order.

The published version is pinned for the whole computation. A correction
uploaded while normalization is running (from either instance) does not affect
the result; the response reports the version that was actually normalized.
The work runs in a worker thread, so uploads and `/variants` queries keep
being served meanwhile.

Response:

```json
{
  "version": 3,
  "chrom": "chr1",
  "count": 4,
  "records": [
    {
      "original": "chr1\t11\t.\tA\tAA\t.\t.\t.",
      "original_pos": 11,
      "original_ref": "A",
      "original_alt": "AA",
      "pos": 1,
      "ref": "C",
      "alt": "CA",
      "key": "chr1:1:C:CA"
    }
  ]
}
```

`original` is the raw VCF line as published; `pos`/`ref`/`alt` are the
normalized 1-based alleles; `key` (`chrom:pos:ref:alt`) is the canonical key
for cross-manifest comparison. Errors: `404` if nothing is published, `400` for
an invalid FASTA, chromosome mismatch, or any reference mismatch, `413` when
the upload is over the size limit.

## Consistency design

- **Atomic publish**: a new version is fully built (sorted, compressed,
  indexed) in a per-instance staging directory, moved into `versions/`, and
  only then the current-version pointer in SQLite is flipped inside a single
  `BEGIN IMMEDIATE` transaction. Readers never see a half-built file set.
- **Optimistic concurrency**: the pointer flip checks `base_version` inside
  the same transaction, so concurrent replacements from both instances are
  serialized and only one succeeds.
- **Old queries finish**: a query pins its version in a `refs` table before
  opening files. A background reaper deletes a superseded version's data and
  index only when no pin references it. Instances heartbeat; pins of crashed
  instances expire and never block reclamation. `/normalize` uses the same
  pinning, so it always reads a consistent (version, data, index) trio even
  while a correction is published, and never blocks the high-traffic paths.
- **Crash recovery**: if the process dies mid-build or mid-publish, the
  staging/orphan directories left behind are swept (with a grace period) and
  the previously published version keeps serving. The pointer flip is the last
  step, so a restart never exposes a partial file set.

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 71 tests
uvicorn --factory app.main:create_app   # single dev instance (uses /data or $VCF_DATA_DIR)
```

Interactive docs (Swagger UI) are served on each instance at `/docs`:
http://localhost:8001/docs and http://localhost:8002/docs. `/` redirects to
`/docs`; the machine-readable schema is at `/openapi.json`.

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`, `VCF_HEARTBEAT_INTERVAL`,
`VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
