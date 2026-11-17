"""FastAPI application: versioned VCF upload + indexed interval queries.

Layout of the shared data directory:

    <data_dir>/meta.sqlite3          metadata, current-version pointer, ref pins
    <data_dir>/versions/<uuid>/      published file sets (variants.vcf.gz + .tbi)
    <data_dir>/staging/<instance>/<uuid>/   in-progress builds

Publishing is atomic: the file set is completed in staging, renamed into
versions/, and only then the current-version pointer is flipped inside one
SQLite transaction. Readers pin a version before opening its files, so an
in-flight query keeps its (version, data, index) trio consistent while a
replacement is published; the reaper frees old files once no pin remains.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import time
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass

import pysam
from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import RedirectResponse
from pydantic import BaseModel, Field

from . import normalize
from .db import Conflict, Store
from .vcfbuild import GZ_NAME, VCFError, validate_and_build

DEFAULT_MAX_UPLOAD = 10 * 1024 * 1024  # 10 MiB


class NormalizedRecord(BaseModel):
    original: str = Field(description="Raw VCF line as published, unmodified")
    original_pos: int = Field(description="POS as published (1-based)")
    original_ref: str = Field(description="REF as published")
    original_alt: str = Field(description="ALT as published")
    pos: int = Field(description="Normalized 1-based POS")
    ref: str = Field(description="Normalized REF")
    alt: str = Field(description="Normalized ALT")
    key: str = Field(description="Canonical key: chrom:pos:ref:alt")


class NormalizeResponse(BaseModel):
    version: int = Field(description="Published version that was normalized (pinned)")
    chrom: str = Field(description="Chromosome shared by the reference and dataset")
    count: int = Field(description="Number of records; one output per original record")
    records: list[NormalizedRecord]


async def _save_upload(file: UploadFile, dest: str, max_bytes: int) -> None:
    """Stream an upload to `dest`, rejecting anything over the byte limit."""
    size = 0
    with open(dest, "wb") as out:
        while True:
            chunk = await file.read(1 << 20)
            if not chunk:
                break
            size += len(chunk)
            if size > max_bytes:
                raise HTTPException(
                    413, f"upload exceeds the {max_bytes}-byte limit"
                )
            out.write(chunk)


@dataclass
class Config:
    data_dir: str
    instance_id: str
    max_upload_bytes: int = DEFAULT_MAX_UPLOAD
    reaper_interval: float = 2.0      # seconds between reaper runs (0 disables)
    heartbeat_interval: float = 2.0   # seconds between heartbeats (0 disables)
    dead_after: float = 15.0          # instance considered dead after this silence
    orphan_grace: float = 60.0        # min age before orphan dirs are swept


def _config_from_env() -> Config:
    e = os.environ
    return Config(
        data_dir=e.get("VCF_DATA_DIR", "/data"),
        instance_id=e.get("VCF_INSTANCE_ID") or uuid.uuid4().hex,
        max_upload_bytes=int(e.get("VCF_MAX_UPLOAD_BYTES", DEFAULT_MAX_UPLOAD)),
        reaper_interval=float(e.get("VCF_REAPER_INTERVAL", "2")),
        heartbeat_interval=float(e.get("VCF_HEARTBEAT_INTERVAL", "2")),
        dead_after=float(e.get("VCF_DEAD_AFTER", "15")),
        orphan_grace=float(e.get("VCF_ORPHAN_GRACE", "60")),
    )


def create_app(config: Config | None = None) -> FastAPI:
    cfg = config or _config_from_env()
    store = Store(os.path.join(cfg.data_dir, "meta.sqlite3"))
    versions_dir = os.path.join(cfg.data_dir, "versions")
    staging_base = os.path.join(cfg.data_dir, "staging")
    staging_root = os.path.join(staging_base, cfg.instance_id)
    os.makedirs(versions_dir, exist_ok=True)
    # Leftover builds of a previous incarnation of this instance can never be
    # published any more; drop them before serving.
    shutil.rmtree(staging_root, ignore_errors=True)
    os.makedirs(staging_root, exist_ok=True)

    def reap_once() -> None:
        now = time.time()
        # 1) DB: drop dead instances' pins, delete reclaimable version rows.
        for _version, rel in store.reap_metadata(now, cfg.dead_after):
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        # 2) FS: sweep version dirs with no DB row (crashed publish/reaper).
        cutoff = now - cfg.orphan_grace
        known = store.all_version_paths()
        for name in os.listdir(versions_dir):
            path = os.path.join(versions_dir, name)
            rel = os.path.relpath(path, cfg.data_dir)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if rel not in known and old:
                shutil.rmtree(path, ignore_errors=True)
        # 3) FS: sweep staging dirs of dead instances.
        live = store.live_instances(now, cfg.dead_after)
        for inst in os.listdir(staging_base):
            path = os.path.join(staging_base, inst)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if inst not in live and old:
                shutil.rmtree(path, ignore_errors=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        store.heartbeat(cfg.instance_id, time.time())
        stopped = asyncio.Event()

        async def heartbeat_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.heartbeat_interval)
                try:
                    store.heartbeat(cfg.instance_id, time.time())
                except Exception:
                    pass

        async def reaper_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.reaper_interval)
                try:
                    await asyncio.to_thread(reap_once)
                except Exception:
                    pass

        tasks = []
        if cfg.heartbeat_interval > 0:
            tasks.append(asyncio.create_task(heartbeat_loop()))
        if cfg.reaper_interval > 0:
            tasks.append(asyncio.create_task(reaper_loop()))
        try:
            yield
        finally:
            stopped.set()
            for t in tasks:
                t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            try:
                store.unregister(cfg.instance_id)
            except Exception:
                pass

    app = FastAPI(
        title="Versioned VCF variant store",
        version="1.0.0",
        lifespan=lifespan,
        description=(
            "Upload single-chromosome biallelic VCFs (max 10 MiB, no sample "
            "columns, REF/ALT in ACGT); records are sorted, BGZF-compressed and "
            "Tabix-indexed. Query 1-based closed intervals; every record whose "
            "reference span [POS, POS+len(REF)-1] intersects the interval is "
            "returned via the index. Uploads are optimistic-concurrency "
            "controlled through `base_version`. `POST /normalize` left-aligns "
            "and trims the published records against an uploaded single-record "
            "ACGT FASTA reference, pinning the published version for the whole "
            "computation. Interactive docs for both operations are served under "
            "`/docs`."
        ),
    )
    app.state.cfg = cfg
    app.state.store = store
    app.state.reap_once = reap_once

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse("/docs")

    @app.post("/vcf", status_code=201)
    async def upload_vcf(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Version the client last saw (0 when empty)"),
    ):
        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = os.path.join(staging, "upload.vcf")
            await _save_upload(file, upload_path, cfg.max_upload_bytes)
            try:
                chrom, n_records = await asyncio.to_thread(
                    validate_and_build, upload_path, staging
                )
            except VCFError as exc:
                raise HTTPException(400, f"invalid VCF: {exc}")
            os.remove(upload_path)

            # File set is complete; move it into place, then flip the pointer.
            final_rel = os.path.join("versions", uuid.uuid4().hex)
            os.rename(staging, os.path.join(cfg.data_dir, final_rel))
            try:
                new_version = await asyncio.to_thread(
                    store.publish, base_version, chrom, final_rel, n_records, time.time()
                )
            except Conflict as c:
                shutil.rmtree(os.path.join(cfg.data_dir, final_rel), ignore_errors=True)
                raise HTTPException(
                    409,
                    f"base_version {base_version} does not match current version "
                    f"{c.current}; re-read and retry",
                )
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
        return {"version": new_version, "chrom": chrom, "n_records": n_records}

    @app.post(
        "/normalize",
        response_model=NormalizeResponse,
        responses={
            400: {"description": "Invalid FASTA, chromosome mismatch, or REF/reference mismatch"},
            404: {"description": "No VCF has been published yet"},
            413: {"description": "Upload exceeds the byte limit"},
        },
    )
    async def normalize_vcf(file: UploadFile = File(...)):
        """Normalize the currently published VCF against an uploaded reference.

        The published version is pinned before reading, so a publish that
        happens while this request is computing does not change the dataset
        being normalized: the response names the version actually used and the
        current manifest is never modified.  The heavy work runs in a worker
        thread, so uploads and interval queries stay responsive.
        """
        pinned = await asyncio.to_thread(store.pin, cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned

        staging = os.path.join(staging_root, uuid.uuid4().hex)
        try:
            os.makedirs(staging)
            fasta_path = os.path.join(staging, "reference.fasta")
            await _save_upload(file, fasta_path, cfg.max_upload_bytes)
            try:
                reference = await asyncio.to_thread(normalize.read_fasta, fasta_path)
            except normalize.FastaError as exc:
                raise HTTPException(400, f"invalid FASTA: {exc}")

            if reference.chrom != vchrom:
                raise HTTPException(
                    400,
                    f"reference chromosome '{reference.chrom}' does not match the "
                    f"published dataset chromosome '{vchrom}'",
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                records = await asyncio.to_thread(
                    normalize.normalize_dataset, gz_path, reference
                )
            except normalize.NormalizationError as exc:
                # Whole request fails; no partial result is returned.
                raise HTTPException(400, f"reference mismatch: {exc}")
            except RuntimeError as exc:
                raise HTTPException(500, f"cannot normalize published VCF: {exc}")
        finally:
            shutil.rmtree(staging, ignore_errors=True)
            await asyncio.to_thread(store.unpin, ref_id)

        return {
            "version": version,
            "chrom": vchrom,
            "count": len(records),
            "records": records,
        }

    @app.get("/variants")
    def get_variants(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                tb = pysam.TabixFile(gz_path)
                try:
                    # 1-based closed [start, end] -> 0-based half-open fetch.
                    lines = tb.fetch(vchrom, start - 1, end)
                    records = []
                    for ln in lines:
                        f = ln.split("\t", 4)
                        pos = int(f[1])
                        # Exact span check on top of the index iterator.
                        if pos <= end and pos + len(f[3]) - 1 >= start:
                            records.append(ln)
                finally:
                    tb.close()
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "version": version,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin(ref_id)

    @app.get("/version")
    def get_version():
        info = store.version_info()
        if info is None:
            return {"published": False, "current_version": 0}
        return {
            "published": True,
            "current_version": info["version"],
            "chrom": info["chrom"],
            "n_records": info["n_records"],
            "created_at": info["created_at"],
        }

    return app


# Run with: uvicorn --factory app.main:create_app
