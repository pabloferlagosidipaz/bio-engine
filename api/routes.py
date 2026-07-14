# SPDX-License-Identifier: MIT
"""
API Routes Definitions
=======================

This module defines the RESTful endpoints for the bio-engine sidecar application.
It exposes functionalities for job lifecycle management, annotation, alignment preview,
and biological database communication (e.g., fetching references, HGVS variants).
"""

import logging
import os
import re
import shutil

from fastapi import APIRouter, BackgroundTasks, Depends, File, UploadFile
from sqlalchemy import Integer, cast, func
from sqlalchemy.orm import Session

from core.cache import cache
from core.database import Base, engine, get_db
from core.exceptions import BioEngineError
from core.paths import resolve_within, safe_filename_join
from core.security import require_trusted_origin
from data.models import (
    AddCommentRequest,
    AddHGVSAlternativesRequest,
    ApprovedVariantResponse,
    ApproveVariantRequest,
    CreateJobRequest,
    HGVSRequest,
    HotspotPoint,
    ImportJobRequest,
    Job,
    JobStatus,
    OCInstallTask,
    OCModule,
    OCStatus,
    OpenCravatConfigRequest,
    ProxyConfigRequest,
    RenameJobRequest,
    ShareJobRequest,
    UpdateJobRequest,
    UpdateVariantStatusRequest,
)
from data.models_db import ApprovedVariant
from services import aligner as aligner_service
from services import opencravat as oc_service
from services import reference as ref_service
from services.job_manager import JobManager
from tasks.worker import annotate_hgvs_background, process_job_background

# Create tables on startup (simple approach for now)
Base.metadata.create_all(bind=engine)

logger = logging.getLogger(__name__)

router = APIRouter(dependencies=[Depends(require_trusted_origin)])
job_manager = JobManager()

# Matches an absolute genomic HGVS/SPDI identifier, e.g. "NC_000019.10:g.38494330C>T"
# or its SPDI form "19:38494329:C:T". Used to check whether a candidate
# already carries genomic (not just transcript-relative) coordinates.
_GENOMIC_HGVS_RE = r"(?:NC_0000\d{2}|NC_012920|(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT|M))[:.](?:g\.)?\d+"
# Same shape as above, but with capture groups for chromosome/position extraction.
#
# The separator before the position is `(?:\.\d+)?:` rather than a bare
# `[:.]`: a versioned accession (e.g. "NC_000017.11:g.43045712T>C") has a
# ".VERSION" segment before the ":g." part, and a bare `[:.]` would match
# that "." first, capturing the version number ("11") as if it were the
# position instead of the real position (43045712). Requiring the
# optional ".VERSION" to be consumed before the mandatory ":" fixes this
# while leaving unversioned HGVS and SPDI forms (which have no "." here)
# unaffected.
_GENOMIC_HGVS_CAPTURE_RE = r"(?:NC_0000(\d{2})|(NC_012920)|(?:chr)?(\d{1,2}|X|Y|MT|M))(?:\.\d+)?:(?:g\.)?(\d+)"


def _resolve_genomic_variant_from_hgvs(job: Job, job_id: str, hgvs_key: str, row: list, col_map: dict) -> dict | None:
    """
    Resolves absolute genomic (chromosome, position) coordinates for a
    variant's HGVS key: checks its cached alternatives first, falling back
    to an on-the-fly Ensembl lookup if none are genomic yet, then parses
    the first genomic HGVS/SPDI candidate found.

    Shared by both variant_key forms handled in update_variant_status
    (patient-scoped row index and grouped genomic position); this logic
    was previously duplicated near-verbatim between the two branches.
    """
    alts = job.hgvs_alternatives.get(hgvs_key, [])
    has_genomic = any(re.match(_GENOMIC_HGVS_RE, a) for a in [hgvs_key] + alts)

    if not has_genomic:
        try:
            from utilities.ensembl_hgvs import EnsemblHGVS
            logger.info(f"Sync: Fetching genomic alternatives for {hgvs_key} on-the-fly...")
            ensembl = EnsemblHGVS(assembly=job.hgvs_config.assembly if job.hgvs_config else "GRCh38")
            results = ensembl.get_equivalents_batch([hgvs_key])
            new_alts = results.get(hgvs_key, [])
            if new_alts:
                alts = list(set(alts + new_alts))
                job_manager.add_job_hgvs_alternatives(job_id, hgvs_key, alts)
        except Exception as e:
            logger.warning(f"On-the-fly HGVS fetch failed for {hgvs_key}: {e}")

    candidates = [hgvs_key] + alts
    logger.info(f"Sync: Checking {len(candidates)} HGVS/SPDI candidates for {hgvs_key}")
    for alt in candidates:
        # Match both NC_0000XX.X:g.XXXX (HGVS) and NC_0000XX.X:XXXX:A:T (SPDI)
        # And chromosome numbers like 19:g.XXXX or 19:XXXX:A:T
        match = re.match(_GENOMIC_HGVS_CAPTURE_RE, alt)
        if match:
            chrom = match.group(1) or (match.group(2) and "MT") or match.group(3)
            if chrom == "M": chrom = "MT"
            if chrom.isdigit(): chrom = str(int(chrom))

            pos = int(match.group(4))
            # SPDI is 0-based, HGVS is 1-based.
            # SPDI usually looks like Ac:Pos:Ref:Alt (3 colons)
            # HGVS looks like Ac:g.PosRef>Alt
            if ":g." not in alt and alt.count(":") >= 2:
                pos += 1 # Normalize SPDI 0-based to 1-based

            target_variant = {
                "chromosome": chrom,
                "position": pos,
                "ref": row[col_map.get("REF", 0)],
                "alt": row[col_map.get("ALT", 0)],
                "gene": row[col_map.get("GENE", 0)] if "GENE" in col_map else job.hgvs_config.gene if job.hgvs_config else None
            }
            logger.info(f"Sync: Resolved to Chr {target_variant['chromosome']} Pos {target_variant['position']} from {alt}")
            return target_variant

    return None


def _resolve_genomic_variant_from_reference(job: Job, row: list, col_map: dict, position: int) -> dict | None:
    """
    Fallback used when HGVS-based resolution didn't produce coordinates: if
    the job's reference is itself a genomic NC_ accession, derive the
    chromosome directly from it (position is already known from the row).

    Shared by both variant_key forms in update_variant_status; previously
    duplicated near-verbatim between the two branches.
    """
    if job.reference.get("type") != "ncbi":
        return None

    ref_val = job.reference.get("value", "")
    match = re.search(r"NC_0000(\d{2})", ref_val)
    if not match:
        return None

    target_variant = {
        "chromosome": str(int(match.group(1))),
        "position": position,
        "ref": row[col_map.get("REF", 0)],
        "alt": row[col_map.get("ALT", 0)],
        "gene": row[col_map.get("GENE", 0)] if "GENE" in col_map else job.hgvs_config.gene if job.hgvs_config else None
    }
    logger.info(f"Sync: Resolved via Job Reference to Chr {target_variant['chromosome']} Pos {target_variant['position']}")
    return target_variant

@router.post("/jobs/{job_id}/comments", response_model=Job)
def add_comment(job_id: str, request: AddCommentRequest):
    """
    Adds a user comment to a specific variant in a job.

    Args:
        job_id (str): The unique identifier of the job.
        request (AddCommentRequest): The details of the comment.

    Returns:
        Job: The updated job object containing the new comment.
    """
    job = job_manager.add_job_comment(
        job_id,
        variant_key=request.variant_key,
        text=request.text,
        author=request.author
    )
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to add comment")
    return job

@router.put("/jobs/{job_id}/variant-status", response_model=Job)
def update_variant_status(job_id: str, request: UpdateVariantStatusRequest, db: Session = Depends(get_db)):
    """
    Updates the status of a specific variant and synchronizes the global hotspots map.
    """
    # 1. Update the job file status
    job = job_manager.update_variant_status(
        job_id,
        variant_key=request.variant_key,
        status=request.status
    )
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to update variant status")

    # 2. Sync with approved_variants DB
    # The variant_key can be "patient_id:index" or just "position" (for grouped variants)
    patient_id, variant_id = request.variant_key.split(':') if ':' in request.variant_key else (None, request.variant_key)

    # Extract variant details from results
    target_variant = None
    if job.results:
        for res in job.results:
            # If patient_id is specified, only look in that patient's results
            if patient_id and res.get("id") != patient_id:
                continue

            alignment_data = res.get("alignment", {})
            v_data = alignment_data.get("variants", {})
            rows, cols = v_data.get("rows", []), v_data.get("columns", [])
            col_map = {col.upper(): i for i, col in enumerate(cols)}

            try:
                # If variant_id is a row index (usually when patient_id is present)
                if patient_id:
                    v_idx = int(variant_id)
                    if v_idx < len(rows):
                        row = rows[v_idx]
                        # Priority: Extract both Chr and Pos from HGVS equivalents (Absolute Genomic Coordinates)
                        hgvs_key = row[col_map["HGVS"]] if "HGVS" in col_map else None
                        if hgvs_key:
                            target_variant = _resolve_genomic_variant_from_hgvs(job, job_id, hgvs_key, row, col_map)

                        # Fallback: if reference is genomic NC_...
                        if not target_variant:
                            target_variant = _resolve_genomic_variant_from_reference(
                                job, row, col_map, int(row[col_map.get("POS", 0)])
                            )

                        if not target_variant:
                            logger.warning(f"Variant {request.variant_key} skipped: Could not resolve genomic coordinates.")
                        break
                else:
                    # If variant_id is a position string (grouped variant)
                    target_pos = int(variant_id)
                    pos_idx = col_map.get("POS", 0)
                    for row in rows:
                        if int(row[pos_idx]) == target_pos:
                            # Extract both Chr and Pos from HGVS (Absolute Genomic Coordinates)
                            hgvs_key = row[col_map["HGVS"]] if "HGVS" in col_map else None
                            if hgvs_key:
                                target_variant = _resolve_genomic_variant_from_hgvs(job, job_id, hgvs_key, row, col_map)

                            # Fallback: if reference is genomic NC_...
                            if not target_variant:
                                target_variant = _resolve_genomic_variant_from_reference(job, row, col_map, target_pos)

                            if target_variant: break

                    if not target_variant:
                        logger.warning(f"Grouped variant at pos {target_pos} in job {job_id} skipped: Could not resolve genomic coordinates.")
                    break
            except Exception as e:
                logger.error(f"Error resolving variant {request.variant_key}: {e}")
                continue

    if target_variant:
        if request.status == "approved":
            # Add to DB if missing
            exists = db.query(ApprovedVariant).filter(
                ApprovedVariant.job_id == job_id,
                ApprovedVariant.chromosome == target_variant["chromosome"],
                ApprovedVariant.position == target_variant["position"]
            ).first()
            if not exists:
                db_variant = ApprovedVariant(
                    chromosome=target_variant["chromosome"],
                    position=target_variant["position"],
                    ref_allele=target_variant["ref"],
                    alt_allele=target_variant["alt"],
                    gene=target_variant["gene"],
                    approved_by="Dashboard",
                    job_id=job_id,
                    patient_id=patient_id,
                    assembly=job.hgvs_config.assembly if job.hgvs_config else "GRCh38"
                )
                db.add(db_variant)
        else:
            # Remove from DB if status changed from approved to anything else
            db.query(ApprovedVariant).filter(
                ApprovedVariant.job_id == job_id,
                ApprovedVariant.chromosome == target_variant["chromosome"],
                ApprovedVariant.position == target_variant["position"]
            ).delete()

        db.commit()

    return job

@router.delete("/jobs/{job_id}/comments/{variant_key}/{comment_id}", response_model=Job)
def delete_comment(job_id: str, variant_key: str, comment_id: str):
    """
    Deletes a specific comment from a variant inside a job.

    Args:
        job_id (str): The unique identifier of the job.
        variant_key (str): The identifier key of the variant.
        comment_id (str): The ID of the comment to delete.

    Returns:
        Job: The updated job object.
    """
    job = job_manager.delete_job_comment(job_id, variant_key, comment_id)
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to delete comment")
    return job

@router.get("/")
def health_check():
    """Returns the basic health status of the engine."""
    return {"status": "online", "engine": "bio-engine"}

@router.post("/config/proxy")
def configure_proxy(request: ProxyConfigRequest):
    """
    Dynamically configures HTTP/HTTPS proxy settings for the engine.
    """
    # 1. Update global settings
    if request.http_proxy is not None:
        from core.config import settings
        settings.http_proxy = request.http_proxy if request.http_proxy != "" else None

    if request.https_proxy is not None:
        from core.config import settings
        settings.https_proxy = request.https_proxy if request.https_proxy != "" else None

    # 2. Update environment and internal state
    from core.proxy_manager import proxy_manager
    proxy_manager.refresh_proxies()

    return {
        "status": "success",
        "http_proxy": os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy"),
        "https_proxy": os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    }

@router.post("/config/opencravat")
def configure_opencravat(request: OpenCravatConfigRequest):
    """
    Dynamically configures the global OpenCRAVAT executable path for the engine.
    """
    from core.config import settings
    if request.oc_path is not None:
        new_path = request.oc_path if request.oc_path != "" else "oc"
        oc_service.resolve_oc_executable(new_path)
        settings.oc_path = new_path
    return {"status": "success", "oc_path": settings.oc_path}

@router.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    """
    Uploads a file to the server's local storage.
    Returns the absolute path to the uploaded file.
    """
    from core.config import settings

    file_path = safe_filename_join(settings.uploads_dir, file.filename)

    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)

        logger.info(f"File uploaded successfully: {file_path}")
        return {"path": file_path}
    except Exception as e:
        logger.error(f"Failed to upload file {file.filename}: {e}")
        raise BioEngineError("Failed to upload file")

@router.get("/check-reference")
def check_reference(id: str):
    """
    Checks if a given reference sequence ID exists and is accessible.

    Args:
        id (str): The NCBI accession or the local identifier to check.
    """
    exists = ref_service.check_ncbi_reference_exists(id)
    return {"exists": exists}

@router.get("/search-reference")
def search_reference(query: str, assembly: str | None = None):
    """
    Searches NCBI nucleotide database for a given query (gene name or accession).
    If assembly is provided, it attempts to prioritize transcripts for that build.
    """
    results = ref_service.search_reference(query, assembly=assembly)
    return {"results": results}

@router.get("/preview-read")
def preview_read(path: str):
    """Returns trace data for a given file path.

    `path` must resolve inside uploads_dir or jobs_dir - the only two
    places read files legitimately live (freshly uploaded, or copied in by
    a job import) - to prevent it being used to read arbitrary files off
    disk (CWE-22).
    """
    from core.config import settings

    safe_path = resolve_within(path, [settings.uploads_dir, settings.jobs_dir])

    try:
        return aligner_service.get_read_preview(safe_path)
    except BioEngineError:
        raise
    except Exception as e:
        logger.error(f"Preview failed for {path}: {e}")
        raise BioEngineError("Preview failed")

@router.post("/create-job", response_model=Job)
def create_job(request: CreateJobRequest):
    """
    Creates a new analysis job based on the requested configuration.

    If an NCBI reference is provided, it attempts to load its FASTA sequence
    before starting the job manager creation phase.

    Args:
        request (CreateJobRequest): Payload with reference, patients, and optional configs.

    Returns:
        Job: Detailed job data including its initial queued status.
    """
    # Resolve reference sequence
    ref_input = request.reference.get("value")
    assembly = request.hgvs_config.assembly if request.hgvs_config else "GRCh38"
    sequence = None
    try:
        if ref_input:
            ref_path = ref_service.load_reference(ref_input, assembly=assembly)
            sequence = ref_service.get_fasta_sequence(ref_path, assembly=assembly)
    except Exception as e:
        logger.error(f"Failed to load reference sequence for {ref_input}: {e}")

    return job_manager.create_job(
        name=request.name,
        reference=request.reference,
        patients=request.patients,
        reference_sequence=sequence,
        app_version=request.app_version,
        config=request.config.model_dump() if request.config else None,
        hgvs_config=request.hgvs_config.model_dump() if request.hgvs_config else (
            {"transcript": ref_input} if request.reference.get("type") == "ncbi" else None
        )
    )

@router.get("/jobs", response_model=list[Job])
def list_jobs():
    """
    Retrieves a list of all jobs currently tracked by the engine.
    """
    return job_manager.list_jobs()

@router.get("/jobs/{job_id}", response_model=Job)
def get_job(job_id: str):
    """
    Fetches the details of a single job.

    Args:
        job_id (str): The ID of the job to retrieve.
    """
    job = job_manager.get_job(job_id)
    if not job:
        raise BioEngineError(f"Job {job_id} not found")
    return job

@router.delete("/jobs/{job_id}")
def delete_job(job_id: str, db: Session = Depends(get_db)):
    """
    Removes a job from the persistent storage and cleans up associated approved variants.
    """
    # 1. Cleanup approved variants from DB
    db.query(ApprovedVariant).filter(ApprovedVariant.job_id == job_id).delete()
    db.commit()

    # 2. Delete job files
    success = job_manager.delete_job(job_id)
    if not success:
        raise BioEngineError(f"Failed to delete job {job_id} (or not found)")
    return {"status": "deleted", "id": job_id}

@router.put("/jobs/{job_id}/rename", response_model=Job)
def rename_job(job_id: str, request: RenameJobRequest):
    """
    Updates the displayed name of a specific job.
    """
    job = job_manager.rename_job(job_id, request.name)
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to rename")
    return job

@router.put("/jobs/{job_id}", response_model=Job)
def update_job(job_id: str, request: UpdateJobRequest):
    """
    Full update of an existing job configuration.

    Args:
        job_id (str): The ID of the job to update.
        request (UpdateJobRequest): The new configuration parameters for the job.
    """
    # Resolve reference sequence
    ref_input = request.reference.get("value")
    assembly = request.hgvs_config.assembly if request.hgvs_config else "GRCh38"
    sequence = None
    try:
        if ref_input:
            ref_path = ref_service.load_reference(ref_input, assembly=assembly)
            sequence = ref_service.get_fasta_sequence(ref_path, assembly=assembly)
    except Exception as e:
        logger.error(f"Failed to load reference sequence for {ref_input}: {e}")

    job = job_manager.update_job(
        job_id,
        name=request.name,
        reference=request.reference,
        patients=request.patients,
        config=request.config.model_dump() if request.config else None,
        hgvs_config=request.hgvs_config.model_dump() if request.hgvs_config else None,
        reference_sequence=sequence
    )
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to update")
    return job

@router.post("/run-job/{job_id}", response_model=Job)
def run_job(job_id: str, background_tasks: BackgroundTasks):
    """
    Starts the analysis for the given job ID in the background.
    Returns the job object immediately.
    """
    job = job_manager.get_job(job_id)
    if not job:
        raise BioEngineError(f"Job {job_id} not found")

    # Update status to RUNNING
    job_manager.update_job_status(job_id, JobStatus.RUNNING)
    job_manager.update_job_progress(job_id, 0, "Initializing job...")

    # Start background task
    background_tasks.add_task(process_job_background, job_id)

    return job_manager.get_job(job_id)

@router.post("/jobs/{job_id}/annotate-hgvs", response_model=Job)
def annotate_hgvs(job_id: str, background_tasks: BackgroundTasks):
    """
    Triggers manual HGVS annotation for the given job.
    """
    job = job_manager.get_job(job_id)
    if not job:
        raise BioEngineError(f"Job {job_id} not found")

    # Update status to RUNNING so frontend polls
    job_manager.update_job_status(job_id, JobStatus.RUNNING)

    # Start background task
    background_tasks.add_task(annotate_hgvs_background, job_id)

    return job_manager.get_job(job_id)

@router.post("/jobs/{job_id}/hgvs-alternatives", response_model=Job)
def add_hgvs_alternatives(job_id: str, request: AddHGVSAlternativesRequest):
    job = job_manager.add_job_hgvs_alternatives(
        job_id,
        principal_hgvs=request.principal_hgvs,
        alternatives=request.alternatives
    )
    if not job:
        raise BioEngineError(f"Job {job_id} not found or failed to add HGVS alternatives")
    return job

@router.post("/tools/hgvs/alternatives", response_model=list[str])
def hgvs_alternatives(request: HGVSRequest):
    """
    Returns full HGVS alternatives for a specific variant.
    """
    try:
        from utilities.ensembl_hgvs import EnsemblHGVS

        ensembl = EnsemblHGVS(assembly=request.assembly)

        # Determine HGVS type (c, g, or n)
        if request.transcript.startswith(('NM_', 'XM_')): h_type = 'c'
        elif request.transcript.startswith(('NC_', 'NG_', 'NW_', 'NT_')): h_type = 'g'
        else: h_type = 'n'

        # Generate primary ID
        primary = EnsemblHGVS.format_hgvs(
            request.transcript,
            h_type,
            request.pos,
            request.ref,
            request.alt
        )

        # Batch lookup (even for one) is faster and more consistent now
        results_map = ensembl.get_equivalents_batch([primary])
        return results_map.get(primary, [primary])

    except Exception as e:
        logger.error(f"HGVS alternatives fetch failed: {e}")
        raise BioEngineError("HGVS alternatives fetch failed")

@router.post("/jobs/{job_id}/share", response_model=dict)
def share_job(job_id: str, request: ShareJobRequest):
    """
    Exports a job to a local directory for sharing.
    """
    try:
        export_path = job_manager.export_job(
            job_id=job_id,
            level=request.level.value,
            target_folder=request.target_folder
        )
        return {"status": "success", "export_path": export_path}
    except Exception as e:
        logger.error(f"Failed to export job {job_id}: {e}")
        raise BioEngineError(f"Failed to export job {job_id}")

@router.post("/jobs/import", response_model=Job)
def import_job(request: ImportJobRequest, db: Session = Depends(get_db)):
    """
    Imports a shared job and synchronizes its approved variants to the global database.
    """
    try:
        job = job_manager.import_job(request.source_folder)

        # Sync approved variants from the imported job
        for variant_key, status in job.variant_statuses.items():
            if status == "approved":
                # Find the variant in results to get coordinates
                # variant_key format is usually "patient_id:variant_id" or similar
                # We need to find the matching entry in job.results
                patient_id, variant_idx_str = variant_key.split(':') if ':' in variant_key else (None, variant_key)

                # Search for the variant in the results
                target_variant = None
                if job.results:
                    for res in job.results:
                        if res.get("id") == patient_id:
                            variants_data = res.get("variants", {})
                            rows = variants_data.get("rows", [])
                            cols = variants_data.get("columns", [])
                            try:
                                v_idx = int(variant_idx_str)
                                if v_idx < len(rows):
                                    row = rows[v_idx]
                                    # Map columns
                                    col_map = {col.upper(): i for i, col in enumerate(cols)}

                                    # Ensure we have a chromosome name
                                    chr_name = row[col_map["CHR"]] if "CHR" in col_map else None
                                    if not chr_name:
                                        continue

                                    target_variant = {
                                        "chromosome": chr_name,
                                        "position": int(row[col_map.get("POS", 0)]),
                                        "ref_allele": row[col_map.get("REF", 0)],
                                        "alt_allele": row[col_map.get("ALT", 0)],
                                        "gene": row[col_map.get("GENE", 0)] if "GENE" in col_map else job.hgvs_config.gene if job.hgvs_config else None
                                    }
                                    break
                            except Exception:
                                continue

                if target_variant:
                    # Check if already exists to avoid duplicates
                    exists = db.query(ApprovedVariant).filter(
                        ApprovedVariant.job_id == job.id,
                        ApprovedVariant.chromosome == target_variant["chromosome"],
                        ApprovedVariant.position == target_variant["position"]
                    ).first()

                    if not exists:
                        db_variant = ApprovedVariant(
                            chromosome=target_variant["chromosome"],
                            position=target_variant["position"],
                            ref_allele=target_variant["ref_allele"],
                            alt_allele=target_variant["alt_allele"],
                            gene=target_variant["gene"],
                            approved_by="Imported",
                            job_id=job.id,
                            patient_id=patient_id,
                            assembly=job.hgvs_config.assembly if job.hgvs_config else "GRCh38"
                        )
                        db.add(db_variant)

        db.commit()
        return job
    except Exception as e:
        logger.error(f"Failed to import job from {request.source_folder}: {e}")
        db.rollback()
        raise BioEngineError("Failed to import job")

@router.post("/variants/approve")
def approve_variant(request: ApproveVariantRequest, db: Session = Depends(get_db)):
    """
    Saves a variant to the approved variants database.
    """
    db_variant = ApprovedVariant(
        chromosome=request.chromosome,
        position=request.position,
        ref_allele=request.ref_allele,
        alt_allele=request.alt_allele,
        gene=request.gene,
        approved_by=request.approved_by,
        job_id=request.job_id,
        patient_id=request.patient_id,
        assembly=request.assembly
    )
    db.add(db_variant)
    db.commit()
    db.refresh(db_variant)
    return {"status": "success", "id": db_variant.id}

@router.get("/variants/hotspots", response_model=list[HotspotPoint])
def get_hotspots(bin_size: int = 1000000, assembly: str = "GRCh38", db: Session = Depends(get_db)):
    """
    Returns aggregated variant counts (hotspots) grouped by chromosome and genomic bins.
    """
    # Group by chromosome and (position / bin_size)
    query = db.query(
        ApprovedVariant.chromosome,
        cast(ApprovedVariant.position / bin_size, Integer).label("bin"),
        func.count(ApprovedVariant.id).label("count")
    ).filter(ApprovedVariant.assembly == assembly).group_by(
        ApprovedVariant.chromosome,
        "bin"
    )

    results = []
    for chromosome, bin_index, count in query.all():
        results.append(HotspotPoint(
            chr=chromosome,
            start=bin_index * bin_size,
            stop=(bin_index + 1) * bin_size,
            count=count
        ))

    return results

@router.post("/cache/flush")
def flush_cache():
    """
    Clears the entire Redis cache.
    """
    success = cache.flush()
    if success:
        return {"status": "success", "message": "Cache cleared successfully"}
    else:
        return {"status": "error", "message": "Failed to clear cache or Redis not available"}, 500


@router.get("/variants/approved", response_model=list[ApprovedVariantResponse])
def get_approved_variants(assembly: str = "GRCh38", db: Session = Depends(get_db)):
    """
    Returns a list of all approved variants.
    """
    return db.query(ApprovedVariant).filter(ApprovedVariant.assembly == assembly).all()


# ---------------------------------------------------------------------------
# OpenCRAVAT management
# ---------------------------------------------------------------------------

@router.get("/opencravat/status", response_model=OCStatus)
def oc_status():
    """
    Returns OpenCRAVAT installation status, version, data directory path,
    and disk space consumed / available.
    """
    from core.config import settings
    return oc_service.get_status(oc_path=settings.oc_path)


@router.get("/opencravat/modules", response_model=list[OCModule])
def oc_list_installed():
    """
    Returns all locally installed OpenCRAVAT modules (annotators, reporters, etc.).
    """
    from core.config import settings
    return oc_service.list_installed(oc_path=settings.oc_path)


@router.get("/opencravat/modules/store", response_model=list[OCModule])
def oc_list_store():
    """
    Returns all modules available in the OpenCRAVAT store, including those
    not yet installed. Requires an internet connection to refresh the store index.
    """
    from core.config import settings
    return oc_service.list_store(oc_path=settings.oc_path)


@router.post("/opencravat/modules/{module_name}/install", response_model=OCInstallTask)
def oc_install_module(module_name: str):
    """
    Starts a background installation of the given OpenCRAVAT module.
    Returns a task object — poll GET /opencravat/tasks/{task_id} for status.
    """
    from core.config import settings
    task_id = oc_service.install_module(module_name, oc_path=settings.oc_path)
    return oc_service.get_task(task_id)


@router.get("/opencravat/tasks", response_model=list[OCInstallTask])
def oc_list_tasks():
    """
    Returns all registered OpenCRAVAT module installation tasks.
    """
    return oc_service.list_tasks()


@router.get("/opencravat/tasks/{task_id}", response_model=OCInstallTask)
def oc_get_task(task_id: str):
    """
    Returns the current status of an install task started via the install endpoint.
    Status values: pending | running | completed | failed
    """
    task = oc_service.get_task(task_id)
    if task is None:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail=f"Task '{task_id}' not found")
    return task


@router.delete("/opencravat/modules/{module_name}")
def oc_uninstall_module(module_name: str):
    """
    Uninstalls the given OpenCRAVAT module and removes its data from disk.
    """
    from core.config import settings
    return oc_service.uninstall_module(module_name, oc_path=settings.oc_path)


@router.get("/version")
def get_version():
    """
    Returns the version of the bio-engine.
    """
    return {"version": "0.8.0"}
