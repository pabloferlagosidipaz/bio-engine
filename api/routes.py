"""
API Routes Definitions
=======================

This module defines the RESTful endpoints for the bio-engine sidecar application.
It exposes functionalities for job lifecycle management, annotation, alignment preview,
and biological database communication (e.g., fetching references, HGVS variants).
"""

import logging
import os

from fastapi import APIRouter, BackgroundTasks

from core.exceptions import BioEngineError
from data.models import (
    AddCommentRequest,
    AddHGVSAlternativesRequest,
    CreateJobRequest,
    HGVSRequest,
    ImportJobRequest,
    Job,
    JobStatus,
    RenameJobRequest,
    ShareJobRequest,
    UpdateJobRequest,
    ProxyConfigRequest,
)
from services import aligner as aligner_service
from services import reference as ref_service
from services.job_manager import JobManager
from tasks.worker import annotate_hgvs_background, process_job_background

logger = logging.getLogger(__name__)

router = APIRouter()
job_manager = JobManager()

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
    if request.http_proxy is not None:
        if request.http_proxy == "":
            os.environ.pop("HTTP_PROXY", None)
            os.environ.pop("http_proxy", None)
        else:
            os.environ["HTTP_PROXY"] = request.http_proxy
            os.environ["http_proxy"] = request.http_proxy
            
    if request.https_proxy is not None:
        if request.https_proxy == "":
            os.environ.pop("HTTPS_PROXY", None)
            os.environ.pop("https_proxy", None)
        else:
            os.environ["HTTPS_PROXY"] = request.https_proxy
            os.environ["https_proxy"] = request.https_proxy
            
    return {"status": "success", "http_proxy": os.environ.get("HTTP_PROXY"), "https_proxy": os.environ.get("HTTPS_PROXY")}

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
def search_reference(query: str):
    """
    Searches NCBI nucleotide database for a given query (gene name or accession).
    """
    results = ref_service.search_reference(query)
    return {"results": results}

@router.get("/preview-read")
def preview_read(path: str):
    """Returns trace data for a given file path."""
    try:
        return aligner_service.get_read_preview(path)
    except Exception as e:
        logger.error(f"Preview failed for {path}: {e}")
        raise BioEngineError(f"Preview failed: {e}")

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
    sequence = None
    try:
        if ref_input:
            ref_path = ref_service.load_reference(ref_input)
            sequence = ref_service.get_fasta_sequence(ref_path)
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
def delete_job(job_id: str):
    """
    Removes a job from the persistent storage.
    
    Args:
        job_id (str): The ID of the job to delete.
    """
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
    sequence = None
    try:
        if ref_input:
            ref_path = ref_service.load_reference(ref_input)
            sequence = ref_service.get_fasta_sequence(ref_path)
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
        raise BioEngineError(f"HGVS alternatives fetch failed: {e}")

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
        raise BioEngineError(f"Failed to export job {job_id}: {e}")

@router.post("/jobs/import", response_model=Job)
def import_job(request: ImportJobRequest):
    """
    Imports a shared job from a local directory back into the engine.
    """
    try:
        job = job_manager.import_job(request.source_folder)
        return job
    except Exception as e:
        logger.error(f"Failed to import job from {request.source_folder}: {e}")
        raise BioEngineError(f"Failed to import job: {e}")
