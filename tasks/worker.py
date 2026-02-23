"""
Background Worker Tasks
=======================

This module contains non-blocking background logic for processing asynchronous
bio-engine jobs. It coordinates parallel alignment and subsequent annotation loops.
"""

import concurrent.futures
import logging
import os
import tempfile
from typing import Any

import orjson

from data.models import AlignmentResponse, Job, JobStatus
from services import aligner as aligner_service
from services import reference as ref_service
from services.job_manager import JobManager

logger = logging.getLogger(__name__)

# Initialize managers
job_manager = JobManager()

def _process_single_read(patient_id: str, read_entry: Any, job: Job, ref_path: str, job_output_base: str, annotator: Any = None) -> dict[str, Any]:
    """
    Submits a single sample read to the alignment service in an isolated context.
    
    Args:
        patient_id (str): The ID of the patient owning this read.
        read_entry (Any): The sequence payload or path containing the trace.
        job (Job): The contextual job.
        ref_path (str): The path to the resolved reference file.
        job_output_base (str): Temporary folder path to store Tracy alignments.
        annotator (Any): Active HGVSAnnotator instance for sequence variant metadata.
        
    Returns:
        dict[str, Any]: A serialized JSON block representing the aligned sequence trace.
    """
    # read_entry is a dict (from JobRead) or string (legacy)
    read_path = read_entry['file'] if isinstance(read_entry, dict) else str(read_entry)

    try:
        # Prepare config for this specific read
        current_config = job.config
        if isinstance(read_entry, dict) and (read_entry.get('trimLeft') is not None or read_entry.get('trimRight') is not None):
            # Create a copy of the config to avoid modifying the global job config
            if current_config:
                # Create a dictionary from the current config
                config_dict = current_config.model_dump()
                if read_entry.get('trimLeft') is not None:
                    config_dict['trimLeft'] = read_entry['trimLeft']
                if read_entry.get('trimRight') is not None:
                    config_dict['trimRight'] = read_entry['trimRight']
                current_config = type(current_config)(**config_dict)

        logger.info(f"Starting alignment for {read_path} (ref: {ref_path})")
        # aligner_service.run_alignment should likely return a Path or str
        alignment_json_path = aligner_service.run_alignment(
            patient_ab1=read_path,
            reference_fasta=ref_path,
            config=current_config,
            hgvs_config=job.hgvs_config,
            alignment_id=f"{patient_id}_{os.path.basename(read_path)}",
            output_base_dir=job_output_base,
            annotator=annotator
        )

        with open(alignment_json_path, "rb") as f:
            data = orjson.loads(f.read())

        # Structure the result using the self-mapping AlignmentResponse model
        # This automatically handles alignment, trace, and all extra fields
        alignment_obj = AlignmentResponse.model_validate(data)

        result_entry = {
            "patientId": patient_id,
            "readPath": read_path,
            "alignment": alignment_obj.model_dump()
        }
        return result_entry

    except Exception as e:
        logger.error(f"Error processing {read_path} for patient {patient_id}: {e}")
        return {
            "patientId": patient_id,
            "readPath": read_path,
            "error": str(e)
        }

def process_job_background(job_id: str):
    """
    Executes an entire Bio-Engine analysis pipeline asynchronously.
    
    1. Loads the specified reference.
    2. Initializes HGVS services via UTA if requested.
    3. Triggers parallel alignment across multiple workers.
    4. Gathers multi-sequence mutations into HGVS formatting.
    5. Dispatches global genomic annotations to Ensembl's VEP API.
    
    Args:
        job_id (str): The identifier of the queued job to run.
    """
    try:
        job = job_manager.get_job(job_id)
        if not job:
            logger.error(f"Background Job {job_id} not found")
            return

        # Resolve Reference
        ref_input = job.reference.get("value") if isinstance(job.reference, dict) else getattr(job.reference, "value", None)
        if not ref_input:
            logger.error(f"Job {job_id} missing reference value")
            job_manager.update_job_status(job_id, JobStatus.FAILED)
            return

        # Use service to load reference
        ref_path = ref_service.load_reference(ref_input)
        job_manager.update_job_progress(job_id, 10, "Reference loaded successfully.")

        # Ensure reference_sequence is populated
        if not job.reference_sequence:
            try:
                seq = ref_service.get_fasta_sequence(ref_path)
                job_manager.update_job_reference_sequence(job_id, seq)
                job.reference_sequence = seq
            except Exception as e:
                logger.error(f"Failed to load reference sequence for job {job_id}: {e}")

        # Extract features (exons, CDS, etc.)
        features = ref_service.get_reference_features(ref_input)
        if features:
            logger.info(f"Extracted {len(features)} features for job {job_id}")
            logger.info("Saving features to job...")
            job_manager.update_job_features(job_id, features)

        job_manager.update_job_progress(job_id, 15, "Features extracted. Starting alignment...")

        results = []

        # Create a temporary directory for this job's output
        job_output_base = tempfile.mkdtemp(prefix=f"ms_analyzer_{job_id}_")
        logger.info(f"Created job output directory: {job_output_base}")

        # Prepare tasks for parallel execution
        max_workers = min(32, (os.cpu_count() or 1) + 4)

        # Initialize HGVS Annotator once
        # Initialize HGVS Annotator once
        annotator = None
        if job.hgvs_config and job.hgvs_config.transcript:
            job_manager.update_job_progress(job_id, 16, "Connecting to HGVS database (UTA)...")
            try:
                import socket

                from utilities.hgvs_utils import HGVSAnnotator, get_uta_connection

                # Attempt connection with timeout
                socket.setdefaulttimeout(10) # 10 seconds timeout
                try:
                    logger.info("Getting shared connection to HGVS UTA...")
                    # use singleton
                    hdp = get_uta_connection()
                    logger.info("HGVS UTA connected.")
                    assembly = job.hgvs_config.assembly if job.hgvs_config.assembly else "GRCh38"
                    # Pass hdp to annotator
                    annotator = HGVSAnnotator(hdp, assembly=assembly)
                    logger.info("Shared HGVS Annotator initialized successfully")
                finally:
                    socket.setdefaulttimeout(None) # Reset timeout

            except Exception as e:
                logger.warning(f"Failed to initialize shared HGVS annotator (network issue?): {e}")
                job_manager.update_job_progress(job_id, 16, "HGVS connection failed. Proceeding without annotation...")
                # Proceed with annotator = None

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_read = {}

            for patient in job.patients:
                patient_id = patient.get("id")
                reads = patient.get("reads", [])

                for read_entry in reads:
                    future = executor.submit(
                        _process_single_read,
                        patient_id,
                        read_entry,
                        job,
                        ref_path,
                        job_output_base,
                        annotator
                    )
                    future_to_read[future] = (patient_id, read_entry)

            total_reads = len(future_to_read)
            completed_reads = 0

            for future in concurrent.futures.as_completed(future_to_read):
                try:
                    res = future.result()
                    results.append(res)

                    read_path_processed = res.get('readPath', 'unknown')
                    logger.info(f"Finished processing read: {read_path_processed}")
                except Exception as exc:
                    pid, re = future_to_read[future]
                    read_p = re['file'] if isinstance(re, dict) else re
                    logger.error(f"Read processing generated an exception: {exc}")
                    results.append({
                        "patientId": pid,
                        "readPath": read_p,
                        "error": str(exc)
                    })

                completed_reads += 1
                # Progress from 15% to 80%
                progress_percent = 15 + int((completed_reads / total_reads) * 65)
                job_manager.update_job_progress(job_id, progress_percent, f"Processed {completed_reads}/{total_reads} reads...")

        # Aggregate HGVS alternatives from all results
        all_alternatives = {}
        for res in results:
            # Check deep inside alignment result for hgvs_alternatives
            # Structure: res -> alignment -> hgvs_alternatives
            if "alignment" in res and isinstance(res["alignment"], dict):
                alts = res["alignment"].pop("hgvs_alternatives", None)
                if alts:
                    all_alternatives.update(alts)

        # Bulk update alternatives to job
        if all_alternatives:
            job_manager.update_job_hgvs_alternatives_bulk(job_id, all_alternatives)

        # Update results in job manager before VEP processing to ensure we have the latest HGVS
        job_manager.update_job_results(job_id, results)

        # Global VEP Annotation (if enabled)
        if job.config and job.config.vep and job.hgvs_config and job.hgvs_config.auto_vep:
            job_manager.update_job_progress(job_id, 80, "Annotating variants with VEP (this may take a while)...")
            try:
                from utilities.vep_utils import VEPAnnotator
                vep_annotator = VEPAnnotator(assembly=job.hgvs_config.assembly or "GRCh38")

                # 1. Collect all unique HGVS from results
                all_hgvs = set()
                for res in results:
                    variants = res.get("alignment", {}).get("variants", {})
                    rows = variants.get("rows", [])
                    header = variants.get("columns", [])
                    if "hgvs" in header:
                        h_idx = header.index("hgvs")
                        for row in rows:
                            if len(row) > h_idx and row[h_idx]:
                                hgvs_val = row[h_idx]
                                best_hgvs = hgvs_val
                                if hgvs_val in all_alternatives:
                                    equivalents = all_alternatives[hgvs_val]
                                    # For VEP, NC_ (genomic chromosome) is the most reliable input.
                                    # If not available, we fall back to NM_ or the original.
                                    nc_alt = next((alt for alt in equivalents if alt.startswith("NC_")), None)
                                    if nc_alt:
                                        best_hgvs = nc_alt
                                    elif not best_hgvs.startswith("NM_"):
                                        nm_alt = next((alt for alt in equivalents if alt.startswith("NM_")), None)
                                        if nm_alt:
                                            best_hgvs = nm_alt
                                all_hgvs.add(best_hgvs)

                # 2. Filter out already cached annotations
                existing_hgvs = set(job.vep_annotations.keys())
                hgvs_to_annotate = list(all_hgvs - existing_hgvs)

                if hgvs_to_annotate:
                    # Filter out HGVS strings that VEP doesn't support (e.g. NG_ genomic refs)
                    # VEP REST API primarily supports Transcript (NM, NR, XM, XR, ENST) or Chromosome (NC)
                    valid_prefixes = ('NM_', 'NR_', 'XM_', 'XR_', 'NC_', 'ENST')
                    valid_hgvs = [h for h in hgvs_to_annotate if h.startswith(valid_prefixes) or ':' not in h]
                    # check for 1:g.123 (chromosome:g.)
                    # Logic: if it starts with NG_ skip it.

                    filtered_hgvs = []
                    for h in hgvs_to_annotate:
                         if h.startswith("NG_"):
                             logger.warning(f"Skipping VEP annotation for {h}: NG_ accession not supported by Ensembl REST")
                             continue
                         filtered_hgvs.append(h)

                    if filtered_hgvs:
                         logger.info(f"Fetching VEP annotations for {len(filtered_hgvs)} variants...")
                         job_manager.update_job_progress(job_id, 85, f"Annotating {len(filtered_hgvs)} variants via Ensembl API...")
                         new_annotations = vep_annotator.get_annotations(filtered_hgvs)
                         if new_annotations:
                             job_manager.update_job_vep_annotations(job_id, new_annotations)
                    else:
                        logger.info("No valid HGVS strings for VEP (filtered all NG_ refs).")
                        job_manager.update_job_progress(job_id, 90, "Skipped VEP (unsupported reference type).")

                else:
                    logger.info("All variants already have VEP annotations (using cache).")
                    job_manager.update_job_progress(job_id, 90, "Using cached VEP annotations.")

            except Exception as e:
                logger.error(f"Global VEP annotation failed for job {job_id}: {e}")
                job_manager.update_job_progress(job_id, 90, f"VEP Annotation failed: {str(e)}")

        # Final status update
        job_manager.update_job_progress(job_id, 100, "Analysis completed successfully.")
        job_manager.update_job_status(job_id, JobStatus.COMPLETED)

    except Exception as e:
        job_manager.update_job_status(job_id, JobStatus.FAILED)
        job_manager.update_job_progress(job_id, 100, f"Analysis failed: {str(e)}")
        logger.error(f"Job {job_id} background task failed: {e}")

def annotate_hgvs_background(job_id: str):
    """
    Independent background task to re-annotate an existing job's results
    with new or alternate HGVS notations via UTA network connection.
    
    Args:
        job_id (str): The target job's identifier.
    """
    try:
        job = job_manager.get_job(job_id)
        if not job or not job.results:
            logger.error(f"Job {job_id} not found or has no results for annotation")
            return

        if not job.hgvs_config or not job.hgvs_config.transcript:
            logger.error(f"Job {job_id} missing HGVS configuration")
            return

        from utilities.hgvs_utils import HGVSAnnotator, get_uta_connection

        # Connect to UTA
        try:
            hdp = get_uta_connection()
        except Exception as e:
            logger.error(f"Failed to connect to UTA: {e}")
            job_manager.update_job_status(job_id, JobStatus.FAILED)
            return

        assembly = job.hgvs_config.assembly if job.hgvs_config.assembly else "GRCh38"
        annotator = HGVSAnnotator(hdp, assembly=assembly)

        updated_results = []
        for res in job.results:
            # Re-process variants to add HGVS
            # We need to reconstruct the data structure expected by add_hgvs_equivalents

            new_res = res.copy()
            alignment_data = new_res.get("alignment", {})

            if alignment_data.get("variants"):
                try:
                    # Direct annotation
                    annotator.annotate_data(alignment_data, job.hgvs_config, full=False)
                except Exception as e:
                    logger.warning(f"HGVS processing failed for one result in job {job_id}: {e}")

            updated_results.append(new_res)

        # Update job
        job_manager.update_job_results(job_id, updated_results)
        job_manager.update_job_status(job_id, JobStatus.COMPLETED)
        logger.info(f"Job {job_id} HGVS annotation completed")

    except Exception as e:
        job_manager.update_job_status(job_id, JobStatus.FAILED)
        logger.error(f"Job {job_id} HGVS annotation failed: {e}")
