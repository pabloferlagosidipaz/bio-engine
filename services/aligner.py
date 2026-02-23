"""
Alignment Service
=================

This module provides functions to run local DNA sequence alignments
using the Tracy CLI wrapper and to extract trace data for previewing.
"""

import logging
import os
import tempfile

from core.exceptions import AlignmentError, FileNotFoundError
from data.models import HGVSConfig, TracyConfig
from utilities.hgvs_utils import HGVSAnnotator
from utilities.tracy_pipeline import TracyPipeline

logger = logging.getLogger(__name__)

def run_alignment(
    patient_ab1: str, 
    reference_fasta: str, 
    config: TracyConfig | None = None, 
    hgvs_config: HGVSConfig | None = None, 
    alignment_id: str = "results", 
    output_base_dir: str | None = None, 
    annotator: HGVSAnnotator | None = None
) -> str:
    """
    Performs local alignment of a patient sequence against a reference.
    
    Args:
        patient_ab1 (str): Path to the patient's trace file (.ab1).
        reference_fasta (str): Path to the reference FASTA file.
        config (TracyConfig | None): Optional configuration for Tracy alignment.
        hgvs_config (HGVSConfig | None): Optional HGVS annotation configuration.
        alignment_id (str): Identifier for the generated output folder.
        output_base_dir (str | None): Base directory for output. Uses a temp dir if None.
        annotator (HGVSAnnotator | None): Instantiated HGVS annotator to use.
        
    Returns:
        str: Path to the generated JSON results parsed from Tracy.
        
    Raises:
        FileNotFoundError: If the input sequence files are missing.
        AlignmentError: If Tracy fails to produce JSON output.
    """
    if not (os.path.exists(reference_fasta) and os.path.exists(patient_ab1)):
        raise FileNotFoundError("Reference or sample file missing.")

    # Determine base directory
    if output_base_dir:
        if not os.path.exists(output_base_dir):
            os.makedirs(output_base_dir)
        temp_base = output_base_dir
    else:
        # Create a temporary directory for this alignment task
        temp_base = tempfile.mkdtemp(prefix="ms_analyzer_")

    output_dir = os.path.join(temp_base, alignment_id)

    tracy = TracyPipeline(output_dir=output_dir)
    json_paths = tracy.process_samples(reference_fasta, [patient_ab1], config=config, hgvs_config=hgvs_config, annotator=annotator)

    if not json_paths:
        raise AlignmentError("Alignment produced no output.")
    return json_paths[0]

def get_read_preview(read_path: str) -> dict:
    """
    Extracts trace metadata and peak data for previewing.
    
    Args:
        read_path (str): Path to the sequence file (e.g. .ab1).
        
    Returns:
        dict: Trace properties formatted as a dictionary.
        
    Raises:
        FileNotFoundError: If the read file doesn't exist.
    """
    if not os.path.exists(read_path):
        raise FileNotFoundError(f"Read file missing: {read_path}")

    tracy = TracyPipeline()
    # We use a temporary directory internally in get_trace_data
    return tracy.get_trace_data(read_path)
