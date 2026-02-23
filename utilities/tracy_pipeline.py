"""
Tracy Alignment Engine Wrapper
==============================

This module wraps the `tracy` CLI tool for sequence alignment, mutation decomposition,
and basecalling. It orchestrates the internal JSON output from `tracy` into structured
application models and delegates to `HGVSAnnotator` for genomic notation integration.
"""

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import orjson

from core.exceptions import AlignmentError
from data.models import HGVSConfig, TracyConfig
from services.reference import ensure_indexed
from utilities.hgvs_utils import HGVSAnnotator, get_uta_connection
from utilities.sequence_utils import get_complement, get_iupac_consensus, get_reverse_complement, two_iupac_consensus

logger = logging.getLogger(__name__)

class TracyPipeline:
    def __init__(self, output_dir: str = "results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.hdp = None

    def _get_hgvs_dataprovider(self):
        if self.hdp is None:
            try:
                self.hdp = get_uta_connection()
            except Exception as e:
                logger.error(f"Failed to connect to UTA: {e}")
        return self.hdp

    def process_samples(self, reference_path: str, sanger_files: list[str], config: TracyConfig | None = None, hgvs_config: HGVSConfig | None = None, annotator: HGVSAnnotator | None = None) -> list[str]:
        """Processes 1 or 2 Sanger sequences using Tracy."""
        config = config or TracyConfig()

        if len(sanger_files) not in [1, 2]:
            raise ValueError("Input array must contain 1 or 2 Sanger sequence files.")

        json_paths = []
        for sanger_file in sanger_files:
            filename = os.path.basename(sanger_file)
            logger.info(f"Processing individual file: {filename}")

            base_name = os.path.splitext(filename)[0]
            output_prefix = str(self.output_dir / base_name)

            json_path = self._run_decompose_safe(reference_path, sanger_file, output_prefix, config, hgvs_config, annotator)
            if json_path:
                json_paths.append(json_path)

        return json_paths

    def _run_decompose_safe(self, ref: str, ab1: str, prefix: str, config: TracyConfig, hgvs_config: HGVSConfig | None = None, annotator: HGVSAnnotator | None = None) -> str | None:
        """Attempts decomposition and handles large reference errors."""
        try:
            return self._execute_decompose(ref, ab1, prefix, config, hgvs_config, annotator)
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)

            if "Reference is larger than 50Kbp" in error_msg:
                logger.warning("Reference > 50Kbp detected. Attempting to auto-index...")
                try:
                    indexed_ref = ensure_indexed(ref)
                    logger.info(f"Reference indexed successfully: {indexed_ref}")
                    return self._execute_decompose(indexed_ref, ab1, prefix, config, hgvs_config, annotator)
                except Exception as index_err:
                    msg = f"Auto-indexing failed: {index_err}"
                    logger.error(msg)
                    raise AlignmentError(msg)
            else:
                logger.error(f"Tracy decompose failed: {error_msg}")
                raise AlignmentError(f"Tracy decompose failed: {error_msg}")

    def _execute_decompose(self, ref: str, ab1: str, prefix: str, config: TracyConfig, hgvs_config: HGVSConfig | None = None, annotator: HGVSAnnotator | None = None) -> str:
        """Executes the tracy decompose command."""
        cmd = [
            "tracy", "decompose",
            "-p", str(config.pratio), "-k", str(config.kmer),
            "-s", str(config.support), "-i", str(config.maxindel),
            "-g", str(config.gapopen), "-e", str(config.gapext),
            "-m", str(config.match), "-n", str(config.mismatch),
            "-t", str(config.trim), "-q", str(config.trimLeft), "-u", str(config.trimRight),
            "-l", str(config.linelimit), "-o", prefix, "-r", ref, ab1
        ]

        if config.annotate:
            cmd.extend(["-a", config.annotate])
        if config.callVariants:
            cmd.append("-v")

        try:
            logger.info(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True, timeout=120) # 2 minute timeout
            logger.info(f"Alignment successful for {prefix}")
        except subprocess.TimeoutExpired:
            logger.error(f"Tracy decompose timed out for {ab1}")
            raise AlignmentError(f"Tracy decompose timed out after 120s for {ab1}")
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.decode() if e.stderr else str(e)
            logger.error(f"Tracy decompose failed with error: {error_msg}")
            raise # Let upstream handle it

        json_file = f"{prefix}.json"
        if os.path.exists(json_file):
            with open(json_file, 'rb') as f:
                data = orjson.loads(f.read())

            # Extract accession from ref path (e.g., /path/to/NG_008866.1.fa -> NG_008866.1)
            ref_basename = os.path.basename(ref)
            ref_ac = ref_basename
            for ext in ['.fasta', '.fa', '.fna', '.gz', '.fsa']:
                if ref_ac.lower().endswith(ext):
                    ref_ac = ref_ac[:-len(ext)]
                    # Handle double extensions like .fa.gz
                    for ext2 in ['.fasta', '.fa', '.fna', '.fsa']:
                        if ref_ac.lower().endswith(ext2):
                            ref_ac = ref_ac[:-len(ext2)]

            normalized_data = self._normalize_tracy_json(data, ref, ab1, hgvs_config, annotator, source_ac=ref_ac)

            with open(json_file, 'wb') as f:
                f.write(orjson.dumps(normalized_data, option=orjson.OPT_INDENT_2))

            return json_file
        raise AlignmentError(f"Alignment succeeded but output JSON not found: {json_file}")

    def _normalize_tracy_json(self, data: dict[str, Any], ref_path: str, ab1_path: str, hgvs_config: HGVSConfig | None = None, annotator: HGVSAnnotator | None = None, source_ac: str | None = None) -> dict[str, Any]:
        """Normalizes and enhances Tracy JSON data."""
        new_data = data.copy()

        is_reverse = data.get("ref1forward", 1) == 0

        if is_reverse:
            new_data["ms_analyzer"] = {"original_orientation": "reverse"}
            self._normalize_reverse_orientation(new_data, data)
        else:
            new_data["ms_analyzer"] = {"original_orientation": "forward"}

        if hgvs_config and hgvs_config.transcript:
            if not annotator:
                hdp = self._get_hgvs_dataprovider()
                if hdp:
                     assembly = hgvs_config.assembly if hgvs_config.assembly else "GRCh38"
                     annotator = HGVSAnnotator(hdp, assembly=assembly)

            self.add_hgvs_equivalents(new_data, hgvs_config, annotator, source_ac=source_ac)

        self._ensure_hgvs_column(new_data)

        # Calculate readSeqRef
        try:
            aligned_data = self.get_aligned_trace_data(ab1_path, ref_path)
            read_seq_ref = aligned_data.get("readSeqRef", "")
            new_data["readSeqRef"] = get_reverse_complement(read_seq_ref) if aligned_data.get("forward", 1) == 0 else read_seq_ref
        except Exception as e:
            logger.error(f"Failed to calculate readSeqRef: {e}")
            new_data["readSeqRef"] = ""

        if new_data.get("primarySeq"):
            new_data["baseCount"] = len(new_data["primarySeq"])

        new_data["consensusAlign"] = self.getConsensusPositioning(new_data)
        new_data["readSeqConsensus"] = self.getReadSeqConsensus(new_data)
        new_data["readSeqConsensusComplementary"] = list(get_complement("".join(new_data["readSeqConsensus"])))

        return new_data

    def _normalize_reverse_orientation(self, new_data: dict[str, Any], original_data: dict[str, Any]):
        """Handles normalization logic when the sequence was aligned in reverse."""
        # Normalize Peaks
        self._normalize_peaks(new_data, original_data)

        # Normalize Sequences
        for key in ["primarySeq", "secondarySeq", "ref1align", "alt1align",
                   "ref2align", "alt2align", "allele1align", "allele2align"]:
            if val := original_data.get(key):
                new_data[key] = get_reverse_complement(val)

        new_data["ref1pos"] = original_data.get("ref1pos", 0)
        new_data["ref1forward"] = 1

        # Normalize Basecalls and Variants
        self._normalize_basecalls_and_variants(new_data, original_data)

    def _ensure_hgvs_column(self, data: dict[str, Any]):
        """Ensures the HGVS column exists in variants table."""
        if variants := data.get("variants"):
            rows = variants.get("rows", [])
            header = variants.get("columns", [])
            if rows and "hgvs" not in header:
                header.append("hgvs")
                variants["columns"] = header
                hgvs_idx = header.index("hgvs")
                for row in rows:
                    if len(row) < len(header):
                        row.append("")
                    else:
                        row[hgvs_idx] = ""

    def getReadSeqConsensus(self, data: dict[str, Any]) -> list[str]:
        primarySeq = data.get("primarySeq", "")
        secondarySeq = data.get("secondarySeq", "")
        return two_iupac_consensus(primarySeq, secondarySeq)

    def getConsensusPositioning(self, data: dict[str, Any]) -> dict[str, Any]:
        alignPositioning1 = self.getAlginPositioning(data, 1)
        alignPositioning2 = self.getAlginPositioning(data, 2)

        # Merge keys
        allKeys = sorted(set(alignPositioning1.keys()) | set(alignPositioning2.keys()))

        consensus_positioning = {}

        for refPos in allKeys:
            item = {
                "refPos": int(refPos),
                "alt1": [], "sangerPos1": [],
                "alt2": [], "sangerPos2": [],
                "cons": []
            }

            if a1 := alignPositioning1.get(refPos):
                item.update(a1)
            if a2 := alignPositioning2.get(refPos):
                item.update(a2)

            if item["alt1"] and item["alt2"]:
                s1 = "".join(item["alt1"])
                s2 = "".join(item["alt2"])
                item["cons"] = get_iupac_consensus(s1, s2)
            consensus_positioning[str(refPos)] = item

        return consensus_positioning

    def getAlginPositioning(self, data: dict[str, Any], allele: int = 1) -> dict[str, Any]:
        """Calculates the alignment positioning relative to reference."""
        align_positioning = {}
        start_pos = data.get(f"ref{allele}pos", 0)
        idx_aligned = 0
        ref_offset = -1

        alt_key = f"alt{allele}"
        sanger_pos_key = f"sangerPos{allele}"

        refalign = data.get(f"ref{allele}align", "")
        altalign = data.get(f"alt{allele}align", "")

        if not refalign: return {}

        gaps_before_first_base = len(refalign) - len(refalign.lstrip("-"))

        is_reverse = data.get("ms_analyzer", {}).get("original_orientation", "forward") == "reverse"

        # Get trim value
        args = data.get("meta", {}).get("arguments", {})
        trim_key = "trimLeft" if not is_reverse else "trimRight"
        sanger_pos = args.get(trim_key, 0)

        achieved_first_base = False

        # Iterating by index to handle both strings simultaneously
        length = len(refalign)

        for idx_aligned in range(length):
            actual_base_ref = refalign[idx_aligned]
            actual_base_alt = altalign[idx_aligned]

            if actual_base_alt != '-':
                sanger_pos += 1

            if actual_base_ref == "-":
                if not achieved_first_base:
                    refPosIdx = start_pos + idx_aligned - gaps_before_first_base
                    refPosKey = str(refPosIdx)
                    if refPosKey not in align_positioning:
                        align_positioning[refPosKey] = {alt_key: [], sanger_pos_key: []}
                    align_positioning[refPosKey][alt_key].append(actual_base_alt)
                    align_positioning[refPosKey][sanger_pos_key].append(sanger_pos)
                else:
                    refPosIdx = start_pos + ref_offset
                    # Append to existing position (insertion relative to ref)
                    refPosKey = str(refPosIdx)
                    if refPosKey in align_positioning:
                        align_positioning[refPosKey][alt_key].append(actual_base_alt)
                        align_positioning[refPosKey][sanger_pos_key].append(sanger_pos)
            else:
                ref_offset += 1
                if not achieved_first_base:
                    achieved_first_base = True

                refPosIdx = start_pos + ref_offset
                align_positioning[str(refPosIdx)] = {
                    alt_key: [actual_base_alt],
                    sanger_pos_key: [sanger_pos]
                }

        return align_positioning

    def add_hgvs_equivalents(self, data: dict[str, Any], hgvs_config: HGVSConfig, annotator: HGVSAnnotator | None = None, source_ac: str | None = None):
        """Adds HGVS equivalents to variants using HGVSAnnotator."""
        try:
            if annotator:
                annotator.annotate_data(data, hgvs_config, full=False, source_ac=source_ac)
                return

            hdp = self._get_hgvs_dataprovider()
            if not hdp:
                return

            assembly = hgvs_config.assembly if hgvs_config.assembly else "GRCh38"
            annotator = HGVSAnnotator(hdp, assembly=assembly)
            annotator.annotate_data(data, hgvs_config, full=False, source_ac=source_ac)

        except Exception as e:
            logger.error(f"HGVS processing failed: {e}")


    def _normalize_peaks(self, new_data: dict[str, Any], data: dict[str, Any]):
        peakA, peakC, peakG, peakT = data.get("peakA", []), data.get("peakC", []), data.get("peakG", []), data.get("peakT", [])
        if peakA:
            new_data["peakA"], new_data["peakT"] = peakT[::-1], peakA[::-1]
            new_data["peakC"], new_data["peakG"] = peakG[::-1], peakC[::-1]
            max_len = len(peakA)
            basecallPos = data.get("basecallPos", [])
            if basecallPos:
                new_data["basecallPos"] = [(max_len - 1 - p) for p in reversed(basecallPos)]

    def _normalize_basecalls_and_variants(self, new_data: dict[str, Any], data: dict[str, Any]):
        basecallPos = data.get("basecallPos", [])
        peakA = data.get("peakA", [])
        max_len = len(peakA) if peakA else 0
        old_basecalls = data.get("basecalls", {})
        new_variants = data.get("variants", {"rows": []})

        if old_basecalls and max_len > 0:
            new_basecalls = {}
            num_calls = len(basecallPos)
            for old_sig_pos_str, call_str in old_basecalls.items():
                old_sig_pos = int(old_sig_pos_str)
                new_sig_pos = max_len - 1 - old_sig_pos
                parts = call_str.split(":")
                if len(parts) < 2: continue

                call_idx = int(parts[0])
                new_call_idx = num_calls - call_idx + 1

                for variant in new_variants.get("rows", []):
                    # Variant structure: [..., pos, ...] - index 10 was used in original code
                    # verification needed: tracy variant array structure.
                    # Assuming index 10 is indeed signal position.
                    if len(variant) > 10 and variant[10] == old_sig_pos:
                        variant[10], variant[9] = new_sig_pos, new_call_idx

                new_bases = [get_complement(b) for b in parts[1].split("|")]
                new_basecalls[str(new_sig_pos)] = f"{new_call_idx}:{'|'.join(new_bases)}"
            new_data["basecalls"] = new_basecalls

        chartConfig = data.get("chartConfig", {})
        if chartConfig and new_variants["rows"]:
            # Recalculate ranges
            # Assuming index 10 is signal pos
            newRanges = [[v[10]-150, v[10]+150] for v in new_variants["rows"] if len(v) > 10]
            if newRanges:
                chartConfig["x"]["axis"]["range"] = [min(r[0] for r in newRanges), max(r[1] for r in newRanges)]
                new_variants["xranges"] = newRanges
        new_data["chartConfig"] = chartConfig
        new_data["variants"] = new_variants

    def get_trace_data(self, ab1_path: str) -> dict[str, Any]:
        """Extracts trace data from an AB1 file using tracy basecall."""
        if not os.path.exists(ab1_path):
            raise FileNotFoundError(f"File not found: {ab1_path}")

        with tempfile.TemporaryDirectory() as temp_dir:
            base_name = os.path.basename(ab1_path)
            output_prefix = os.path.join(temp_dir, os.path.splitext(base_name)[0])
            json_file = f"{output_prefix}.json"

            cmd = ["tracy", "basecall", ab1_path, "-f", "json", "-o", json_file]

            try:
                subprocess.run(cmd, check=True, capture_output=True)

                if os.path.exists(json_file):
                    with open(json_file, 'rb') as f:
                        data = orjson.loads(f.read())

                    return {
                        "traceA": data.get("peakA", []),
                        "traceC": data.get("peakC", []),
                        "traceG": data.get("peakG", []),
                        "traceT": data.get("peakT", []),
                        "peakLocations": data.get("basecallPos", []),
                        "basecalls": data.get("basecalls", {})
                    }
                else:
                     raise AlignmentError("Tracy basecall produced no JSON output.")

            except subprocess.CalledProcessError as e:
                err = e.stderr.decode() if e.stderr else str(e)
                logger.error(f"Tracy basecall failed: {err}")
                raise AlignmentError(f"Tracy basecall failed: {err}")

    def get_aligned_trace_data(self, ab1_path: str, reference_path: str) -> dict[str, Any]:
        """
        Extracts trace data and alignment from an AB1 file using tracy align.
        Calculates readSeqRef which holds the reference nucleotide for each Sanger position.
        """
        if not os.path.exists(ab1_path):
            raise FileNotFoundError(f"AB1 file not found: {ab1_path}")
        if not os.path.exists(reference_path):
            raise FileNotFoundError(f"Reference file not found: {reference_path}")

        with tempfile.TemporaryDirectory() as temp_dir:
            base_name = os.path.basename(ab1_path)
            output_prefix = os.path.join(temp_dir, os.path.splitext(base_name)[0])
            json_file = f"{output_prefix}.json"

            # Run tracy align
            cmd = ["tracy", "align", "-r", reference_path, ab1_path, "-o", output_prefix]

            try:
                subprocess.run(cmd, check=True, capture_output=True)

                if os.path.exists(json_file):
                    with open(json_file, 'rb') as f:
                        data = orjson.loads(f.read())

                    # Extract alignment strings (try both standard and numbered keys)
                    refalign = data.get("refalign")
                    altalign = data.get("altalign")

                    read_seq_ref = self._calculate_read_seq_ref(refalign, altalign)
                    data["readSeqRef"] = read_seq_ref

                    return data
                else:
                     raise AlignmentError("Tracy align produced no JSON output.")

            except subprocess.CalledProcessError as e:
                err = e.stderr.decode() if e.stderr else str(e)
                logger.error(f"Tracy align failed: {err}")
                raise AlignmentError(f"Tracy align failed: {err}")

    def _calculate_read_seq_ref(self, refalign: str, altalign: str) -> str:
        """
        Calculates reference sequence corresponding to the read.
        Iterate through refalign and altalign:
          - If alt char is '-', skip (gap in read).
          - If ref char is '-', append '-' (gap in ref / insertion in read).
          - Otherwise, append ref char.
        """
        if not refalign or not altalign or len(refalign) != len(altalign):
            return ""

        read_seq_ref = []
        for r_char, a_char in zip(refalign, altalign):
            if a_char == '-':
                continue
            read_seq_ref.append(r_char)

        return "".join(read_seq_ref)
