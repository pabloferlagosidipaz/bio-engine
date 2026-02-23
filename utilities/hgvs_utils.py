"""
HGVS Utilities
==============

This module provides the `HGVSAnnotator` class and connection singletons for 
interfacing with the Universal Transcript Archive (UTA) via the `hgvs` package.
It supports deep variant notation mapping (e.g., genomic to transcript).
"""

import logging
import threading

import hgvs.assemblymapper
import hgvs.dataproviders.interface
import hgvs.edit
import hgvs.location
import hgvs.parser
import hgvs.posedit
import hgvs.sequencevariant

from data.models import HGVSConfig

logger = logging.getLogger(__name__)

_shared_hdp = None

def get_uta_connection():
    """Returns a shared HGVS UTA connection (Singleton)."""
    global _shared_hdp
    if _shared_hdp is None:
        try:
            logger.info("Initializing shared HGVS UTA connection...")
            import hgvs.dataproviders.uta
            _shared_hdp = hgvs.dataproviders.uta.connect()
            logger.info("Shared HGVS UTA connection established.")
        except Exception as e:
            logger.error(f"Failed to connect to HGVS UTA: {e}")
            raise
    return _shared_hdp

class HGVSAnnotator:
    """
    Stateful orchestrator for mapping genetic variants utilizing `hgvs` mappers.
    
    Resolves mapping between genomic (NC_), coding (NM_), and arbitrary custom
    identifiers. Encapsulates an internal lock to ensure thread safety across
    background workers mutating common data frames.
    """
    def __init__(self, hdp: hgvs.dataproviders.interface.Interface | None = None, assembly: str = "GRCh38"):
        if hdp is None:
             hdp = get_uta_connection()

        self.hdp = hdp
        self.assembly = assembly
        self._lock = threading.Lock()
        try:
            self.am = hgvs.assemblymapper.AssemblyMapper(
                hdp,
                assembly_name=assembly,
                alt_aln_method='splign',
                replace_reference=True
            )
            self.hp = hgvs.parser.Parser()
        except Exception as e:
            logger.error(f"Failed to initialize HGVS mappers: {e}")
            raise

    @staticmethod
    def get_hgvs_type(accession: str) -> str:
        if accession.startswith(('NM_', 'XM_')):
            return 'c'
        if accession.startswith(('NR_', 'XR_')):
            return 'n'
        if accession.startswith(('NC_', 'NG_', 'NW_', 'NT_')):
            return 'g'
        return 'n' # fall back to non-coding transcript

    def annotate_data(self, data: dict, hgvs_config: HGVSConfig, full: bool = False, col_name: str = "hgvs", mode: str = "primary", source_ac: str | None = None):
        """
        Adds HGVS equivalents to variants rows in the provided data dict.
        mode: 'primary' (default), 'full', or 'genomic'
        source_ac: The accession used for 'pos' in the data (if different from target transcript)
        """
        with self._lock:
            variants_rows = data.get("variants", {}).get("rows", [])
            header = data.get("variants", {}).get("columns", [])

            if col_name not in header:
                header.append(col_name)
                data["variants"]["columns"] = header

            try:
                pos_idx = header.index("pos")
                ref_idx = header.index("ref")
                alt_idx = header.index("alt")
                hgvs_idx = header.index(col_name)
            except ValueError:
                logger.warning("Required variant columns (pos, ref, alt) not found in header.")
                return

            # Add hgvs_alternatives to data if not present
            if "hgvs_alternatives" not in data:
                data["hgvs_alternatives"] = {}

            target_ac = hgvs_config.transcript


            for i, row in enumerate(variants_rows):
                try:
                    pos = int(row[pos_idx])
                    ref = row[ref_idx]
                    alt = row[alt_idx]

                    if not source_ac:
                         current_source = target_ac if target_ac else None
                    else:
                         current_source = source_ac

                    if not current_source:
                        continue # Cannot map without a starting point

                    source_type = self.get_hgvs_type(current_source)

                    equivalents = self.find_equivalents(current_source, source_type, pos, ref, alt)

                    selected_hgvs = ""

                    candidates_nm = []
                    candidates_xm = []
                    exact_match = None

                    for eq in equivalents:
                        if target_ac and eq.startswith(target_ac + ":"):
                            exact_match = eq
                        elif eq.startswith("NM_"):
                            candidates_nm.append(eq)
                        elif eq.startswith("XM_"):
                            candidates_xm.append(eq)

                    if candidates_nm:
                        selected_hgvs = candidates_nm[0] # Pick first NM
                    elif exact_match:
                        selected_hgvs = exact_match
                    elif candidates_xm:
                        selected_hgvs = candidates_xm[0]
                    elif equivalents:
                        selected_hgvs = equivalents[0]

                    # Store selected
                    if len(row) < len(header):
                        row.append(selected_hgvs)
                    else:
                        row[hgvs_idx] = selected_hgvs

                    # Store alternatives (excluding selected)
                    if selected_hgvs:
                        equivalents.sort(key=lambda x: (0 if x.startswith("NM") else 1 if x.startswith("NP") else 2 if x.startswith("NC") else 3, x))
                        if equivalents:
                            data["hgvs_alternatives"][selected_hgvs] = equivalents

                except Exception as ex:
                    logger.error(f"Error processing variant row {i}: {ex}")
                    if len(row) < len(header):
                        row.append("")

    def find_equivalents(self, ac: str, ref_type: str, pos: int, ref: str, alt: str) -> list[str]:
        """
        Explores all equivalent HGVS variants using breadth-first search (BFS).
        
        Args:
            ac (str): Accession string (e.g. 'NM_00123.4').
            ref_type (str): Type of reference ('c', 'g', 'n', 'p').
            pos (int): Origin position for the variant.
            ref (str): Reference nucleotide/amino acid.
            alt (str): Alternate sequence.
            
        Returns:
            list[str]: Sorted list of unique mapped HGVS string variations.
        """
        try:
            edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
            posedit = hgvs.posedit.PosEdit(
                pos=hgvs.location.Interval(
                    start=hgvs.location.SimplePosition(pos),
                    end=hgvs.location.SimplePosition(pos + len(ref) - 1)
                ),
                edit=edit
            )
            var_primary = hgvs.sequencevariant.SequenceVariant(ac=ac, type=ref_type, posedit=posedit)

            all_vars = {str(var_primary): var_primary}
            to_process = [var_primary]
            processed_acs = {ac}

            depth = 0
            while to_process and depth < 2:  # Reduced depth to 2 to prevent excessive network calls
                next_to_process = []
                for var in to_process:
                    v_type = var.type
                    v_ac = var.ac

                    # 1. Map to equivalents via mapping options (genomic <-> transcript)
                    try:
                        options = self.hdp.get_tx_mapping_options(v_ac)
                        if options:
                            for opt in options[:10]: # Limit to first 10 options to avoid huge fan-out
                                target_ac = opt[1]
                                if target_ac in processed_acs: continue

                                try:
                                    target_type = self.get_hgvs_type(target_ac)
                                    v_mapped = None

                                    if v_type == 'c' and target_type == 'g':
                                        v_mapped = self.am.c_to_g(var)
                                    elif v_type == 'n' and target_type == 'g':
                                        v_mapped = self.am.n_to_g(var)
                                    elif v_type == 'g':
                                        if target_type == 'c': v_mapped = self.am.g_to_c(var, target_ac)
                                        elif target_type == 'n': v_mapped = self.am.g_to_n(var, target_ac)

                                    if v_mapped:
                                        v_str = str(v_mapped)
                                        if v_str not in all_vars:
                                            all_vars[v_str] = v_mapped
                                            next_to_process.append(v_mapped)
                                        processed_acs.add(v_mapped.ac)
                                except: pass
                                processed_acs.add(target_ac)
                    except: pass

                    # 2. If genomic, find all overlapping transcripts
                    if v_type == 'g':
                        try:
                            interval = var.posedit.pos
                            transcripts = self.hdp.get_tx_for_region(v_ac, 'splign', interval.start.base, interval.end.base)
                            if transcripts:
                                # Prioritize NM_ transcripts first to find the best coding matches quickly
                                sorted_transcripts = sorted(transcripts, key=lambda x: 0 if x[0].startswith("NM_") else 1)

                                for tx_row in sorted_transcripts[:5]: # Limit to top 5 overlapping transcripts
                                    tx_ac = tx_row[0]
                                    if tx_ac in processed_acs: continue

                                    try:
                                        tx_type = self.get_hgvs_type(tx_ac)
                                        v_tx = None
                                        if tx_type == 'c': v_tx = self.am.g_to_c(var, tx_ac)
                                        elif tx_type == 'n': v_tx = self.am.g_to_n(var, tx_ac)

                                        if v_tx:
                                            v_str = str(v_tx)
                                            if v_str not in all_vars:
                                                all_vars[v_str] = v_tx
                                                next_to_process.append(v_tx)
                                            processed_acs.add(v_tx.ac)
                                    except: pass
                                    processed_acs.add(tx_ac)
                        except: pass

                    # 3. If coding transcript, try to find protein
                    if v_type == 'c':
                        try:
                            var_p = self.am.c_to_p(var)
                            if var_p:
                                v_str = str(var_p)
                                if v_str not in all_vars:
                                    all_vars[v_str] = var_p
                                processed_acs.add(var_p.ac)
                        except: pass

                to_process = next_to_process
                depth += 1

            return sorted(list(all_vars.keys()))

        except Exception as e:
            logger.error(f"Failed to find equivalents for {ac}:{pos}{ref}>{alt}: {e}")
            return []

