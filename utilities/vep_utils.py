import logging
import re
import time
from datetime import datetime
from typing import Any

import httpx

# Configure traceability and logging levels
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class VEPAnnotator:
    """
    Advanced genomic variant annotation engine interacting with Ensembl's REST API.
    Implements a dual pipeline that automatically detects RefSeq nomenclatures 
    (NM_, NP_, NC_, NG_) and standardizes them using the Variant Recoder middleware
    before executing VEP consequence predictions.
    """

    def __init__(self, assembly: str = "GRCh38", timeout: float = 30.0):
        """
        Initializes the network client and routes traffic to the appropriate cluster based on the assembly.
        """
        self.assembly = assembly.upper()
        # DNS switching to support retrospective clinical genomic data on GRCh37
        if self.assembly == "GRCh37":
            self.base_url = "https://grch37.rest.ensembl.org"
        else:
            self.base_url = "https://rest.ensembl.org"

        self.client = httpx.Client(timeout=timeout)
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

        # Optimized lexical pattern for intercepting RefSeq prefixes
        self.refseq_pattern = re.compile(r'^(NM_|NP_|NC_|NG_)', re.IGNORECASE)

    def __del__(self):
        """
        Class destructor ensuring graceful closure of underlying HTTP sockets to prevent 
        memory leaks and port exhaustion during massive operations.
        """
        try:
            self.client.close()
        except Exception:
            pass

    def _handle_rate_limit(self, response: httpx.Response) -> bool:
        """
        Self-protection mechanism against network throttling (Rate Limiting).
        Analyzes HTTP 429 codes and obeys the 'Retry-After' header dictated by EMBL-EBI.
        """
        if response.status_code == 429:
            retry_after = int(response.headers.get("Retry-After", 1))
            logging.warning(f"Restricción de tráfico detectada (HTTP 429). Suspendiendo ejecución por {retry_after} segundos...")
            time.sleep(retry_after)
            return True
        return False

    def recode_variants_batch(self, variants: list[str]) -> dict[str, str]:
        """
        The Translator Middleware. Receives protein (p.) or RefSeq transcript HGVS strings 
        and requests Ensembl to recode them to absolute genomic coordinates (HGVSg).
        """
        if not variants:
            return {}

        # Use the Variant Recoder bulk (POST) endpoint for humans
        endpoint = f"{self.base_url}/variant_recoder/human"

        # Limit the response to extract only the stable genomic coordinate
        payload = {
            "ids": variants,
            "fields": "hgvsg,id"
        }

        translation_map = {}
        retries = 3

        for attempt in range(retries):
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)

                if self._handle_rate_limit(response):
                    continue

                response.raise_for_status()
                data = response.json()

                # Iterate over the resulting JSON array of decoded variants
                for item in data:
                    # Ensembl returns a dictionary per ID, indexed by the alternate allele
                    for allele_info in item.values():
                        if not isinstance(allele_info, dict):
                            continue

                        original_input = allele_info.get("input")
                        if not original_input:
                            continue

                        # Extract direct genomic mapping
                        hgvsg_list = allele_info.get("hgvsg", [])
                        if hgvsg_list and isinstance(hgvsg_list, list) and len(hgvsg_list) > 0:
                            # Prioritize the first provided genomic coordinate
                            translation_map[original_input] = hgvsg_list[0]
                        else:
                            # Fallback: retain original string if mapping fails
                            translation_map[original_input] = original_input

                        # Move to next item after processing current definition
                        break
                break

            except httpx.HTTPError as e:
                logging.error(f"Network failure during genomic recoding (Attempt {attempt + 1}): {e}")
                if attempt == retries - 1:
                    logging.warning("Retry limit exceeded. Injecting original strings to predictor.")
                    translation_map = {v: v for v in variants}
                time.sleep(2)

        return translation_map

    def annotate_hgvs_batch(self, hgvs_variants: list[str]) -> list[dict[str, Any]]:
        """
        Main orchestrator of the VEP pipeline. Triage, recoding, dispatch to the master predictor, 
        and assembly of the complex payload with implicit support for RefSeq models.
        """
        if not hgvs_variants:
            return []

        # Triage Phase: Binary classification based on ontological prefix
        variants_to_recode = []
        clean_variants = []

        for variant in hgvs_variants:
            if self.refseq_pattern.match(variant):
                variants_to_recode.append(variant)
            else:
                clean_variants.append(variant)

        # Recoding Phase: Middleware invocation
        recoded_map = {}
        if variants_to_recode:
            logging.info(f"Routing {len(variants_to_recode)} ambiguous/RefSeq sequences via Variant Recoder...")
            recoded_map = self.recode_variants_batch(variants_to_recode)

        # Synthesize universal operational block for the VEP engine
        final_vep_payload = clean_variants + list(recoded_map.values())
        final_vep_payload = list(set(final_vep_payload))  # Deduplicate for network efficiency

        endpoint = f"{self.base_url}/vep/human/hgvs"

        # Advanced Configuration Payload
        payload = {
            "hgvs_notations": final_vep_payload,
            "ignore_invalid": 1,      # Algorithmic firewall against full batch corruption
            "refseq": 1,              # Required direct resolution for RefSeq alignments
            "mane": 1,                # Aggressive filtering towards high-confidence canonical transcripts
            "sift": "b",              # Activate predictive SIFT plugin (score + class)
            "polyphen": "b",          # Activate structural PolyPhen-2 plugin (score + class)
            "symbol": 1,              # Expose taxonomic gene name
            "transcript_version": 1   # Strict sequence patch linkage
        }

        results = []
        retries = 3

        for attempt in range(retries):
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)

                if self._handle_rate_limit(response):
                    continue

                # Structural degradation fallback: Convert batch to serial requests if highly complex payload crashes (HTTP 400).
                if response.status_code == 400 and len(final_vep_payload) > 1:
                    logging.warning("Batch collapsed in VEP (HTTP 400). Degrading to serialized individual pipeline...")
                    return self._fallback_individual_requests(final_vep_payload, payload)

                response.raise_for_status()
                results = response.json()
                break

            except httpx.HTTPError as e:
                logging.error(f"Direct connection failure with VEP Engine (Attempt {attempt + 1}): {e}")
                if attempt == retries - 1:
                    return []
                time.sleep(2)

        # Topological Reconstruction:
        # VEP operated natively on genomic coordinates (NC_); However, users requested annotations 
        # using RefSeq identifiers (NP_ or NM_). Overwrite result labels to maintain transparent logic end-to-end.
        inverse_recoded_map = {v: k for k, v in recoded_map.items()}

        for record in results:
            input_variant = record.get("input")
            if input_variant in inverse_recoded_map:
                # Restore the user requested nomenclature
                record["original_input"] = inverse_recoded_map[input_variant]
            else:
                record["original_input"] = input_variant

        return results

    def _fallback_individual_requests(self, variants: list[str], base_payload: dict) -> list[dict[str, Any]]:
        """
        Data survival strategy. Decomposes a massive failed batch into atomic operations 
        to isolate the pathological variant and save the rest of the dataset.
        """
        endpoint = f"{self.base_url}/vep/human/hgvs"
        results = []

        for variant in variants:
            payload = base_payload.copy()
            payload["hgvs_notations"] = [variant]
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)
                if self._handle_rate_limit(response):
                    response = self.client.post(endpoint, json=payload, headers=self.headers)

                if response.status_code == 200:
                    results.extend(response.json())
            except Exception as e:
                logging.debug(f"Irresolvable variant isolated ({variant}): {e}")
                continue

        return results

    def _rank_impact(self, impact: str) -> int:
        """
        Hierarchical heuristic classifier to prioritize transcriptional damage.
        Quantifies the clinical severity of the variant based on VEP ontological consequences.
        """
        ranking = {
            "HIGH": 4,       # Truncations, stop codon loss, severe splicing alteration
            "MODERATE": 3,   # Missense mutations, in-frame deletions
            "LOW": 2,        # Synonymous mutations marginally altering translation
            "MODIFIER": 1    # Intronic, UTR variants or intergenic regions
        }
        return ranking.get(impact.upper(), 0)

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        """
        Final semantic processor. Ingests the raw list of variants, invokes the prediction cascade, 
        evaluates the multiplicity of affected transcripts, and distills the result into the most 
        critical and representative biological consequence for each original input.
        """
        raw_results = self.annotate_hgvs_batch(hgvs_variants)
        structured_annotations = {}

        for record in raw_results:
            # Unbreakable link via user's original identifier
            identifier = record.get("original_input", record.get("input"))

            # Taxonomic template for final variant report
            annotation = {
                "gene_symbol": "",
                "consequence": "",
                "impact": "",
                "hgvs_c": "",
                "hgvs_p": "",
                "sift": "",
                "polyphen": "",
                "clin_sig": [],
                "phenotype": [],
                "vep_raw": record,
                "retrieved_at": datetime.now().isoformat()
            }

            transcripts = record.get("transcript_consequences", [])
            most_severe_transcript = None
            highest_impact_score = -1

            # Transcript prioritization algorithm: scan all affected isoforms
            for transcript in transcripts:
                impact_str = transcript.get("impact", "")
                current_score = self._rank_impact(impact_str)

                # Clinical Modifier: Promote MANE Select universal transcripts
                is_mane = 'mane_select' in transcript

                # Elite selection of most severe consequence
                if current_score > highest_impact_score or (current_score == highest_impact_score and is_mane):
                    highest_impact_score = current_score
                    most_severe_transcript = transcript

            if most_severe_transcript:
                # Extraction and aesthetic mapping for human analyst
                annotation["gene_symbol"] = most_severe_transcript.get("gene_symbol", "")
                annotation["impact"] = most_severe_transcript.get("impact", "")

                consequence_terms = most_severe_transcript.get("consequence_terms", [])
                annotation["consequence"] = ", ".join(consequence_terms)

                annotation["hgvs_c"] = most_severe_transcript.get("hgvsc", "")
                annotation["hgvs_p"] = most_severe_transcript.get("hgvsp", "")

                # Consolidation of in-silico pathogenicity matrices
                if "sift_prediction" in most_severe_transcript and "sift_score" in most_severe_transcript:
                    annotation["sift"] = f"{most_severe_transcript['sift_prediction']} ({most_severe_transcript['sift_score']})"

                if "polyphen_prediction" in most_severe_transcript and "polyphen_score" in most_severe_transcript:
                    annotation["polyphen"] = f"{most_severe_transcript['polyphen_prediction']} ({most_severe_transcript['polyphen_score']})"

            # Extract clinical significance and phenotypes from colocated variants
            colocated_variants = record.get("colocated_variants", [])
            clin_sigs = set()
            phenotypes = []

            for cv in colocated_variants:
                if "clin_sig" in cv:
                    for sig in cv["clin_sig"]:
                        clin_sigs.add(sig.replace("_", " "))

                if "phenotype_or_disease" in cv and cv["phenotype_or_disease"] == 1:
                    # Depending on VEP config, specific disease names might be in `phenotype` or `var_synonyms`
                    # but typically if phenotype_or_disease is 1, it's structurally relevant.
                    if "id" in cv and not cv["id"].startswith("COS"):  # Filter out generic COSMIC structural IDs
                        phenotypes.append(cv["id"])
                    
                    # Some variants hold 'pubmed' arrays which are highly relevant for clinical
                    if "pubmed" in cv:
                         for pmid in cv["pubmed"]:
                             # Avoid massive arrays
                             if len(phenotypes) < 5:
                                 phenotypes.append(f"PMID:{pmid}")

            annotation["clin_sig"] = sorted(list(clin_sigs))
            # Deduplicate and cap phenotypes to prevent UI bloat
            annotation["phenotype"] = sorted(list(set(phenotypes)))[:10]

            structured_annotations[identifier] = annotation

        return structured_annotations
