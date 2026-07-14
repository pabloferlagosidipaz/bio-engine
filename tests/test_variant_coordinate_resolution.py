# SPDX-License-Identifier: MIT
"""
Characterization tests for Phase 6: extracting the shared HGVS/SPDI genomic
coordinate resolver out of update_variant_status's two branches (patient-
scoped row index vs grouped genomic position) must not change resolution
behavior. These tests pin down the pre-refactor behavior of
_resolve_genomic_variant_from_hgvs and _resolve_genomic_variant_from_reference
directly, plus an end-to-end check that both call sites in the route still
produce identical results for equivalent inputs.
"""
import os
import sys
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from core.database import Base, get_db
from data.models import HGVSConfig, Job, JobStatus
from main import app


def _make_job(reference=None, hgvs_alternatives=None, hgvs_config=None, **overrides):
    defaults = dict(
        id="job-1",
        name="Test Job",
        status=JobStatus.COMPLETED,
        created_at="2026-01-01T00:00:00",
        updated_at="2026-01-01T00:00:00",
        reference=reference or {"type": "file", "value": "ref.fasta"},
        patients=[],
        hgvs_alternatives=hgvs_alternatives or {},
        hgvs_config=hgvs_config,
    )
    defaults.update(overrides)
    return Job(**defaults)


class TestResolveGenomicVariantFromHgvs(unittest.TestCase):
    def setUp(self):
        from api.routes import _resolve_genomic_variant_from_hgvs
        self.resolve = _resolve_genomic_variant_from_hgvs
        self.row = ["T", "C", "BRCA1"]
        self.col_map = {"REF": 0, "ALT": 1, "GENE": 2}

    def test_already_genomic_hgvs_resolves_without_ensembl_call(self):
        job = _make_job()
        with patch("utilities.ensembl_hgvs.EnsemblHGVS") as mock_cls:
            result = self.resolve(job, "job-1", "NC_000017.11:g.43045712T>C", self.row, self.col_map)
            mock_cls.assert_not_called()

        self.assertEqual(result, {
            "chromosome": "17", "position": 43045712,
            "ref": "T", "alt": "C", "gene": "BRCA1",
        })

    def test_versioned_hgvs_form_resolves_position_correctly(self):
        """
        Regression guard for a bug found while writing these Phase 6
        characterization tests: a versioned accession like "NC_000017.11"
        used to have its version number ("11") captured as the position
        instead of the real position (43045712), because the shared
        regex's separator matched the "." before the version rather than
        requiring the ":" immediately before "g." (fixed in Phase 7 - see
        the comment on _GENOMIC_HGVS_CAPTURE_RE in api/routes.py).
        """
        job = _make_job()
        result = self.resolve(job, "job-1", "NC_000017.11:g.43045712T>C", self.row, self.col_map)
        self.assertEqual(result["position"], 43045712)

    def test_unversioned_hgvs_form_resolves_position_correctly(self):
        job = _make_job()
        result = self.resolve(job, "job-1", "NC_000017:g.43045712T>C", self.row, self.col_map)
        self.assertEqual(result["position"], 43045712)

    def test_multi_digit_version_number_does_not_leak_into_position(self):
        job = _make_job()
        result = self.resolve(job, "job-1", "NC_000017.123:g.43045712T>C", self.row, self.col_map)
        self.assertEqual(result["position"], 43045712)

    def test_spdi_form_normalizes_zero_based_to_one_based(self):
        job = _make_job()
        result = self.resolve(job, "job-1", "17:43045711:T:C", self.row, self.col_map)
        self.assertEqual(result["chromosome"], "17")
        self.assertEqual(result["position"], 43045712)  # 43045711 + 1

    def test_mitochondrial_accession_resolves_to_mt(self):
        job = _make_job()
        result = self.resolve(job, "job-1", "NC_012920.1:g.100A>T", self.row, self.col_map)
        self.assertEqual(result["chromosome"], "MT")
        self.assertEqual(result["position"], 100)

    def test_chr_m_prefix_normalizes_to_mt(self):
        job = _make_job()
        # SPDI form (3 colons, no ":g.") is 0-based, so 100 -> 101.
        result = self.resolve(job, "job-1", "chrM:100:A:T", self.row, self.col_map)
        self.assertEqual(result["chromosome"], "MT")
        self.assertEqual(result["position"], 101)

    def test_zero_padded_chromosome_is_normalized(self):
        job = _make_job()
        result = self.resolve(job, "job-1", "NC_000001.11:g.100A>T", self.row, self.col_map)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 100)

    def test_uses_cached_genomic_alternative_without_ensembl_call(self):
        job = _make_job(hgvs_alternatives={
            "NM_007294.4:c.5379A>C": ["NC_000017.11:g.43045712T>C"]
        })
        with patch("utilities.ensembl_hgvs.EnsemblHGVS") as mock_cls:
            result = self.resolve(job, "job-1", "NM_007294.4:c.5379A>C", self.row, self.col_map)
            mock_cls.assert_not_called()
        self.assertEqual(result["chromosome"], "17")
        self.assertEqual(result["position"], 43045712)

    def test_falls_back_to_ensembl_when_no_genomic_alternative_cached(self):
        job = _make_job(hgvs_config=HGVSConfig(assembly="GRCh38"))
        mock_ensembl = MagicMock()
        mock_ensembl.get_equivalents_batch.return_value = {
            "NM_007294.4:c.5379A>C": ["NC_000017.11:g.43045712T>C"]
        }
        with patch("utilities.ensembl_hgvs.EnsemblHGVS", return_value=mock_ensembl), \
             patch("api.routes.job_manager") as mock_job_manager:
            result = self.resolve(job, "job-1", "NM_007294.4:c.5379A>C", self.row, self.col_map)
            mock_job_manager.add_job_hgvs_alternatives.assert_called_once_with(
                "job-1", "NM_007294.4:c.5379A>C", ["NC_000017.11:g.43045712T>C"]
            )
        self.assertEqual(result["chromosome"], "17")
        self.assertEqual(result["position"], 43045712)

    def test_ensembl_failure_is_swallowed_and_returns_none(self):
        job = _make_job(hgvs_config=HGVSConfig(assembly="GRCh38"))
        with patch("utilities.ensembl_hgvs.EnsemblHGVS", side_effect=RuntimeError("network down")):
            result = self.resolve(job, "job-1", "NM_007294.4:c.5379A>C", self.row, self.col_map)
        self.assertIsNone(result)

    def test_no_matching_candidate_returns_none(self):
        job = _make_job(hgvs_alternatives={"NM_007294.4:c.5379A>C": ["not-a-genomic-id"]})
        result = self.resolve(job, "job-1", "NM_007294.4:c.5379A>C", self.row, self.col_map)
        self.assertIsNone(result)


class TestResolveGenomicVariantFromReference(unittest.TestCase):
    def setUp(self):
        from api.routes import _resolve_genomic_variant_from_reference
        self.resolve = _resolve_genomic_variant_from_reference
        self.row = ["T", "C", "BRCA1"]
        self.col_map = {"REF": 0, "ALT": 1, "GENE": 2}

    def test_ncbi_genomic_reference_resolves_chromosome(self):
        job = _make_job(reference={"type": "ncbi", "value": "NC_000017.11"})
        result = self.resolve(job, self.row, self.col_map, 43045712)
        self.assertEqual(result, {
            "chromosome": "17", "position": 43045712,
            "ref": "T", "alt": "C", "gene": "BRCA1",
        })

    def test_non_ncbi_reference_returns_none(self):
        job = _make_job(reference={"type": "file", "value": "ref.fasta"})
        result = self.resolve(job, self.row, self.col_map, 100)
        self.assertIsNone(result)

    def test_ncbi_transcript_reference_without_nc_prefix_returns_none(self):
        # e.g. an NG_ or NM_ reference - not a genomic NC_0000XX accession.
        job = _make_job(reference={"type": "ncbi", "value": "NG_008866.1"})
        result = self.resolve(job, self.row, self.col_map, 100)
        self.assertIsNone(result)


class TestUpdateVariantStatusEndToEnd(unittest.TestCase):
    """
    Exercises both call sites (patient-scoped row index and grouped genomic
    position) through the real route, confirming both branches still reach
    the shared resolver and produce a persisted ApprovedVariant.
    """

    def setUp(self):
        # Sandbox the DB so these tests never touch the real user
        # bio-engine.db - override the get_db dependency with an isolated
        # in-memory SQLite database for the duration of the test.
        test_engine = create_engine(
            "sqlite:///:memory:",
            connect_args={"check_same_thread": False},
            poolclass=StaticPool,
        )
        Base.metadata.create_all(bind=test_engine)
        TestSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)

        def _override_get_db():
            db = TestSessionLocal()
            try:
                yield db
            finally:
                db.close()

        app.dependency_overrides[get_db] = _override_get_db
        self.addCleanup(app.dependency_overrides.pop, get_db, None)
        self.client = TestClient(app)

    def _job_with_results(self, variant_statuses):
        return _make_job(
            reference={"type": "ncbi", "value": "NC_000017.11"},
            results=[{
                "id": "patient-1",
                "alignment": {
                    "variants": {
                        "columns": ["POS", "REF", "ALT", "GENE", "HGVS"],
                        "rows": [
                            [43045712, "T", "C", "BRCA1", "NC_000017.11:g.43045712T>C"],
                        ],
                    }
                },
            }],
            variant_statuses=variant_statuses,
        )

    def test_patient_scoped_variant_key_persists_approved_variant(self):
        job = self._job_with_results({"patient-1:0": "approved"})
        with patch("api.routes.job_manager.update_variant_status", return_value=job):
            resp = self.client.put(
                "/jobs/job-1/variant-status",
                json={"variant_key": "patient-1:0", "status": "approved"},
            )
        self.assertEqual(resp.status_code, 200)

        resp = self.client.get("/variants/approved", params={"assembly": "GRCh38"})
        self.assertEqual(resp.status_code, 200)
        matches = [v for v in resp.json() if v["job_id"] == "job-1" and v["chromosome"] == "17" and v["position"] == 43045712]
        self.assertEqual(len(matches), 1)

    def test_grouped_variant_key_persists_approved_variant(self):
        job = self._job_with_results({"43045712": "approved"})
        with patch("api.routes.job_manager.update_variant_status", return_value=job):
            resp = self.client.put(
                "/jobs/job-1/variant-status",
                json={"variant_key": "43045712", "status": "approved"},
            )
        self.assertEqual(resp.status_code, 200)

        resp = self.client.get("/variants/approved", params={"assembly": "GRCh38"})
        self.assertEqual(resp.status_code, 200)
        matches = [v for v in resp.json() if v["job_id"] == "job-1" and v["chromosome"] == "17" and v["position"] == 43045712]
        self.assertEqual(len(matches), 1)


if __name__ == "__main__":
    unittest.main()
