"""
Job Manager Service
===================

This module provides the `JobManager` class which handles the complete lifecycle
of analysis jobs. It is responsible for creating, updating, retrieving, and 
saving job state using JSON files in the local filesystem.
"""

import logging
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any

import orjson

from core.persistence import PersistenceManager
from data.models import Comment, HGVSConfig, Job, JobStatus, TracyConfig

logger = logging.getLogger(__name__)

class JobManager:
    """
    Manages persistence and state mutations for Bio-Engine jobs.
    
    Reads and writes `Job` models to disk as JSON files, ensuring that status
    updates, annotations, and analysis results are correctly stored.
    """
    def __init__(self, persistence: PersistenceManager | None = None):
        self.persistence = persistence or PersistenceManager()
        self.jobs_dir = Path(self.persistence.get_jobs_dir())
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Initialized JobManager with jobs directory: {self.jobs_dir}")

    def _get_job_path(self, job_id: str) -> Path:
        return self.jobs_dir / f"{job_id}.json"

    def create_job(self, name: str, reference: dict, patients: list[dict], reference_sequence: str | None = None, app_version: str | None = None, config: dict | None = None, hgvs_config: dict | None = None) -> Job:
        """
        Creates a new primary job record and initializes its status to CREATED.
        
        Args:
            name (str): Label for the job.
            reference (dict): Reference definition.
            patients (list[dict]): Patient definitions including sequence files.
            reference_sequence (str | None): Full nucleotide sequence of the reference.
            app_version (str | None): UI client software version.
            config (dict | None): Alignment Tracy configuration.
            hgvs_config (dict | None): HGVS annotation configuration.
            
        Returns:
            Job: The newly created job.
        """
        job_id = str(uuid.uuid4())
        timestamp = datetime.now().isoformat()

        # We construct the model directly to validate inputs
        # But we need to handle the dicts passed from API

        # Construct Job object
        job = Job(
            id=job_id,
            name=name,
            status=JobStatus.CREATED,
            created_at=timestamp,
            updated_at=timestamp,
            reference=reference,
            patients=patients,
            reference_sequence=reference_sequence,
            app_version=app_version,
            config=TracyConfig(**config) if config else None,
            hgvs_config=HGVSConfig(**hgvs_config) if hgvs_config else None
        )

        self._save_job(job)

        logger.info(f"Created job {job_id} ({name})")
        return job

    def get_job(self, job_id: str) -> Job | None:
        """
        Retrieves a job by ID from the filesystem.
        """
        job_path = self._get_job_path(job_id)
        if not job_path.exists():
            return None

        try:
            with open(job_path, "rb") as f:
                data = orjson.loads(f.read())
            return Job(**data)
        except orjson.JSONDecodeError as e:
            file_size = job_path.stat().st_size
            logger.error(f"Failed to parse job {job_id} at {job_path} (size: {file_size} bytes): {e}")
            return None
        except Exception as e:
            logger.error(f"Failed to load job {job_id} from {job_path}: {e}")
            return None

    def update_job_status(self, job_id: str, status: JobStatus):
        # Optimistic update: ideally we'd lock the file
        job = self.get_job(job_id)
        if job:
            job.status = status
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} status to {status.value}")

    def update_job_progress(self, job_id: str, progress: int, message: str | None = None):
        job = self.get_job(job_id)
        if job:
            job.progress = progress
            if message:
                job.status_message = message
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} progress to {progress}%: {message}")

    def update_job(self, job_id: str, name: str, reference: dict, patients: list[dict], config: dict | None = None, hgvs_config: dict | None = None, reference_sequence: str | None = None) -> Job | None:
        """
        Updates the primary configuration properties of an existing job.
        """
        job = self.get_job(job_id)
        if job:
            job.name = name
            job.reference = reference
            job.patients = patients
            if reference_sequence is not None:
                job.reference_sequence = reference_sequence
            if config is not None:
                job.config = TracyConfig(**config) if isinstance(config, dict) else config
            if hgvs_config is not None:
                job.hgvs_config = HGVSConfig(**hgvs_config) if isinstance(hgvs_config, dict) else hgvs_config
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} data")
            return job
        return None

    def update_job_results(self, job_id: str, results: list[dict]) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.results = results
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} results")
            return job
        return None
    def update_job_vep_annotations(self, job_id: str, annotations: dict[str, Any]) -> Job | None:
        job = self.get_job(job_id)
        if job:
            # annotations is a dict keyed by HGVS
            # We convert raw dicts to VariantAnnotation objects if needed (Pydantic handles this usually)
            job.vep_annotations.update(annotations)
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated VEP annotations for job {job_id}")
            return job
        return None

    def update_job_features(self, job_id: str, features: list[dict]) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.features = features
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} features")
            return job
        return None

    def update_job_reference_sequence(self, job_id: str, sequence: str) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.reference_sequence = sequence
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated job {job_id} reference sequence")
            return job
        return None

    def add_job_comment(self, job_id: str, variant_key: str, text: str, author: str) -> Job | None:
        job = self.get_job(job_id)
        if job:
            comment_id = str(uuid.uuid4())
            timestamp = datetime.now().isoformat()
            new_comment = Comment(
                id=comment_id,
                text=text,
                author=author,
                created_at=timestamp
            )

            if variant_key not in job.comments:
                job.comments[variant_key] = []

            job.comments[variant_key].append(new_comment)
            job.updated_at = timestamp
            self._save_job(job)
            logger.info(f"Added comment {comment_id} to job {job_id} for variant {variant_key}")
            return job
        return None

    def add_job_hgvs_alternatives(self, job_id: str, principal_hgvs: str, alternatives: list[str]) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.hgvs_alternatives[principal_hgvs] = alternatives
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated HGVS alternatives for {principal_hgvs} in job {job_id}")
            return job
        return None

    def update_job_hgvs_alternatives_bulk(self, job_id: str, alternatives: dict[str, list[str]]) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.hgvs_alternatives.update(alternatives)
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Updated HGVS alternatives in bulk for job {job_id} ({len(alternatives)} items)")
            return job
        return None

    def delete_job_comment(self, job_id: str, variant_key: str, comment_id: str) -> Job | None:
        job = self.get_job(job_id)
        if job:
            if variant_key in job.comments:
                original_count = len(job.comments[variant_key])
                job.comments[variant_key] = [c for c in job.comments[variant_key] if c.id != comment_id]

                if len(job.comments[variant_key]) < original_count:
                    job.updated_at = datetime.now().isoformat()
                    self._save_job(job)
                    logger.info(f"Deleted comment {comment_id} from job {job_id}")
                    return job
        return None


    def delete_job(self, job_id: str) -> bool:
        job_path = self._get_job_path(job_id)
        if job_path.exists():
            try:
                job_path.unlink()
                logger.info(f"Deleted job {job_id}")
                return True
            except OSError as e:
                logger.error(f"Error deleting job {job_id}: {e}")
                return False
        return False

    def rename_job(self, job_id: str, new_name: str) -> Job | None:
        job = self.get_job(job_id)
        if job:
            job.name = new_name
            job.updated_at = datetime.now().isoformat()
            self._save_job(job)
            logger.info(f"Renamed job {job_id} to {new_name}")
            return job
        return None

    def _save_job(self, job: Job):
        """
        Atomically saves the Job object to a JSON file.
        """
        job_path = self._get_job_path(job.id)
        # Use atomic write via temp file
        temp_path = job_path.with_suffix(".tmp")
        try:
            with open(temp_path, "wb") as f:
                f.write(orjson.dumps(job.model_dump(), option=orjson.OPT_INDENT_2))

            # Atomic move
            temp_path.replace(job_path)

        except Exception as e:
            logger.error(f"Failed to save job {job.id}: {e}")
            if temp_path.exists():
                temp_path.unlink()
            raise

    def list_jobs(self) -> list[Job]:
        """
        Scans the jobs directory and returns a unified list of all tracked jobs.
        Optimizes by deleting the heavy 'results' payload in memory before returning.
        """
        jobs = []
        if not self.jobs_dir.exists():
            return jobs
        for file_path in self.jobs_dir.glob("*.json"):
            try:
                # To check for performance: maybe reading only header?
                # orjson loads is fast.
                with open(file_path, "rb") as f:
                    data = orjson.loads(f.read())
                if 'results' in data:
                    del data['results']

                jobs.append(Job(**data))
            except Exception as e:
                logger.warning(f"Skipping invalid job file {file_path}: {e}")

        jobs.sort(key=lambda x: x.created_at, reverse=True)
        return jobs
