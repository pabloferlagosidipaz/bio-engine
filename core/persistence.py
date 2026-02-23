import logging
"""
Persistence Management Module
=============================

This module provides the `PersistenceManager` class, responsible for 
resolving and creating the necessary local directories (logs, jobs, cache) 
based on the host operating system conventions.
"""

import os
import pathlib
import sys

logger = logging.getLogger(__name__)

class PersistenceManager:
    """
    Manages local storage paths for the application, ensuring that 
    required directories exist across different operating systems.
    """
    def __init__(self, app_name: str = "ms-analyzer"):
        self.app_name = app_name
        self.base_dir = self._get_user_data_dir()
        self._ensure_dirs()

    def _get_user_data_dir(self) -> pathlib.Path:
        home = pathlib.Path.home()
        if sys.platform == "win32":
            return home / "AppData" / "Roaming" / self.app_name
        elif sys.platform == "darwin":
            return home / "Library" / "Application Support" / self.app_name
        else:  # Linux/Unix
            return home / ".local" / "share" / self.app_name

    def _ensure_dirs(self):
        dirs = [
            self.get_jobs_dir(),
            self.get_logs_dir(),
            self.get_cache_dir()
        ]
        for d in dirs:
            try:
                os.makedirs(d, exist_ok=True)
            except OSError as e:
                logger.error(f"Could not create directory {d}: {e}")
                # Fallback to local 'data' directory if permission denied
                fallback = pathlib.Path(__file__).parent.parent / "data" / os.path.basename(d)
                os.makedirs(fallback, exist_ok=True)
                logger.warning(f"Falling back to {fallback}")

    def get_jobs_dir(self) -> str:
        return str(self.base_dir / "jobs")

    def get_logs_dir(self) -> str:
        return str(self.base_dir / "logs")

    def get_cache_dir(self) -> str:
        return str(self.base_dir / "ncbi_cache")

    def get_log_file(self, filename: str = "bio-engine.log") -> str:
        return os.path.join(self.get_logs_dir(), filename)
