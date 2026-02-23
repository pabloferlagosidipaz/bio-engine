import os

from pydantic_settings import BaseSettings, SettingsConfigDict

from core.persistence import PersistenceManager


class Settings(BaseSettings):
    app_name: str = "Bio-Engine Sidecar"

    # Network
    host: str = "127.0.0.1"
    port: int = 8000

    # External Services
    entrez_email: str = "example@example.com"

    # Paths (Defaults handled by PersistenceManager if not set)
    # We use a factory to delay PersistenceManager instantiation if possible,
    # but for simplicity we can instantiate it here as it just resolves paths.
    _persistence: PersistenceManager = PersistenceManager()

    cache_dir: str = os.getenv("BIO_CACHE_DIR", _persistence.get_cache_dir())
    jobs_dir: str = _persistence.get_jobs_dir()
    logs_dir: str = _persistence.get_logs_dir()

    # Logging
    log_level: str = "INFO"
    json_logs: bool = True

    model_config = SettingsConfigDict(env_prefix="BIO_")

settings = Settings()

