import logging
import sys

import orjson

from core.config import settings


class ORJSONFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        log_record = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
        }

        if record.exc_info:
            log_record["exception"] = self.formatException(record.exc_info)

        # Merge extra attributes
        if hasattr(record, "payload"):
            log_record.update(record.payload)

        return orjson.dumps(log_record).decode("utf-8")

def setup_logging():
    handlers = []

    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    if settings.json_logs:
        console_handler.setFormatter(ORJSONFormatter())
    else:
        console_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
    handlers.append(console_handler)

    # File Handler
    log_file = f"{settings.logs_dir}/bio-engine.log"
    file_handler = logging.FileHandler(log_file)
    if settings.json_logs:
        file_handler.setFormatter(ORJSONFormatter())
    else:
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
    handlers.append(file_handler)

    # Root Logger Configuration
    logging.basicConfig(
        level=getattr(logging, settings.log_level.upper(), logging.INFO),
        handlers=handlers,
        force=True
    )

    # Silence uvicorn access logs if we are strictly controlling logs,
    # but usually we want them. We can just let them propagate.
