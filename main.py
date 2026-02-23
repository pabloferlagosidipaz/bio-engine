"""
Bio-Engine Sidecar Entry Point
==============================

This module is the entry point for the PS Analyzer Bio-Engine sidecar.
It initializes the FastAPI application, sets up logging, CORS, and
global exception handlers, and configures Uvicorn for serving the API.
"""

import logging
import multiprocessing
import os
import signal
import sys

import uvicorn
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from api.routes import router
from core.config import settings
from core.exceptions import BioEngineError
from core.logging import setup_logging

# Handle PyInstaller one-file mode: Add temp dir to PATH so tracy is found
if getattr(sys, 'frozen', False):
    bundle_dir = sys._MEIPASS
    os.environ["PATH"] = bundle_dir + os.pathsep + os.environ["PATH"]

# Initialize Logging
setup_logging()
logger = logging.getLogger(__name__)

if getattr(sys, 'frozen', False):
    logger.info(f"Running in frozen mode. Added {sys._MEIPASS} to PATH. Tracy should be available.")


app = FastAPI(title="Bio-Engine Sidecar")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.exception_handler(BioEngineError)
async def bio_engine_exception_handler(request: Request, exc: BioEngineError):
    """
    Global exception handler for BioEngineError.
    
    Transforms internal application exceptions into a standardized JSON response 
    with a 500 status code, exposing the exception type, message, and context.
    """
    return JSONResponse(
        status_code=500,
        content={
            "type": type(exc).__name__,
            "message": exc.message,
            "context": exc.context
        }
    )

# Include the API router
app.include_router(router)

def signal_handler(sig, frame):
    """
    Handles OS signals for graceful shutdown.
    
    Args:
        sig: The signal number.
        frame: The current stack frame.
    """
    logger.info(f"Received signal {sig}, shutting down...")
    sys.exit(0)

if __name__ == "__main__":
    multiprocessing.freeze_support()

    # Register signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # By default, uvicorn uses its own logging config.
    # To include uvicorn logs in our file, we pass log_config=None to use our root logger setup.
    try:
        uvicorn.run(app, host=settings.host, port=settings.port, log_config=None)
    except Exception as e:
        logger.error(f"Failed to start server: {e}")
        sys.exit(1)
