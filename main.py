# SPDX-License-Identifier: MIT
"""
Bio-Engine Sidecar Entry Point
==============================

This module is the entry point for the PS Analyzer Bio-Engine sidecar.
It initializes the FastAPI application, sets up logging, CORS, and
global exception handlers, and configures Uvicorn for serving the API.
"""

from core.logging import setup_logging

setup_logging()

import argparse
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
from core.security import TRUSTED_ORIGIN_REGEX

logger = logging.getLogger(__name__)

if getattr(sys, 'frozen', False):
    bundle_dir = sys._MEIPASS
    exe_dir = os.path.dirname(sys.executable)
    logger.info("Running in frozen mode.")
    logger.info(f"Bundle dir (_MEIPASS): {bundle_dir}")
    logger.info(f"Executable dir: {exe_dir}")

    # Add bundle_dir to PATH (for internal PyInstaller files)
    os.environ["PATH"] = bundle_dir + os.pathsep + os.environ["PATH"]

    # Try to find tracy and bgzip in the executable's directory
    # They might have the ps-analyzer- prefix and the target triple
    target_triple = "x86_64-pc-windows-msvc" if sys.platform == "win32" else "x86_64-unknown-linux-gnu"

    for tool_name, setting_attr in [("tracy", "tracy_path"), ("bgzip", "bgzip_path"), ("samtools", "samtools_path")]:
        potential_names = [
            f"ps-analyzer-{tool_name}-{target_triple}.exe" if sys.platform == "win32" else f"ps-analyzer-{tool_name}-{target_triple}",
            f"ps-analyzer-{tool_name}.exe" if sys.platform == "win32" else f"ps-analyzer-{tool_name}",
            f"{tool_name}.exe" if sys.platform == "win32" else tool_name
        ]

        for name in potential_names:
            path = os.path.join(exe_dir, name)
            if os.path.exists(path):
                logger.info(f"Auto-detected {tool_name} sidecar at: {path}")
                setattr(settings, setting_attr, path)
                break

app = FastAPI(title="Bio-Engine Sidecar")

@app.on_event("startup")
def startup_event():
    """
    FastAPI startup event to configure and initialize OpenCRAVAT if available.
    """
    try:
        from services.opencravat import init_opencravat
        init_opencravat(oc_path=settings.oc_path)
    except Exception as e:
        logger.error(f"Failed to initialize OpenCRAVAT on startup: {e}")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost",
        "http://127.0.0.1",
        "tauri://localhost"
    ],
    allow_origin_regex=TRUSTED_ORIGIN_REGEX,
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
        status_code=exc.status_code,
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

    parser = argparse.ArgumentParser(description="Bio-Engine Sidecar")
    parser.add_argument("--tracy-path", type=str, help="Path to the tracy binary")
    parser.add_argument("--samtools-path", type=str, help="Path to the samtools binary")
    parser.add_argument("--bgzip-path", type=str, help="Path to the bgzip binary")
    parser.add_argument("--resource-path", type=str, help="Path to the application resources (where DLLs are located)")
    parser.add_argument("--data-dir", type=str, help="Path to the application data directory (logs, jobs, cache)")
    parser.add_argument("--host", type=str, help="Host to bind the FastAPI server")
    parser.add_argument("--port", type=int, help="Port to bind the FastAPI server")
    args, unknown = parser.parse_known_args()

    def clean_win_path(p):
        if p and p.startswith('\\\\?\\'):
            return p[4:]
        return p

    if args.resource_path:
        resource_path = clean_win_path(args.resource_path)
        logger.info(f"Adding resource path to PATH: {resource_path}")

        # We search for potential DLL locations relative to the resource_path
        # Best Practice: Cover multiple common Tauri 2 directory patterns
        potential_dll_dirs = [
             resource_path,                                    # Root
             os.path.join(resource_path, "binaries"),          # binaries folder
             os.path.join(resource_path, "resources"),         # resources folder
             os.path.join(resource_path, "resources", "binaries") # resources/binaries folder
        ]

        for p in potential_dll_dirs:
             if os.path.isdir(p):
                 if p not in os.environ["PATH"]:
                     os.environ["PATH"] = p + os.pathsep + os.environ["PATH"]
                     logger.info(f"Added to PATH: {p}")
             else:
                 logger.debug(f"Path does not exist, skipping: {p}")

    # Also add the directory of the currently running executable to PATH
    # This ensures sidecar-relative DLLs are found.
    # We use both sys.executable and sys._MEIPASS (if frozen)
    if getattr(sys, 'frozen', False):
        exe_dir = os.path.dirname(sys.executable)
        if exe_dir not in os.environ["PATH"]:
            os.environ["PATH"] = exe_dir + os.pathsep + os.environ["PATH"]
            logger.info(f"Added executable dir to PATH: {exe_dir}")

        bundle_dir = getattr(sys, '_MEIPASS', None)
        if bundle_dir and bundle_dir not in os.environ["PATH"]:
             os.environ["PATH"] = bundle_dir + os.pathsep + os.environ["PATH"]
             logger.info(f"Added bundle dir to PATH: {bundle_dir}")

    # FINAL CHECK: Look for where zlib1.dll actually is in the app directory
    # (Aggressive discovery fallback)
    if sys.platform == "win32" and args.resource_path:
        found_zlib = False
        for root, _dirs, files in os.walk(clean_win_path(args.resource_path)):
            if "zlib1.dll" in files:
                if root not in os.environ["PATH"]:
                    os.environ["PATH"] = root + os.pathsep + os.environ["PATH"]
                    logger.info(f"Discovered and added DLL folder to PATH: {root}")
                    found_zlib = True
                    break
        if not found_zlib:
            logger.warning("Could not find zlib1.dll in resource path. Sidecars might fail.")

    # Log the final path on Windows for debugging
    if sys.platform == "win32":
        logger.info(f"Final environment PATH contains {len(os.environ['PATH'].split(os.pathsep))} entries.")
        # Only log full PATH in debug mode to avoid log pollution
        logger.debug(f"Full effective PATH: {os.environ['PATH']}")

    logger.info(f"Received CLI arguments: {sys.argv}")
    logger.info(f"Unknown arguments: {unknown}")
    logger.info(f"Environment BIO_SAMTOOLS_PATH: {os.getenv('BIO_SAMTOOLS_PATH')}")
    logger.info(f"Environment BIO_BGZIP_PATH: {os.getenv('BIO_BGZIP_PATH')}")
    logger.info(f"Environment TRACY_PATH: {os.getenv('TRACY_PATH')}")

    if args.tracy_path:
        logger.info(f"Overriding tracy_path from CLI: {args.tracy_path}")
        settings.tracy_path = args.tracy_path

    if args.samtools_path:
        logger.info(f"Overriding samtools_path from CLI: {args.samtools_path}")
        settings.samtools_path = args.samtools_path

    if args.bgzip_path:
        logger.info(f"Overriding bgzip_path from CLI: {args.bgzip_path}")
        settings.bgzip_path = args.bgzip_path

    if args.data_dir:
        data_dir = clean_win_path(args.data_dir)
        logger.info(f"Overriding data directory from CLI: {data_dir}")
        settings.logs_dir = os.path.join(data_dir, "logs")
        settings.jobs_dir = os.path.join(data_dir, "jobs")
        settings.cache_dir = os.path.join(data_dir, "ncbi_cache")
        settings.uploads_dir = os.path.join(data_dir, "uploads")

        # Ensure directories exist
        os.makedirs(settings.logs_dir, exist_ok=True)
        os.makedirs(settings.jobs_dir, exist_ok=True)
        os.makedirs(settings.cache_dir, exist_ok=True)
        os.makedirs(settings.uploads_dir, exist_ok=True)

        # Re-initialize logging with new path
        setup_logging()
        logger.info(f"Logging re-initialized at: {settings.logs_dir}")

    if args.host:
        logger.info(f"Overriding host from CLI: {args.host}")
        settings.host = args.host

    if args.port:
        logger.info(f"Overriding port from CLI: {args.port}")
        settings.port = args.port

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
