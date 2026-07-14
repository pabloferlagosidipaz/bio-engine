# SPDX-License-Identifier: MIT
"""
OpenCRAVAT management service.

Wraps the `oc` CLI for module listing, install/uninstall, disk usage reporting,
and ad-hoc variant annotation. Install tasks run in background threads so the
API returns immediately; poll /opencravat/tasks/{task_id} for progress.
"""

import logging
import os
import re
import shutil
import subprocess
import threading
import uuid
from pathlib import Path
from typing import Any

from core.exceptions import BioEngineError

logger = logging.getLogger(__name__)

# In-memory install/uninstall task registry. Survives for the process lifetime.
_tasks: dict[str, dict[str, Any]] = {}


def init_opencravat(oc_path: str = "oc"):
    """
    Initializes OpenCRAVAT by configuring the modules directory to a persistent
    location inside the /app/storage folder and installing base modules if needed.
    """
    # Programmatically opt out of OpenCRAVAT metrics/funding email prompts to prevent hangs in CLI execution
    try:
        import yaml
        rc, out, err = _run("config", "system", oc_path=oc_path)
        if rc == 0:
            match = re.search(r"Configuration file path:\s*(.*)", out)
            if match:
                config_path = match.group(1).strip()
                if os.path.exists(config_path):
                    with open(config_path) as f:
                        config_data = yaml.safe_load(f) or {}
                    if not config_data.get("user_email_opt_out") or config_data.get("save_metrics") is not False:
                        config_data["user_email"] = "no-email"
                        config_data["user_email_opt_out"] = True
                        config_data["save_metrics"] = False
                        with open(config_path, "w") as f:
                            yaml.safe_dump(config_data, f)
                        logger.info(f"Successfully opted out of OpenCRAVAT metrics in {config_path}")
    except Exception as e:
        logger.warning(f"Could not programmatically opt-out of OpenCRAVAT metrics: {e}")

    if os.path.exists("/app/storage"):
        target_dir = "/app/storage/open-cravat"
        os.makedirs(target_dir, exist_ok=True)

        logger.info(f"Setting OpenCRAVAT modules directory to persistent storage: {target_dir}")
        rc, out, err = _run("config", "md", target_dir, oc_path=oc_path)
        if rc != 0:
            logger.error(f"Failed to configure OpenCRAVAT modules directory: {err}")
            return

        rc, out, err = _run("module", "ls", oc_path=oc_path)
        if rc == 0:
            if "vcf-converter" not in out:
                logger.info("Initializing OpenCRAVAT base modules (this may take a few seconds)...")
                rc_base, out_base, err_base = _run("module", "install-base", oc_path=oc_path, timeout=300)
                if rc_base != 0:
                    logger.error(f"Failed to run oc module install-base: {err_base}")
                else:
                    logger.info("OpenCRAVAT base modules installed successfully.")
            else:
                logger.info("OpenCRAVAT base modules are already initialized.")
        else:
            logger.warning(f"Could not check OpenCRAVAT modules status: {err}")

_SIZE_PATTERN = re.compile(r"([\d.]+)\s*(B|KB|MB|GB|TB)", re.IGNORECASE)
_SIZE_UNITS = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}


def resolve_oc_executable(oc_path: str) -> str:
    """
    Validates that `oc_path` resolves to a real, executable binary before it
    is trusted as argv[0] of a subprocess call (see CWE for uncontrolled
    executable path). Bare names (e.g. "oc") are resolved via PATH, same as
    subprocess.run would do; explicit paths must point at an existing file.

    Returns the resolved path unchanged (callers keep using whatever form -
    bare name or absolute path - they were given; this only validates it).
    Raises BioEngineError (400) if nothing executable is found.
    """
    resolved = shutil.which(oc_path)
    if not resolved:
        raise BioEngineError(
            f"'{oc_path}' does not resolve to an executable binary",
            context="oc_path_validation",
            status_code=400,
        )
    return oc_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _run(
    *args: str,
    oc_path: str = "oc",
    timeout: int = 60,
    stdin_input: str | None = None,
) -> tuple[int, str, str]:
    try:
        result = subprocess.run(
            [oc_path, *args],
            capture_output=True,
            text=True,
            timeout=timeout,
            input=stdin_input,
        )
        return result.returncode, result.stdout.strip(), result.stderr.strip()
    except FileNotFoundError:
        return (
            -1,
            "",
            f"OpenCRAVAT binary not found at '{oc_path}'. "
            "Install with: pip install open-cravat && oc module install-base",
        )
    except subprocess.TimeoutExpired:
        return -1, "", f"Command timed out after {timeout}s"


def _parse_size_bytes(size_str: str) -> int | None:
    m = _SIZE_PATTERN.search(size_str)
    if not m:
        return None
    return int(float(m.group(1)) * _SIZE_UNITS.get(m.group(2).upper(), 1))


def _dir_size_bytes(path: str) -> int:
    total = 0
    for f in Path(path).rglob("*"):
        if f.is_file():
            try:
                total += f.stat().st_size
            except OSError:
                pass
    return total


def _get_data_dir(oc_path: str) -> str | None:
    rc, out, _ = _run("config", "md", oc_path=oc_path)
    if rc == 0 and out:
        return out
    default = Path.home() / ".open-cravat"
    return str(default) if default.exists() else None


def _parse_module_table_fallback(stdout: str, default_installed: bool | None = None) -> list[dict[str, Any]]:
    known_types = {"annotator", "converter", "postaggregator", "reporter", "mapper", "webviewerwidget"}
    modules: list[dict[str, Any]] = []

    for line in stdout.splitlines():
        line = line.strip()
        if not line or re.match(r"^[-=\s]*$", line):
            continue
        if re.match(r"(?i)^\s*(name|module)\s", line):
            continue

        parts = line.split()
        if not parts:
            continue

        name = parts[0]
        mod: dict[str, Any] = {
            "name": name,
            "version": None,
            "type": None,
            "size_bytes": None,
            "installed": default_installed,
            "title": None,
        }

        for p in parts[1:]:
            if re.match(r"^\d[\d.]{1,}", p):
                mod["version"] = p
                break

        for p in parts[1:]:
            if p.lower() in known_types:
                mod["type"] = p.lower()
                break

        for p in parts[1:]:
            if p.lower() == "yes":
                mod["installed"] = True
                break
            if p.lower() == "no":
                mod["installed"] = False
                break

        for i, p in enumerate(parts[:-1]):
            if re.match(r"^[\d.]+$", p) and parts[i + 1].upper() in _SIZE_UNITS:
                mod["size_bytes"] = _parse_size_bytes(f"{p} {parts[i + 1]}")
                break

        modules.append(mod)

    return modules


def _parse_module_table(stdout: str, default_installed: bool | None = None) -> list[dict[str, Any]]:
    """
    Parse the tabular output of `oc module ls` or `oc module ls -a`.
    """
    header = None
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"(?i)^\s*(name|module)\s", line):
            header = line
            break

    if not header:
        return _parse_module_table_fallback(stdout, default_installed)

    lower_header = header.lower()
    idx_name = lower_header.find("name")
    idx_title = lower_header.find("title")
    idx_type = lower_header.find("type")
    idx_installed = lower_header.find("installed")

    idx_version = lower_header.find("version")
    if idx_version == -1:
        idx_version = lower_header.find("store ver")
    if idx_version == -1:
        idx_version = lower_header.find("local ver")

    idx_size = lower_header.find("size")

    cols = []
    if idx_name != -1: cols.append(("name", idx_name))
    if idx_title != -1: cols.append(("title", idx_title))
    if idx_type != -1: cols.append(("type", idx_type))
    if idx_installed != -1: cols.append(("installed", idx_installed))
    if idx_version != -1: cols.append(("version", idx_version))
    if idx_size != -1: cols.append(("size", idx_size))

    cols.sort(key=lambda x: x[1])

    modules: list[dict[str, Any]] = []

    for line in stdout.splitlines():
        if not line.strip() or re.match(r"^[-=\s]*$", line):
            continue
        if re.match(r"(?i)^\s*(name|module)\s", line):
            continue

        mod: dict[str, Any] = {
            "name": "",
            "version": None,
            "type": None,
            "size_bytes": None,
            "installed": default_installed,
            "title": None,
        }

        for i, (col_name, start_idx) in enumerate(cols):
            end_idx = cols[i+1][1] if i + 1 < len(cols) else len(line)
            val = line[start_idx:end_idx].strip()

            if col_name == "name":
                mod["name"] = val
            elif col_name == "title":
                mod["title"] = val if val else None
            elif col_name == "type":
                mod["type"] = val if val else None
            elif col_name == "installed":
                if val.lower() in ("yes", "installed"):
                    mod["installed"] = True
                elif val.lower() == "no":
                    mod["installed"] = False
            elif col_name == "version":
                version_parts = val.split()
                mod["version"] = version_parts[0] if version_parts else None
            elif col_name == "size":
                mod["size_bytes"] = _parse_size_bytes(val)

        if not mod["name"]:
            continue

        if mod["installed"] is None and default_installed is not None:
            mod["installed"] = default_installed
        elif mod["installed"] is None and default_installed is None:
            idx_local = lower_header.find("local ver")
            if idx_local != -1:
                idx_local_data = lower_header.find("local data ver")
                end_local = idx_local_data if idx_local_data != -1 else lower_header.find("size")
                local_val = line[idx_local:end_local].strip()
                if local_val:
                    mod["installed"] = True
                    mod["version"] = local_val.split()[0] if local_val.split() else mod["version"]
                else:
                    mod["installed"] = False
            else:
                mod["installed"] = False

        modules.append(mod)

    return modules


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_status(oc_path: str = "oc") -> dict[str, Any]:
    rc, version_out, err = _run("version", oc_path=oc_path)
    if rc != 0:
        return {"installed": False, "error": err}

    data_dir = _get_data_dir(oc_path)
    disk_used_bytes: int | None = None
    disk_free_bytes: int | None = None

    if data_dir and Path(data_dir).exists():
        disk_used_bytes = _dir_size_bytes(data_dir)
        disk_free_bytes = shutil.disk_usage(data_dir).free

    return {
        "installed": True,
        "version": version_out,
        "data_dir": data_dir,
        "disk_used_bytes": disk_used_bytes,
        "disk_free_bytes": disk_free_bytes,
        "error": None,
    }


def list_installed(oc_path: str = "oc") -> list[dict[str, Any]]:
    rc, stdout, stderr = _run("module", "ls", oc_path=oc_path, timeout=30)
    if rc != 0:
        logger.error(f"oc module ls failed: {stderr}")
        return []
    return _parse_module_table(stdout, default_installed=True)


def list_store(oc_path: str = "oc") -> list[dict[str, Any]]:
    """List all modules available in the OC store (installed + available)."""
    rc, stdout, stderr = _run("module", "ls", "-a", oc_path=oc_path, timeout=60)
    if rc != 0:
        logger.error(f"oc module ls -a failed: {stderr}")
        return []
    return _parse_module_table(stdout)


def get_task(task_id: str) -> dict[str, Any] | None:
    return _tasks.get(task_id)


def _do_install(task_id: str, module_name: str, oc_path: str) -> None:
    _tasks[task_id]["status"] = "running"
    rc, _, stderr = _run(
        "module", "install", "-y", module_name,
        oc_path=oc_path,
        timeout=3600,   # large installs (gnomAD ~60 GB) can take a long time
        stdin_input="y\n",
    )
    if rc == 0:
        _tasks[task_id]["status"] = "completed"
        logger.info(f"OpenCRAVAT module '{module_name}' installed successfully")
    else:
        _tasks[task_id]["status"] = "failed"
        _tasks[task_id]["error"] = stderr
        logger.error(f"OpenCRAVAT install '{module_name}' failed: {stderr[:300]}")


def install_module(module_name: str, oc_path: str = "oc") -> str:
    """Start a background install. Returns task_id to poll for status."""
    task_id = uuid.uuid4().hex[:8]
    _tasks[task_id] = {
        "task_id": task_id,
        "module": module_name,
        "status": "pending",
        "error": None,
    }
    threading.Thread(
        target=_do_install, args=(task_id, module_name, oc_path), daemon=True
    ).start()
    return task_id


def uninstall_module(module_name: str, oc_path: str = "oc") -> dict[str, Any]:
    rc, _, stderr = _run("module", "uninstall", "-y", module_name, oc_path=oc_path, timeout=120)
    if rc == 0:
        logger.info(f"OpenCRAVAT module '{module_name}' uninstalled")
        return {"success": True, "error": None}
    return {"success": False, "error": stderr}


def list_tasks() -> list[dict[str, Any]]:
    """Returns all registered installation tasks."""
    return list(_tasks.values())
