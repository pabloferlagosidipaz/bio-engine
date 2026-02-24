# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all, collect_submodules
import sys
import os
import site

# Get the environment prefix dynamically
conda_prefix = sys.prefix

# Get the actual site-packages directory
site_packages = site.getsitepackages()
print(f"INFO: sys.prefix = {sys.prefix}")
print(f"INFO: site-packages = {site_packages}")

datas = []
binaries = []

# Tracy will be handled by Tauri, not bundled here

# Optionally include these shared libraries only if they exist (Linux/Conda specific fixes)
for lib in ['libcrypto.so.3', 'libssl.so.3']:
    lib_path = os.path.join(conda_prefix, 'lib', lib)
    if os.path.exists(lib_path):
        binaries.append((lib_path, '.'))

# Verify that critical packages can be imported before proceeding
import importlib
for pkg in ['uvicorn', 'fastapi', 'h11', 'starlette']:
    try:
        mod = importlib.import_module(pkg)
        print(f"INFO: {pkg} found at {mod.__file__}")
    except ImportError as e:
        print(f"ERROR: {pkg} NOT FOUND: {e}")

hiddenimports = [
    'fastapi', 
    'uvicorn', 
    'uvicorn.logging',
    'uvicorn.loops',
    'uvicorn.loops.auto',
    'uvicorn.protocols',
    'uvicorn.protocols.http',
    'uvicorn.protocols.http.auto',
    'uvicorn.protocols.http.httptools_impl', 
    'uvicorn.protocols.http.h11_impl', 
    'uvicorn.protocols.websockets',
    'uvicorn.protocols.websockets.auto',
    'uvicorn.protocols.websockets.websockets_impl', 
    'uvicorn.lifespan',
    'uvicorn.lifespan.on',
    'uvicorn.lifespan.off',
    'uvicorn.config',
    'uvicorn.main',
    'uvicorn.server',
    'uvicorn._types',
    'h11',
    'h11._connection',
    'h11._events',
    'h11._state',
    'h11._readers',
    'h11._writers',
    'h11._util',
    'h11._abnf',
    'h11._headers',
    'h11._receivebuffer',
    'click',
    'anyio',
    'anyio._backends',
    'anyio._backends._asyncio',
    'sniffio',
    'starlette',
    'starlette.routing',
    'starlette.responses',
    'starlette.requests',
    'starlette.middleware',
    'starlette.middleware.cors',
    'psycopg2',
    'pkg_resources'
]

# Collect all submodules explicitly for critical packages
for pkg in ['uvicorn', 'fastapi', 'starlette', 'h11', 'anyio']:
    try:
        submods = collect_submodules(pkg)
        hiddenimports += submods
        print(f"INFO: Collected {len(submods)} submodules for {pkg}")
    except Exception as e:
        print(f"WARNING: Could not collect submodules for {pkg}: {e}")

# Collect all dependencies (data, binaries, and hidden imports)
packages_to_collect = [
    'fastapi', 'uvicorn', 'h11', 'click', 'anyio', 'sniffio', 'starlette',
    'Bio', 'hgvs', 'ometa', 'parsley', 'terml', 
    'psycopg2', 'bioutils', 'pkg_resources'
]

for package in packages_to_collect:
    try:
        tmp_ret = collect_all(package)
        datas += tmp_ret[0]
        binaries += tmp_ret[1]
        hiddenimports += tmp_ret[2]
        print(f"INFO: collect_all({package}) -> {len(tmp_ret[0])} datas, {len(tmp_ret[1])} binaries, {len(tmp_ret[2])} hiddenimports")
    except Exception as e:
        print(f"WARNING: Could not collect {package}: {e}")

# Build search paths: include site-packages so PyInstaller can find all installed modules
pathex = ['.']
pathex += site_packages
pathex.append(os.path.join(conda_prefix, 'lib'))

a = Analysis(
    ['main.py'],
    pathex=pathex,
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='bio-engine',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
