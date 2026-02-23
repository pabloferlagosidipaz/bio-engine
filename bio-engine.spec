# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all
import sys
import os

# Get the environment prefix dynamically
conda_prefix = sys.prefix

datas = []
binaries = []

import shutil

# Try to find Tracy in the PATH or Conda environment
tracy_path = shutil.which('tracy')
if not tracy_path and os.path.exists(os.path.join(conda_prefix, 'bin', 'tracy')):
    tracy_path = os.path.join(conda_prefix, 'bin', 'tracy')

if tracy_path:
    print(f"Bundling Tracy binary from: {tracy_path}")
    binaries.append((tracy_path, '.'))
else:
    print("WARNING: tracy binary not found. It will not be bundled.")

# Optionally include these shared libraries only if they exist (Linux/Conda specific fixes)
for lib in ['libcrypto.so.3', 'libssl.so.3']:
    lib_path = os.path.join(conda_prefix, 'lib', lib)
    if os.path.exists(lib_path):
        binaries.append((lib_path, '.'))

hiddenimports = [
    'fastapi', 
    'uvicorn', 
    'uvicorn.logging',
    'uvicorn.loops',
    'uvicorn.loops.auto',
    'uvicorn.protocols.http.auto',
    'uvicorn.protocols.http.httptools_impl', 
    'uvicorn.protocols.http.h11_impl', 
    'uvicorn.protocols.websockets.auto',
    'uvicorn.protocols.websockets.websockets_impl', 
    'uvicorn.lifespan.on',
    'psycopg2',
    'pkg_resources'
]

# Collect all dependencies
packages_to_collect = [
    'fastapi', 'uvicorn', 'Bio', 'hgvs', 'ometa', 'parsley', 'terml', 
    'psycopg2', 'bioutils', 'pkg_resources'
]

for package in packages_to_collect:
    try:
        tmp_ret = collect_all(package)
        datas += tmp_ret[0]
        binaries += tmp_ret[1]
        hiddenimports += tmp_ret[2]
    except Exception as e:
        print(f"WARNING: Could not collect {package}: {e}")

a = Analysis(
    ['main.py'],
    pathex=['.', os.path.join(conda_prefix, 'lib')],
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
