#!/usr/bin/env python
"""
Prepare thermodynamic database files for aqequil wheel builds.

This script:
1. Downloads the latest database files from WORM-db repository:
   - data0.wrm (thermodynamic database source)
   - wrm_data_latest.csv (main thermodynamic data)
   - elements.csv (element properties)
   - solid_solutions.csv (solid solution parameters)
   - wrm_data_logK.csv (equilibrium constants)
   - wrm_data_logK_S.csv (equilibrium constants for S species)
2. Runs EQPT to compile data0.wrm into data1.wrm
3. Copies all files to the aqequil/databases directory

This script should be run AFTER compile_eq3_6.py (requires EQPT executable).
"""

import os
import sys
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.request import urlopen


WORM_DB_BASE_URL = "https://raw.githubusercontent.com/worm-portal/WORM-db/master"

# Database files to download from WORM-db
DATABASE_FILES = {
    "data0.wrm": f"{WORM_DB_BASE_URL}/data0.wrm",
    "wrm_data_latest.csv": f"{WORM_DB_BASE_URL}/wrm_data_latest.csv",
    "elements.csv": f"{WORM_DB_BASE_URL}/elements.csv",
    "solid_solutions.csv": f"{WORM_DB_BASE_URL}/solid_solutions.csv",
    "wrm_data_logK.csv": f"{WORM_DB_BASE_URL}/wrm_data_logK.csv",
    "wrm_data_logK_S.csv": f"{WORM_DB_BASE_URL}/wrm_data_logK_S.csv",
    "speciation_groups_WORM.txt": f"{WORM_DB_BASE_URL}/speciation_groups_WORM.txt",
}


def download_file(filename, url, dest_path):
    """Download a file from a URL."""
    print(f"Downloading {filename} from {url}...")

    try:
        with urlopen(url) as response:
            content = response.read()

        with open(dest_path, 'wb') as f:
            f.write(content)

        size_kb = len(content) / 1024
        print(f"[OK] Downloaded {filename} ({size_kb:.1f} KB)")
        return True

    except Exception as e:
        print(f"[ERROR] Failed to download {filename}: {e}", file=sys.stderr)
        return False


def download_all_databases(dest_dir):
    """Download all database files from WORM-db."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    success = True
    downloaded_files = {}

    for filename, url in DATABASE_FILES.items():
        dest_path = dest_dir / filename
        if download_file(filename, url, dest_path):
            downloaded_files[filename] = dest_path
        else:
            success = False

    return success, downloaded_files


def find_eqpt():
    """Find the EQPT executable."""
    script_dir = Path(__file__).parent
    bin_dir = script_dir / "aqequil" / "bin"

    if sys.platform == "win32":
        eqpt_path = bin_dir / "eqpt.exe"
    else:
        eqpt_path = bin_dir / "eqpt"

    if eqpt_path.exists():
        return str(eqpt_path)

    # Check if EQPT is in PATH
    eqpt_in_path = shutil.which("eqpt")
    if eqpt_in_path:
        return eqpt_in_path

    return None


def run_eqpt(data0_path, work_dir):
    """Run EQPT to convert data0 to data1."""
    eqpt_path = find_eqpt()

    if not eqpt_path:
        print("[ERROR] EQPT executable not found. Run compile_eq3_6.py first.", file=sys.stderr)
        return False

    print(f"Running EQPT on {data0_path}...")
    print(f"Using EQPT at: {eqpt_path}")

    try:
        # Run EQPT with the data0 file path as argument
        # EQPT output is decoded explicitly rather than with text=True. text=True
        # uses the locale codec and raises UnicodeDecodeError on any byte it
        # cannot handle, which would be reported as a failure to run EQPT even
        # though EQPT itself ran fine. backslashreplace keeps undecodable bytes
        # visible as escapes and keeps the result printable on a non-UTF-8
        # console.
        result = subprocess.run(
            [eqpt_path, str(data0_path)],
            cwd=work_dir,
            capture_output=True,
            encoding='utf-8',
            errors='backslashreplace',
            timeout=120
        )

        # EQPT creates output files in the working directory
        # Check for data1 output (could be "data1" or "data0.d1")
        data1_path = Path(work_dir) / "data1"
        data0_d1_path = Path(work_dir) / "data0.d1"

        if data1_path.exists():
            return str(data1_path)
        elif data0_d1_path.exists():
            # Rename to data1
            final_path = Path(work_dir) / "data1"
            data0_d1_path.rename(final_path)
            return str(final_path)
        else:
            print("[ERROR] EQPT did not produce data1 output", file=sys.stderr)
            if result.stdout:
                print(f"stdout: {result.stdout}", file=sys.stderr)
            if result.stderr:
                print(f"stderr: {result.stderr}", file=sys.stderr)
            return None

    except subprocess.TimeoutExpired:
        print("[ERROR] EQPT timed out", file=sys.stderr)
        return None
    except Exception as e:
        print(f"[ERROR] Failed to run EQPT: {e}", file=sys.stderr)
        return None


def prepare_databases():
    """Main function to prepare database files."""
    script_dir = Path(__file__).parent
    databases_dir = script_dir / "aqequil" / "databases"
    test_data_dir = script_dir / "aqequil" / "test_data"

    # Ensure directories exist
    databases_dir.mkdir(parents=True, exist_ok=True)
    test_data_dir.mkdir(parents=True, exist_ok=True)

    # Work in a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Download all database files
        print("Downloading database files from WORM-db...")
        success, downloaded_files = download_all_databases(tmpdir)

        if not success:
            print("[ERROR] Failed to download some database files", file=sys.stderr)
            return False

        # Run EQPT to create data1 from data0
        data0_tmp = downloaded_files.get("data0.wrm")
        if data0_tmp:
            print("\nCompiling data0.wrm to data1.wrm with EQPT...")
            data1_tmp = run_eqpt(data0_tmp, tmpdir)

            if not data1_tmp:
                print("[ERROR] Failed to create data1.wrm", file=sys.stderr)
                return False

            data1_tmp = Path(data1_tmp)
            size_kb = data1_tmp.stat().st_size / 1024
            print(f"[OK] Created data1.wrm ({size_kb:.1f} KB)")

            # Add data1.wrm to the files to copy
            downloaded_files["data1.wrm"] = data1_tmp

        # Copy all files to databases directory
        print("\nCopying files to databases directory...")
        for filename, src_path in downloaded_files.items():
            dest_path = databases_dir / filename
            shutil.copy(src_path, dest_path)
            print(f"[OK] Copied {filename} to {dest_path}")

        # Also copy data0.wrm to test_data directory
        if "data0.wrm" in downloaded_files:
            test_data0_dest = test_data_dir / "data0.wrm"
            shutil.copy(downloaded_files["data0.wrm"], test_data0_dest)
            print(f"[OK] Copied data0.wrm to {test_data0_dest}")

    return True


def main():
    """Main entry point."""
    # EQPT output is relayed to these streams, and a UnicodeEncodeError while
    # printing it would hide the reason the database preparation failed
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")

    print("=" * 60)
    print("Preparing thermodynamic databases for aqequil")
    print("=" * 60)
    print(f"Platform: {sys.platform}")
    print(f"Python: {sys.version.split()[0]}")

    try:
        success = prepare_databases()
        if success:
            print("=" * 60)
            print("[OK] Database preparation successful!")
            print("=" * 60)
            return 0
        else:
            print("=" * 60)
            print("[ERROR] Database preparation failed!")
            print("=" * 60)
            return 1
    except Exception as e:
        print("=" * 60, file=sys.stderr)
        print(f"[ERROR] {e}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
