#!/usr/bin/env python
"""
Compile EQ3/6 executables for the current platform.

This script is called during the wheel build process to compile the EQ3/6
Fortran executables into platform-specific binaries.
"""

import os
import sys
import subprocess
import shutil
import platform
from pathlib import Path


def find_gfortran():
    """Find gfortran compiler on the system."""
    # Check if explicit gfortran path is provided via environment variable
    explicit_path = os.environ.get('GFORTRAN_PATH')
    if explicit_path and os.path.isfile(explicit_path):
        print(f"Using explicit gfortran path from GFORTRAN_PATH: {explicit_path}")
        return explicit_path

    # On Windows, check common MinGW locations first
    if sys.platform == "win32":
        mingw_paths = [
            r'C:\msys64\mingw64\bin\gfortran.exe',
            r'C:\msys64\ucrt64\bin\gfortran.exe',
            r'C:\mingw64\bin\gfortran.exe',
        ]
        for path in mingw_paths:
            if os.path.isfile(path):
                print(f"Found gfortran at: {path}")
                return path

    # Try common gfortran locations
    candidates = ['gfortran', 'gfortran-11', 'gfortran-12', 'gfortran-13', 'gfortran-14']

    for compiler in candidates:
        compiler_path = shutil.which(compiler)
        if compiler_path:
            return compiler_path

    raise RuntimeError(
        "gfortran compiler not found. Please install gfortran:\n"
        "  Ubuntu/Debian: sudo apt-get install gfortran\n"
        "  macOS: brew install gcc\n"
        "  Windows: Install MinGW-w64 with gfortran support"
    )


def compile_eq3_6():
    """Compile EQ3/6 source code to executables."""
    # Determine paths
    script_dir = Path(__file__).parent
    source_dir = script_dir / "aqequil" / "eq3_6_src"
    output_dir = script_dir / "aqequil" / "bin"

    # Verify source directory exists
    if not source_dir.exists():
        raise FileNotFoundError(f"EQ3/6 source not found: {source_dir}")

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find gfortran
    gfortran = find_gfortran()
    print(f"Using compiler: {gfortran}")

    # Get the path to gfortran directory for make
    gfortran_dir = Path(gfortran).parent
    print(f"Gfortran directory: {gfortran_dir}")

    # Build EQ3/6 using Make
    print(f"Building EQ3/6 from {source_dir}")

    # Set environment variables for compilation
    env = os.environ.copy()

    # Ensure gfortran and make are in PATH
    # Use correct path separator for platform (: for Unix, ; for Windows)
    path_sep = ';' if sys.platform == 'win32' else ':'

    # Add gfortran directory to PATH
    if str(gfortran_dir) not in env.get('PATH', ''):
        env['PATH'] = f"{gfortran_dir}{path_sep}{env.get('PATH', '')}"

    # On Windows, also add MSYS2 usr/bin for make
    if sys.platform == 'win32':
        msys2_usr_bin = r'C:\msys64\usr\bin'
        if msys2_usr_bin not in env.get('PATH', ''):
            env['PATH'] = f"{msys2_usr_bin}{path_sep}{env.get('PATH', '')}"

    # Set FC to explicitly use gfortran
    env['FC'] = gfortran

    # Find make command
    if sys.platform == 'win32':
        # On Windows, use the full path to make.exe in MSYS2
        make_cmd = r'C:\msys64\usr\bin\make.exe'
        if not os.path.isfile(make_cmd):
            raise RuntimeError("make not found. Please install with: pacman -S make")
    else:
        make_cmd = 'make'

    # Platform-specific FFLAGS
    if sys.platform == "darwin":
        # macOS: detect architecture
        archflags = os.environ.get('ARCHFLAGS', '')
        print(f"ARCHFLAGS environment variable: '{archflags}'")

        if 'x86_64' in archflags:
            env['FFLAGS'] = '-O3 -arch x86_64'
            print("Building for macOS x86_64")
        elif 'arm64' in archflags:
            env['FFLAGS'] = '-O3 -arch arm64'
            print("Building for macOS ARM64")
        else:
            # Fallback to machine architecture
            machine = platform.machine()
            arch = 'arm64' if machine == 'arm64' else 'x86_64'
            env['FFLAGS'] = f'-O3 -arch {arch}'
            print(f"Building for macOS {arch} (detected)")
    elif sys.platform == "win32":
        # Windows: static linking to avoid runtime dependencies
        # Add -fcheck=bounds -fbacktrace for debugging array out-of-bounds
        env['FFLAGS'] = '-O0 -g -static -fcheck=bounds -fbacktrace'
        print("Building for Windows with static linking and bounds checking")
    else:
        # Linux
        env['FFLAGS'] = '-O3'
        print("Building for Linux")

    try:
        # Clean any previous builds
        print("Cleaning previous builds...")

        # Remove old binaries from source directory
        source_bin_dir = source_dir / "bin"
        if source_bin_dir.exists():
            shutil.rmtree(source_bin_dir)
            print(f"Removed {source_bin_dir}")

        # Remove old binaries from output directory (preserve __init__.py)
        if output_dir.exists():
            # Save __init__.py if it exists
            init_file = output_dir / "__init__.py"
            init_content = None
            if init_file.exists():
                init_content = init_file.read_text(encoding="utf-8")

            # Remove all files in the directory
            for item in output_dir.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
            print(f"Cleaned {output_dir}")

            # Restore __init__.py
            if init_content is not None:
                init_file.write_text(init_content, encoding="utf-8")
                print(f"Preserved {init_file}")
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

        # Run make clean
        subprocess.run(
            [make_cmd, 'clean'],
            cwd=source_dir,
            env=env,
            check=False,  # Don't fail if clean fails
            capture_output=True
        )

        # Build all executables
        print("Compiling EQ3/6 executables...")
        # Compiler output is decoded explicitly rather than with text=True.
        # text=True uses the locale codec and raises UnicodeDecodeError on any
        # byte it cannot handle, which hides the build log behind a decode
        # error. backslashreplace keeps undecodable bytes visible as escapes
        # and keeps the result printable on a non-UTF-8 console.
        result = subprocess.run(
            [make_cmd, 'all'],
            cwd=source_dir,
            env=env,
            check=True,
            capture_output=True,
            encoding='utf-8',
            errors='backslashreplace'
        )

        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        # Copy executables to output directory
        source_bin_dir = source_dir / "bin"
        if not source_bin_dir.exists():
            raise FileNotFoundError(f"Build failed: {source_bin_dir} not found")

        print(f"Copying executables to {output_dir}")

        # Define expected executables
        if sys.platform == "win32":
            executables = ['eq3nr.exe', 'eq6.exe', 'eqpt.exe', 'xcon3.exe', 'xcon6.exe']
        else:
            executables = ['eq3nr', 'eq6', 'eqpt', 'xcon3', 'xcon6']

        # Copy each executable
        for exe in executables:
            src = source_bin_dir / exe
            dst = output_dir / exe
            if src.exists():
                shutil.copy2(src, dst)
                # Make executable on Unix-like systems
                if sys.platform != "win32":
                    dst.chmod(0o755)
                size_kb = dst.stat().st_size / 1024
                print(f"  [OK] {exe} ({size_kb:.1f} KB)")
            else:
                print(f"  [WARNING] {exe} not found, skipping")

        # Verify at least the main executables exist
        eq3nr = output_dir / ('eq3nr.exe' if sys.platform == "win32" else 'eq3nr')
        eq6 = output_dir / ('eq6.exe' if sys.platform == "win32" else 'eq6')

        if not eq3nr.exists() or not eq6.exists():
            raise RuntimeError("Build failed: eq3nr or eq6 not found")

        print(f"[OK] Successfully compiled EQ3/6 executables")
        return True

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Compilation failed!", file=sys.stderr)
        print(f"Exit code: {e.returncode}", file=sys.stderr)
        if e.stdout:
            print(f"stdout: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"stderr: {e.stderr}", file=sys.stderr)
        raise RuntimeError(f"Failed to compile EQ3/6: {e}")


def main():
    """Main entry point."""
    # the compiler log is relayed to these streams, and a UnicodeEncodeError
    # while printing it would hide the result of the build
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")

    print("=" * 60)
    print("Compiling EQ3/6 executables for aqequil")
    print("=" * 60)
    print(f"Platform: {sys.platform}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Machine: {platform.machine()}")

    try:
        compile_eq3_6()
        print("=" * 60)
        print("Compilation successful!")
        print("=" * 60)
        return 0
    except Exception as e:
        print("=" * 60, file=sys.stderr)
        print(f"ERROR: {e}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
