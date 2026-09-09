"""
Setup script for aqequil.

This setup.py is needed to:
1. Force the creation of platform-specific wheels since we include compiled
   Fortran executables (EQ3/6 binaries).
2. Compile the Fortran executables when building from source distribution.
"""

import os
import sys
import subprocess
from setuptools import setup
from setuptools.dist import Distribution
from setuptools.command.build import build


class BinaryDistribution(Distribution):
    """
    Distribution which always forces a binary package with platform name.

    This ensures that wheels are built as platform-specific (e.g.,
    aqequil-0.41.0-py3-none-linux_x86_64.whl) rather than pure Python
    wheels (py3-none-any.whl).
    """

    def has_ext_modules(self):
        """
        Tell setuptools this distribution has extension modules.
        This forces platform-specific wheel names.
        """
        return True


class BuildWithFortran(build):
    """
    Custom build command that compiles Fortran executables before building.

    This ensures that when installing from source distribution (sdist),
    the EQ3/6 Fortran executables are compiled automatically.
    """

    def run(self):
        # Get the directory containing setup.py
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        compile_script = os.path.join(setup_dir, 'compile_eq3_6.py')
        bin_dir = os.path.join(setup_dir, 'aqequil', 'bin')

        # Check if executables already exist (e.g., when building wheel with cibuildwheel)
        if sys.platform == 'win32':
            eq3nr_path = os.path.join(bin_dir, 'eq3nr.exe')
        else:
            eq3nr_path = os.path.join(bin_dir, 'eq3nr')

        executables_exist = os.path.isfile(eq3nr_path)

        if not executables_exist and os.path.isfile(compile_script):
            print("=" * 60)
            print("Compiling EQ3/6 Fortran executables...")
            print("=" * 60)
            # The compile script relays the compiler log, which can contain
            # bytes the locale codec cannot handle. Pin both sides to UTF-8 so
            # neither the child's print nor the decode here can fail and mask
            # the real build result.
            child_env = os.environ.copy()
            child_env['PYTHONIOENCODING'] = 'utf-8:backslashreplace'

            try:
                result = subprocess.run(
                    [sys.executable, compile_script],
                    cwd=setup_dir,
                    env=child_env,
                    check=True,
                    capture_output=True,
                    encoding='utf-8',
                    errors='backslashreplace'
                )
                if result.stdout:
                    print(result.stdout)
                print("=" * 60)
                print("Fortran compilation successful!")
                print("=" * 60)
            except subprocess.CalledProcessError as e:
                print("=" * 60, file=sys.stderr)
                print("WARNING: Failed to compile Fortran executables!", file=sys.stderr)
                print(f"Error: {e}", file=sys.stderr)
                if e.stdout:
                    print(f"stdout: {e.stdout}", file=sys.stderr)
                if e.stderr:
                    print(f"stderr: {e.stderr}", file=sys.stderr)
                print("=" * 60, file=sys.stderr)
                print("The package will still install, but you'll need to:", file=sys.stderr)
                print("  1. Install gfortran (apt-get install gfortran / brew install gcc)", file=sys.stderr)
                print("  2. Run: python -m aqequil.compile_eq3_6", file=sys.stderr)
                print("  Or set the EQ36CO environment variable to point to EQ3/6 executables.", file=sys.stderr)
            except FileNotFoundError:
                print("=" * 60, file=sys.stderr)
                print("WARNING: gfortran not found - skipping Fortran compilation.", file=sys.stderr)
                print("The package will still install, but you'll need to:", file=sys.stderr)
                print("  1. Install gfortran (apt-get install gfortran / brew install gcc)", file=sys.stderr)
                print("  2. Run: python -m aqequil.compile_eq3_6", file=sys.stderr)
                print("  Or set the EQ36CO environment variable to point to EQ3/6 executables.", file=sys.stderr)
                print("=" * 60, file=sys.stderr)
        elif executables_exist:
            print("EQ3/6 executables already exist, skipping compilation.")

        # Continue with the standard build process
        build.run(self)


if __name__ == "__main__":
    setup(
        distclass=BinaryDistribution,
        cmdclass={
            'build': BuildWithFortran,
        },
    )
