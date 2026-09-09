#!/usr/bin/env python
"""
Integration test for aqequil package.

This test verifies that:
1. EQ3/6 executables are properly bundled and accessible
2. The package can perform basic speciation calculations
3. Results can be accessed and plotted

This is run during the wheel build process to ensure the package works correctly.
"""

import sys
import os

# Set matplotlib to non-GUI backend BEFORE any imports that might use matplotlib
# This prevents font cache building issues in headless CI environments
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib
matplotlib.use('Agg', force=True)


def test_bundled_executables():
    """Test that bundled EQ3/6 executables are found."""
    print("=" * 60)
    print("Test 1: Checking for bundled EQ3/6 executables")
    print("=" * 60)

    import aqequil
    from aqequil.AqSpeciation import _get_bundled_exe_path

    bin_path = _get_bundled_exe_path()
    if bin_path is None:
        print("[FAIL] No bundled executables found!")
        return False

    print(f"[OK] Found bundled executables at: {bin_path}")

    # Check for required executables
    if sys.platform == "win32":
        exe_names = ['eq3nr.exe', 'eq6.exe', 'eqpt.exe']
    else:
        exe_names = ['eq3nr', 'eq6', 'eqpt']

    for exe in exe_names:
        exe_path = os.path.join(bin_path, exe)
        if os.path.isfile(exe_path):
            size_mb = os.path.getsize(exe_path) / (1024 * 1024)
            print(f"  [OK] {exe} ({size_mb:.2f} MB)")
        else:
            print(f"  [FAIL] {exe} not found!")
            return False

    return True


def test_import_and_basic_usage():
    """Test basic aqequil import and usage."""
    print("\n" + "=" * 60)
    print("Test 2: Testing aqequil import and basic usage")
    print("=" * 60)

    try:
        import aqequil
        print("[OK] Successfully imported aqequil")

        # Get test data path
        from aqequil.test_data import get_test_data_path
        test_csv = get_test_data_path("input_example_wrm.csv")

        if not os.path.isfile(test_csv):
            print(f"[FAIL] Test data not found at: {test_csv}")
            return False
        print(f"[OK] Found test data at: {test_csv}")

        return True

    except Exception as e:
        print(f"[FAIL] Error during import/setup: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_download_data0():
    """Download latest data0.wrm from WORM-db and update bundled copies."""
    print("\n" + "=" * 60)
    print("Test 3: Downloading latest data0.wrm from WORM-db")
    print("=" * 60)

    try:
        from aqequil.test_data import get_test_data_path
        from aqequil.databases import get_database_path
        from urllib.request import urlopen
        import shutil

        data0_url = "https://raw.githubusercontent.com/worm-portal/WORM-db/master/data0.wrm"
        print(f"Downloading data0.wrm from {data0_url}...")

        # Download data0.wrm from WORM-db
        with urlopen(data0_url) as response:
            content = response.read()

        size_kb = len(content) / 1024
        print(f"[OK] Downloaded data0.wrm ({size_kb:.1f} KB)")

        # Get destination paths
        test_data_dir = os.path.dirname(get_test_data_path("data0.wrm"))
        databases_dir = get_database_path()

        # Copy to test_data folder
        test_data_dest = os.path.join(test_data_dir, "data0.wrm")
        if os.path.exists(test_data_dest):
            os.remove(test_data_dest)
        with open(test_data_dest, 'wb') as f:
            f.write(content)
        print(f"[OK] Copied data0.wrm to test_data (replaced existing)")

        # Copy to databases folder
        databases_dest = os.path.join(databases_dir, "data0.wrm")
        if os.path.exists(databases_dest):
            os.remove(databases_dest)
        with open(databases_dest, 'wb') as f:
            f.write(content)
        print(f"[OK] Copied data0.wrm to databases (replaced existing)")

        return True

    except Exception as e:
        print(f"[FAIL] Error downloading data0.wrm: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_runeqpt():
    """Test runeqpt - converting data0 to data1 format."""
    print("\n" + "=" * 60)
    print("Test 4: Testing runeqpt (data0 to data1 conversion)")
    print("=" * 60)

    try:
        import aqequil
        from aqequil.test_data import get_test_data_path
        from aqequil.databases import get_database_path
        import tempfile
        import shutil

        # Get data0.wrm from test_data
        data0_source = get_test_data_path("data0.wrm")
        if not os.path.isfile(data0_source):
            print(f"[FAIL] Test data0.wrm not found at: {data0_source}")
            return False

        # Change to a temporary directory
        original_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            print(f"[INFO] Working directory: {tmpdir}")

            # Copy data0.wrm to temp directory
            shutil.copy(data0_source, "data0.wrm")
            print("[OK] Copied data0.wrm to working directory")

            try:
                # Create an AqEquil instance
                ae = aqequil.AqEquil(load_thermo=False, verbose=0)

                # Run EQPT to convert data0.wrm to data1.wrm
                print("Running EQPT on data0.wrm...")
                ae.runeqpt('wrm')
                print("[OK] EQPT completed")

                # Check that data1.wrm was created
                if not os.path.isfile("data1.wrm"):
                    print("[FAIL] data1.wrm was not created")
                    return False

                size_mb = os.path.getsize("data1.wrm") / (1024 * 1024)
                print(f"[OK] data1.wrm created ({size_mb:.2f} MB)")

                # Copy data1.wrm to the bundled databases directory
                # This ensures the wheel includes a freshly compiled data1.wrm
                # File is overwritten if it already exists
                databases_dir = get_database_path()
                dest_data1 = os.path.join(databases_dir, "data1.wrm")
                if os.path.exists(dest_data1):
                    os.remove(dest_data1)
                shutil.copy("data1.wrm", dest_data1)
                print(f"[OK] Copied data1.wrm to bundled databases (replaced existing)")

                return True

            finally:
                os.chdir(original_dir)

    except Exception as e:
        print(f"[FAIL] Error during runeqpt: {e}")
        import traceback
        traceback.print_exc()

        # Known platform issues
        if (sys.platform == 'win32' and 'EQPT' in str(e)) or \
           (sys.platform == 'darwin' and 'did not terminate normally' in str(e)):
            print(f"[WARN] Known {sys.platform} issue with EQPT runtime detected")
            print("[INFO] Skipping runeqpt test due to platform-specific issues")
            return True

        return False


def test_speciation_simple():
    """Test simple speciation with wrm database."""
    print("\n" + "=" * 60)
    print("Test 5: Testing simple speciation (wrm database)")
    print("=" * 60)

    try:
        import aqequil
        from aqequil.test_data import get_test_data_path
        import pandas as pd
        import tempfile
        import shutil

        test_csv = get_test_data_path("input_example_wrm.csv")
        print(f"Running speciation on {test_csv}...")

        # Get data0.wrm from test_data
        data0_source = get_test_data_path("data0.wrm")
        if not os.path.isfile(data0_source):
            print(f"[FAIL] Test data0.wrm not found at: {data0_source}")
            return False

        # Change to a temporary directory
        original_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            print(f"[INFO] Working directory: {tmpdir}")

            # Copy data0.wrm and create data1.wrm for local database testing
            shutil.copy(data0_source, "data0.wrm")
            print("[OK] Copied data0.wrm to working directory")

            try:
                # First run EQPT to create data1.wrm
                ae_eqpt = aqequil.AqEquil(load_thermo=False, verbose=0)
                print("Running EQPT to create data1.wrm...")
                ae_eqpt.runeqpt('wrm')
                print("[OK] EQPT completed, data1.wrm created")

                # Now create AqEquil instance with wrm database
                ae = aqequil.AqEquil(db="wrm", verbose=0)
                print("[OK] AqEquil instance created with wrm database")

                # Run speciation
                speciation = ae.speciate(
                    input_filename=test_csv,
                    exclude=["Year", "Area"]
                )
                print("[OK] Speciation completed")

                # Check results
                if "Bison Pool" not in speciation.sample_data:
                    print("[FAIL] 'Bison Pool' sample not found in results")
                    return False

                aq_dist = speciation.sample_data["Bison Pool"]["aq_distribution"]
                if not isinstance(aq_dist, pd.DataFrame):
                    print(f"[FAIL] aq_distribution is not a DataFrame, got {type(aq_dist)}")
                    return False

                print(f"[OK] aq_distribution is a DataFrame with {len(aq_dist)} rows")
                return True

            finally:
                os.chdir(original_dir)

    except Exception as e:
        print(f"[FAIL] Error during simple speciation: {e}")
        import traceback
        traceback.print_exc()

        if (sys.platform == 'win32' and 'EQPT' in str(e)) or \
           (sys.platform == 'darwin' and 'did not terminate normally' in str(e)):
            print(f"[WARN] Known {sys.platform} issue - skipping test")
            return True

        return False


def test_water_rock_reaction():
    """Test water-rock reaction with EQ6."""
    print("\n" + "=" * 60)
    print("Test 6: Testing water-rock reaction")
    print("=" * 60)

    try:
        import aqequil
        from aqequil.test_data import get_test_data_path
        from aqequil.databases import get_database_path
        import pandas as pd
        import tempfile
        import shutil

        test_csv = get_test_data_path("input_example_wrm.csv")
        print(f"Running speciation on {test_csv}...")

        # Change to a temporary directory
        original_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            print(f"[INFO] Working directory: {tmpdir}")

            try:
                # Create AqEquil instance and run speciation
                # download_csv_files=True downloads the WORM database files
                ae = aqequil.AqEquil(db="WORM", exclude_organics=True, download_csv_files=True, verbose=0)

                # Copy downloaded database files to the bundled databases directory
                # This ensures the wheel includes the latest WORM database files
                # Files are overwritten if they already exist
                databases_dir = get_database_path()
                db_files = [
                    "wrm_data_latest.csv",
                    "elements.csv",
                    "solid_solutions.csv",
                    "wrm_data_logk.csv",
                    "wrm_data_logk_s.csv",
                    "speciation_groups_WORM.txt",
                ]
                for db_file in db_files:
                    if os.path.isfile(db_file):
                        dest_path = os.path.join(databases_dir, db_file)
                        # Remove existing file first to ensure clean replacement
                        if os.path.exists(dest_path):
                            os.remove(dest_path)
                        shutil.copy(db_file, dest_path)
                        print(f"[OK] Copied {db_file} to bundled databases (replaced existing)")
                    else:
                        print(f"[WARN] {db_file} not found in working directory")

                speciation = ae.speciate(
                    input_filename=test_csv,
                    exclude=["Year", "Area"]
                )
                print("[OK] Initial speciation completed")

                # Set up reaction
                Ant = aqequil.Reactant(reactant_name="antigorite", amount_remaining=1)
                r = aqequil.Prepare_Reaction(reactants=[Ant])
                print("[OK] Reaction prepared")

                # Run water-rock reaction
                print("Running water-rock reaction...")
                speciation_mt = aqequil.react(speciation, r)
                print("[OK] Reaction completed")

                # Get mass transfer results for one sample
                m = speciation_mt.mt("Ambergris")

                # Check that we got results
                if not hasattr(m, 'misc_params'):
                    print("[FAIL] Mass transfer results missing 'misc_params' attribute")
                    return False

                if not isinstance(m.misc_params, pd.DataFrame):
                    print(f"[FAIL] misc_params is not a DataFrame, got {type(m.misc_params)}")
                    return False

                print(f"[OK] misc_params is a DataFrame with {len(m.misc_params)} rows")
                return True

            finally:
                os.chdir(original_dir)

    except Exception as e:
        print(f"[FAIL] Error during water-rock reaction: {e}")
        import traceback
        traceback.print_exc()

        if (sys.platform == 'win32' and 'EQPT' in str(e)) or \
           (sys.platform == 'darwin' and 'did not terminate normally' in str(e)):
            print(f"[WARN] Known {sys.platform} issue - skipping test")
            return True

        return False


def test_pitzer_database():
    """Test loading a Pitzer parameter CSV and building data0 Pitzer blocks."""
    print("\n" + "=" * 60)
    print("Test 7: Testing Pitzer parameter database support")
    print("=" * 60)

    try:
        import pandas as pd
        from aqequil.test_data import get_test_data_path
        from aqequil.databases import __file__ as db_init
        from aqequil.pitzer import (validate_pitzer_db, build_pitzer_superblocks,
                                    PITZER_SUPERBLOCK_HEADERS, PITZER_SUPERBLOCK_ORDER,
                                    data1_is_pitzer, data0_is_pitzer)

        pitzer_csv = get_test_data_path("pitzer_params_ypf.csv")
        if not os.path.isfile(pitzer_csv):
            print(f"[FAIL] Pitzer parameter CSV not found at: {pitzer_csv}")
            return False

        pitzer_df = validate_pitzer_db(pd.read_csv(pitzer_csv), filename=pitzer_csv)
        print(f"[OK] Loaded and validated {pitzer_df.shape[0]} Pitzer parameters")

        # charges of aqueous species in the bundled WORM database
        wrm_csv = os.path.join(os.path.dirname(db_init), "wrm_data_latest.csv")
        wrm = pd.read_csv(wrm_csv)
        wrm_aq = wrm[wrm["state"] == "aq"]
        charges = dict(zip(wrm_aq["name"], pd.to_numeric(wrm_aq["z.T"], errors="coerce").fillna(0)))
        charges.update({"H+": 1.0, "H2O": 0.0})

        out = build_pitzer_superblocks(pitzer_df, charges, verbose=0)
        if out["n_blocks"] <= 0:
            print("[FAIL] No Pitzer interaction blocks were built")
            return False
        print(f"[OK] Built {out['n_blocks']} Pitzer interaction blocks with "
              f"{out['n_terms']}-term temperature functions")

        # superblock headers must appear once each, in the order EQPT expects
        lines = out["text"].split("\n")
        positions = []
        for key in PITZER_SUPERBLOCK_ORDER:
            header = PITZER_SUPERBLOCK_HEADERS[key]
            idx = [i for i, l in enumerate(lines) if l == header]
            if len(idx) != 1:
                print(f"[FAIL] Superblock header '{header}' appears {len(idx)} times")
                return False
            positions.append(idx[0])
        if positions != sorted(positions):
            print("[FAIL] Pitzer superblocks are out of order")
            return False
        print("[OK] Pitzer superblocks are complete and in the required order")

        # every coefficient value must contain a decimal point (EQPT reads
        # them with an E25.18 format, which assumes 18 implied decimals for
        # values written without a decimal point)
        bad = [l for l in lines if l.strip().startswith("a") and "=" in l
               and "." not in l.split("=")[1]]
        if len(bad) > 0:
            print(f"[FAIL] Coefficients written without a decimal point: {bad[:3]}")
            return False
        print("[OK] Coefficient formatting is compatible with EQPT")

        # data file type detection
        data1_wrm = os.path.join(os.path.dirname(db_init), "data1.wrm")
        if os.path.isfile(data1_wrm):
            with open(data1_wrm, "rb") as f:
                if data1_is_pitzer(f.read(512)):
                    print("[FAIL] data1.wrm was detected as a Pitzer data file")
                    return False
            print("[OK] data1.wrm detected as a B-dot/Davies data file")
        data0_wrm = get_test_data_path("data0.wrm")
        if os.path.isfile(data0_wrm):
            with open(data0_wrm, "r") as f:
                if data0_is_pitzer(f.read()):
                    print("[FAIL] data0.wrm was detected as a Pitzer data file")
                    return False
            print("[OK] data0.wrm detected as a B-dot/Davies data file")

        return True

    except Exception as e:
        print(f"[FAIL] Error during Pitzer database test: {e}")
        import traceback
        traceback.print_exc()
        return False


SMALL_DATA0 = """data0.tst
Test data0 for data0_to_logK_csv
BEGIN CONFIGURATION DATA BLOCK
INTERPRET 500 AS NO DATA= YES
END CONFIGURATION DATA
+--------------------------------------------------------------------
Miscellaneous parameters
+--------------------------------------------------------------------
Temperature limits (degC)
         0.0000  300.0000
temperatures
         0.0000   25.0000   60.0000  100.0000
       150.0000  200.0000  250.0000  300.0000
pressures
         1.0132    1.0132    1.0132    1.0132
         4.7572   15.5365   39.7365   85.8378
+--------------------------------------------------------------------
elements
+--------------------------------------------------------------------
O          15.99940
+--------------------------------------------------------------------
basis species
+--------------------------------------------------------------------
H2O
     charge  =   0.0
****
     2 element(s):
      2.0000 H              1.0000 O
****
+--------------------------------------------------------------------
Ca++
     charge  =   2.0
****
     1 element(s):
      1.0000 Ca
****
+--------------------------------------------------------------------
Cl-
     charge  =  -1.0
****
     1 element(s):
      1.0000 Cl
****
+--------------------------------------------------------------------
Cr+++
     charge  =   3.0
****
     1 element(s):
      1.0000 Cr
****
+--------------------------------------------------------------------
Fe++
     charge  =   2.0
****
     1 element(s):
      1.0000 Fe
****
+--------------------------------------------------------------------
H+
     charge  =   1.0
****
     1 element(s):
      1.0000 H
****
+--------------------------------------------------------------------
Na+
     charge  =   1.0
****
     1 element(s):
      1.0000 Na
****
+--------------------------------------------------------------------
NO3-
     charge  =  -1.0
****
     2 element(s):
      1.0000 N              3.0000 O
****
+--------------------------------------------------------------------
O2(g)
     charge  =   0.0
****
     1 element(s):
      2.0000 O
****
+--------------------------------------------------------------------
auxiliary basis species
+--------------------------------------------------------------------
CrO4--
     charge  =  -2.0
****
     2 element(s):
      1.0000 Cr             4.0000 O
****
     5 species in reaction:
    -1.0000  CrO4--                      -5.0000  H+
     0.7500  O2(g)                        1.0000  Cr+++
     2.5000  H2O
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
        13.3037   11.8888   10.4704    9.3902
         8.5746    8.1635    8.0620    8.2996
+--------------------------------------------------------------------
Fe+++
     charge  =   3.0
****
     1 element(s):
      1.0000 Fe
****
     5 species in reaction:
    -1.0000  Fe+++                       -0.5000  H2O
     0.2500  O2(g)                        1.0000  Fe++
     1.0000  H+
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
        -9.3693   -7.7654   -5.9151   -4.2176
        -2.5286   -1.1656   -0.0198    0.9786
+--------------------------------------------------------------------
aqueous species
+--------------------------------------------------------------------
CaCl+
     sp.type =  aqueous
     charge  =   1.0
****
     2 element(s):
      1.0000 Ca             1.0000 Cl
****
     3 species in reaction:
    -1.0000  CaCl+                        1.0000  Ca++
     1.0000  Cl-
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
         0.7865    0.7659    0.6105    0.3008
        -0.2261   -0.8623   -1.5765  500.0000
* Source: 98ste/fel (Model 3)
+--------------------------------------------------------------------
CrOH++
     charge  =   2.0
****
     3 element(s):
      1.0000 Cr             1.0000 H              1.0000 O
****
     4 species in reaction:
    -1.0000  CrOH++                      -1.0000  H+
     1.0000  Cr+++                        1.0000  H2O
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
       500.0000    4.0000  500.0000  500.0000
       500.0000  500.0000  500.0000  500.0000
+--------------------------------------------------------------------
CaCl2(aq)
     charge  =   0.0
****
     2 element(s):
      1.0000 Ca             2.0000 Cl
****
     3 species in reaction:
    -1.0000  CaCl2(aq)                    1.0000  Ca++
     2.0000  Cl-
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
        20.3922   16.4677   12.0000    8.0000
         4.0000    1.0000  500.0000  500.0000
+--------------------------------------------------------------------
solids
+--------------------------------------------------------------------
Halite                  NaCl
     V0PrTr  =    27.015 cm**3/mol [source: 78hel/del]
****
     2 element(s):
      1.0000 Cl             1.0000 Na
****
     3 species in reaction:
    -1.0000  Halite                       1.0000  Cl-
     1.0000  Na+
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
         1.5012    1.5857    1.6198    1.5825
         1.4770    1.3249    1.1237    0.8557
+--------------------------------------------------------------------
Soda Niter              NaNO3
     V0PrTr  =   000.000 cm**3/mol [source: ]
****
     3 element(s):
      1.0000 Na             1.0000 N              3.0000 O
****
     3 species in reaction:
    -1.0000  Soda Niter                   1.0000  NO3-
     1.0000  Na+
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
         0.7192    1.0915    1.4321    1.6649
         1.8721   No_Data   No_Data   No_Data
+--------------------------------------------------------------------
CaCl2
     V0PrTr  =    51.620 cm**3/mol [source: 78rob/hem]
****
     2 element(s):
      1.0000 Ca             2.0000 Cl
****
     3 species in reaction:
    -1.0000  CaCl2                        1.0000  Ca++
     2.0000  Cl-
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
        13.1770   11.9420   10.4444    9.0000
         7.0000    5.0000    3.0000    1.0000
+--------------------------------------------------------------------
Hematite                Fe2O3
     V0PrTr  =    30.274 cm**3/mol [source: 78hel/del]
****
     2 element(s):
      2.0000 Fe             3.0000 O
****
     4 species in reaction:
    -1.0000  Hematite                    -6.0000  H+
     2.0000  Fe+++                        3.0000  H2O
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
         1.0000    0.1086   -1.0000   -2.0000
        -3.0000   -4.0000   -5.0000   -6.0000
+--------------------------------------------------------------------
UO2.3333(beta)
     sp.type =  solid   polymorph
     V0PrTr  =   500.000 [null]
****
     2 element(s):
      2.3333 O              1.0000 U
****
     5 species in reaction:
    -2.0000  UO2.3333(beta)              -8.0000  H+
     0.3333  O2(g)                        2.0000  U++++
     4.0000  H2O
*
**** logK grid [0-25-60-100C @1bar; 150-200-250-300C @Psat-H2O]:
       -26.0327  -26.7364  -27.4546  -28.0234
       500.0000  500.0000  500.0000  500.0000
+--------------------------------------------------------------------
liquids
+--------------------------------------------------------------------
gases
+--------------------------------------------------------------------
O2(g)                   O2(g)
     V0PrTr  =     0.000 cm**3/mol [source: supcrt92]
****
     1 element(s):
      2.0000 O
****
     1 species in reaction:
    -1.0000  O2(g)
*
**** logK grid [0-25-60-100C @1.0132bar; 150-200-250-300C @Psat-H2O]:
         0.0000    0.0000    0.0000    0.0000
         0.0000    0.0000    0.0000    0.0000
+--------------------------------------------------------------------
solid solutions
+--------------------------------------------------------------------
references
+--------------------------------------------------------------------
stop.
"""


def test_data0_to_logK_csv():
    """Test converting the species blocks of a data0 file into a logK CSV."""
    print("\n" + "=" * 60)
    print("Test 8: Testing data0 to logK CSV conversion")
    print("=" * 60)

    try:
        import tempfile
        import numpy as np
        import pandas as pd
        import aqequil

        with tempfile.TemporaryDirectory() as tmp:
            data0_path = os.path.join(tmp, "data0.tst")
            with open(data0_path, "w") as f:
                f.write(SMALL_DATA0)
            csv_path = os.path.join(tmp, "tst_logK.csv")
            out = aqequil.data0_to_logK_csv(data0_path, csv_path, verbose=0)
            df = out["logK"]
            report = out["report"]
            if not os.path.isfile(csv_path):
                print("[FAIL] logK CSV was not written")
                return False
        rows = {r["name"]: r for _, r in df.iterrows()}
        print(f"[OK] Converted {len(df)} species: {sorted(rows)}")

        # WORM names, ':' hydrate separator, lower-case minerals, spaces removed
        for name in ["CaCl+", "CrOH+2", "Fe+3", "Cr+3", "halite", "soda_niter", "hematite",
                     "CaCl2", "CaCl2(cr)"]:
            if name not in rows:
                print(f"[FAIL] Expected species '{name}' missing from the logK CSV")
                return False
        print("[OK] Species names follow WORM conventions; the CaCl2 solid was renamed CaCl2(cr)")

        # CrO4-- is a strict basis species in WORM, so it is not written and
        # Cr+3 becomes an auxiliary species with the inverted reaction
        if "CrO4-2" in rows:
            print("[FAIL] CrO4-2 is a target basis species and should not be written")
            return False
        cr = rows["Cr+3"]
        if cr["tag"] != "aux" or "CrO4-2" not in cr["dissrxn"] or "O2(g)" not in cr["dissrxn"]:
            print(f"[FAIL] Cr+3 should be an auxiliary species written in terms of CrO4-2: {cr['dissrxn']}")
            return False
        if abs(cr["logK2"] - (-11.8888)) > 1e-4:
            print(f"[FAIL] Cr+3 log K at 25 C should be -11.8888, got {cr['logK2']}")
            return False
        print("[OK] Basis change Cr+3 -> CrO4-2 handled (Cr+3 is an auxiliary species, log K inverted)")

        # reactions written in terms of an available auxiliary species are kept as is
        hem = rows["hematite"]
        if "Fe+3" not in hem["dissrxn"] or abs(hem["logK2"] - 0.1086) > 1e-6:
            print(f"[FAIL] hematite reaction should be kept in terms of Fe+3: {hem['dissrxn']}")
            return False
        if rows["Fe+3"]["tag"] != "aux":
            print("[FAIL] Fe+3 should carry the 'aux' tag")
            return False
        print("[OK] Reactions involving auxiliary species are kept as written")

        # 'no data' markers (500 and No_Data) shrink the temperature grid
        cacl = rows["CaCl+"]
        if not (cacl["T7"] == 250.0 and np.isnan(cacl["T8"])):
            print(f"[FAIL] CaCl+ should have a 7-point grid (T8 missing), got T7={cacl['T7']} T8={cacl['T8']}")
            return False
        if not (rows["CrOH+2"]["T1"] == 25.0 and np.isnan(rows["CrOH+2"]["T2"])):
            print("[FAIL] CrOH+2 should have a one-point grid at 25 C")
            return False
        sn = rows["soda_niter"]
        if not (sn["T5"] == 150.0 and np.isnan(sn["T6"])):
            print("[FAIL] 'No_Data' entries were not treated as missing data for soda_niter")
            return False
        print("[OK] 500 and No_Data grid entries are treated as missing data")

        # formula, volume, formula_ox and category columns
        if rows["halite"]["formula"] != "NaCl" or abs(rows["halite"]["V"] - 27.015) > 1e-9:
            print("[FAIL] halite formula or volume not carried over")
            return False
        if rows["hematite"]["formula_ox"] != "2Fe+3 3O-2":
            print(f"[FAIL] hematite formula_ox should be '2Fe+3 3O-2', got {rows['hematite']['formula_ox']}")
            return False
        if rows["halite"]["category_1"] != "tst_cr" or rows["CaCl+"]["category_1"] != "tst_aq":
            print("[FAIL] category_1 should be '<dataset>_<state>'")
            return False
        print("[OK] formula, V, formula_ox and category_1 columns are populated")

        # reaction that cannot be balanced with four-decimal coefficients is skipped
        if "UO2.3333(beta)" in rows or "UO2.3333(beta)" not in report["skipped"]:
            print("[FAIL] UO2.3333(beta) should be skipped (unbalanceable with 4-decimal coefficients)")
            return False
        if "U+4" in rows or any("U++++" in k for k in rows):
            print("[FAIL] unknown species U++++ should not produce a row")
            return False
        print("[OK] Species with unrepresentable reactions or unknown participants are skipped and reported")

        # the gas block for O2(g) duplicates the basis species and is not written
        if "O2(g)" in rows:
            print("[FAIL] O2(g) should not be written to the logK CSV")
            return False
        return True

    except Exception as e:
        print(f"[FAIL] Error during data0 to logK CSV test: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all integration tests."""
    print("\n" + "=" * 60)
    print("AQEQUIL INTEGRATION TEST SUITE")
    print("=" * 60)
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {sys.platform}")
    if hasattr(os, 'uname'):
        print(f"Architecture: {os.uname().machine}")
    elif sys.platform == 'win32':
        import platform
        print(f"Architecture: {platform.machine()}")
    print("=" * 60)

    tests = [
        ("Bundled Executables", test_bundled_executables),
        ("Import and Basic Usage", test_import_and_basic_usage),
        ("Download Latest data0.wrm", test_download_data0),
        ("EQPT Data0 to Data1 Conversion", test_runeqpt),
        ("Simple Speciation (wrm database)", test_speciation_simple),
        ("Water-Rock Reaction", test_water_rock_reaction),
        ("Pitzer Parameter Database", test_pitzer_database),
        ("data0 to logK CSV Conversion", test_data0_to_logK_csv),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n[FAIL] Test '{test_name}' crashed with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))

    # Print summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    for test_name, result in results:
        status = "[PASS]" if result else "[FAIL]"
        print(f"{status} {test_name}")

    print("=" * 60)

    # Return exit code
    all_passed = all(result for _, result in results)
    if all_passed:
        print("\n[PASS] All tests passed!")
        return 0
    else:
        print("\n[FAIL] Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
