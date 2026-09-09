FIXED_SPECIES = ["H2O", "H+", "O2(g)", "water", "e-", "OH-", "O2", "H2O(g)"] # ["H2O", "water", "e-", "H+"]
WORM_DB_BASE_URL = "https://raw.githubusercontent.com/worm-portal/WORM-db/master"
# CSV databases that accompany data0 and data1 files, keyed by letter code.
# Functions like plot_reaction_paths() need information in these CSVs that is
# not present in a data0 or data1 file.
WORM_COMPANION_CSVS = {"wrm":["wrm_data_latest.csv", "wrm_data.csv"]}
WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT = {
        'name':'str', 'abbrv':'str', 'formula':'str',
        'state':'str', 'ref1':'str', 'ref2':'str',
        'date': 'str', 'E_units':'str',
        'G':'float', 'H':'float', 'S':'float',
        'Cp':'float', 'V':'float', 'a1.a':'float',
        'a2.b':'float', 'a3.c':'float', 'a4.d':'float',
        'c1.e':'float', 'c2.f':'float',
        'omega.lambda':'float', 'z.T':'float',
        'azero':'float', 'neutral_ion_type':'float',
        'dissrxn':'str', 'tag':'str',
        'formula_ox':'str', 'category_1':'str',
    }

import os
import re
import sys
import shutil
import copy
import collections
import dill
import math
import itertools
from fractions import Fraction
import functools

from ipywidgets import IntProgress
import time

from urllib.request import urlopen
from io import StringIO
import warnings
import subprocess
import pandas as pd
import numpy as np
from chemparse import parse_formula

from IPython.display import display, HTML

import periodictable
from collections import Counter

import pychnosz
from ._HKF_cgl import OBIGT2eos, calc_logK
from .check_TP_grid import check_TP_grid
from .preprocess_for_EQ3 import preprocess, write_3i_file
from .output_3o_mine import main_3o_mine
from .create_data0 import create_data0 as create_data0_py
from .fill_data0_header import fill_data0_head
from .pitzer import validate_pitzer_db, pitzer_param_keys, format_pitzer_param_key, data0_is_pitzer, data1_is_pitzer, PitzerDatabaseError


# matplotlib for static plots
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotly.express as px
import plotly.io as pio

from plotly.subplots import make_subplots
import plotly.graph_objects as go

from wormutils import Error_Handler, chemlabel, format_equation, check_balance, format_coeff, get_colors, isnotebook, assign_worm_db_col_dtypes, import_package_file

from .speciation_groups import SpeciationGroups

# Import new Python modules for dissociation reaction processing
from .process_dissrxns import process_dissrxns
from .redox_suppression import parse_pseudoelement_name

def load(filename, messages=True, hide_traceback=True):
    
    """
    Load a speciation file.

    Parameters
    ----------
    filename : str
        Name of the speciation file.

    messages : bool, default True
        Print messages produced by this function?

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this function?
        When True, error messages handled by this class will be short and to
        the point.

    Returns
    ----------
    An object of class `Speciation`.
    """

    err_handler = Error_Handler(clean=hide_traceback)
    
    if len(filename) <= 12:
        print("Attempting to load "+str(filename)+".speciation ...")
        filename = filename+".speciation"
    
    if 'speciation' in filename[-11:]:
        if os.path.exists(filename) and os.path.isfile(filename):
            pass
        else:
            err = "Cannot locate input file {}/{}".format(os.getcwd(), filename)
            err_handler.raise_exception(err)
    else:
        err = ("Input file {}".format(filename) + " "
            "must be in {} format.".format(ext_dict[ext]))
        err_handler.raise_exception(err)
    
    if os.path.getsize(filename) > 0:
        with open(filename, 'rb') as handle:
            speciation = dill.load(handle)
            if messages:
                print("Loaded '{}'".format(filename))
            return speciation
    else:
        msg = "Cannot open " + str(filename) + " because the file is empty."
        err_handler.raise_exception(msg)


def _get_duplicates(array):
    """
    Return a list of duplicate elements present in a list
    https://stackoverflow.com/questions/46554866/efficiently-finding-duplicates-in-a-list
    """
    c = Counter(array)
    return [k for k in c if c[k] > 1]


def _all_equal(iterable):
    # check that all elements of a list are equal
    g = itertools.groupby(iterable)
    return next(g, True) and not next(g, False)


def _is_url(s):
    """
    Is this string a web address?
    """
    if not isinstance(s, str):
        return False
    return s[0:8].lower() == "https://" or s[0:7].lower() == "http://" or s[0:4].lower() == "www."


def _get_bundled_exe_path():
    """
    Get the path to bundled EQ3/6 executables.

    Returns the path to the aqequil/bin directory containing pre-compiled
    EQ3/6 executables. Returns None if the directory doesn't exist or is empty.
    """
    try:
        # Get the directory where this module is located
        module_dir = os.path.dirname(os.path.abspath(__file__))
        bin_dir = os.path.join(module_dir, 'bin')

        # Check if bin directory exists and contains executables
        if os.path.isdir(bin_dir):
            # Check for at least eq3nr executable (with or without .exe extension)
            eq3nr = os.path.join(bin_dir, 'eq3nr')
            eq3nr_exe = os.path.join(bin_dir, 'eq3nr.exe')

            if os.path.isfile(eq3nr) or os.path.isfile(eq3nr_exe):
                return bin_dir

        return None
    except Exception:
        return None


def _get_bundled_data_path():
    """
    Get the path to bundled EQ3/6 data files (e.g., data1.wrm).

    Returns the path to the aqequil/databases directory containing bundled
    thermodynamic database files. Returns None if the directory doesn't exist
    or doesn't contain data1 files.
    """
    try:
        # Get the directory where this module is located
        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(module_dir, 'databases')

        # Check if databases directory exists and contains data1 files
        if os.path.isdir(data_dir):
            # Check for at least data1.wrm
            data1_wrm = os.path.join(data_dir, 'data1.wrm')

            if os.path.isfile(data1_wrm):
                return data_dir

        return None
    except Exception:
        return None

class AqEquil(object):

    """
    Class containing functions to speciate aqueous water chemistry data using
    existing or custom thermodynamic datasets.
    
    Parameters
    ----------
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored. 
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.
    
    db : str, pd.DataFrame, or list, default "WORM"
        Determines which thermodynamic database is used in the speciation
        calculation. There are several options available:

        - "WORM" will load the default WORM thermodynamic database,
        solid solution database, and logK database. These files are retrieved
        from https://github.com/worm-portal/WORM-db to ensure they are
        up-to-date.
        - Three letter file extension for the desired data1 database, e.g.,
        "wrm". This will use a data1 file with this file extension, e.g.,
        "data1.wrm" located in the path stored in the 'EQ36DA' environment
        variable used by EQ3NR.
        - The name of a data0 file located in the current working directory,
        e.g., "data0.wrm". This data0 file will be compiled by EQPT
        automatically during the speciation calculation.
        - The name of a CSV file containing thermodynamic data located in
        the current working directory, e.g., "wrm_data.csv". The CSV file
        will be used to generate a data0 file for each sample (using
        additional arguments from `db_args` if desired).
        - The URL of a data0 file, e.g.,
        "https://raw.githubusercontent.com/worm-portal/WORM-db/master/data0.wrm"
        - The URL of a CSV file containing thermodynamic data, e.g.,
        "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_latest.csv"

        If a data0 or data1 file has a WORM CSV database that accompanies it,
        e.g., "wrm_data_latest.csv" for "data1.wrm", then that CSV is loaded as
        well. The CSV is not used for EQ3/6 calculations, but functions like
        plot_reaction_paths() require the additional information it contains. A
        copy of the CSV in the current working directory is used if there is
        one; otherwise it is retrieved from
        https://github.com/worm-portal/WORM-db

        A list of two or more databases can also be given, in which case the
        databases are concatenated into one database, e.g.,
        `db=["WORM", "my_extra_species.csv"]`. Only WORM-style CSV databases
        can be concatenated, so each entry in the list must be "WORM", the
        filepath of a CSV file, the URL of a CSV file, or a Pandas DataFrame.
        Three letter codes (e.g., "wrm"), data0 files, and data1 files cannot
        be concatenated because they are not CSV databases. If a species is
        defined in more than one of the databases, the entry in the last
        database that defines it is the one that is used, so databases later
        in the list override earlier ones.

    solid_solutions : str, pd.DataFrame, or list
        Filepath of a CSV file containing parameters for solid solutions, e.g.,
        "my_solid_solutions.csv". If `db` is set to "WORM" and `solid_solutions`
        is not defined, then parameters for solid solutions will be retrieved
        from "Solid_solutions.csv" at https://github.com/worm-portal/WORM-db
        A list of CSV filepaths, URLs, and/or Pandas DataFrames can be given to
        concatenate multiple solid solution databases. Solid solutions defined
        in more than one database are taken from the last database in the list
        that defines them.

    logK : str, pd.DataFrame, or list
        Filepath of a CSV file containing equilibrium constants for chemical
        species, e.g., "my_logK_entries.csv". If `db` is set to "WORM" and `logK`
        is not defined, then equilibrium constants will be retrieved from
        "wrm_data_logK.csv" at https://github.com/worm-portal/WORM-db
        A species in the logK CSV that also exists in the main CSV
        thermodynamic database replaces the main database entry (strict basis
        species excepted), so a logK CSV can be used to override log K values
        of minerals and aqueous species. Such a CSV can be generated from a
        Geochemist's Workbench Pitzer dataset with `aqequil.gwb_tdat_to_csv`.
        A list of CSV filepaths, URLs, and/or Pandas DataFrames can be given to
        concatenate multiple logK databases. Species defined in more than one
        logK database are taken from the last database in the list that defines
        them.

    logK_S : str, pd.DataFrame, or list
        Filepath of a CSV file containing equilibrium constants for chemical
        species, e.g., "my_logK_S_entries.csv". If `db` is set to "WORM" and `logK_S`
        is not defined, then equilibrium constants will be retrieved from
        "wrm_data_logK_S.csv" at https://github.com/worm-portal/WORM-db
        A list of CSV filepaths, URLs, and/or Pandas DataFrames can be given to
        concatenate multiple logK_S databases. Species defined in more than one
        logK_S database are taken from the last database in the list that
        defines them.

    pitzer : str, pd.DataFrame, or list, optional
        Filepath or URL of a CSV file containing Pitzer ion-interaction
        parameters, e.g., "my_pitzer_params.csv". The CSV has one row per
        parameter with the columns species1, species2, species3, param
        (beta0, beta1, beta2, cphi, theta, lambda, psi, zeta, or mu), a1
        through a4 (temperature function coefficients), and optionally alpha,
        ref1, ref2, date, and note. See `aqequil.pitzer` for a full
        description of the format. When a Pitzer database is loaded alongside
        a CSV thermodynamic database, data0 files can be generated for the
        Pitzer activity model (see the `activity_model` parameter of
        `speciate`). This parameter has no effect if the thermodynamic
        database is a data0 or data1 file, because those files already
        determine which activity model must be used.
        A list of CSV filepaths, URLs, and/or Pandas DataFrames can be given to
        concatenate multiple Pitzer parameter databases. A parameter defined in
        more than one database is taken from the last database in the list that
        defines it.
        Two practical notes: EQPT checks every possible ion pair and triplet
        for Pitzer parameters, so compiling a Pitzer data0 file built from a
        large CSV database (hundreds of aqueous species) can take a few
        minutes and several gigabytes of memory per data0 file. Reducing the
        database with `exclude_organics` or `exclude_category` speeds this up
        considerably. Also, Pitzer parameters already account for short-range
        ion interactions, so aqueous ion pairs in the thermodynamic database
        (e.g., NaCl, KCl) usually should be suppressed with the `suppress`
        parameter of `speciate` when Pitzer's equations are used.

    logK_extrapolate : str, default "none"
        What method should be used to extrapolate equilibrium constants in the
        logK database (defined by parameter `logK`) as a function of
        temperature? Can be either "none", "flat", "poly", or "linear".
    
    download_csv_files : bool, default False
        Download copies of database CSV files to your current working directory?

    exclude_organics : bool, default False
        Exclude organic molecules from thermodynamic database? Organic species
        are excluded from the main thermodynamic database CSV and the
        equilibrium constant (logK) CSV database. This parameter has no effect
        if the thermodynamic database is a data0 or data1 file.
        Requires that the databases have a 'category_1' column that designates
        organic molecules.
        The purpose of this parameter is to quickly toggle:
        `exclude_category={'category_1':["organic_aq", "organic_cr"]}`

    exclude_organics_except : list of str, optional
        Suppress the formation of all organic molecules except for the ones
        listed and their protonation states and complexes. For example:
        `exclude_organics_except=['acetic-acid', 'formic-acid'],`
        This will suppress all organics except acetic acid, formic acid,
        acetate, formate, and all formate and acetate complexes.
    
    exclude_category : dict
        Exclude species from thermodynamic databases based on column values.
        For instance,
        `exclude_category={'category_1':["organic_aq", "organic_cr"]}`
        will exclude all species that have "organic_aq" or "organic_cr" in
        the column "category_1".
        Species are excluded from the main thermodynamic database CSV and the
        equilibrium constant (logK) CSV database. This parameter has no effect
        if the thermodynamic database is a data0 or data1 file.
        
    suppress_redox : list of str, default []
        Suppress equilibrium between oxidation states of listed elements
        (Cl, H, and O cannot be included).
        
    input_template : str, default "none"
        Can be either "strict", "basis", "all", or "none" (default). If any
        option other than "none" is chosen, a sample input file template CSV
        file customized to this thermodynamic dataset called
        "sample_input_template.csv" will be generated in the current directory.
        This template can be populated with water sample data to be speciated by
        the `speciate` function. The "strict" option is highly recommended for
        most users. This is because strict basis species speciate into auxiliary
        and non-basis species, but not the other way around.
        Columns in the template include 'Sample', 'Temperature', 'logfO2', and
        others, depending on the chosen option. If "strict", columns for strict
        basis species will be included. If "basis", columns for both strict and
        auxiliary basis species will be included. If "all", then columns for all
        aqueous species will be included.
        
    water_model : str, default "SUPCRT92"
        Desired water model. Can be either "SUPCRT92", "IAPWS95", or "DEW".
        These models are described here: http://chnosz.net/manual/water.html
        SUPCRT92 is valid from 0.01 °C upward. IAPWS95 extrapolates into the
        supercooled liquid region, so it allows calculations at sub-zero
        temperatures (down to about -30 °C at 1 bar), e.g., for cold brines
        with the FREZCHEM or ColdChem Pitzer datasets. DEW is experimental
        and not yet fully supported.
        
    exceed_Ttr : bool, default True
        Calculate Gibbs energies of mineral phases and other species
        beyond their transition temperatures?
        
    verbose : int, 0, 1, or 2, default 1
        Level determining how many messages are returned during a
        calculation. 2 for all messages, 1 for errors or warnings only,
        0 for silent.

    load_thermo : bool, default True
        Load thermodynamic database(s) when instantiating this class?

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this class?
        When True, error messages handled by this class will be short and to
        the point.
    
    Attributes
    ----------
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored.
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.
        
    df_input_processed : pd.DataFrame
        Pandas dataframe containing user-supplied sample chemistry data that has
        been processed by `speciate`.
    
    half_cell_reactions : pd.DataFrame
        Pandas dataframe containing half cell reactions that can be combined
        into redox reactions for calculating chemical affinity and energy supply
        values during speciation.
        
    """

    def __init__(self,
                 eq36da=None,
                 eq36co=None,
                 db="WORM",
                 elements=None,
                 solid_solutions=None,
                 logK=None,
                 logK_S=None,
                 pitzer=None,
                 mineral_db="default",
                 logK_extrapolate="none",
                 download_csv_files=False,
                 exclude_organics=False,
                 exclude_organics_except=None,
                 exclude_category=None,
                 suppress_redox=None,
                 input_template="none",
                 water_model="SUPCRT92",
                 exceed_Ttr=True,
                 verbose=1,
                 load_thermo=True,
                 hide_traceback=True):

        if not isinstance(eq36da, str):
            # Try environment variable first
            eq36da = os.environ.get('EQ36DA')
            # If not set, try to use bundled data files
            if eq36da is None:
                bundled_data_path = _get_bundled_data_path()
                if bundled_data_path is not None:
                    eq36da = bundled_data_path
            # Final fallback: current working directory
            if eq36da is None:
                eq36da = os.getcwd()

        if not isinstance(eq36co, str):
            # Try environment variable first
            eq36co = os.environ.get('EQ36CO')
            # If not set, try to use bundled executables
            if eq36co is None:
                bundled_path = _get_bundled_exe_path()
                if bundled_path is not None:
                    eq36co = bundled_path

        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        self.water_model = water_model

        if not isinstance(suppress_redox, list):
            suppress_redox = []

        self.half_cell_reactions = import_package_file(__name__, 'half_cell_reactions.csv')
        self._half_cell_reactions_original_copy = copy.deepcopy(self.half_cell_reactions)
        
        self.verbose = verbose
        self.hide_traceback = hide_traceback
        self.err_handler = Error_Handler(clean=self.hide_traceback)

        self.raw_3_input_dict = {}
        self.raw_3_output_dict = {}
        self.raw_3_pickup_dict_bottom = {}
        self.raw_3_pickup_dict_top = {}
        
        self.batch_T = []
        self.batch_P = []
        
        self.logK_models = {}

        if "'" in os.getcwd():
            self.err_handler.raise_exception("The current working directory "
                    "({})".format(str(os.getcwd()))+" "
                    "has a parent folder with an invalid name that "
                    "contains an apostrophe ' . Please remove the apostrophe "
                    "and then retry.")
        elif "*" in os.getcwd():
            self.err_handler.raise_exception("The current working directory "
                    "({})".format(str(os.getcwd()))+" "
                    "has a parent folder with an invalid name that "
                    "contains an asterisk * . Please remove the asterisk "
                    "and then retry.")
        else:
            pass
        
        if load_thermo:
            
            # attributes to add to AqEquil class
            self.db = db
            self.elements = elements
            self.solid_solutions = solid_solutions
            self.mineral_db = mineral_db

            if exclude_category is None:
                exclude_category = dict()

            if isinstance(exclude_organics_except, str):
                exclude_organics_except = [exclude_organics_except]
            
            if exclude_organics and not isinstance(exclude_organics_except, list):
                if not isinstance(exclude_category.get("category_1"), list):
                    exclude_category["category_1"] = ["organic_aq", "organic_cr"]
                else:
                    if "organic_aq" not in exclude_category["category_1"]:
                        exclude_category["category_1"] = exclude_category["category_1"].append("organic_aq")
                    if "organic_cr" not in exclude_category["category_1"]:
                        exclude_category["category_1"] = exclude_category["category_1"].append("organic_cr")
            elif isinstance(exclude_organics_except, list):
                # exclude all organics except what is specified by the user in exclude_organics_except
                # and species in associated speciation group
                
                sp = Speciation({})
                sp._make_speciation_group_dict()
                reactant_dict_scalar = copy.deepcopy(sp.reactant_dict_scalar)
                sp_dict_unpacked = copy.deepcopy(sp.speciation_group_dict_unpacked)
                sp_dict_groups = copy.deepcopy(sp.sp_dict_groups)
                self.reactant_dict_scalar = reactant_dict_scalar
                self.sp_dict_unpacked = sp_dict_unpacked
                self.sp_dict_groups = sp_dict_groups
                del sp

                exclude_organics_except_orig = copy.copy(exclude_organics_except)
                for org in exclude_organics_except_orig:
                    for group in self.sp_dict_groups:
                        if org in self.sp_dict_groups[group]:
                            exclude_organics_except += self.sp_dict_groups[group]
                exclude_organics_except = list(set(exclude_organics_except))

            if mineral_db != "default":
                exclude_category["category_1"] = exclude_category["category_1"].append("inorganic_cr")
            
            self.exclude_category = exclude_category
            self.exclude_organics_except = exclude_organics_except
            self.logK = logK
            self.logK_S = logK_S
            self.pitzer = pitzer
            self.logK_extrapolate = logK_extrapolate
            self.download_csv_files = download_csv_files
            self.suppress_redox = suppress_redox
            self.exceed_Ttr = exceed_Ttr
            self.input_template = input_template
            
            self.thermo = AqEquil.Thermodata(AqEquil_instance=self) # outer instance passed to inner instance

            self.data1 = self.thermo.data1

        elif download_csv_files:
            # Download database files without loading them into memory
            self._download_database_files(verbose=verbose)


    def _download_database_files(self, verbose=1):
        """
        Download WORM database files to the current working directory.

        This method downloads the database files without loading them into memory,
        useful for setting up a local copy of the databases.
        """
        from urllib.request import urlopen

        files_to_download = {
            "wrm_data_latest.csv": f"{WORM_DB_BASE_URL}/wrm_data_latest.csv",
            "elements.csv": f"{WORM_DB_BASE_URL}/elements.csv",
            "solid_solutions.csv": f"{WORM_DB_BASE_URL}/solid_solutions.csv",
            "wrm_data_logK.csv": f"{WORM_DB_BASE_URL}/wrm_data_logK.csv",
            "wrm_data_logK_S.csv": f"{WORM_DB_BASE_URL}/wrm_data_logK_S.csv",
            "speciation_groups_WORM.txt": f"{WORM_DB_BASE_URL}/speciation_groups_WORM.txt",
            "data0.wrm": f"{WORM_DB_BASE_URL}/data0.wrm",
        }

        if verbose > 0:
            print("Downloading WORM database files to current working directory...")

        for filename, url in files_to_download.items():
            try:
                if verbose > 0:
                    print(f"Downloading {filename}...")
                with urlopen(url) as response:
                    content = response.read()
                with open(filename, 'wb') as f:
                    f.write(content)
                if verbose > 0:
                    size_kb = len(content) / 1024
                    print(f"  [OK] {filename} ({size_kb:.1f} KB)")
            except Exception as e:
                if verbose > 0:
                    print(f"  [ERROR] Failed to download {filename}: {e}")

        if verbose > 0:
            print("Database file download complete.")


    def _report_3o_6o_errors(self, lines, samplename):

        recording = False
        start_index = 0
        end_index = -1
        error_list = []
        normal_exit = False
        for i,line in enumerate(lines):
            if "* Error" in line:
                recording = True
                start_index = i
                
            if (recording and line == "") or (recording and i == len(lines)-1):
                end_index = i
                recording = False
                error_lines = lines[start_index:end_index+1]
                error_lines = "\n".join(error_lines)
                error_list.append(error_lines)

            if "Normal exit" in line:
                normal_exit = True
                
        if not normal_exit and len(error_list)==0:
            error_list.append("\n * Error - (EQ3/6) The calculation did not terminate normally.")
                
        error_lines = "\n".join(error_list)
        
        if len(error_lines) > 0:
            if self.verbose > 0:
                print("\nThe sample '"+samplename+"' experienced errors during the reaction:")
                print(error_lines+"\n")
            return True
        else:
            return False
            
    
        
    def __file_exists(self, filename, ext='.csv'):
        """
        Check that a file exists and that it has the correct extension.
        Returns True if so, raises exception if not.
        """
        
        ext_dict = {
            ".csv" : "comma separated values (.csv)",
            ".txt" : "standard text (.txt)",
            ".rds" : "R Data (.rds)",
        }

        if ext in filename[-4:]:
            
            if os.path.exists(filename) and os.path.isfile(filename):
                return True
            else:
                err = "Cannot locate input file {}/{}".format(os.getcwd(), filename)
                self.err_handler.raise_exception(err)
        else:
            err = ("Input file {}".format(filename) + " "
                "must be in {} format.".format(ext_dict[ext]))
            self.err_handler.raise_exception(err)
        
        return False


    @staticmethod
    def _ensure_newline_ending(filename):
        """Ensure CSV file ends with newline to prevent warnings in R"""
        if os.path.getsize(filename) > 0:  # Check file isn't empty
            with open(filename, 'rb') as f:
                f.seek(-1, 2)  # Go to last byte
                if f.read(1) != b'\n':
                    with open(filename, 'a') as f:
                        f.write('\n')

    
    def _check_sample_input_file(self, input_filename, exclude, db,
                                       dynamic_db, charge_balance_on,
                                       suppress_missing,
                                       redox_suppression, redox_flag):
        """
        Check for problems in sample input file.
        """
        
        # does the input file exist? Is it a CSV?
        if self.__file_exists(input_filename):
            
            self._ensure_newline_ending(input_filename)
            
            df_in = pd.read_csv(input_filename, header=None) # no headers for now so colname dupes can be checked
        else:
            self.err_handler.raise_exception("_check_sample_input() error!")
        
        # are there any samples?
        if df_in.shape[0] <= 2:
            err_no_samples = ("The file {}".format(input_filename) + " "
                "must contain at least three rows: the "
                "first for column names, the second for column subheaders, "
                "followed by one or more rows for sample data.")
            self.err_handler.raise_exception(err_no_samples)
        
        err_list = [] # for appending errors found in the sample input file
        
        # get header list
        col_list = list(df_in.iloc[0, 1:])
        
        # are there blank headers?
        if True in [isinstance(x, float) and x != x for x in col_list]:
            # isinstance(x, float) and x != x is a typesafe way to check for nan
            err_blank_header = ("One or more columns in the sample input "
                "file have blank headers. These might be empty columns. "
                "Only the first column may have a blank header. Remove any "
                "empty columns and/or give each header a name.")
            self.err_handler.raise_exception(err_blank_header)
        
        # are there duplicate headers?
        dupe_cols = list(set([x for x in col_list if col_list.count(x) > 1]))
        if len(dupe_cols) > 0:
            err_dupe_cols = ("Duplicate column names are not allowed. "
                "Duplicate column names were found for:\n"
                "{}".format(str(dupe_cols)))
            err_list.append(err_dupe_cols)
        
        df_in.columns = df_in.iloc[0] # set column names
        df_in = df_in.drop(df_in.index[0], axis=0) # drop column name row
        df_in_headercheck = copy.deepcopy(df_in.iloc[:,1:]) # drop first column. Deepcopy slice because drop() doesn't work well with unnamed columns.
        
        # drop excluded headers
        for exc in exclude:
            if exc == df_in.columns[0]: # skip if 'sample' column is excluded
                continue
            try:
                df_in_headercheck = df_in_headercheck.drop(exc, axis=1) # drop excluded columns
            except:
                err_bad_exclude = (
                        "Could not exclude the header '{}'".format(exc)+". "
                        "This header could not be found in "
                        "{}".format(input_filename)+"")
                err_list.append(err_bad_exclude)
        
        # get row list
        row_list = list(df_in.iloc[1:, 0])
        
        # are there blank rows?
        if True in [isinstance(x, float) and x != x for x in row_list]:
            # isinstance(x, float) and x != x is a typesafe way to check for nan
            err_blank_row = ("One or more rows in the sample input "
                "file have blank sample names. These might be empty rows. "
                "Remove any empty rows and/or give each sample a name. Sample "
                "names go in the first column.")
            self.err_handler.raise_exception(err_blank_row)
            
        # are there duplicate rows?
        dupe_rows = list(set([x for x in row_list if row_list.count(x) > 1]))
        if len(dupe_rows) > 0:
            err_dupe_rows = ("Duplicate sample names are not allowed. "
                "Duplicate sample names were found for:\n"
                "{}".format(str(dupe_rows)))
            err_list.append(err_dupe_rows)
        
        # does the first column contain non-string values?
        non_string_names = [n for n in list(df_in.iloc[1:, 0])
                            if not isinstance(n, str)]
        if len(non_string_names) > 0:
            self.err_handler.raise_exception(
                "The first column of the input file should contain sample "
                "names (text), but non-text values were found: "
                "{}. This can happen when a CSV file is saved with a "
                "numeric index column. Remove the extra index column and "
                "try again.".format(non_string_names))

        # are there any leading or trailing spaces in sample names?
        invalid_sample_names = [n for n in list(df_in.iloc[1:, 0])
                                if str(n[0])==" " or str(n[-1])==" "]
        if len(invalid_sample_names) > 0:
            err_sample_leading_trailing_spaces = ("The following sample names "
                "have leading or trailing spaces. Remove spaces and try again: "
                "{}".format(invalid_sample_names))
            self.err_handler.raise_exception(err_sample_leading_trailing_spaces)
        
        # are column names valid entries in the database?
        if self.thermo.custom_data0:
            if "data0" in db:
                data_path = db
            else:
                data_path = "data0." + db
        elif self.thermo.dynamic_db:
            data_path = self.thermo.thermo_db_filename
        else:
            data_path = os.path.join(self.eq36da, f"data0.{db}")
        
        if self.thermo.thermo_db_type == "data0" and self.thermo.thermo_db_source == "URL":
            data_path = "data0." + self.thermo.data0_lettercode
        
        # (data_path names every database that was concatenated when db is a
        # list, so there is no single file to look for in that case)
        if not (os.path.exists(data_path) or os.path.isfile(data_path)) and self.thermo.thermo_db_source=="file" and not isinstance(db, list):
            warn_no_data0 = ("Warning: Could not locate {}.".format(data_path) + " "
                "Unable to determine if column headers included in "
                "{} ".format(input_filename) + "match entries for species "
                "in the requested thermodynamic database '{}'.".format(db))
            if self.verbose > 0:
                print(warn_no_data0)
            
        if self.thermo.thermo_db_type == "data0":
            data0_lines = self.thermo.thermo_db.split("\n")

            # Species listed in the bdot parameters section (B-dot/Davies
            # data0 files only; Pitzer data0 files do not have this section).
            db_species = []
            start_index = None
            end_index = None
            for i,s in enumerate(data0_lines):
                if start_index is not None and "+---" in s:
                    end_index = i
                    break
                if '*  species name' in s:
                    start_index = i+1
            if start_index is not None and end_index is not None:
                db_species = [l.split()[0] for l in data0_lines[start_index:end_index] if len(l.split()) > 0]

            # Also extract species from basis and gas sections of data0,
            # since gas basis species (e.g., S2(g)) appear there but not
            # in the bdot parameters section.
            section_headers = {"basis species", "auxiliary basis species",
                               "aqueous species", "gases"}
            current_section = None
            for i, s in enumerate(data0_lines):
                stripped = s.strip()
                if stripped.startswith("+---"):
                    if i + 1 < len(data0_lines):
                        next_line = data0_lines[i + 1].strip()
                        if next_line in section_headers:
                            current_section = next_line
                        elif next_line.startswith("+---") or next_line in {
                            "elements", "solids", "liquids",
                            "solid solutions", "references"}:
                            current_section = None
                        elif current_section is not None:
                            species_name = next_line.split()[0]
                            if species_name not in db_species:
                                db_species.append(species_name)
        elif self.thermo.thermo_db_type == "CSV":
            df_OBIGT = self.thermo.thermo_db
            db_species = list(df_OBIGT["name"])

        if charge_balance_on == 'pH':
            err_charge_balance_on_pH = ("To balance charge on pH, use "
                "charge_balance_on='H+'")
            err_list.append(err_charge_balance_on_pH)
        elif charge_balance_on in ['Temperature', 'logfO2']:
            err_charge_balance_invalid_type = ("Cannot balance charge "
                "on {}.".format(charge_balance_on))
            err_list.append(err_charge_balance_invalid_type)
        elif charge_balance_on == 'column':
            # Check if 'charge_balance_on' column exists - if not, we'll handle it later
            # by defaulting to 'none' (no charge balancing) with a warning message
            pass
        elif charge_balance_on != "none" and charge_balance_on not in list(set(df_in_headercheck.columns)):
            err_charge_balance_invalid_sp = ("The species chosen for charge balance"
                " '{}'".format(charge_balance_on)+""
                " was not found among the headers of the sample input file.")
            err_list.append(err_charge_balance_invalid_sp)
            
        if self.thermo.thermo_db_type in ["data0", "CSV"]:
            
            for species in list(dict.fromkeys(df_in_headercheck.columns)):
                if species not in db_species and species not in ['Temperature', 'logfO2', 'pH', 'Pressure', 'Eh', 'pe']+FIXED_SPECIES:
                    err_species_not_in_db = ("The column '{}'".format(species) + " in the input file "
                        "was not found in the currently loaded chemical species from {}".format(data_path) + ". "
                        "If the column '"+format(species)+"' contains data that should not be "
                        "included in the speciation calculation (such as sample metadata), add the "
                        "column name to the 'exclude' parameter. Otherwise, the "
                        "species might not be loaded because it is not in the thermodynamic database, "
                        "or it was excluded via the 'exclude_category', 'exclude_organics', or 'exclude_organics_except' "
                        "parameters when the database was loaded. Species may also be excluded when "
                        "element redox states are isolated with the parameter 'suppress_redox' and "
                        "appropriate basis species cannot be found to represent them.")
                    err_list.append(err_species_not_in_db)
                elif species == 'pH':
                    err_species_pH = ("Please rename the 'pH' column in "
                        "the sample input file to 'H+' with the subheader "
                        "unit 'pH'.")
                    err_list.append(err_species_pH)

        # are subheader units valid?
        subheaders = df_in_headercheck.iloc[0,]
        valid_subheaders = ["degC", "ppm", "ppb", "Suppressed", "Molality",
                            "Molarity", "mg/L", "mg/kg.sol", "Alk., eq/kg.H2O",
                            "Alk., eq/L", "Alk., eq/kg.sol", "Alk., mg/L CaCO3",
                            "Alk., mg/L HCO3-", "Log activity", "Log act combo",
                            "Log mean act", "pX", "pH", "pe", "pHCl", "pmH", "pmX",
                            "Hetero. equil.", "Homo. equil.", "Make non-basis",
                            "logfO2", "Mineral", "bar", "volts"]
        for i, subheader in enumerate(subheaders):
            if subheader not in valid_subheaders:
                err_valid_sub = ("The subheader '{}'".format(subheader) + " "
                    "for the column '{}'".format(df_in_headercheck.columns[i]) + " "
                    "is not recognized. Valid subheaders are {}".format(str(valid_subheaders)) + ". "
                    "If the column {}".format(df_in_headercheck.columns[i]) + " "
                    "contains data that is not meant for the "
                    "speciation calculation, add the column name "
                    "to the 'exclude' argument. Try help(AqEquil.AqEquil.speciate) "
                    "for more information about 'exclude'.")
                err_list.append(err_valid_sub)
            
        # is a 'Temperature' column present?
        if "Temperature" not in df_in_headercheck.columns and "Temperature" not in exclude:
            err_temp = ("The column 'Temperature' was not found in the input file. "
                "Please include a column with 'Temperature' in the first row, "
                "'degC' in the second row, and a temperature value for each "
                "sample in degrees Celsius.")
            err_list.append(err_temp)


        # handle automatic selection of redox flag
        if redox_flag == "auto":
            redox_selection_made = False
            if "logfO2" in df_in_headercheck.columns:
                if df_in_headercheck["logfO2"][1] == "logfO2":
                    redox_flag = "logfO2"
                    redox_selection_made = True
            if "O2" in df_in_headercheck.columns:
                redox_flag = "O2"
                redox_selection_made = True
            if "Eh" in df_in_headercheck.columns:
                if df_in_headercheck["Eh"][1] == "volts":
                    redox_flag = "Eh"
                    redox_selection_made = True
            if "pe" in df_in_headercheck.columns:
                if df_in_headercheck["pe"][1] == "pe":
                    redox_flag = "pe"
                    redox_selection_made = True
            if "O2(g)" in df_in_headercheck.columns:
                if df_in_headercheck["O2(g)"][1] == "Hetero. equil.":
                    redox_flag = "O2(g)"
                    redox_selection_made = True
            if redox_selection_made:
                if redox_flag == "O2":
                    redox_flag = "logfO2"
                if self.verbose > 0:
                    print("The input file column '"+redox_flag+"' will be used to "
                          "set sample redox state. If a another column is desired, "
                          "set it manually using the redox_flag parameter.")
            else:
                if self.verbose > 0:
                    print("A column to set sample redox state could not be "
                          "automatically detected in the input file. Valid columns "
                          "include 'logfO2' (with subheader 'logfO2'), 'O2' "
                          "(with a valid concentration subheader), 'Eh' "
                          "(with subheader 'volts'), 'pe' (with subheader 'pe'), and "
                          "'O2(g)' (with subheader 'Hetero. equil.'). A "
                          "default logfO2 value will be used to set sample redox "
                          "state. It is also possible to set sample redox state "
                          "based on a valid auxiliary basis species/basis species "
                          "combo (e.g., Fe+3/Fe+2) by setting the redox_flag "
                          "parameter to 'redox aux' and the redox_aux parameter to "
                          "the name of the desired auxiliary basis species.")
                redox_flag = "logfO2"
        
        # is the 'O2(g)' redox flag set up correctly?
        if redox_flag == "O2(g)" or redox_flag == -3:
            if "O2(g)" not in df_in_headercheck.columns:
                err_O2g = ("The redox flag 'O2(g)' was selected but there is no "
                           "O2(g) column in the input file.")
                err_list.append(err_O2g)
            else:
                if df_in_headercheck["O2(g)"][1] != "Hetero. equil.":
                    err_O2g = ("The redox flag 'O2(g)' was selected but the "
                               "O2(g) column in the input file is set to not set to "
                               "heterogeneous equilibrium with minerals, as is "
                               "required. You can do this by changing the subheader "
                               "of the O2(g) column from '"+df_in_headercheck["O2(g)"][1]+"' "
                               "to 'Hetero. equil.' and then specifying the name "
                               "of a mineral on each sample row that should be used "
                               "to set redox state (e.g., hematite)."
                              )
                    err_list.append(err_O2g)

        # raise an exception that summarizes all errors found
        if len(err_list) > 0:
            errs = "\n\n*".join(err_list)
            errs = ("The input file {}".format(input_filename)+" encountered"
                " errors:\n\n*" + errs)
            self.err_handler.raise_exception(errs)
        
        # warn about "suppress_redox" in db_args if "Hetero. equil." among subheaders.
        # Redox suppression won't work for an element constrained by heterogeneous equilibrium
        if redox_suppression and "Hetero. equil." in list(subheaders):
            if self.verbose > 0:
                print("Warning: 'suppress_redox' does not currently work with the "
                      "heterogeneous equilibrium option if the mineral or gas "
                      "contains a redox-suppressed element.")
        
        sample_temps = [float(t) for t in list(df_in["Temperature"])[1:]]
        if "Pressure" in df_in.columns:
            sample_press = [float(p) if p.lower() != 'psat' else 'psat' for p in list(df_in["Pressure"])[1:]]
        else:
            sample_press = ['psat']*len(sample_temps)
        
        
        return sample_temps, sample_press, redox_flag
        

    def __move_eqpt_extra_output(self):
        """
        Moves all EQPT output and data0 into the eqpt_files folder
        """
        
        self.__mk_check_del_directory("eqpt_files")
        if os.path.exists("eqpt_log.txt") and os.path.isfile("eqpt_log.txt"):
            shutil.move("eqpt_log.txt", "eqpt_files/eqpt_log.txt")
        if os.path.exists("data1f.txt") and os.path.isfile("data1f.txt"):
            shutil.move("data1f.txt", "eqpt_files/data1f.txt")
        if os.path.exists("slist.txt") and os.path.isfile("slist.txt"):
            shutil.move("slist.txt", "eqpt_files/slist.txt")

            
    def runeqpt(self, db, dynamic_db=False):
        
        """
        Convert a data0 into a data1 file with EQPT.
        
        Parameters
        ----------
        db : str
            Three letter code of database.
        """

        if os.path.exists("data0."+db) and os.path.isfile("data0."+db):
            pass
        else:
            self.err_handler.raise_exception(" ".join(["Error: could not locate custom database",
                            "data0.{} in {}.".format(db, os.getcwd())]))

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            os.remove("data1."+db)

        self.__move_eqpt_extra_output()

        # Build path to data0 file using os.path.join for cross-platform compatibility
        data0_path = os.path.join(os.getcwd(), f"data0.{db}")

        try:
            # Run EQPT with cross-platform execution method
            self.__run_eq_executable('eqpt', [data0_path], cwd=os.getcwd())
        except:
            self.err_handler.raise_exception(
                "Error: EQPT failed to run on {}.".format("data0."+db))

        if os.path.exists("data1") and os.path.isfile("data1"):
            os.rename("data1", "data1."+db)
        if os.path.exists("data0.d1") and os.path.isfile("data0.d1"):
            os.rename("data0.d1", "data1."+db)
        if os.path.exists("data0.po") and os.path.isfile("data0.po"):
            os.rename("data0.po", "eqpt_log.txt")
        if os.path.exists("data0.d1f") and os.path.isfile("data0.d1f"):
            os.rename("data0.d1f", "data1f.txt")
        if os.path.exists("data0.s") and os.path.isfile("data0.s"):
            os.rename("data0.s", "slist.txt")

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            if self.verbose > 0:
                if not dynamic_db:
                    print("Successfully created a data1."+db+" from data0."+db)
        else:
            if dynamic_db:
                msg = ("EQPT has encounted a problem processing the database "
                       "for this sample. Check eqpt_log.txt for details.")
            else:
                msg = ("EQPT could not create data1."+db+" from "
                       "data0."+db+". Check eqpt_log.txt for details.")
            self.err_handler.raise_exception(msg)
        
        self.__move_eqpt_extra_output()

    
    def runeq3(self,
               filename_3i,
               db,
               samplename=None,
               path_3i="",
               path_3o="",
               path_3p="",
               data1_path="",
               dynamic_db_name=None):
        
        """
        Call EQ3 on a .3i input file.
        
        Parameters
        ----------
        filename_3i : str
            Name of 3i input file.
        
        db : str
            Three letter code of database.
        
        path_3i : path str, default current working directory
            Path of .3i input files.
            
        path_3o : path str, default current working directory
            Path of .3o output files.
        
        path_3p : path str, default current working directory
            Path of .3p pickup files.
        
        data1_path : str, default None
            File path of data1 file.
            
        dynamic_db_name : str
            Database name to be printed if dynamic databases are being used.
            This parameter is for internal use.
        """

        # get current working dir
        cwd = os.getcwd()
        cwdd = cwd + "/"
        
        if samplename == None:
            samplename = filename_3i[:-3]
        
        if self.verbose > 0 and dynamic_db_name == None:
            print('Using ' + db + ' to speciate ' + samplename)
        elif self.verbose > 0 and isinstance(dynamic_db_name, str):
            print('Using ' + dynamic_db_name + ' to speciate ' + samplename)

        # Build paths using os.path.join for cross-platform compatibility
        data1_file = os.path.join(data1_path, f"data1.{db}")

        # If data1 file doesn't exist in the specified path, check current directory as fallback
        if not os.path.exists(data1_file) or not os.path.isfile(data1_file):
            cwd_data1_file = os.path.join(cwd, f"data1.{db}")
            if os.path.exists(cwd_data1_file) and os.path.isfile(cwd_data1_file):
                data1_file = cwd_data1_file

        work_dir = os.path.join(cwd, path_3i)

        input_file_path = os.path.join(work_dir, filename_3i)
        if not os.path.exists(input_file_path):
            self.err_handler.raise_exception(
                "Could not find EQ3 input file '{}'.".format(input_file_path))

        # Run EQ3NR with cross-platform execution method
        # Pass just the filename (not full path) since we're setting cwd to work_dir
        self.__run_eq_executable('eq3nr', [data1_file, filename_3i], cwd=work_dir)
        
        filename_3o = filename_3i[:-1] + 'o'
        filename_3p = filename_3i[:-1] + 'p'

        # The new eq36 build truncates names, e.g., MLS.Source.3i creates MLS.3o
        # Correct for this here:
        input_dir = os.path.join(cwd, path_3i)
        output_dir = os.path.join(cwd, path_3o)
        pickup_dir = os.path.join(cwd, path_3p)

        files_3o = [file for file in os.listdir(input_dir) if ".3o" in file]
        files_3p = [file for file in os.listdir(input_dir) if ".3p" in file]

        if len(files_3o) == 0:
            self.err_handler.raise_exception(
                "EQ3 failed to produce .3o output for '{}'.".format(filename_3i))
        elif len(files_3o) == 1:
            file_3o = files_3o[0]
            try:
                # move output
                shutil.move(os.path.join(input_dir, file_3o), os.path.join(output_dir, filename_3o))
            except:
                self.err_handler.raise_exception("Error: could not move", os.path.join(path_3i, file_3o), "to", os.path.join(path_3o, filename_3o))
        else:
            # multiple 3o output files are present in the directory
            # this might happen when using runeq3() by itself in a directory with 3o files
            pass

        if len(files_3p) == 0:
            self.err_handler.raise_exception(
                "EQ3 failed to produce .3p pickup file for '{}'.".format(filename_3i))
        elif len(files_3p) == 1:
            file_3p = files_3p[0]
            try:
                # move output
                shutil.move(os.path.join(input_dir, file_3p), os.path.join(pickup_dir, filename_3p))
            except:
                self.err_handler.raise_exception("Error: could not move", os.path.join(path_3i, file_3p), "to", os.path.join(path_3p, filename_3p))
        else:
            # multiple 3p output files are present in the directory
            # this might happen when using runeq3() by itself in a directory with 3p files
            pass


    def runeq6(self,
               filename_6i,
               db,
               samplename=None,
               path_6i="",
               data1_path=None,
               dynamic_db_name=None):
        
        """
        Call EQ6 on a .6i input file.
        
        Parameters
        ----------
        filename_6i : str
            Name of 6i input file.
        
        db : str
            Three letter code of database.
        
        samplename : str
            The name of the sample, used to announce which sample is being run.
        
        path_6i : path str, default current working directory
            Path of directory containing .6i input files.
            
        data1_path : path str, default current working directory
            Path of directory where the data1 thermodynamic database file is
            stored. The data1 file will be called from this location to
            perform the speciation. The data1 file must be named
            data1.xyz, where xyz matches `db`, the three letter code of your
            chosen database.
            
        dynamic_db_name : str
            Database name to be printed if dynamic databases are being used.
            This parameter is for internal use.
        """

        if data1_path == None:
            data1_path = self.eq36da
        
        # get current working dir
        cwd = os.getcwd()
        cwdd = cwd + "/"
        
        if samplename == None:
            samplename = filename_6i[:-3]
        
        if self.verbose > 0 and dynamic_db_name == None:
            print('Using ' + db + ' to react ' + samplename)
        elif self.verbose > 0 and isinstance(dynamic_db_name, str):
            print('Using ' + dynamic_db_name + ' to react ' + samplename)

        # Build paths using os.path.join for cross-platform compatibility
        data1_file = os.path.join(data1_path, f"data1.{db}")

        # If data1 file doesn't exist in the specified path, check current directory as fallback
        if not os.path.exists(data1_file) or not os.path.isfile(data1_file):
            cwd_data1_file = os.path.join(cwd, f"data1.{db}")
            if os.path.exists(cwd_data1_file) and os.path.isfile(cwd_data1_file):
                data1_file = cwd_data1_file

        work_dir = os.path.join(cwd, path_6i)

        input_file_path = os.path.join(work_dir, filename_6i)
        if not os.path.exists(input_file_path):
            self.err_handler.raise_exception(
                "Could not find EQ6 input file '{}'.".format(input_file_path))

        # Run EQ6 with cross-platform execution method
        # Pass just the filename (not full path) since we're setting cwd to work_dir
        self.__run_eq_executable('eq6', [data1_file, filename_6i], cwd=work_dir)
        
                
    def __mk_check_del_directory(self, path):
        
        """
        Checks for the dir being created. If it is already present, delete it
        before recreating it.
        """
        
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            shutil.rmtree(path)
            os.makedirs(path)

    
    def __run_script_and_wait(self, args):

        """
        Runs shell commands.

        DEPRECATED: This method is kept for backward compatibility but should not be used.
        Use __run_eq_executable() instead for cross-platform support.
        """

        # DEVNULL and STDOUT needed to suppress all warnings
        subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True).wait()


    def __run_eq_executable(self, exe_name, args, cwd=None):
        """
        Run an EQ3/6 executable in a cross-platform way.

        Parameters
        ----------
        exe_name : str
            Name of the executable (without .exe extension)
        args : list
            List of command-line arguments to pass to the executable
        cwd : str, optional
            Working directory to run the command in
        """
        # Check that eq36co is set
        if self.eq36co is None:
            err_msg = (
                "EQ3/6 executables not found. This can happen when:\n"
                "  1. Installing from source distribution without gfortran installed\n"
                "  2. The EQ36CO environment variable is not set\n\n"
                "To fix this:\n"
                "  - Install gfortran:\n"
                "      Ubuntu/Debian: sudo apt-get install gfortran make\n"
                "      macOS: brew install gcc\n"
                "      Windows: Install MinGW-w64 with gfortran\n"
                "  - Then reinstall aqequil: pip install --force-reinstall aqequil\n\n"
                "Or set the EQ36CO environment variable to point to a directory\n"
                "containing the EQ3/6 executables (eq3nr, eq6, eqpt)."
            )
            self.err_handler.raise_exception(err_msg)

        if sys.platform == "win32":
            # On Windows, run the MinGW executables directly
            exe_path = os.path.join(self.eq36co, f"{exe_name}.exe")
            cmd = [exe_path] + args

            # MinGW executables need DLLs from the MinGW bin directory
            # Add it to PATH for this subprocess
            env = os.environ.copy()
            mingw_bin = r'C:\msys64\mingw64\bin'
            env['PATH'] = mingw_bin + os.pathsep + env.get('PATH', '')
        else:
            # On Linux/Mac, run directly
            exe_path = os.path.join(self.eq36co, exe_name)
            cmd = [exe_path] + args
            env = None  # Use parent environment

        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        if result.returncode != 0:
            stderr_text = result.stdout.decode(errors="replace").strip()
            msg = "{} failed with return code {}.".format(exe_name, result.returncode)
            if stderr_text:
                msg += "\nOutput:\n" + stderr_text
            self.err_handler.raise_exception(msg)

            
    def _delete_rxn_folders(self):
        
        """
        Deletes folders storing raw EQ3/6 input and output.
        """
        
        if os.path.exists('rxn_3i') and os.path.isdir('rxn_3i'):
            shutil.rmtree('rxn_3i')
        if os.path.exists('rxn_3o') and os.path.isdir('rxn_3o'):
            shutil.rmtree('rxn_3o')
        if os.path.exists('rxn_3p') and os.path.isdir('rxn_3p'):
            shutil.rmtree('rxn_3p')
        if os.path.exists('rxn_6i') and os.path.isdir('rxn_6i'):
            shutil.rmtree('rxn_6i')
        if os.path.exists('rxn_6o') and os.path.isdir('rxn_6o'):
            shutil.rmtree('rxn_6o')
        if os.path.exists('rxn_6p') and os.path.isdir('rxn_6p'):
            shutil.rmtree('rxn_6p')
        if os.path.exists('eqpt_files') and os.path.isdir('eqpt_files'):
            shutil.rmtree('eqpt_files')
        if os.path.exists('rxn_data0') and os.path.isdir('rxn_data0'):
            shutil.rmtree('rxn_data0')


    @staticmethod
    def __f(x, poly_coeffs):
        # return values from a polynomial fit
        value = 0
        for i in range(0,len(poly_coeffs)):
            value += poly_coeffs[i]*x**i
        return value


    def __plot_TP_grid_polyfit(self, xvals, yvals, poly_coeffs_1, poly_coeffs_2,
                               res=500, width=600, height=300):

#         print("R COEFFS")
#         print(poly_coeffs_1)
#         print(poly_coeffs_2)

        if (len(xvals) % 2) == 0:
            # if xvals has an even length
            n_mid1 = math.floor(len(xvals)/2)-1
            n_mid2 = n_mid1+1
        else:
            # if xvals has an odd length
            n_mid1 = math.floor(len(xvals)/2)
            n_mid2 = n_mid1+1
        
        
        f1_x = np.linspace(xvals[0], xvals[n_mid1], num=res)
        f2_x = np.linspace(xvals[n_mid1], xvals[len(xvals)-1], num=res)
        f1_y = [self.__f(x, poly_coeffs_1) for x in f1_x]
        f2_y = [self.__f(x, poly_coeffs_2) for x in f2_x]

        fig = go.Figure()

        fig.add_trace(go.Scatter(x=f1_x, y=f1_y,
                            mode='lines',
                            name='f1'))
        fig.add_trace(go.Scatter(x=f2_x, y=f2_y,
                            mode='lines',
                            name='f2'))
        fig.add_trace(go.Scatter(x=xvals, y=yvals,
                            mode='markers',
                            name='TP points'))
        
        fig.update_layout(legend_title=None,
                          title={'text':"TP grid polyfit"}, autosize=False,
                          width=width, height=height,
                          margin={"t": 40}, xaxis={'fixedrange':True},
                          yaxis={'fixedrange':True}, template="simple_white")

        fig['layout']['xaxis']['title']='Temperature, °C'
        fig['layout']['yaxis']['title']='Pressure, bar'
            
        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': [],
                  'toImageButtonOptions': {
                                           'format': 'png', # one of png, svg, jpeg, webp
                                           'filename': "TP_grid_fit",
                                           'height': height,
                                           'width': width,
                                           'scale': 1,
                                           },
                 }

        fig.show(config=config)

    
    def _interpolate_logK(self, T, logK_grid, T_grid, logK_extrapolate="none"):

        logK_grid_trunc = [t for t in logK_grid if not math.isnan(t)]
        grid_len = len(logK_grid_trunc)
        logK_grid = logK_grid_trunc
        T_grid = T_grid[0:grid_len]

        if logK_extrapolate=="none" and (T > max(T_grid) or T < min(T_grid)):
            return np.nan, "no fit"
        elif logK_extrapolate=="no fit":
            return np.nan, "no fit"
        
        # # turns off poor polyfit warning
        # # TODO: restore polyfit warning setting afterward
        # warnings.simplefilter('ignore', np.RankWarning)
        
        if len(T_grid) >= 4:
            if (len(T_grid) % 2) == 0:
                # if T_grid has an even length
                n_mid1 = math.floor(len(T_grid)/2)-1
                n_mid2 = n_mid1+1
            else:
                # if T_grid has an odd length
                n_mid1 = math.floor(len(T_grid)/2)
                n_mid2 = n_mid1+1

            poly_coeffs_1 = np.polyfit(T_grid[:n_mid2], logK_grid[:n_mid2], len(T_grid[:n_mid2])-1)
            poly_coeffs_2 = np.polyfit(T_grid[n_mid1:], logK_grid[n_mid1:], len(T_grid[n_mid1:])-1)

            model_1 = np.poly1d(poly_coeffs_1)
            model_2 = np.poly1d(poly_coeffs_2)

            if T >= T_grid[0] and T <= T_grid[n_mid1]:
                logK = model_1(T)
                model = "model 1"
            elif T > T_grid[n_mid1] and T <= T_grid[-1]:
                logK = model_2(T)
                model = "model 2"
            else:
                # dependent on extrapolation option
                if logK_extrapolate=="none":
                    logK = np.nan
                    model = "no fit"
                elif logK_extrapolate=="poly":
                    if T < T_grid[0]:
                        logK = model_1(T)
                        model = "model 1"
                    elif T > T_grid[-1]:
                        logK = model_2(T)
                        model = "model 2"
                    else:
                        logK = np.nan
                        model = "no fit"
                elif logK_extrapolate=="linear":
                    poly_coeffs_1 = np.polyfit(T_grid[0:2], logK_grid[0:2], 1)
                    linear_model_1 = np.poly1d(poly_coeffs_1)
                    poly_coeffs_2 = np.polyfit(T_grid[-2:], logK_grid[-2:], 1)
                    linear_model_2 = np.poly1d(poly_coeffs_2)
                    if T < T_grid[0]:
                        logK = linear_model_1(T)
                        model = "linear model 1"
                    elif T > T_grid[-1]:
                        logK = linear_model_2(T)
                        model = "linear model 2"
                    else:
                        logK = np.nan
                        model = "no fit"
                elif logK_extrapolate=="flat":
                    if T < T_grid[0]:
                        logK = logK_grid[0]
                        model = "flat extrap. 1"
                    elif T > T_grid[-1]:
                        logK = logK_grid[-1]
                        model = "flat extrap. 2"
                    else:
                        logK = np.nan
                        model = "no fit"
        elif len(T_grid) >= 2:
            poly_coeffs = np.polyfit(T_grid, logK_grid, len(T_grid)-1)
            model_fit = np.poly1d(poly_coeffs)
            if T >= T_grid[0] and T <= T_grid[-1]:
                logK = model_fit(T)
                model = "model 1"
            else:
                # dependent on extrapolation option
                if logK_extrapolate=="none":
                    logK = np.nan
                    model = "no fit"
                elif logK_extrapolate in ["poly", "linear"]:
                    logK = model_fit(T)
                    model = "model 1"
                elif logK_extrapolate=="flat":
                    if T < T_grid[0]:
                        logK = logK_grid[0]
                        model = "flat extrap."
                    elif T > T_grid[-1]:
                        logK = logK_grid[-1]
                        model = "flat extrap."
                    else:
                        logK = np.nan
                        model = "no fit"
        else:
            # only one T_grid value
            if T == T_grid[0]:
                logK = logK_grid[0]
                model = "single point extrap."
            else:
                # dependent on extrapolation option
                if logK_extrapolate=="none":
                    logK = np.nan
                    model = "no fit"
                elif logK_extrapolate!="none":
                    logK = logK_grid[0]
                    model = "single point extrap."
            
#         ### TEST
#         from matplotlib import pyplot as plt
#         plt.plot(T_grid, logK_grid, 'o')
#         T_m1 = np.linspace(min(T_grid[:n_mid2]), max(T_grid[:n_mid2]), 100)
#         T_m2 = np.linspace(min(T_grid[n_mid1:]), max(T_grid[n_mid1:]), 100)
#         plt.plot(T_m1, model_1(T_m1))
#         plt.plot(T_m2, model_2(T_m2))
#         ###
        
        return logK, model


    def _resolve_activity_model(self, activity_model):
        """
        Choose the aqueous activity coefficient model to use for a
        calculation, checking that it is compatible with the loaded
        thermodynamic database.

        EQ3/6 requires the activity model (the iopg(1) switch of an input
        file) to match the format of the data file: a Pitzer data0/data1 file
        must be used with Pitzer's equations, and a B-dot/Davies (SEDH)
        data0/data1 file cannot be used with Pitzer's equations. A CSV
        database can generate either kind of data0 file, but a Pitzer
        data0 file requires a loaded Pitzer parameter database.

        Parameters
        ----------
        activity_model : str
            "auto", "b-dot", "davies", or "pitzer".

        Returns
        -------
        str
            "b-dot", "davies", or "pitzer".
        """

        valid_models = ["auto", "b-dot", "davies", "pitzer"]
        if activity_model not in valid_models:
            self.err_handler.raise_exception("Unrecognized activity_model "
                "'{}'. Valid options are {}.".format(activity_model, valid_models))

        db_model = self.thermo.db_activity_model
        db_name = str(self.thermo.thermo_db_filename)

        if db_model == "pitzer":
            if activity_model in ["auto", "pitzer"]:
                return "pitzer"
            self.err_handler.raise_exception("The thermodynamic database "
                "{} is formatted for Pitzer's equations and cannot be used "
                "with the '{}' activity model. Use activity_model='pitzer' "
                "(or 'auto'), or load a B-dot/Davies database.".format(db_name, activity_model))

        elif db_model == "sedh":
            if activity_model == "pitzer":
                self.err_handler.raise_exception("The thermodynamic database "
                    "{} is formatted for the B-dot and Davies equations and "
                    "does not contain Pitzer interaction parameters, so it "
                    "cannot be used with activity_model='pitzer'. To use "
                    "Pitzer's equations, load a Pitzer data0 or data1 file "
                    "(e.g., db='data0.ypf'), or load a CSV thermodynamic "
                    "database together with a Pitzer parameter CSV "
                    "(e.g., AqEquil(db='wrm_data.csv', pitzer='pitzer_params.csv')).".format(db_name))
            return "b-dot" if activity_model == "auto" else activity_model

        else:
            # CSV or DataFrame database: a data0 is generated for each sample
            if activity_model == "auto":
                return "pitzer" if self.thermo.pitzer_active else "b-dot"
            if activity_model == "pitzer" and not self.thermo.pitzer_active:
                self.err_handler.raise_exception("activity_model='pitzer' requires "
                    "a Pitzer parameter database, but none was loaded. Load one "
                    "with the 'pitzer' parameter when creating the AqEquil object, "
                    "e.g., AqEquil(db='wrm_data.csv', pitzer='pitzer_params.csv').")
            return activity_model


    def speciate(self,
                 input_filename,
                 logK_extrapolate=None,
                 activity_model="auto",
                 redox_flag="auto",
                 redox_aux="Fe+3",
                 default_logfO2=-6,
                 exclude=[],
                 suppress=[],
                 alter_options=[],
                 charge_balance_on="column",
                 suppress_missing=True,
                 blanks_are_0=False,
                 strict_minimum_pressure=True,
                 aq_scale=1,
                 verbose=1,
                 report_filename=None,
                 get_aq_dist=True,
                 aq_dist_type="log_activity",
                 get_mass_contribution=True,
                 mass_contribution_other=True,
                 get_mineral_sat=True,
                 mineral_sat_type="affinity",
                 get_redox=True,
                 redox_type="Eh",
                 get_ion_activity_ratios=True,
                 get_fugacity=True,
                 get_basis_totals=True,
                 get_solid_solutions=True,
                 mineral_reactant_energy_supplies=False,
                 get_charge_balance=True,
                 batch_3o_filename=None,
                 delete_generated_folders=False,
                 db_args={}):
        
        """
        Calculate the equilibrium distribution of chemical species in solution.
        Additionally, calculate chemical affinities and energy supplies for
        user-specified reactions.
        
        Parameters
        ----------
        input_filename : str
            User-supplied utf8-encoded comma separated value (csv) file
            containing sample data intended for speciation. The file must
            follow this format:
            
            - the first row is a header row that must contain the names of the
              species to be included in the speciation calculation. There
              cannot be duplicate headers.
            - the second row must contain subheaders for each species in the
              header row. These subheaders must be taken from the following:
              
                    degC
                    ppm
                    ppb
                    Suppressed
                    Molality
                    Molarity
                    mg/L
                    mg/kg.sol
                    Alk., eq/kg.H2O
                    Alk., eq/L
                    Alk., eq/kg.sol
                    Alk., mg/L CaCO3
                    Alk., mg/L HCO3-
                    Log activity
                    Log act combo
                    Log mean act
                    pX
                    pH
                    pHCl
                    pmH
                    pmX
                    Hetero. equil.
                    Homo. equil.
                    Make non-basis
                    
            - 'Temperature' must be included as a header, with 'degC' as its
              subheader.
            - The first column must contain sample names. There cannot be
              duplicate sample names.
        
        activity_model : str, default "auto"
            Aqueous species activity coefficient model to use for speciation.
            Can be "auto", "b-dot", "davies", or "pitzer".

            EQ3/6 requires that the activity model match the format of the
            thermodynamic data file. A data0 or data1 file is formatted either
            for the B-dot and Davies equations (e.g., data0.wrm) or for
            Pitzer's equations (e.g., data0.ypf), and cannot be used with the
            other kind of model. A CSV thermodynamic database can produce
            either kind of data0 file; the "pitzer" model requires that a
            Pitzer parameter database was loaded with the `pitzer` parameter
            of `AqEquil`.

            "auto" selects "pitzer" when the loaded database is a Pitzer
            data0 or data1 file or when a Pitzer parameter database was loaded
            alongside a CSV database, and "b-dot" otherwise.
        
        redox_flag : str, default "auto"
            Determines which column in the sample input file sets the overall
            redox state of the samples. Options for redox_flag include 'O2(g)',
            'pe', 'Eh', 'logfO2', and 'redox aux'. The code will search your
            sample spreadsheet file (see `filename`) for a column corresponding
            to the option you chose:

            * 'auto' to automatically decide which redox variable to use based
              on the user's input file.
            * 'O2(g)' set to heterogeneous equilibrium with a mineral
            * 'pe' with subheader pe
            * 'Eh' with subheader volts
            * 'logfO2' with subheader logfO2
            * 'redox aux' will search for a column corresponding to the
              auxilliary basis species selected to form a redox couple with its
              linked strict basis species (see `redox_aux`). For example, the
              redox couple Fe+2/Fe+3 would require a column named Fe+3
            
            If an appropriate header or redox data cannot be found to define
            redox state, `default_logfO2` is used to set sample logfO2.
            
            There is a special case where dissolved oxygen can be used to impose
            sample redox state if `redox_flag` is set to logfO2 and a column named
            logfO2 does not appear in your sample spreadsheet. If there is a
            column corresponding to dissolved oxygen measurements, logfO2 is
            calculated from the equilibrium reaction O2(aq) = O2(g) at the
            temperature and pressure of the sample using the revised Helgeson-
            Kirkham-Flowers (HKF) equation of state (JC Tanger IV and HC
            Helgeson, Am. J. Sci., 1988, 288, 19).
        
        redox_aux : default "Fe+3", optional
            Ignored unless `redox_flag` equals 1. Name of the auxilliary species
            whose reaction links it to a basis species (or another auxilliary
            species) such that they form a redox couple that controls sample
            fO2. For instance, Fe+3 is linked to Fe+2 in many supporting data
            files, so selecting `redox_flag` = 1 and `redox_aux` = "Fe+3" will
            set sample fO2 based on the Fe+2/Fe+3 redox couple.
        
        default_logfO2 : float, default -6
            Default value for sample logfO2 in case redox data cannot be found
            in the user-supplied sample spreadsheet.
        
        exclude : list of str, default []
            Names of columns in the user-supplied sample spreadsheet that should
            not be considered aqueous species. Useful for excluding columns
            containing sample metatadata, such as "Year" and "Location".
            
        suppress : list of str, default []
            Names of chemical species that will be prevented from forming in the
            speciation calculation.
        
        alter_options : list, default []
            A list of lists, e.g.,
            [["CaOH+", "Suppress"], ["CaCl+", "AugmentLogK", -1]]
            The first element of each interior list is the name of a species.
            The second element is an option to alter the species, and can be:
            
            - Suppress : suppress the formation of the species. (See also:
            `suppress`).
            - Replace : replace the species' log K value with a desired value.
            - AugmentLogK : augment the value of the species' log K.
            - AugmentG : augment the Gibbs free energy of the species by a
            desired value, in kcal/mol.
            
            The third element is a numeric value corresponding to the chosen
            option. A third element is not required for Suppress.
            
        charge_balance_on : str, default "column"
            Specifies how to handle electrical charge balancing in the speciation
            calculation. Options include:

            - "column" (default): Read charge balance constraints from a
              'charge_balance_on' column in the input file, allowing different
              constraints for each sample. If no such column exists, defaults
              to "none" (no charge balancing) with a warning message.
            - "none": Do not balance electrical charge between cations and anions.
            - Species name (e.g., "H+"): The activity of the specified species
              will be allowed to change until charge balance is obtained. For
              example, charge_balance_on = "H+" will calculate what pH a sample
              must have to have zero net charge.
        
        suppress_missing : bool, default True
            Suppress the formation of an aqueous species if it is missing a
            value in the user-supplied sample spreadsheet?

        blanks_are_0 : bool, default False
            Assume all blank values in the water chemistry input file are 0?
            
        strict_minimum_pressure : bool, default True
            Ensure that the minimum pressure in the speciation calculation does
            not go below the minimum pressure in the TP grid of the data0 file?
        
        aq_scale : float, default 1
            Scale factor for the mass of the aqueous phase. By default, the
            aqueous phase is 1 kg of solvent.
        
        verbose : int, 0, 1, or 2, default 1
            Level determining how many messages are returned during a
            calculation. 2 for all messages, 1 for errors or warnings only,
            0 for silent.
            
        report_filename : str, optional
            Name of the comma separated values (csv) report file generated when
            the calculation is complete. If this argument is not defined, a
            report file is not generated.
            
        get_aq_dist : bool, default True
            Calculate distributions of aqueous species?
        
        aq_dist_type : str, default "log_activity"
            Desired units of measurement for reported distributions of aqueous
            species. Can be "molality", "log_molality", "log_gamma", or
            "log_activity". Ignored if `get_aq_dist` is False.
        
        get_mass_contribution : bool, default True
            Calculate basis species contributions to mass balance of aqueous
            species?
        
        mass_contribution_other : bool, default True
            Include an "other" species for the sake of summing percents of basis
            species contributions to 100%? Ignored if `get_mass_contribution` is
            False.
        
        get_mineral_sat : bool, default True
            Calculate saturation states of pure solids?
        
        mineral_sat_type : str, default "affinity"
            Desired units of measurement for reported saturation states of pure
            solids. Can be "logQoverK" or "affinity". Ignored if
            `get_mineral_sat` is False.
        
        get_redox : bool, default True
            Calculate potentials of redox couples?
            
        redox_type : str, default "Eh"
            Desired units of measurement for reported redox potentials. Can be
            "Eh", "pe", "logfO2", or "Ah". Ignored if `get_redox` is False.
        
        get_ion_activity_ratios : bool, default True
            Calculate ion/H+ activity ratios and neutral species activities?
        
        get_fugacity : bool, default True
            Calculate gas fugacities?

        get_basis_totals : bool, default True
            Report total compositions of basis aqueous species?

        get_solid_solutions : bool, default True
            Permit the calculation of solid solutions and include them in the
            speciation report?

        mineral_reactant_energy_supplies : bool, default False
            Report energy supplies for reactions with mineral reactants? This
            option is False by default because mineral reactants are considered
            to be unlimited. As a result, energy supplies from reactions with
            reactant minerals tend to be artificially high, especially in
            systems where the reactant minerals are unstable.
        
        rxn_filename : str, optional
            Name of .txt file containing reactions used to calculate affinities
            and energy supplies. Ignored if `get_affinity_energy` is False.

        get_charge_balance : bool, default True
            Calculate charge balance and ionic strength?
        
        batch_3o_filename : str, optional
            Name of rds (R object) file exported after the speciation
            calculation? No file will be generated if this argument is not
            defined.
            
        delete_generated_folders : bool, default False
            Delete the 'rxn_3i', 'rxn_3o', 'rxn_3p', and 'eqpt_files' folders
            containing raw EQ3NR input, output, pickup, and EQPT files once the
            speciation calculation is complete?
           
        db_args : dict, default {}
            Dictionary of arguments to modify how the thermodynamic database is
            processed. Only used when `db` points to thermodynamic data in a CSV
            file. Ignored if `db` points to a data0 file (because a data0 file
            is already ready for a speciation calculation). Options for
            `db_args` are passed to the `create_data0` function, so refer to
            `create_data0` for more information about what options are possible.
            
            - Example of `db_args` where organics are excluded and redox is
            suppressed for Fe and S:
            
            db_args = {
               "exclude_category":{"category_1":["organic_aq"]},
               "suppress_redox":["Fe", "S"],
            }
            
        
        Returns
        -------
        speciation : object of class Speciation
            Contains the results of the speciation calculation.
        
        """

        self.batch_T = []
        self.batch_P = []
        
        self.verbose = verbose

        db = self.thermo.db
            
        if self.thermo.thermo_db_type in ["CSV", "Pandas DataFrame"]:
            db_args["db"] = "dyn"
            
        dynamic_db = self.thermo.dynamic_db
        data0_lettercode = self.thermo.data0_lettercode # needs to be this way
        
        if (self.thermo.thermo_db_type == "data0" or self.thermo.thermo_db_type == "data1") and len(db_args) > 0:
            if self.verbose > 0:
                print("Warning: Ignoring db_args because a premade data0 or data1 file is being used: '" + db + "'")
            
        redox_suppression = False
        if "suppress_redox" in db_args.keys() and self.thermo.thermo_db_type != "data0" and self.thermo.thermo_db_type != "data1":
            if len(db_args["suppress_redox"]) > 0:
                redox_suppression = True

        # Exclude the 'charge_balance_on' column from being treated as a chemical species
        # (if it exists in the input file). This allows users to have this column in their CSV
        # for use with charge_balance_on='column', but still use uniform charge balancing
        # (e.g., charge_balance_on='H+') without errors.
        exclude_modified = list(exclude)  # Create a copy to avoid modifying the original
        # Check if the column exists before trying to exclude it
        df_check = pd.read_csv(input_filename, header=0, nrows=0)
        if 'charge_balance_on' in df_check.columns and 'charge_balance_on' not in exclude_modified:
            exclude_modified.append('charge_balance_on')

        # choose the activity model that is compatible with the loaded database
        activity_model = self._resolve_activity_model(activity_model)

        # check input sample file for errors
        sample_temps, sample_press, redox_flag = self._check_sample_input_file(
                                      input_filename, exclude_modified, db,
                                      dynamic_db, charge_balance_on, suppress_missing,
                                      redox_suppression, redox_flag)

        if aq_dist_type not in ["molality", "log_molality", "log_gamma", "log_activity"]:
            self.err_handler.raise_exception("Unrecognized aq_dist_type. Valid "
                "options are 'molality', 'log_molality', 'log_gamma', 'log_activity'")
        if mineral_sat_type not in ["logQoverK", "affinity"]:
            self.err_handler.raise_exception("Unrecognized mineral_sat_type. Valid "
                "options are 'logQoverK' or 'affinity'")
        if redox_type not in ["Eh", "pe", "logfO2", "Ah"]:
            self.err_handler.raise_exception("Unrecognized redox_type. Valid "
                "options are 'Eh', 'pe', 'logfO2', or 'Ah'")

        if redox_flag == "O2(g)" or redox_flag == "O2" or redox_flag == -3:
            redox_flag = -3
        elif redox_flag == "pe" or redox_flag == -2:
            redox_flag = -2
        elif redox_flag == "Eh" or redox_flag == -1:
            redox_flag = -1
        elif redox_flag == "logfO2" or redox_flag == 0:
            redox_flag = 0
        elif redox_flag == "redox aux" or redox_flag == 1:
            redox_flag = 1
        else:
            self.err_handler.raise_exception("Unrecognized redox flag. Valid "
                    "options are 'auto', 'O2(g)', 'pe', 'Eh', 'logfO2', and "
                    "'redox aux'")
        
        # handle batch_3o naming
        if batch_3o_filename != None:
            if ".rds" in batch_3o_filename[-4:]:
                batch_3o_filename = batch_3o_filename
            else:
                batch_3o_filename = "batch_3o_{}.pkl".format(data0_lettercode)
        # else: keep batch_3o_filename as None

        # reset logK_models whenever speciate() is called
        # (prevents errors when speciations are run back-to-back)
        self.logK_models = {}

        
        
        # dynamic data0 creation per sample
        if dynamic_db:
            db_args["fill_data0"] = False
            db_args["dynamic_db"] = True
            db_args["verbose"] = self.verbose
            db_args["activity_model"] = activity_model
            db_args["dynamic_db_sample_temps"] = sample_temps
            db_args["dynamic_db_sample_press"] = sample_press
            
            if logK_extrapolate != None:
                db_args["logK_extrapolate"] = logK_extrapolate
            elif self.thermo.logK_active:
                db_args["logK_extrapolate"] = self.thermo.logK_extrapolate
                logK_extrapolate = self.thermo.logK_extrapolate
            else:
                logK_extrapolate = "none"
                    
            if self.verbose > 0:
                print("Getting", self.thermo.thermo_db_filename, "ready. This will take a moment...")
                
            thermo_df, data0_file_lines, grid_temps, grid_press, data0_lettercode, water_model, P1, plot_poly_fit = self.create_data0(**db_args)

        if self.thermo.custom_data0 and not dynamic_db:
            self.__mk_check_del_directory('eqpt_files')
            if self.thermo.thermo_db_type != "data1":
                self.runeqpt(data0_lettercode)
            
            if os.path.exists("data1."+data0_lettercode) and os.path.isfile("data1."+data0_lettercode):
                try:
                    # store contents of data1 file in AqEquil object
                    with open("data1."+data0_lettercode, mode='rb') as data1:
                        self.data1["all_samples"] = data1.read()
                    # move or copy data1
                    if self.thermo.thermo_db_type != "data1":
                        shutil.move("data1."+data0_lettercode, "eqpt_files/data1."+data0_lettercode)
                    else:
                        shutil.copyfile("data1."+data0_lettercode, "eqpt_files/data1."+data0_lettercode)
                        
                except:
                    if self.verbose > 0:
                        print('Error: Could not move', "data1."+data0_lettercode, "to eqpt_files")
            
            data1_path = os.path.join(os.getcwd(), "eqpt_files")  # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA

            data0_path = "data0." + data0_lettercode
            
        elif dynamic_db:
            self.__mk_check_del_directory('eqpt_files')

        else:
            data0_path = os.path.join(self.eq36da, f"data0.{data0_lettercode}")

        # gather information from data0 file and perform checks
        if not dynamic_db:
            if os.path.exists(data0_path) and os.path.isfile(data0_path):
                with open(data0_path) as data0:
                    data0_lines = data0.readlines()
                    data0_lines = [line.rstrip()+"\n" for line in data0_lines]
                    start_index = [i+1 for i, s in enumerate(data0_lines) if s == 'temperatures\n']
                    # the T-P grid is followed by 'debye huckel a (adh)' in
                    # B-dot/Davies data0 files and 'debye huckel aphi' in
                    # Pitzer data0 files
                    end_index = [i for i, s in enumerate(data0_lines) if s.strip().lower() in ['debye huckel a (adh)', 'debye huckel aphi']]
                    if len(start_index) == 0 or len(end_index) == 0:
                        self.err_handler.raise_exception("Could not find the temperature and "
                            "pressure grid in {}".format(data0_path))

                    db_grids_unformatted = [i.split("pressures")[0] for i in data0_lines[start_index[0]:end_index[0]]]
                    db_grids = [" ".join(i.split()) for i in db_grids_unformatted if i != '']
                    grid_temps = db_grids[0] + " " + db_grids[1]
                    grid_press = db_grids[2] + " " + db_grids[3]
                    grid_temps = grid_temps.split(" ")
                    grid_press = grid_press.split(" ")

                    try:
                        n_TP_points = data0_lines[2].split("points: ")[1] # extract number of TP points from the third line of data0 file
                        n_TP_points = n_TP_points.replace("\n", "")
                        n_TP_points = int(n_TP_points)
                    except:
                        n_TP_points = 8
                    if n_TP_points == 1:
                        grid_temps = grid_temps[0]
                        grid_press = grid_press[0]

                    try:
                        water_model = data0_lines[1].split("model: ")[1] # extract water model from the second line of data0 file
                        water_model = water_model.replace("\n", "")
                    except:
                        water_model = "SUPCRT92"
    #                     print("Water model could not be referenced from {}".format(data0_path)+""
    #                           ". Defaulting to SUPCRT92 water model...")


                    if(water_model not in ["SUPCRT92", "IAPWS95", "DEW"]):
                        water_model = "SUPCRT92" # the default for EQ3/6
                        print("Water model given in {}".format(data0_path)+" was not "
                              "recognized. Defaulting to SUPCRT92 water model...")
                    
            else: # if a data0 file can't be found, assume default water model, 0-350 C and PSAT
                water_model = "SUPCRT92"
                grid_temps = ["0.0100", "50.0000", "100.0000", "150.0000",
                             "200.0000", "250.0000", "300.0000", "350.0000"]
                grid_press = ["1.0000", "1.0000", "1.0132", "4.7572",
                              "15.5366", "39.7366", "85.8379", "165.2114"]
                
            grid_press_numeric = [float(n) for n in grid_press]
   
            list_tp = check_TP_grid(grid_temps=grid_temps,
                                    grid_press=grid_press,
                                    P1=min(grid_press_numeric),
                                    water_model=water_model,
                                    check_for_errors=False,
                                    verbose=self.verbose)

            grid_temps = list_tp["grid_temps"]
            grid_press = list_tp["grid_press"]
            poly_coeffs_1 = list_tp["poly_coeffs_1"]
            poly_coeffs_2 = list_tp["poly_coeffs_2"]
            
            
        else:
            grid_temps = None
            grid_press = None
            poly_coeffs_1 = None
            poly_coeffs_2 = None

        # handle Alter/Suppress options
        # e.g. [["CaCl+", "AugmentLogK", -1], ["CaOH+", "Suppress"]]
        alter_options_dict = {}
        if len(alter_options) > 0:
            for ao in alter_options:
                if not isinstance(ao, list):
                    err = ("alter_options must be a list of lists, e.g.,\n"
                          "[['CaCl+', 'AugmentLogK', -1], ['CaOH+', 'Suppress']]"
                          "\nor\n[['CaHCO3+', 'Suppress']]")
                    self.err_handler.raise_exception(err)
                key = ao[0]
                if ao[1] == "Suppress" and len(ao) == 2:
                    ao += ["0"]
                alter_options_dict[key] = list(ao[1:])
        alter_options = alter_options_dict
        
        input_dir = "rxn_3i"
        output_dir = "rxn_3o"
        pickup_dir = "rxn_3p"

        # If using column-based charge balancing, read the charge_balance_on values
        # before preprocessing (since they will be excluded during preprocessing)
        charge_balance_dict = {}
        if charge_balance_on == 'column':
            df_temp = pd.read_csv(input_filename, header=0)
            if 'charge_balance_on' in df_temp.columns:
                sample_names = df_temp.iloc[1:, 0].tolist()  # First column contains sample names, skip subheader row
                cb_values = df_temp['charge_balance_on'].iloc[1:].tolist()  # Skip subheader row
                for i, sample_name in enumerate(sample_names):
                    charge_balance_dict[sample_name] = cb_values[i]
            else:
                # No charge_balance_on column found, default to no charge balancing
                if self.verbose > 0:
                    print("No 'charge_balance_on' column found in input file. "
                          "Defaulting to no charge balancing for all samples.")
                charge_balance_on = 'none'

        # preprocess for EQ3
        input_processed = preprocess(input_filename=input_filename,
                                               exclude=exclude_modified,
                                               grid_temps=grid_temps,
                                               grid_press=grid_press,
                                               strict_minimum_pressure=strict_minimum_pressure,
                                               dynamic_db=dynamic_db,
                                               poly_coeffs_1=poly_coeffs_1,
                                               poly_coeffs_2=poly_coeffs_2,
                                               water_model=water_model,
                                               verbose=self.verbose)

        self.df_input_processed = input_processed["df"]

        if blanks_are_0:
            self.df_input_processed = self.df_input_processed.fillna(1E-18)
        
        self.__mk_check_del_directory('rxn_3i')
        self.__mk_check_del_directory('rxn_3o')
        self.__mk_check_del_directory('rxn_3p')
        if dynamic_db:
            self.__mk_check_del_directory('rxn_data0')

        # Has the user been warned about redox column during write_3i_file()?
        # Prevents repeated warnings.
        warned_about_redox_column = False
        
        self.batch_T = input_processed["temp_degC"]
        self.batch_P = input_processed["pressure_bar"]
        
        # create and run a 3i file for each sample
        for sample_row_index in range(0, self.df_input_processed.shape[0]):
            
            temp_degC = input_processed["temp_degC"][sample_row_index]
            pressure_bar = input_processed["pressure_bar"][sample_row_index]

            df = self.df_input_processed.iloc[[sample_row_index]] # double brackets to keep as df row instead of series
            
            samplename = str(df.index[0])

            # handle dynamic data0 creation
            if dynamic_db:

                self.__fill_data0(thermo_df=thermo_df,
                                  data0_file_lines=copy.deepcopy(data0_file_lines),
                                  grid_temps=[temp_degC],
                                  grid_press=[pressure_bar],
                                  db=data0_lettercode,
                                  water_model=water_model,
                                  activity_model=activity_model,
                                  P1=P1,
                                  plot_poly_fit=plot_poly_fit,
                                  logK_extrapolate=logK_extrapolate,
                                  dynamic_db=dynamic_db,
                                  verbose=verbose)

                if self.thermo.thermo_db_type != "data1":
                    self.runeqpt(data0_lettercode, dynamic_db=True)
                
                if os.path.exists("data1."+data0_lettercode) and os.path.isfile("data1."+data0_lettercode):
                    # store contents of data1 file in AqEquil object
                    with open("data1."+data0_lettercode, mode='rb') as data1:
                        self.data1[samplename] = data1.read()
                    try:
                        # move data1
                        shutil.move("data1."+data0_lettercode, "eqpt_files/data1."+data0_lettercode)
                    except:
                        if self.verbose > 0:
                            print('Error: Could not move', "data1."+data0_lettercode, "to eqpt_files")

                data1_path = os.path.join(os.getcwd(), "eqpt_files")  # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA

                data0_path = "data0." + data0_lettercode

            else:
                pressure_bar = input_processed["pressure_bar"][sample_row_index]
                data1_path = self.thermo.eq36da
            
            # allowed aq block species are left after any category exclusion in db_args
            allowed_aq_block_species = ["all"]
            if dynamic_db:
                allowed_aq_block_species = list(thermo_df["name"]) + FIXED_SPECIES

            # Handle per-sample charge balance constraint
            if charge_balance_on == 'column':
                # Get charge balance constraint from the dictionary for this sample
                samplename_for_cb = self.df_input_processed.iloc[sample_row_index, self.df_input_processed.columns.get_loc("Sample")]
                sample_charge_balance = charge_balance_dict.get(samplename_for_cb, 'none')
                if pd.isna(sample_charge_balance) or sample_charge_balance == '':
                    sample_charge_balance = 'none'
                else:
                    sample_charge_balance = str(sample_charge_balance)
            else:
                sample_charge_balance = charge_balance_on

            warned_about_redox_column = write_3i_file(df=df,
                               temp_degC=temp_degC,
                               pressure_bar=pressure_bar,
                               minimum_pressure=input_processed["minimum_pressure"],
                               strict_minimum_pressure=strict_minimum_pressure,
                               pressure_override=dynamic_db,
                               suppress_missing=suppress_missing,
                               exclude=input_processed["exclude"],
                               allowed_aq_block_species=allowed_aq_block_species,
                               charge_balance_on=sample_charge_balance,
                               suppress=suppress,
                               alter_options=alter_options,
                               aq_scale=aq_scale,
                               get_solid_solutions=get_solid_solutions,
                               input_dir=input_dir,
                               redox_flag=redox_flag,
                               redox_aux=redox_aux,
                               default_logfO2=default_logfO2,
                               water_model=water_model,
                               warned_about_redox_column=warned_about_redox_column,
                               activity_model=activity_model,
                               verbose=self.verbose)

            # run EQ3 on each 3i file
            samplename = self.df_input_processed.iloc[sample_row_index, self.df_input_processed.columns.get_loc("Sample")]
            filename_3i = self.df_input_processed.index[sample_row_index]+".3i"
            filename_3o = filename_3i[:-1] + 'o'
            filename_3p = filename_3i[:-1] + 'p'
            
            if dynamic_db:
                dynamic_db_name = self.thermo.thermo_db_filename
            else:
                dynamic_db_name = None
            
            self.runeq3(filename_3i=filename_3i,
                        db=data0_lettercode,
                        samplename=samplename,
                        path_3i=input_dir,
                        path_3o=output_dir,
                        path_3p=pickup_dir,
                        data1_path=data1_path,
                        dynamic_db_name=dynamic_db_name)

            # store input, output, and pickup as dicts in AqEquil object
            try:
                with open(input_dir + "/" + filename_3i, "r") as f:
                    lines=f.readlines()
                self.raw_3_input_dict[samplename] = lines
            except:
                pass
            try:
                with open(output_dir + "/" + filename_3o, "r") as f:
                    lines = [line.rstrip() for line in f.readlines()]
                self.raw_3_output_dict[samplename] = lines
                EQ3_errors_found = self._report_3o_6o_errors(lines, samplename)
            except:
                pass
            try:
                with open(pickup_dir + "/" + filename_3p, "r") as f:
                    lines=f.readlines()
                    
                # capture everything after "start of the bottom half" of 3p
                top_half = []
                bottom_half = []
                capture_it = False
                for line in lines:
                    if "Start of the bottom half of the input file" in line:
                        capture_it = True
                    if capture_it:
                        bottom_half.append(line)
                    else:
                        top_half.append(line)
                        
                self.raw_3_pickup_dict_top[samplename] = top_half # top half of the 3p file, including header for mixing calcs
                self.raw_3_pickup_dict_bottom[samplename] = bottom_half # the bottom half
                
            except:
                pass
            
            if dynamic_db:
                shutil.move("data0.dyn", "rxn_data0/"+filename_3i[0:-3]+"_data0.dat")

        if self.thermo.custom_data0:
            # delete straggling data1 files generated after running eq3
            if os.path.exists("data1") and os.path.isfile("data1"):
                os.remove("data1")
        
        files_3o = [file+".3o" for file in self.df_input_processed.index]
        
        df_input_processed_names = list(self.df_input_processed.columns)

        batch_3o = main_3o_mine(
            files_3o=files_3o,
            get_aq_dist=get_aq_dist,
            get_mass_contribution=get_mass_contribution,
            get_mineral_sat=get_mineral_sat,
            get_redox=get_redox,
            get_charge_balance=get_charge_balance,
            get_ion_activity_ratios=get_ion_activity_ratios,
            get_fugacity=get_fugacity,
            get_basis_totals=get_basis_totals,
            get_solid_solutions=get_solid_solutions,
            mass_contribution_other=mass_contribution_other,
            csv_filename=report_filename,
            aq_dist_type=aq_dist_type,
            mineral_sat_type=mineral_sat_type,
            redox_type=redox_type,
            input_filename=input_filename,
            input_pressures=input_processed["pressure_bar"],
            batch_3o_filename=batch_3o_filename,
            df_input_processed=self.df_input_processed,
            df_input_processed_names=df_input_processed_names,
            verbose=self.verbose,
        )

        if len(batch_3o) == 0:
            self.err_handler.raise_exception("Could not compile a speciation report. This is "
                            "likely because errors occurred during "
                            "the speciation calculation.")
            return
        
        if get_mass_contribution:
            mass_contribution = batch_3o['mass_contribution']

        df_report = batch_3o['report']
        report_divs = batch_3o['report_divs']

        input_cols = report_divs['input']
        df_input = df_report[input_cols].copy()
        
        # add a pressure column to df_input
        df_input["Pressure_bar"] = pd.Series(dtype='float')
        sample_data = batch_3o['sample_data']
        for sample_name, sample in sample_data.items():
            df_input.loc[sample['name'], "Pressure_bar"] = float(sample['pressure'])
        report_divs['input'] = input_cols + ["Pressure_bar"]
            
        # handle headers and subheaders of input section
        headers = [col.split("_")[0] for col in list(df_input.columns)]
        headers = ["pH" if header == "H+" else header for header in headers]
        headers = [header+"_(input)" if header not in ["Temperature", "logfO2", "Pressure"]+exclude else header for header in headers]
        report_divs['input'] = headers # modify headers in the 'input' section of report_divs
        subheaders = [subheader[1] if len(subheader) > 1 else "" for subheader in [
            col.split("_") for col in list(df_input.columns)]]
        multicolumns = pd.MultiIndex.from_arrays(
            [headers, subheaders], names=['Sample', ''])
        
        df_input.columns = multicolumns

        df_join = df_input

        if get_aq_dist:
            aq_distribution_cols = report_divs['aq_distribution']
            df_aq_distribution = df_report[aq_distribution_cols]
            df_aq_distribution = df_aq_distribution.apply(pd.to_numeric, errors='coerce')
            
            # create a pH column from H+
            df_aq_distribution["pH"] = np.nan # pH values are assigned when sample data is assembled later
            
            # handle headers of aq_distribution section
            headers = df_aq_distribution.columns
            subheaders = [aq_dist_type]*(len(headers)-1) # -1 because the last column will have subheader pH (see next line)
            subheaders = subheaders + ["pH"]
            
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_aq_distribution.columns = multicolumns
            
            # ensure final pH column is included in report_divs aq_distribution section
            report_divs["aq_distribution"] = list(headers)
            
            df_join = df_join.join(df_aq_distribution)

        if get_mineral_sat:
            mineral_sat_cols = report_divs['mineral_sat']
            mineral_sat_cols = [col for col in mineral_sat_cols if col != "df"] # TO DO: why is df appearing in mineral sat cols and redox sections?
            df_mineral_sat = df_report[mineral_sat_cols]
            df_mineral_sat = df_mineral_sat.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_mineral_sat section
            if mineral_sat_type == "affinity":
                mineral_sat_unit = "affinity_kcal"
            elif mineral_sat_type == "logQoverK":
                mineral_sat_unit = "logQ/K"
            else:
                self.err_handler.raise_exception(
                    "mineral_sat_type must be either 'affinity' or 'logQoverK'")

            headers = df_mineral_sat.columns
            subheaders = [mineral_sat_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_mineral_sat.columns = multicolumns
            df_join = df_join.join(df_mineral_sat)

        if get_redox:
            redox_cols = report_divs['redox']
            redox_cols = [col for col in redox_cols if col != "df"] # TO DO: why is df appearing in mineral sat cols and redox sections?
            df_redox = df_report[redox_cols]
            df_redox = df_redox.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_redox section
            if redox_type == "Eh":
                redox_unit = "Eh_volts"
            elif redox_type == "pe":
                redox_unit = "pe"
            elif redox_type == "logfO2":
                redox_unit = "logfO2"
            elif redox_type == "Ah":
                redox_unit = "Ah_kcal"
            else:
                self.err_handler.raise_exception(
                    "redox_type must be either 'Eh', 'pe', 'logfO2', or 'Ah'")

            headers = df_redox.columns
            subheaders = [redox_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_redox.columns = multicolumns
            df_join = df_join.join(df_redox)

        if get_charge_balance:
            charge_balance_cols = report_divs['charge_balance']
            df_charge_balance = df_report[charge_balance_cols]
            df_charge_balance = df_charge_balance.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_charge_balance section
            headers = df_charge_balance.columns
            subheaders = ["%"]*2 + ['eq/kg.H2O', 'molality'] + \
                ['eq/kg.H2O']*4 + ['molality']
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_charge_balance.columns = multicolumns
            df_join = df_join.join(df_charge_balance)
            
        if get_ion_activity_ratios:
            if 'ion_activity_ratios' in list(report_divs.keys()):
                ion_activity_ratio_cols = report_divs['ion_activity_ratios']

                df_ion_activity_ratios = df_report[ion_activity_ratio_cols]
                df_ion_activity_ratios = df_ion_activity_ratios.apply(pd.to_numeric, errors='coerce')

                # handle headers of df_ion_activity_ratios section
                headers = df_ion_activity_ratios.columns
                #headers = [header+"_Log_ion-H+_activity_ratio" for header in headers]
                subheaders = ["Log ion-H+ activity ratio"]*len(headers)
                multicolumns = pd.MultiIndex.from_arrays(
                    [headers, subheaders], names=['Sample', ''])
                df_ion_activity_ratios.columns = multicolumns
                df_join = df_join.join(df_ion_activity_ratios)
            
        if get_fugacity:
            fugacity_cols = report_divs['fugacity']
            fugacity_cols = [col for col in fugacity_cols if col != "df"] # TO DO: why is df appearing in fugacity, mineral sat cols and redox sections?
            df_fugacity = df_report[fugacity_cols]
            df_fugacity = df_fugacity.apply(pd.to_numeric, errors='coerce')
            
            # handle headers of fugacity section
            headers = df_fugacity.columns
            subheaders = ["log_fugacity"]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_fugacity.columns = multicolumns
            df_join = df_join.join(df_fugacity)

        if get_basis_totals:
            sc_cols = report_divs['basis_totals']
            df_sc = df_report[sc_cols]
            df_sc = df_sc.apply(pd.to_numeric, errors='coerce')
            
            # handle headers of basis_totals section
            headers = df_sc.columns
            subheaders = ["molality"]*(len(headers))
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_sc.columns = multicolumns
            df_join = df_join.join(df_sc)

        out_dict = {'sample_data': {},
                    'report': df_join,
                    'input': df_input,
                    'report_divs': report_divs}
        
        if get_mass_contribution:
            out_dict['mass_contribution'] = mass_contribution

        sample_data = batch_3o['sample_data']

        # assemble sample data
        for sample_name, sample in sample_data.items():
            dict_sample_data = {
                "filename": sample['filename'],
                "name": sample['name'],
                "temperature": float(sample['temperature']),
                "pressure": float(sample['pressure']),
                "logact_H2O": float(sample['logact_H2O']),
                "H2O_density": float(sample['H2O_density']),
                "H2O_molality": float(sample['H2O_molality']),
                "H2O_log_molality": float(sample['H2O_log_molality']),
                "ionic_strength": float(sample['ionic_strength']),
                }

            if get_aq_dist:
                sample_aq_dist = sample['aq_distribution']
                sample_aq_dist = sample_aq_dist.apply(pd.to_numeric, errors='coerce')

                sample_pH = -sample_aq_dist.loc["H+", "log_activity"]
                out_dict["report"].loc[sample['name'], "pH"] = sample_pH

                dict_sample_data.update({"aq_distribution": sample_aq_dist})

            if get_mass_contribution:
                sample_mass_contribution = mass_contribution[mass_contribution["sample"] == sample['name']]
                dict_sample_data.update(
                    {"mass_contribution": sample_mass_contribution})

            if get_mineral_sat:
                mineral_sat_df = sample['mineral_sat']
                mineral_sat_df = mineral_sat_df.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"mineral_sat": mineral_sat_df})
                # replace sample mineral_sat entry with None if there is no mineral saturation data.
                if(len(dict_sample_data['mineral_sat'].index) == 1 and dict_sample_data['mineral_sat'].index[0] == 'None'):
                    dict_sample_data['mineral_sat'] = None

            if get_redox:
                redox_df = sample['redox']
                redox_df = redox_df.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"redox": redox_df})

            if get_charge_balance:
                dict_sample_data.update({"charge_balance": df_charge_balance.loc[sample['name'], :]})

            if get_ion_activity_ratios:
                if 'ion_activity_ratios' in sample:
                    dict_sample_data.update({"ion_activity_ratios": sample['ion_activity_ratios']})
                else:
                    dict_sample_data['ion_activity_ratios'] = None

            if get_fugacity:
                fugacity_df = sample['fugacity']
                fugacity_df = fugacity_df.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"fugacity": fugacity_df})
                # replace sample fugacity entry with None if there is no fugacity data.
                if(len(dict_sample_data['fugacity'].index) == 1 and dict_sample_data['fugacity'].index[0] == 'None'):
                    dict_sample_data['fugacity'] = None
                else:
                    dict_sample_data["fugacity"]["fugacity"] = 10**dict_sample_data["fugacity"]["log_fugacity"]

            if get_basis_totals:
                sc_dist = sample['basis_totals']
                sc_dist = sc_dist.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"basis_totals": sc_dist})

            if get_solid_solutions:
                sample_solid_solutions = sample.get('solid_solutions')

                if sample_solid_solutions is not None:

                    ss_df_list = []
                    for ss_name, ss_data in sample_solid_solutions.items():
                        df_ss_ideal = ss_data.get("ideal solution", pd.DataFrame())
                        df_ss_mineral = ss_data.get("mineral", pd.DataFrame())
                        if df_ss_mineral.empty:
                            continue
                        if not df_ss_ideal.empty and 'component' in df_ss_ideal.columns:
                            df_merged = pd.merge(df_ss_mineral, df_ss_ideal, left_on='mineral', right_on='component', how='left')
                            del df_merged['component']
                        else:
                            df_merged = df_ss_mineral.copy()
                        df_merged.insert(0, 'solid solution', ss_name)
                        ss_df_list.append(df_merged)

                    dict_sample_data.update(
                        {"solid_solutions": pd.concat(ss_df_list)})

            out_dict["sample_data"].update({sample_name: dict_sample_data})

        out_dict.update({"batch_3o": batch_3o})
        
        out_dict.update({"water_model":water_model, "grid_temps":grid_temps, "grid_press":grid_press,
                         "activity_model":activity_model})
        
        speciation = Speciation(out_dict, hide_traceback=self.hide_traceback)

        speciation.half_cell_reactions = self.half_cell_reactions
        speciation._half_cell_reactions_original_copy = self._half_cell_reactions_original_copy

        if report_filename != None:
            # Prepare report for CSV export with proper header structure
            import csv

            report_df = out_dict["report"].copy()

            # Store the multiindex column structure
            col_level_0 = [col[0] for col in report_df.columns]
            col_level_1 = [col[1] for col in report_df.columns]

            # Create a temporary dataframe with proper structure for CSV
            # First, reset the index to make it a column
            report_df_export = report_df.reset_index()

            # Rename the index column to 'Sample'
            report_df_export.rename(columns={'index': 'Sample'}, inplace=True)

            # Create header rows
            header_row_0 = ['Sample'] + col_level_0
            header_row_1 = [''] + col_level_1

            # Write the CSV with proper formatting
            filename = report_filename if ".csv" in report_filename[-4:] else report_filename + ".csv"

            # Write headers using csv.writer for proper escaping
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header_row_0)
                writer.writerow(header_row_1)

            # Append the data without headers
            report_df_export.to_csv(filename, mode='a', header=False, index=False)

        if delete_generated_folders:
            self._delete_rxn_folders()
            try:
                # delete straggler data1 file
                os.remove("data1")
            except:
                pass
        
        if self.verbose > 0:
            print("Finished!")
        
        speciation.raw_3_input_dict = self.raw_3_input_dict
        speciation.raw_3_output_dict = self.raw_3_output_dict
        speciation.raw_3_pickup_dict_top = self.raw_3_pickup_dict_top
        speciation.raw_3_pickup_dict_bottom = self.raw_3_pickup_dict_bottom
        speciation.raw_6_input_dict = {}
        speciation.raw_6_output_dict = {}
        speciation.raw_6_pickup_dict = {}
        speciation.raw_6_pickup_dict_top = {}
        speciation.thermo = self.thermo
        speciation.data1 = self.data1
        speciation.verbose = self.verbose
        
        speciation.logK_models = self.logK_models
        speciation.batch_T = self.batch_T
        speciation.batch_P = self.batch_P
        
        return speciation
    

    @staticmethod
    def __s_d(x, k):
        # specify how many decimals are printed
        # e.g. 12.433 becomes "12.4330" if k=4
        kstr = '{:.'+str(k)+'f}'
        return kstr.format(round(x, k)).strip()
    
    
    def __fill_data0(self, thermo_df, data0_file_lines, grid_temps, grid_press, db,
                   water_model, activity_model, P1, plot_poly_fit, logK_extrapolate,
                   dynamic_db, verbose):

        list_tp = check_TP_grid(grid_temps=grid_temps,
                                 grid_press=grid_press,
                                 P1=P1,
                                 water_model=water_model,
                                 check_for_errors=True,
                                 verbose=self.verbose)

        grid_temps = list_tp["grid_temps"]
        grid_press = list_tp["grid_press"]
        
        if plot_poly_fit and len(grid_temps) >= 8:
            self.__plot_TP_grid_polyfit(xvals=grid_temps,
                                        yvals=grid_press,
                                        poly_coeffs_1=list_tp["poly_coeffs_1"],
                                        poly_coeffs_2=list_tp["poly_coeffs_2"],
                                        res=500)

        # calculate logK at each T and P for every species
        out_dfs = []
        for i,Tc in enumerate(grid_temps):
            out_dfs.append(calc_logK(thermo_df, Tc=Tc, P=grid_press[i], TP_i=i, water_model=water_model))
        
        dissrxn_logK_dict = {'name': out_dfs[0]["name"],
                             'logK_0': out_dfs[0]["dissrxn_logK_0"]}
        
        if len(grid_temps) >= 8:
            for i in range(1, len(grid_temps)):
                dissrxn_logK_dict['logK_'+str(i)] = out_dfs[i]["dissrxn_logK_"+str(i)]
                
        if len(grid_temps) == 1:
            dissrxn_logK_dict['logK_1'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_2'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_3'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_4'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_5'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_6'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
            dissrxn_logK_dict['logK_7'] = float('nan') #out_dfs[0]["dissrxn_logK_0"]
    
        dissrxn_logK = pd.DataFrame(dissrxn_logK_dict)
        
        # remove duplicate rows (e.g., for mineral polymorphs)
        dissrxn_logK = dissrxn_logK.drop_duplicates("name")

        # handle free logK values
        free_logK_names = []
        if "logK1" in thermo_df.columns:
            
            free_logK_df = thermo_df.dropna(subset=['logK1'])
            free_logK_names = list(free_logK_df["name"])
    
            sp_dupes = []
            for i,sp in enumerate(free_logK_names):
                logK_grid = list(free_logK_df[["logK1", "logK2", "logK3",
                                               "logK4", "logK5", "logK6",
                                               "logK7", "logK8"]].iloc[i]) # logK at T and P in datasheet
                
                T_grid = list(free_logK_df[["T1", "T2", "T3",
                                            "T4", "T5", "T6",
                                            "T7", "T8"]].iloc[i]) # T for free logK grid
                
                P_grid = list(free_logK_df[["P1", "P2", "P3",
                                            "P4", "P5", "P6",
                                            "P7", "P8"]].iloc[i]) # P for free logK grid
                
                
                for ii,T in enumerate(grid_temps):
                    
                    logK, model = self._interpolate_logK(T, logK_grid, T_grid, logK_extrapolate)
                    
                    dissrxn_logK.loc[(dissrxn_logK.name == sp), "logK_"+str(ii)] = logK
                
                
                self.logK_models[sp] = {"logK_grid":logK_grid,
                                        "T_grid":T_grid,
                                        "P_grid":P_grid,
                                        "logK_extrapolate":logK_extrapolate,
                                        "type":"free logK values",
                                        }
        
                # check that there aren't duplicates between OBIGT-style datasheet and
                # the 'free logK' datasheet
                if sp in dissrxn_logK["name"]:
                    sp_errs.append(sp)
                    
            if len(sp_dupes) > 0:
                msg = ("The following species are duplicated between the "
                       "thermodynamic datafiles used: " + ",".join(sp_errs))
                self.err_handler.raise_exception(msg)

        # calculate and process logK values of species in the OBIGT-style datasheet
        for idx in range(0, dissrxn_logK.shape[0]):

            name = dissrxn_logK.iloc[idx, dissrxn_logK.columns.get_loc('name')]

            # format the logK reaction block of this species' data0 entry
            logK_grid = list(dissrxn_logK.iloc[idx, 1:len(dissrxn_logK)+1])
            
            if not dynamic_db:
                if name not in self.logK_models.keys() and name not in free_logK_names:
                    self.logK_models[name] = {"logK_grid":logK_grid,
                                  "T_grid":grid_temps,
                                  "P_grid":grid_press,
                                  "logK_extrapolate":logK_extrapolate,
                                  "type":"calculated logK values",
                                  }

            elif dynamic_db:
                if name not in self.logK_models.keys() and name not in free_logK_names:
                    self.logK_models[name] = {"logK_grid":[logK_grid[0]],
                                              "T_grid":[grid_temps[0]],
                                              "P_grid":[grid_press[0]],
                                              "logK_extrapolate":"no fit",
                                              "type":"calculated logK values",
                                              }
                    
                elif name not in free_logK_names:
                    self.logK_models[name]["logK_grid"] += [logK_grid[0]]
                    self.logK_models[name]["T_grid"] += [grid_temps[0]]
                    self.logK_models[name]["P_grid"] += [grid_press[0]]

            # filter out strict basis species
            # TODO: do this by species tag, not just whether it has a logK grid of all 0s
            if len(set(logK_grid)) == 1:
                if set(logK_grid) == set([0]):
                    continue

            # loop through logK values and format for data0
            logK_list = []
            for i in range(0, len(logK_grid)):
                logK_val = self.__s_d(logK_grid[i], 4)
                
                # conditional formatting based on position
                if (i+1) == 1 or (i+1) % 7 == 0: # first entry of a line
                    max_length = 11
                    end_char = ""
                elif (i+1) % 6 == 0 and (i+1) != len(logK_grid): # last entry of a line
                    max_length = 6
                    end_char = "\n"
                else:
                    max_length = 6
                    end_char = ""

                # get decimal position and format spaces accordingly
                decimal_position = logK_val.find(".")
                logK_val = "".join([" "]*(max_length-decimal_position)) + logK_val + end_char
                # append to logk list
                logK_list.append(logK_val)

            logK_list = "".join(logK_list)
            
            # todo: make this more robust to catch any potential logK_grid skips
            if "logK_grid_"+name in data0_file_lines:
                data0_file_lines[data0_file_lines.index("logK_grid_"+name)] = logK_list
        
        # handle data0 header section using Python version
        data0_file_lines = fill_data0_head(data0_template=data0_file_lines,
                                           db=db,
                                           grid_temps=grid_temps,
                                           grid_press=grid_press,
                                           water_model=water_model,
                                           activity_model=activity_model)

        with open("data0."+db, 'w') as f:
            for item in data0_file_lines:
                f.write("%s" % item)

    
    def plot_logK_fit(self, name, plot_out=False, res=200, internal=True,
                      logK_extrapolate=None, T_vals=[], plot_association=False):
        """
        Plot the fit of dissociation reaction logK values used in the
        speciation.

        Parameters
        ----------
        name : str
            Name of the chemical species.
        
        plot_out : bool, default False
            Return a Plotly figure object? If False, a figure is simply shown.
            If True, the function returns a Plotly figure object and does
            not show the plot.
        
        res : int
            Resolution of the fit line. Higher resolutions will be smoother.
            
        internal : bool, default True
            Reuse calculated fits if they already exist?
        
        logK_extrapolate : str, optional
            Option for extrapolating logK values in the plot. Possible values
            for this parameter include 'poly', 'linear', 'flat', or 'none'.
            This is for planning and visualization only and does not affect
            results in `speciate()` or `create_data0()`. Those functions have
            their own parameters for setting logK extrapolation options.
        
        T_vals : list, optional
            Option for visualizing how the fit of logK values will be
            used to estimate the logK values at the temperatures specified in
            the list given to this parameter. This is useful for visualizing
            logK extrapolation options defined by `logK_extrapolate`.

        plot_association : bool, default False
            Plot species association constants instead of dissociation
            constants?
        
        Returns
        ----------
        fig : a Plotly figure object
            Returned if `plot_out` is True.

        """

        if internal and len(self.logK_models.keys()) > 0:
            # use internally calculated logK models already stored...
            if name not in self.logK_models.keys():
                if name not in list(self.thermo.df_rejected_species["name"]):
                    msg = "The chemical species " + str(name) + " is not recognized."
                    self.err_handler.raise_exception(msg)
                else:
                    reject_reason = list(self.thermo.df_rejected_species.loc[self.thermo.df_rejected_species['name'] == name, 'reason for rejection'])[0]
                    
                    msg = ("The chemical species " + str(name) + " cannot be "
                           "plotted because it was rejected from the "
                           "speciation:\n" + str(reject_reason))
                    self.err_handler.raise_exception(msg)

            logK_grid = self.logK_models[name]["logK_grid"]
            T_grid = self.logK_models[name]["T_grid"]
            P_grid = self.logK_models[name]["P_grid"]
        
        else:
            # load logK models from Thermodata class's logK_db
            df_logK = self.thermo.logK_db
            
            i = list(df_logK["name"]).index(name)
            
            logK_grid = list(df_logK[["logK1", "logK2", "logK3",
                                      "logK4", "logK5", "logK6",
                                      "logK7", "logK8"]].iloc[i]) # logK at T and P in datasheet

            T_grid = list(df_logK[["T1", "T2", "T3",
                                   "T4", "T5", "T6",
                                   "T7", "T8"]].iloc[i]) # T for free logK grid

            P_grid = list(df_logK[["P1", "P2", "P3",
                                   "P4", "P5", "P6",
                                   "P7", "P8"]].iloc[i]) # P for free logK grid
            
            if not isinstance(logK_extrapolate, str):
                logK_extrapolate = self.thermo.logK_extrapolate
        
        if not isinstance(logK_extrapolate, str):
            logK_extrapolate = self.logK_models[name]["logK_extrapolate"]
        
        if len(T_vals) == 0:
            grid_temps = self.batch_T
        else:
            grid_temps = T_vals
        
        grid_press = self.batch_P
        
        T_grid = [t for t in T_grid if not pd.isna(t)]
        P_grid = [p for p in P_grid if not pd.isna(p)]
        logK_grid = [k for k in logK_grid if not pd.isna(k)]

        if len(logK_grid) == 0 or len(T_grid) == 0 or len(P_grid) == 0:
            self.err_handler.raise_exception("This species has no valid logK "
                    "values, temperature values, or pressure values. It is "
                    "possible that this is a strict basis species with no "
                    "dissociation reaction for which to calculate logK values.")

        if plot_association:
            logK_grid = [-y for y in logK_grid]
        
        fig = px.scatter(x=T_grid, y=logK_grid)
        
        if len(grid_temps) > 0:
            if min(grid_temps) <= min(T_grid):
                plot_T_min = min(grid_temps)
            else:
                plot_T_min = min(T_grid)
            if max(grid_temps) >= max(T_grid):
                plot_T_max = max(grid_temps)
            else:
                plot_T_max = max(T_grid)
        else:
            plot_T_min = min(T_grid)
            plot_T_max = max(T_grid)
        
        plot_temps = np.linspace(plot_T_min, plot_T_max, res)

        pred_logK = []
        pred_model = []
        for t in plot_temps:
            logK, model = self._interpolate_logK(t, logK_grid, T_grid, logK_extrapolate)
            pred_logK.append(logK)
            pred_model.append(model)
        
        df_plot = pd.DataFrame({"T":plot_temps, "logK":pred_logK, "model":pred_model})
        
        if logK_extrapolate != "no fit":
            fig = px.line(df_plot, x='T', y='logK', color='model', title=name, template="simple_white")
        else:
            fig = px.line(x=[0], y=[0], title=name, template="simple_white") # dummy figure
            
        fig.update_traces(hovertemplate="T = %{x} °C<br>Predicted logK = %{y}<extra></extra>")
        fig.update_layout(xaxis_range=[min(plot_temps) - 0.15*(max(plot_temps) - min(plot_temps)),
                                       max(plot_temps) + 0.15*(max(plot_temps) - min(plot_temps))],
                          xaxis_title="T,°C", yaxis_title="logK")
        
        logK_label = "fitted logK value(s)"
        annotation = ""
        
        if len(grid_temps) > 0:
            for i,gt in enumerate(grid_temps):
                # make vertical lines representing batch temperatures

                if i==0:
                    showlegend=True
                else:
                    showlegend=False

                if isinstance(grid_press, str):
                    ht_samples= "T = "+str(gt) + " °C<br>P = PSAT<extra></extra>"
                else:
                    if len(grid_press) > 0:
                        ht_samples= "T = "+str(gt) + " °C<br>P = " + str(grid_press[i]) + " bar(s)<extra></extra>"
                    else:
                        ht_samples= "T = "+str(gt) + " °C<extra></extra>"
                        
                if len(T_grid) > 1:
                    
                    if logK_extrapolate == "none" and (gt > max(T_grid) or gt < min(T_grid)):
                        viz_logK = max(logK_grid)
                    else:
                        viz_logK, _ = self._interpolate_logK(gt, logK_grid, T_grid, logK_extrapolate)
                    
                    vline_y_vals = [min(logK_grid)-0.15*(max(logK_grid)-min(logK_grid)), viz_logK]
                    
                    
                if logK_extrapolate == "no fit":
                    vline_y_vals = [min(logK_grid)-0.15*(max(logK_grid)-min(logK_grid)), logK_grid[i]]
                    logK_label = "calculated dissociation reaction LogK value(s)"
                    annotation = ("LogK values are calculated from<br>ΔG of dissociation into basis species"
                                  "<br>at the T and P of the speciated samples<br>and do not require a fit.")

                if _all_equal(logK_grid):
                    # if a flat horizontal logK fit line...
                    # then fix the y-axis range to prevent zoomed-in steppy wierdness
                    fig.update_layout(yaxis_range=[logK_grid[0]-1,logK_grid[0]+1])
                    vline_y_vals = [logK_grid[0]-1, logK_grid[0]]

                
                fig.add_trace(
                    go.Scatter(x=[gt, gt],
                               y=vline_y_vals,
                               mode="lines",
                               line=dict(color='rgba(255, 0, 0, 0.75)', width=3, dash="dot"),
                               legendgroup='batch temperatures',
                               name='batch temperatures',
                               showlegend=showlegend,
                               hovertemplate=ht_samples,
                              ),
                )
        
        # add fitted logK points
        fig.add_trace(go.Scatter(x=T_grid, y=logK_grid, name=logK_label,
                                 mode='markers', marker=dict(color="black"),
                                 text = P_grid,
                                 hovertemplate="T = %{x} °C<br>P = %{text} bar(s)<br>logK = %{y}<extra></extra>",
                                 ),
                      )
        
        fig.add_annotation(x=0, y=0, xref="paper", yref="paper", align='left',
                           text=annotation, bgcolor="rgba(255, 255, 255, 0.5)",
                           showarrow=False)
        
        if plot_out:
            return fig
        else:
            fig.show()

        
    def __get_i_of_valid_free_logK_sp(self, free_logK_df, grid_temps,
                                      grid_press, dynamic_db,
                                      logK_extrapolate, db_sp_names):
            """
            Check for species in the free logK database with pressure values that
            are permitted in the context of grid_temps and grid_press, then
            return their indices.
            """
        
            if not isinstance(grid_press, list):
                # "Psat" to ["psat"]
                grid_press_list = [grid_press.lower()]
            else:
                grid_press_list = grid_press

            valid_sp_i = []
            rejected_sp_i_dict = {}
            
            for i,sp in enumerate(list(free_logK_df["name"])):
                
                sp_temps_grid = [free_logK_df.iloc[i]["T"+str(ii)] for ii in range(1,9) if not math.isnan(free_logK_df.iloc[i]["T"+str(ii)])]
                sp_press_grid_init = [float(free_logK_df.iloc[i]["P"+str(ii)]) if free_logK_df.iloc[i]["P"+str(ii)] not in ["Psat", "psat"] else 'psat' for ii in range(1,9)]

                sp_grid_len = len(sp_temps_grid)
                
                sp_press_grid = []
                for p in sp_press_grid_init:
                    if isinstance(p, str):
                        sp_press_grid.append(p)
                    elif not math.isnan(p):
                        sp_press_grid.append(p)
                        
                if dynamic_db:
                    sp_press_grid_numeric = [pychnosz.water("Psat", T=273.15+T, messages=False)+0.0001 if P=="Psat" or P=="psat" or P=="PSAT" else P for P,T in zip(sp_press_grid, sp_temps_grid) ]
                    grid_press_list_numeric = [pychnosz.water("Psat", T=273.15+T, messages=False)+0.0001 if P=="Psat" or P=="psat" or P=="PSAT" else P for P,T in zip(grid_press_list, grid_temps) ]

                    sp_temps_grid_psatted = []
                    for ii,T in enumerate(sp_temps_grid):
                        if T < 100 and sp_press_grid_numeric[ii] == 1:
                            sp_temps_grid_psatted.append(25)
                        else:
                            sp_temps_grid_psatted.append(T)
                    
                    grid_temps_psatted = []
                    for ii,T in enumerate(grid_temps):
                        if T < 100 and grid_press_list_numeric[ii] == 1:
                            grid_temps_psatted.append(25)
                        else:
                            grid_temps_psatted.append(T)
                    
                    sp_TP = [(T,P) for T,P in zip(sp_temps_grid_psatted, sp_press_grid_numeric)]
                    sample_TP = [(T,P) for T,P in zip(grid_temps_psatted, grid_press_list_numeric)]

                
                if sp_press_grid == grid_press_list and sp_temps_grid == grid_temps:
                    # If pressures and temperature grid exactly matches that of the sp...
                    # need to test this!
                    valid_sp_i.append(i)
                elif dynamic_db and all(P in sp_TP for P in sample_TP):
                    # If every sample T,P point matches any of the T,P points in the species grid...
                    valid_sp_i.append(i)
                elif (logK_extrapolate != "none" or (min(grid_temps) >= min(sp_temps_grid) and max(grid_temps) <= max(sp_temps_grid))) and _all_equal(sp_press_grid + grid_press_list):
                    # If all grid temperatures are within minimum and maximum file temperatures,
                    # and all pressures in file for the sp are equal, and all grid pressures match
                    # file pressure, then the species is valid
                    valid_sp_i.append(i)
                else:
                    # species is invalid. Define reasons.
                    
                    reject_reason_list = []

                    if min(grid_temps) < min(sp_temps_grid) and _all_equal(sp_press_grid + grid_press_list) and logK_extrapolate == "none":

                        min_sp = str(min(sp_temps_grid))
                        min_grid = str(min(grid_temps))
                        if dynamic_db:
                            reject_reason_list.append("Minimum temperature in this batch of samples is "+min_grid+"°C, which is below the minimum applicability temperature of this species is "+min_sp+"°C.")
                        else:
                            reject_reason_list.append("Minimum temperature in this data0 file is "+min_grid+"°C, which is below the minimum applicability temperature of this species is "+min_sp+"°C.")
                    
                    if max(grid_temps) > max(sp_temps_grid) and _all_equal(sp_press_grid + grid_press_list) and logK_extrapolate == "none":
                        max_sp = str(max(sp_temps_grid))
                        max_grid = str(max(grid_temps))
                        if dynamic_db:
                            reject_reason_list.append("Maximum temperature in this batch of samples is "+max_grid+"°C, which is above the maximum applicability temperature of this species is "+max_sp+"°C.")
                        else:
                            reject_reason_list.append("Maximum temperature in this data0 file is "+max_grid+"°C, which is above the maximum applicability temperature of this species is "+max_sp+"°C.")
                    
                    if not _all_equal(sp_press_grid + grid_press_list):
                        if dynamic_db:
                            reject_reason_list.append("Mismatch between pressures of samples in this batch and the applicable pressures for "+str(sp)+" given in the logK thermodynamic database.")
                        else:
                            reject_reason_list.append("Mismatch between desired pressure grid of data0 file and the applicable pressures for "+str(sp)+" given in the logK thermodynamic database.")
                    
                    if len(reject_reason_list) == 0:
                        reject_reason_list.append("Unknown")
                    
                    rejected_sp_i_dict[i] = "\n".join(reject_reason_list)

            # loop through valid species and reject them if their dissociation reactions
            # contain species that have been rejected.

            valid_sp_i = list(dict.fromkeys(valid_sp_i))
            while True:
                valid_sp_i_before = copy.deepcopy(valid_sp_i)
                valid_sp_i, rejected_sp_i_dict = self._check_valid_free_logK_sp_dissrxn(valid_sp_i, rejected_sp_i_dict, free_logK_df, db_sp_names)
                if valid_sp_i_before == valid_sp_i:
                    break
            
            reject_indices = list(rejected_sp_i_dict.keys())
            reject_names = list(free_logK_df.iloc[reject_indices]["name"])
            reject_states = list(free_logK_df.iloc[reject_indices]["state"])
            reject_reasons =list(rejected_sp_i_dict.values())
            
            for i,n in enumerate(reject_names):
                self.thermo._reject_species(name=n, reason=reject_reasons[i])
        
            return valid_sp_i
            
            
    def _check_valid_free_logK_sp_dissrxn(self, valid_sp_i, rejected_sp_i_dict, free_logK_df, db_sp_names):
        
        valid_sp_names = list(free_logK_df.iloc[valid_sp_i]["name"])
        
        rejected_sp_names = list(free_logK_df.iloc[list(rejected_sp_i_dict.keys())]["name"])

        for i in valid_sp_i:
            dissrxn_i = free_logK_df.iloc[i]["dissrxn"]
            dissrxn_sp = dissrxn_i.split(" ")[1::2] # get species names from dissrxn
            dissrxn_sp = dissrxn_sp[1:] # ignore the species itself
            
            for sp in dissrxn_sp:
                if sp in rejected_sp_names and sp not in valid_sp_names and sp not in db_sp_names:
                    valid_sp_i.remove(i)
                    rejected_sp_i_dict[i] = "Dissociation reaction contains the species " + sp + ", which has been rejected."
                    return valid_sp_i, rejected_sp_i_dict
                    
        return valid_sp_i, rejected_sp_i_dict
            
        
    def create_data0(self,
                     db,
                     filename_ss=None,
                     activity_model="auto",
                     exceed_Ttr=True,
                     grid_temps=[0.0100, 50.0000, 100.0000, 150.0000,
                                 200.0000, 250.0000, 300.0000, 350.0000],
                     grid_press="Psat",
                     P1=True,
                     plot_poly_fit=False,
                     logK_extrapolate="none",
                     fill_data0=True,
                     dynamic_db=False,
                     dynamic_db_sample_temps=[],
                     dynamic_db_sample_press=[],
                     verbose=1):
        """
        Create a data0 file from a custom thermodynamic dataset.
        
        Parameters
        ----------
        db : str
            Desired three letter code of data0 output.
            
        filename_ss : str, optional
            Name of file containing solid solution parameters.

        activity_model : str, default "auto"
            Aqueous species activity coefficient model the data0 file is
            meant for. Can be "auto", "b-dot", "davies", or "pitzer". The
            "b-dot" and "davies" options both produce a data0 file with B-dot
            (hard core diameter) parameters. The "pitzer" option produces a
            data0 file with Pitzer interaction parameter blocks and requires
            that a Pitzer parameter database was loaded with the `pitzer`
            parameter of `AqEquil`. "auto" chooses "pitzer" if a Pitzer
            parameter database is loaded and "b-dot" otherwise.

        grid_temps : list of eight float, default [0.0100, 50.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000]
            Eight temperature values that make up the T-P grid.
        
        grid_press : list of float, default "Psat"
            Eight pressure values that make up the T-P grid. "Psat" for
            calculations along the liquid-vapor saturation curve.
        
        P1 : bool, default True,
            Use pressure of 1 bar below 100 degrees C instead of calculated
            values of Psat? Ignored if `grid_press` is not "Psat".
        
        plot_poly_fit : bool, default False
            Plot the polynomial fit of the temperature pressure grid?
        
        dynamic_db : bool, default False
            Are data0 files being created dynamically? If unsure, use False.
            Used by `speciate` to display valid messages.
        
        verbose : int, 0, 1, or 2, default 1
            Level determining how many messages are returned during a
            calculation. 2 for all messages, 1 for errors or warnings only,
            0 for silent.
        """

        thermo_df = self.thermo.thermo_db
        db_logK = self.thermo.logK_db
        water_model = self.thermo.water_model

        self.verbose = verbose

        activity_model = self._resolve_activity_model(activity_model)
        pitzer_df = self.thermo.pitzer_db if activity_model == "pitzer" else None
        
        self.batch_T = grid_temps
        self.batch_P = grid_press
        
        if not dynamic_db:
            if self.verbose >= 1:
                print("Creating data0.{}...".format(db), flush=True)
        
        # if len(grid_temps) not in [1, 8]:
        #     self.err_handler.raise_exception("'grid_temps' must have either one or eight values.")
        # if isinstance(grid_press, list):
        #     if len(grid_press) not in [1, 8]:
        #         self.err_handler.raise_exception("'grid_press' must have either one or eight values.")
        
        if sum([T >= 10000 for T in grid_temps]):
            self.err_handler.raise_exception("Grid temperatures must be below 10k °C.")
        
        if isinstance(grid_press, list):
            if sum([P >= 10000 for P in grid_press]) and water_model != "DEW":
                self.err_handler.raise_exception("Grid pressures must be below 10 kilobars.")
                
        if water_model == "SUPCRT92":
            min_T = 0
            max_T = 2250
            min_P = 0
            max_P = 30000
        elif water_model == "IAPWS95":
            min_T = 0
            max_T = 1000
            min_P = 0
            max_P = 10000
        elif water_model == "DEW":
            min_T = 0
            max_T = 1000
            min_P = 1000
            max_P = 60000
        else:
            self.err_handler.raise_exception("The water model '{}' ".format(water_model)+"is not "
                "recognized. Try 'SUPCRT92', 'IAPWS95', or 'DEW'.")
        
        # check that T and P are above minimum values
        if sum([T <= min_T for T in grid_temps]):
            print("WARNING: one or more temperatures in 'grid_temps' is below "
                  "or equal to {} °C".format(min_T)+" and is outside the valid "
                  "temperature range for the {} water model.".format(water_model))
        if isinstance(grid_press, list):
            if sum([P < min_P for P in grid_press]):
                print("WARNING: one or more pressures in 'grid_press' is below "
                      "{} bar".format(min_P)+", the minimum valid "
                      "pressure for the {} water model.".format(water_model))
        
        # check that T and P are below maximum values
        if sum([T > max_T for T in grid_temps]):
            print("WARNING: one or more temperatures in 'grid_temps' is above "
                  "{} °C".format(max_T)+", the maximum valid "
                  "temperature for the {} water model.".format(water_model))
        if isinstance(grid_press, list):
            if sum([P > max_P for P in grid_press]):
                print("WARNING: one or more pressures in 'grid_press' is above "
                      "{} bar".format(max_P)+", the maximum valid "
                      "pressure for the {} water model.".format(water_model))
        
        if water_model == "IAPWS95":
            if self.verbose > 1:
                print("Using the IAPWS95 water model. Water properties (and therefore "
                      "calculations) are available for supercooled water down to "
                      "about -30 °C at 1 bar.")
        elif water_model != "SUPCRT92":
            print("WARNING: water models other than SUPCRT92 and IAPWS95 are not yet fully supported.")
        
        # reset logK_models whenever create_data0() is called
        # (prevents errors when create_data0() functions are run back-to-back)
        self.logK_models = {}
        
        # interpolate logK values from "free logK" datasheet at T and P
        if isinstance(db_logK, pd.DataFrame):

            if len(dynamic_db_sample_temps) > 0:
                grid_or_sample_temps = dynamic_db_sample_temps
            else:
                grid_or_sample_temps = grid_temps
                
            if len(dynamic_db_sample_press) > 0:
                grid_or_sample_press = dynamic_db_sample_press
            else:
                grid_or_sample_press = grid_press
            
            free_logK_df = assign_worm_db_col_dtypes(self.thermo.logK_db)

            valid_i = self.__get_i_of_valid_free_logK_sp(
                free_logK_df,
                grid_or_sample_temps,
                grid_or_sample_press,
                dynamic_db,
                logK_extrapolate,
                db_sp_names=thermo_df["name"],
                )
            free_logK_df_valid = copy.deepcopy(free_logK_df.iloc[valid_i])
            thermo_df = pd.concat([thermo_df, free_logK_df_valid], ignore_index=True)
            
            thermo_df = assign_worm_db_col_dtypes(thermo_df)
        
        if self.thermo.solid_solutions_active:
            solid_solution_df = self.thermo.solid_solution_db
        else:
            solid_solution_df = None

        template = import_package_file(__name__, 'data0.min')

        out_list = self.thermo.out_list

        # Extract dissrxns, basis_pref, and redox_elem_states from out_list
        # out_list is an R list with 'dissrxns', 'basis_pref', and 'redox_elem_states' components
        if hasattr(out_list, 'rx2'):
            # If out_list is still an R object, pass the R objects directly
            dissrxns = out_list.rx2("dissrxns")
            basis_pref = out_list.rx2("basis_pref")
            redox_elem_states = out_list.rx2("redox_elem_states") if "redox_elem_states" in list(out_list.names) else {}
        else:
            # Already a Python dict
            dissrxns = out_list.get("dissrxns", {})
            basis_pref = out_list.get("basis_pref", [])
            redox_elem_states = out_list.get("redox_elem_states", {})

        # assemble data0 file using Python version
        data0_file_lines = create_data0_py(
                          thermo_df=thermo_df,
                          element_df=None,
                          solid_solution_df=solid_solution_df,
                          db=db,
                          water_model=water_model,
                          template=template,
                          dissrxns=dissrxns,
                          basis_pref=basis_pref,
                          exceed_Ttr=exceed_Ttr,
                          fixed_species=FIXED_SPECIES,
                          redox_elem_states=redox_elem_states,
                          activity_model=activity_model,
                          pitzer_df=pitzer_df,
                          verbose=self.verbose,
                          )

        data0_file_lines = data0_file_lines.split("\n")
        
        if fill_data0:

            # begin TP-dependent processes
            self.__fill_data0(thermo_df=thermo_df,
                              data0_file_lines=copy.deepcopy(data0_file_lines),
                              grid_temps=grid_temps,
                              grid_press=grid_press,
                              db=db,
                              water_model=water_model,
                              activity_model=activity_model,
                              P1=P1,
                              plot_poly_fit=plot_poly_fit,
                              logK_extrapolate=logK_extrapolate,
                              dynamic_db=dynamic_db,
                              verbose=self.verbose)
    
        else:
            return thermo_df, data0_file_lines, grid_temps, grid_press, db, water_model, P1, plot_poly_fit

        if self.verbose > 0:
            print("Finished creating data0.{}.".format(db))

    
    class Thermodata(object):
        """
        Metaclass to store and load thermodynamic databases.
        Inherits attributes from its outer class, AqEquil.
        
        """

        def __init__(self, AqEquil_instance):

            self.AqEquil_instance = AqEquil_instance

            # Initialize CHNOSZ thermodynamic system if not already initialized
            # This must happen BEFORE loading custom element databases
            thermo_sys = pychnosz.thermo()
            if not thermo_sys.is_initialized():
                thermo_sys.reset(messages=False)

            # attributes to add to AqEquil class
            self.db = self.AqEquil_instance.db
            self.elements = self.AqEquil_instance.elements
            solid_solutions = self.AqEquil_instance.solid_solutions
            self.exclude_category = self.AqEquil_instance.exclude_category
            self.exclude_organics_except = self.AqEquil_instance.exclude_organics_except
            self.water_model = self.AqEquil_instance.water_model
            self.elements = self.AqEquil_instance.elements
            logK = self.AqEquil_instance.logK
            logK_S = self.AqEquil_instance.logK_S
            pitzer = self.AqEquil_instance.pitzer
            mineral_db = self.AqEquil_instance.mineral_db
            download_csv_files = self.AqEquil_instance.download_csv_files
            suppress_redox = self.AqEquil_instance.suppress_redox
            exceed_Ttr = self.AqEquil_instance.exceed_Ttr
            input_template = self.AqEquil_instance.input_template
            verbose = self.AqEquil_instance.verbose
            
            self.hide_traceback = self.AqEquil_instance.hide_traceback
            self.err_handler = Error_Handler(clean=self.hide_traceback)

            self.eq36da = self.AqEquil_instance.eq36da
            self.eq36co = self.AqEquil_instance.eq36co

            self.df_rejected_species = pd.DataFrame({'database name':[],
                                                     'database index':[],
                                                     "name":[],
                                                     "state":[],
                                                     "reason for rejection":[]})
            
            # active thermo db attributes
            self.thermo_db = None
            self.thermo_db_type = None
            self.thermo_db_source = None
            self.thermo_db_filename = None
            self.custom_data0 = None
            self.data0_lettercode = None
            self.dynamic_db = None
            self.db_csv_name = None

            # data1 attributes
            self.data1 = {}

            # data0 attributes
            self.data0_db = None
            self.data0_db_type = None
            self.data0_db_source = None
            self.data0_db_filename = None

            # csv attributes
            self.csv_db = None
            self.csv_db_type = None
            self.csv_db_source = None
            self.csv_db_filename = None

            # element attributes
            self.element_db = None
            self.element_db_source = None
            self.element_db_filename = None
            self.element_active = None

            # solid solution attributes
            self.solid_solutions_active = False
            self.solid_solution_db = None
            self.solid_solution_db_source = None
            self.solid_solution_db_filename = None

            # logK attributes
            self.logK_active = False
            self.logK_extrapolate = self.AqEquil_instance.logK_extrapolate
            self.logK_db = None
            self.logK_db_source = None
            self.logK_db_filename = None

            # logK attributes
            self.logK_S_active = False
            self.logK_S_db = None
            self.logK_S_db_source = None
            self.logK_S_db_filename = None

            # Pitzer parameter attributes
            self.pitzer_active = False
            self.pitzer_db = None
            self.pitzer_db_source = None
            self.pitzer_db_filename = None

            # activity model that the active database is formatted for:
            # "pitzer" or "sedh" (B-dot/Davies) for data0 and data1 files,
            # None for CSV databases (which can generate either kind of data0)
            self.db_activity_model = None

            self.verbose=verbose

            # "WORM", by itself or as an entry in a list of databases, brings
            # along the WORM element, solid solution, logK, and logK_S databases
            db_entries = self.db if isinstance(self.db, list) else [self.db]
            worm_requested = "WORM" in [d for d in db_entries if isinstance(d, str)]

            if isinstance(self.db, list):
                # the main thermodynamic database can be a list of WORM-style
                # CSV databases that get concatenated into one database
                if len(self.db) == 0:
                    self.err_handler.raise_exception("An empty list was given "
                            "for 'db'. At least one thermodynamic database is "
                            "required.")

                db_list = self._resolve_db_list(self.db)

                if len(db_list) > 1:
                    if self.verbose > 0:
                        print("Loading and concatenating {} thermodynamic databases...".format(len(db_list)))
                    self.db = db_list
                else:
                    # a list of one database is just a single database
                    if self.verbose > 0:
                        if worm_requested:
                            print("Loading Water-Organic-Rock-Microbe (WORM) thermodynamic databases...")
                        else:
                            print("Loading a user-supplied thermodynamic database...")
                    self.db = db_list[0]

                self._set_active_db(db=self.db, download_csv_files=download_csv_files)

            elif isinstance(self.db, str):
                if worm_requested:
                    if self.verbose > 0:
                        print("Loading Water-Organic-Rock-Microbe (WORM) thermodynamic databases...")
                    self.db = WORM_DB_BASE_URL+"/wrm_data_latest.csv"
                else:
                    if self.verbose > 0:
                        print("Loading a user-supplied thermodynamic database...")
                self._set_active_db(db=self.db, download_csv_files=download_csv_files)

            elif isinstance(self.db, pd.DataFrame):
                if self.verbose > 0:
                    print("Loading a user-supplied Pandas DataFrame thermodynamic database...")
                self._set_active_db(db=self.db, download_csv_files=download_csv_files)

            else:
                self.err_handler.raise_exception("The thermodynamic database "
                        "'db' must be a string, a Pandas DataFrame, or a list "
                        "of WORM-style CSV databases to concatenate, not a "
                        "{}.".format(type(self.db).__name__))

            if worm_requested:
                if self.elements is None:
                    self._load_elements(WORM_DB_BASE_URL+"/elements.csv", download_csv_files=download_csv_files)
                if solid_solutions is None:
                    self._load_solid_solutions(WORM_DB_BASE_URL+"/solid_solutions.csv", download_csv_files=download_csv_files)
                if logK is None:
                    self._load_logK(WORM_DB_BASE_URL+"/wrm_data_logK.csv", download_csv_files=download_csv_files)
                if logK_S is None:
                    self._load_logK_S(WORM_DB_BASE_URL+"/wrm_data_logK_S.csv", download_csv_files=download_csv_files)
                if download_csv_files:
                    self._download_txt_file(WORM_DB_BASE_URL+"/speciation_groups_WORM.txt")

            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                self._validate_thermodynamic_database()

            # elements must be loaded if thermo_db_type is a CSV
            if self.elements is not None:
                self._load_elements(self.elements, download_csv_files=download_csv_files)
            if not self.element_active and self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                self._load_elements(WORM_DB_BASE_URL+"/elements.csv", download_csv_files=download_csv_files)

            if solid_solutions is not None:
                self._load_solid_solutions(solid_solutions, download_csv_files=download_csv_files)

            if logK is not None:
                self._load_logK(logK, download_csv_files=download_csv_files)

            # must be loaded after the logK database
            if logK_S is not None:
                self._load_logK_S(logK_S, download_csv_files=download_csv_files)

            if pitzer is not None:
                self._load_pitzer(pitzer, download_csv_files=download_csv_files)

            # which activity model does the active database call for?
            self.db_activity_model = self._detect_db_activity_model()
            if self.verbose > 0 and self.db_activity_model == "pitzer":
                print(self.thermo_db_filename, "is formatted for Pitzer's equations. "
                      "The 'pitzer' activity model will be used for speciation.")

            if self.logK_active:
                self._apply_logK_overrides()
                self.thermo_db = pd.concat([self.thermo_db, self.logK_db], ignore_index=True)

            if mineral_db != "WORM" and self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                # load a different mineral database
                pass
            
            # process dissociation reactions
            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                self._suppress_redox_and_generate_dissrxns(
                    suppress_redox=suppress_redox,
                    exceed_Ttr=exceed_Ttr)
            elif len(suppress_redox) > 0 and self.verbose > 0:
                print("Warning: redox suppression option is not recognized if a data0 or data1 database is used.")
                
            if self.logK_active:
                self.logK_db = self.thermo_db[~self.thermo_db["logK1"].isnull()]
                self.thermo_db = self.thermo_db[self.thermo_db["logK1"].isnull()]
                
            # generate input file template
            # (after species have been excluded)
            if input_template != "none":
                if input_template == 'strict':
                    template_names = list(self.thermo_db[self.thermo_db["tag"]=="basis"]["name"])
                elif input_template == 'basis':
                    template_names = list(self.thermo_db[self.thermo_db["tag"].isin(["basis", "aux"])]["name"])
                elif input_template == 'all':
                    template_names = list(self.thermo_db[self.thermo_db["state"]=="aq"]["name"])

                template_names = sorted(template_names)
                input_template = pd.DataFrame({"Sample":["id"], "H+":["pH"], "Temperature":["degC"], "logfO2":["logfO2"]})
                input_template_2 = pd.DataFrame({name:["Molality"] for name in template_names})
                input_template = pd.concat([input_template, input_template_2], axis=1)

                input_template.to_csv("sample_input_template.csv", index=False)


        def _validate_thermodynamic_database(self):

            # check that the df has all required columns
            req_cols = list(WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT.keys())
            missing_cols = [col for col in req_cols if col not in list(self.thermo_db.columns)]
            if len(missing_cols) > 0:
                self.err_handler.raise_exception("The following required columns are not present "
                        "in the loaded thermodynamic database: " + str(missing_cols))
            self.thermo_db = self.thermo_db.astype(WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT)

            # check that all aq and gas species have unique names
            dupe_list_aq = _get_duplicates(self.thermo_db[self.thermo_db["state"] == "aq"]["name"])
            dupe_list_gas = _get_duplicates(self.thermo_db[self.thermo_db["state"] == "gas"]["name"])
            dupe_list_liq = _get_duplicates(self.thermo_db[self.thermo_db["state"] == "liq"]["name"])
            if len(dupe_list_aq) > 0:
                self.err_handler.raise_exception("The loaded thermodynamic database "
                        "contains duplicate entries for aqueous species: " + str(dupe_list_aq))
            elif len(dupe_list_gas) > 0:
                self.err_handler.raise_exception("The loaded thermodynamic database "
                        "contains duplicate entries for gaseous species: " + str(dupe_list_gas))
            elif len(dupe_list_liq) > 0:
                self.err_handler.raise_exception("The loaded thermodynamic database "
                        "contains duplicate entries for liquid species: " + str(dupe_list_liq))

        
        def _reject_species(self, name, reason):
            
            dbs_to_search = ["csv_db", "logK_db", "logK_S_db"]
            
            db_filename = None
            for db_name in dbs_to_search:
                if isinstance(self.__getattribute__(db_name), pd.DataFrame):
                    if name in list(self.__getattribute__(db_name)["name"]):
                        idx_list = [i for i,n in enumerate(self.__getattribute__(db_name)["name"]) if n==name]
                        state_list = [self.__getattribute__(db_name)["state"].iloc[idx] for i,idx in enumerate(idx_list)]
                        for i,idx in enumerate(idx_list):
                            if name not in list(self.df_rejected_species["name"]) or i not in list(self.df_rejected_species.loc[self.df_rejected_species["name"]==name, "database index"]):
                                d = pd.DataFrame({'database name': [self.__getattribute__(db_name+"_filename")], 'database index': [int(idx)], 'name': [name], 'state': [state_list[i]], 'reason for rejection': [reason]})
                                self.df_rejected_species = pd.concat([self.df_rejected_species, d], ignore_index=True)
                        break
                    
                    
        def _remove_missing_G_species(self):
            # remove species that are missing a gibbs free energy value.
            # handle minerals first. Reject any that have missing G in any polymorph.
            # Berman minerals are be allowed through despite having missing Gibbs.
            mineral_name_reject = list(set(self.csv_db[(self.csv_db["G"].isnull()) & (self.csv_db['state'].str.contains('cr')) & (self.csv_db['model'].str.contains('CGL'))]["name"]))
            
            idx = list(self.csv_db[self.csv_db["name"].isin(mineral_name_reject)].index)
            names = self.csv_db["name"].loc[idx]

            for name in names:
                self._reject_species(name=name, reason="missing Gibbs free energy for at least one polymorph")
            
            self.csv_db = self.csv_db[~self.csv_db["name"].isin(mineral_name_reject)]
    
            # TODO: other states besides minerals
                
        def _resolve_db_list(self, db_list):
            """
            Check a list of thermodynamic databases that are meant to be
            concatenated and expand every "WORM" entry into the URL of the WORM
            thermodynamic database CSV. Only WORM-style CSV databases can be
            concatenated, so each entry must be "WORM", the filepath of a CSV
            file, the URL of a CSV file, or a Pandas DataFrame.
            """

            resolved = []
            for entry in db_list:
                if isinstance(entry, pd.DataFrame):
                    resolved.append(entry)
                elif isinstance(entry, str) and entry == "WORM":
                    resolved.append(WORM_DB_BASE_URL+"/wrm_data_latest.csv")
                elif isinstance(entry, str) and entry[-4:].lower() == ".csv":
                    resolved.append(entry)
                else:
                    self.err_handler.raise_exception("The thermodynamic "
                        "database '{}'".format(str(entry))+" cannot be "
                        "concatenated with other databases because it is not a "
                        "WORM-style CSV database. Each entry in a list given to "
                        "'db' must be:"
                        "\n - 'WORM'"
                        "\n - the name of a CSV file in your working directory. e.g., 'wrm_data.csv'"
                        "\n - a URL directing to a CSV file. e.g.,"
                        "\n\t 'https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_latest.csv'"
                        "\n - a Pandas DataFrame containing a WORM-style thermodynamic database"
                        "\nThree letter codes (e.g., 'wrm'), data0 files, and "
                        "data1 files cannot be concatenated because they are "
                        "not CSV databases.")

            return resolved


        def _read_csv_dbs(self, db, description, download_csv_files=False,
                          key="name", key_label=str, validate=None):
            """
            Read one or more WORM-style CSV databases and concatenate them into
            a single Pandas DataFrame. `db` can be the filepath of a CSV file,
            the URL of a CSV file, a Pandas DataFrame, or a list of any
            combination of these.

            Entries that are defined in more than one database are taken from
            the last database in the list that defines them; matching entries
            in earlier databases are dropped. Entries are matched with the
            column named by `key`, or, if `key` is a function, with the keys
            that `key` returns for each row of a database. `key_label` formats
            a key for display in messages.

            `validate` is an optional function that is applied to each database
            as it is read, e.g., to normalize columns before entries are
            matched.

            Returns a (filename, database, source) tuple describing all the
            databases that were read.
            """

            db_list = db if isinstance(db, list) else [db]

            if len(db_list) == 0:
                self.err_handler.raise_exception("An empty list was given for "
                        "the "+description+". At least one database is required.")

            filenames = []
            dfs = []
            sources = []
            for entry in db_list:
                if isinstance(entry, pd.DataFrame):
                    filenames.append("User-supplied Pandas DataFrame")
                    dfs.append(copy.deepcopy(entry))
                    sources.append("user-supplied")

                elif _is_url(entry):
                    # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK.csv"
                    filename, df = self.__df_from_url(entry, download_csv_files=download_csv_files)
                    filenames.append(filename)
                    dfs.append(df)
                    sources.append("URL")

                elif isinstance(entry, str):
                    # e.g., "wrm_data_logK.csv"
                    if os.path.isfile(entry):
                        filenames.append(entry)
                        dfs.append(pd.read_csv(entry))
                        sources.append("file")
                    else:
                        self.err_handler.raise_exception("Could not locate the "
                                "CSV file '"+entry+"' given for the "+description+".")

                else:
                    self.err_handler.raise_exception("The "+description+" must "
                            "be the name of a CSV file, the URL of a CSV file, "
                            "a Pandas DataFrame, or a list of any combination "
                            "of these. '{}'".format(str(entry))+" is not a "
                            "valid option.")

            if validate is not None:
                dfs = [validate(df, filenames[i]) for i,df in enumerate(dfs)]

            if len(dfs) == 1:
                return filenames[0], dfs[0], sources[0]

            dfs = self._drop_overridden_entries(dfs, filenames, key=key,
                                                key_label=key_label,
                                                description=description)

            filename = " + ".join(filenames)
            source = sources[0] if len(set(sources)) == 1 else "multiple sources"

            if self.verbose > 0:
                print("Concatenated {}".format(len(dfs)), description+"s:", filename)

            return filename, pd.concat(dfs, ignore_index=True), source


        def _drop_overridden_entries(self, dfs, filenames, key, description,
                                     key_label=str):
            """
            Drop entries of databases that are being concatenated if they are
            also defined in a database later in the list. All rows that share a
            key are dropped together so that, for instance, every polymorph of
            an overridden mineral is replaced.
            """

            if callable(key):
                keys = [list(key(df)) for df in dfs]
            else:
                keys = [list(df[key]) if key in df.columns else [] for df in dfs]

            out = [None]*len(dfs)
            later_keys = set()
            for i in reversed(range(len(dfs))):
                overridden = [k for k in dict.fromkeys(keys[i]) if k in later_keys]

                if len(overridden) > 0:
                    if self.verbose > 0:
                        labels = [key_label(k) for k in overridden]
                        if len(labels) <= 20:
                            print("Entries for", "; ".join(labels), "in", filenames[i],
                                  "are replaced by entries in a", description,
                                  "that comes later in the list.")
                        else:
                            print(len(labels), "entries in", filenames[i], "are replaced "
                                  "by entries in a", description, "that comes later in "
                                  "the list, including:", "; ".join(labels[:10]) + "; ...")
                    out[i] = dfs[i][[k not in later_keys for k in keys[i]]]
                else:
                    out[i] = dfs[i]

                later_keys.update(keys[i])

            return out


        def _set_active_db(self, db=None, download_csv_files=False):
            """
            Set the main active thermodynamic database to a Pandas DataFrame,
            CSV file, data0 file, a data1 file on the server, a local file,
            or from a URL address.
            """

            if isinstance(db, list):
                # a list of WORM-style CSV databases to concatenate
                self._load_csv(db, download_csv_files=download_csv_files)

                self.thermo_db = self.csv_db
                self.thermo_db_filename = self.csv_db_filename
                self.thermo_db_type = "CSV"
                self.thermo_db_source = self.csv_db_source
                self.dynamic_db = True
                self.custom_data0 = False
                self.data0_lettercode = None

            elif isinstance(db, pd.DataFrame):
                self.thermo_db = copy.deepcopy(self.db)
                self.db = "custom"
                self.thermo_db_type = "Pandas DataFrame"
                self.thermo_db_source = "user-supplied"
                self.thermo_db_filename = "User-supplied Pandas DataFrame thermodynamic database"
                self.dynamic_db = True
                self.custom_data0 = False
                self.data0_lettercode = None

            elif isinstance(db, str):
                if len(db) == 3:
                    # e.g., "wrm"

                    self.data0_lettercode = db
                    self.dynamic_db = False
                    
                    # search for a data1 file in the eq36da directory
                    data1_path = os.path.join(self.eq36da, f"data1.{db}")
                    if os.path.exists(data1_path) and os.path.isfile(data1_path):
                        self.thermo_db = None
                        self.thermo_db_type = "data1"
                        self.thermo_db_source = "file"
                        self.thermo_db_filename = "data1."+db

                        # store contents of data1 file in AqEquil object
                        with open(data1_path, mode='rb') as data1_file:
                            self.data1["all_samples"] = data1_file.read()
    
                    elif os.path.exists("data0." + db) and os.path.isfile("data0." + db):
                        
                        if self.verbose > 0:
                            print("data1." + db + " was not found in the EQ36DA directory "
                                  "but a data0."+db+" was found in the current working "
                                  "directory. Using it...")
    
                        self._load_data0("data0." + db, source="file")
    
                        self.thermo_db = self.data0_db
                        self.thermo_db_filename = self.data0_db_filename
                        self.thermo_db_type = "data0"
                        self.thermo_db_source = "file"
                        self.custom_data0 = True
                        self.data0_lettercode = db[-3:].lower()
                        self.eq36da = os.path.join(os.getcwd(), "eqpt_files")

                    elif os.path.exists("data1." + db) and os.path.isfile("data1." + db):

                        if self.verbose > 0:
                            print("data1." + db + " was not found in the EQ36DA directory "
                                  "but a data1."+db+" was found in the current working "
                                  "directory. Using it...")

                        self.custom_data0 = True
                        self.thermo_db = None
                        self.eq36da = os.path.join(os.getcwd(), "eqpt_files")
    
                        # search for a data1 locally
    
                        # store contents of data1 file in AqEquil object
                        with open("data1." + db, mode='rb') as data1_file:
                            self.data1["all_samples"] = data1_file.read()
                            self.thermo_db_type = "data1"
                            self.thermo_db_source = "file"
                            self.thermo_db_filename = "data1."+db
    
                    else:
                        msg = ("Could not locate a 'data1."+db+"' file in the EQ36DA "
                              "directory, nor a 'data0."+db+"' or 'data1."+db+"' file in "
                              "the current working directory.")
                        self.err_handler.raise_exception(msg)
    
                elif "data0." in db[-9:].lower() and db[-4:].lower() != ".csv" and (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
                    # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/data0.wrm"
    
                    self._load_data0(db, source="URL")
    
                    self.thermo_db = self.data0_db
                    self.thermo_db_filename = self.data0_db_filename
                    self.thermo_db_type = "data0"
                    self.thermo_db_source = "URL"
                    self.custom_data0 = True
                    self.data0_lettercode = db[-3:]
                    self.dynamic_db = False
    
                elif db[0:-4].lower() == "data0" and not (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
                    # e.g., "data0.wrm"
    
                    self._load_data0(db, source="file")
    
                    self.thermo_db = self.data0_db
                    self.thermo_db_filename = self.data0_db_filename
                    self.thermo_db_type = "data0"
                    self.thermo_db_source = "file"
                    self.custom_data0 = True
                    self.data0_lettercode = db[-3:].lower()
                    self.dynamic_db = False
    
                elif db[-4:].lower() == ".csv" and not (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
                    # e.g., "wrm_data.csv"
    
                    self._load_csv(db)

                    self.thermo_db = self.csv_db
                    self.thermo_db_filename = self.csv_db_filename
                    self.thermo_db_type = "CSV"
                    self.thermo_db_source = "file"
                    self.dynamic_db = True
                    self.custom_data0 = False
                    self.data0_lettercode = None
    
                elif db[-4:].lower() == ".csv" and (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
                    # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv
                    
                    
                    self._load_csv(db, download_csv_files=download_csv_files)
    
                    self.thermo_db = self.csv_db
                    self.thermo_db_filename = self.csv_db_filename
                    self.thermo_db_type = "CSV"
                    self.thermo_db_source = "URL"
                    self.dynamic_db = True
                    self.custom_data0 = False
                    self.data0_lettercode = None
                    
                else:
                    self.err_handler.raise_exception("Unrecognized thermodynamic "
                        "database '{}'".format(db)+" specified for db. A database can specified as:"
                        "\n - a three letter code designating a data0 file. e.g., db='wrm'"
                        "\n - a data0 file in your working directory. e.g., db='data0.wrm'"
                        "\n - a csv file in your working directory. e.g., db='wrm_data.csv'"
                        "\n - a URL directing to a data0 file. e.g.,"
                        "\n\t db='https://raw.githubusercontent.com/worm-portal/WORM-db/master/data0.wrm'"
                        "\n\t (note the data0 file in the URL must have 'data0.' followed by a three letter code)"
                        "\n - a URL directing to a valid csv file. e.g.,"
                        "\n\t db='https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_latest.csv'")

                self.db = db
            
            if self.verbose > 0:
                print(self.thermo_db_filename, "is now set as the active thermodynamic database.")

                if self.thermo_db_filename in ['data0.wrm', 'data1.wrm']:
                    print("This database is meant for rapid calculations between 0 and 350 °C at water saturation pressure.")
                elif self.thermo_db_filename == "wrm_data.csv":
                    print("This database is meant for calculations between 0 and 1000 °C and up to 5 kb pressure.")

            # a data0 or data1 file does not contain everything that functions
            # like plot_reaction_paths() need, so load its companion CSV too
            if self.thermo_db_type in ["data0", "data1"] and not isinstance(self.csv_db, pd.DataFrame):
                self._load_companion_csv(download_csv_files=download_csv_files)

            


        def __df_from_url(self, url, download_csv_files=False):
            """
            Get a filename and dataframe from a URL pointing to a CSV file.
            Falls back to bundled databases if the URL cannot be reached.
            """
            filename = url.split("/")[-1]

            try:
                # Download from URL and decode as UTF-8 text.
                if download_csv_files and self.verbose > 0:
                    print(f"Downloading {filename}...")

                with urlopen(url) as webpage:
                    content = webpage.read().decode('utf-8')

                if download_csv_files:
                    with open(filename, 'w', encoding='utf-8') as output:
                        output.write(content)
                    if self.verbose > 0:
                        size_kb = len(content.encode('utf-8')) / 1024
                        print(f"  [OK] {filename} ({size_kb:.1f} KB)")

                return filename, pd.read_csv(StringIO(content), sep=",")

            except Exception as e:
                # Fallback to bundled databases if URL cannot be reached
                bundled_path = _get_bundled_data_path()
                if bundled_path is not None:
                    bundled_file = os.path.join(bundled_path, filename)
                    if os.path.isfile(bundled_file):
                        if self.verbose > 0:
                            print(f"  [FALLBACK] Cannot reach URL ({e}), using bundled database: {filename}")

                        # If download_csv_files is True, copy bundled file to cwd
                        if download_csv_files:
                            shutil.copy(bundled_file, filename)
                            if self.verbose > 0:
                                size_kb = os.path.getsize(bundled_file) / 1024
                                print(f"  [OK] Copied bundled {filename} ({size_kb:.1f} KB)")

                        return filename, pd.read_csv(bundled_file, sep=",")

                # If no bundled fallback available, raise the original error
                self.err_handler.raise_exception("The webpage "+str(url)+" cannot"
                        " be reached at this time and no bundled database is available.")


        def __str_from_url(self, url):
            """
            Get a filename and contents from a URL pointing to a txt file.
            """

            filename = url.split("/")[-1]

            if self.verbose > 0:
                print(f"Downloading {filename}...")

            try:
                # Download from URL and decode as UTF-8 text.
                with urlopen(url) as webpage:
                    txt_content = webpage.read().decode('utf-8')
            except:
                self.err_handler.raise_exception("The webpage "+str(url)+" cannot"
                        " be reached at this time.")

            with open(filename, 'w', encoding='utf-8') as output:
                output.write(txt_content)

            if self.verbose > 0:
                size_kb = len(txt_content.encode('utf-8')) / 1024
                print(f"  [OK] {filename} ({size_kb:.1f} KB)")

            return filename, txt_content


        def _download_txt_file(self, url):
            """
            Download a text file from a URL to the current working directory.
            """
            filename = url.split("/")[-1]

            if self.verbose > 0:
                print(f"Downloading {filename}...")

            try:
                with urlopen(url) as webpage:
                    content = webpage.read().decode('utf-8')

                with open(filename, 'w', encoding='utf-8') as output:
                    output.write(content)

                if self.verbose > 0:
                    size_kb = len(content.encode('utf-8')) / 1024
                    print(f"  [OK] {filename} ({size_kb:.1f} KB)")

            except Exception as e:
                if self.verbose > 0:
                    print(f"  [ERROR] Failed to download {filename}: {e}")


        def _load_elements(self, db, download_csv_files=False):
            """
            Load one or more element database CSV files from files, URLs, or
            Pandas DataFrames. Multiple databases are concatenated.
            """

            self.element_db_filename, self.element_db, self.element_db_source = \
                self._read_csv_dbs(db, description="element database",
                                   download_csv_files=download_csv_files,
                                   key="element")

            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                if self.verbose > 0:
                    print("Element database", self.element_db_filename, "is active.")
                self.element_active = True

                # Update CHNOSZ's global thermo object with the element database
                # This is needed for functions like entropy() and mass() to work with custom elements
                if self.element_db is not None:
                    pychnosz.thermo(element=self.element_db)
            else:
                if self.verbose > 0:
                    print("Element database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")


        def _load_solid_solutions(self, db, download_csv_files=False):
            """
            Load one or more solid solution database CSV files from files, URLs,
            or Pandas DataFrames. Multiple databases are concatenated.
            """

            self.solid_solution_db_filename, self.solid_solution_db, self.solid_solution_db_source = \
                self._read_csv_dbs(db, description="solid solution database",
                                   download_csv_files=download_csv_files,
                                   key="name")

            if self.thermo_db_type == "CSV":
                if self.verbose > 0:
                    print("Solid solution database", self.solid_solution_db_filename, "is active.")
                self.solid_solutions_active = True
            else:
                if self.verbose > 0:
                    print("Solid solution database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")


        def _load_logK(self, db, download_csv_files=False):
            """
            Load one or more logK database CSV files from files, URLs, or Pandas
            DataFrames. Multiple databases are concatenated.
            """

            self.logK_db_filename, self.logK_db, self.logK_db_source = \
                self._read_csv_dbs(db, description="logK database",
                                   download_csv_files=download_csv_files,
                                   key="name")

            if self.thermo_db_type == "CSV":
                if self.verbose > 0:
                    print("LogK database", self.logK_db_filename, "is active.")
                self.logK_active = True
            else:
                if self.verbose > 0:
                    print("LogK database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")

            self.logK_db = self._exclude_category(df=self.logK_db, df_name=self.logK_db_filename)


        def _load_pitzer(self, db, download_csv_files=False):
            """
            Load one or more Pitzer parameter database CSV files from files,
            URLs, or Pandas DataFrames. Multiple databases are concatenated.
            See `aqequil.pitzer` for the format.
            """

            def validate(df, filename):
                try:
                    return validate_pitzer_db(df, filename=filename)
                except PitzerDatabaseError as e:
                    self.err_handler.raise_exception(str(e))

            # each database is validated before they are concatenated so that
            # species names and parameter names are comparable between them
            self.pitzer_db_filename, self.pitzer_db, self.pitzer_db_source = \
                self._read_csv_dbs(db, description="Pitzer parameter database",
                                   download_csv_files=download_csv_files,
                                   key=pitzer_param_keys,
                                   key_label=format_pitzer_param_key,
                                   validate=validate)

            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                if self.verbose > 0:
                    print("Pitzer parameter database", self.pitzer_db_filename, "is active.")
                self.pitzer_active = True
            else:
                if self.verbose > 0:
                    print("Pitzer parameter database is not active because the active "
                          "thermodynamic database is a", self.thermo_db_type, "and not a CSV. "
                          "A", self.thermo_db_type, "file already determines which activity "
                          "model is used.")


        def _detect_db_activity_model(self):
            """
            Determine which aqueous activity coefficient model the active
            thermodynamic database is formatted for. Returns "pitzer" for
            Pitzer data0/data1 files, "sedh" for B-dot/Davies data0/data1
            files, and None for CSV databases, which can be used to generate
            data0 files for either kind of model.
            """
            if self.thermo_db_type == "data0" and isinstance(self.thermo_db, str):
                return "pitzer" if data0_is_pitzer(self.thermo_db) else "sedh"
            elif self.thermo_db_type == "data1":
                data1_bytes = self.data1.get("all_samples", b"")
                return "pitzer" if data1_is_pitzer(data1_bytes) else "sedh"
            return None


        def _load_logK_S(self, db, download_csv_files=False):
            """
            Load one or more logK_S database CSV files from files, URLs, or
            Pandas DataFrames. Multiple databases are concatenated.
            """

            self.logK_S_db_filename, self.logK_S_db, self.logK_S_db_source = \
                self._read_csv_dbs(db, description="logK_S database",
                                   download_csv_files=download_csv_files,
                                   key="name")

            if self.logK_active and self.element_active:
                if self.verbose > 0:
                    print("LogK_S database", self.logK_S_db_filename, "is active.")
                self.logK_S_active = True
            else:
                if self.verbose > 0:
                    print("LogK_S database is not active because there is no active logK database.")

            self.logK_S_db = self._exclude_category(df=self.logK_S_db, df_name=self.logK_S_db_filename)

            if self.logK_S_active:

                for i,sp in enumerate(self.logK_S_db["name"]):

                    i = self.logK_S_db.index[i]
                    
                    logK_25C = float(self.logK_S_db["logK_25"][i])

                    if "logK_25_IS" in self.logK_S_db.columns:
                        IS_ref = float(self.logK_S_db["logK_25_IS"][i])
                    else:
                        IS_ref = 0
                    
                    T_list = self.logK_S_db["T_vals"][i].split(" ")
                    T_list = [float(T) for T in T_list]

                    Delta_S = float(self.logK_S_db["DeltaS"][i])


                    if IS_ref > 0:
                        # extrapolate to ionic strength 0
                        
                        # collect azero values for metal, ligand, and complex
                        cation_name = self.logK_S_db["cation_name"][i]
    
                        if self.thermo_db["name"].isin([cation_name]).any():
                            metal_azero = list(self.thermo_db[self.thermo_db["name"] == cation_name]["azero"])[0]
                            metal_charge = float(list(self.thermo_db[self.thermo_db["name"] == cation_name]["z.T"])[0])
    #                         elif:
    #                             # todo: get metal azero and charge from other databases, e.g., logK db
    #                             pass
                        elif cation_name == "H+":
                            metal_azero = 4
                            metal_charge = 1
                        else:
                            # todo: throw error
                            pass
                        
                        ligand_name = self.logK_S_db["ligand_name"][i]
                        ligand_basis = self.logK_S_db["ligand_basis"][i]
    
    
                        if isinstance(self.logK_S_db["ligand_charge"][i], float):
                            ligand_charge = float(self.logK_S_db["ligand_charge"][i])
                        elif self.thermo_db["name"].isin([ligand_name]).any():
                            ligand_charge = float(list(self.thermo_db[self.thermo_db["name"] == ligand_name]["z.T"])[0])
                        elif self.thermo_db["name"].isin([ligand_basis]).any():
                            ligand_charge = float(list(self.thermo_db[self.thermo_db["name"] == ligand_basis]["z.T"])[0])
                        else:
                            formula_dict = parse_formula(self.logK_S_db["ligand_formula"])
                            if "+" in list(formula_dict.keys()):
                                ligand_charge = float(formula_dict["+"])
                            elif "-" in list(formula_dict.keys()):
                                ligand_charge = float(formula_dict["-"])
                            else:
                                ligand_charge = 0
                            # self.err_handler.raise_exception("The ligand '"+ligand_name+"' was not found in the "
                            #         "currently loaded databases. Is this ligand present in the database? Was this "
                            #         "ligand excluded (e.g., exclude_organics=True)?")
    
                        if isinstance(self.logK_S_db["ligand_azero"][i], float):
                            ligand_azero = float(self.logK_S_db["ligand_azero"][i])
                        elif self.thermo_db["name"].isin([ligand_name]).any():
                            ligand_azero = float(list(self.thermo_db[self.thermo_db["name"] == ligand_name]["azero"])[0])
                        elif self.thermo_db["name"].isin([ligand_basis]).any():
                            ligand_azero = float(list(self.thermo_db[self.thermo_db["name"] == ligand_basis]["azero"])[0])
                        else:
                            # todo: elif ligand_name in logK database names, get azero from there...
                            # or maybe this is not necessary if logK is merged with thermo_db at this point
    
                            ligand_azero = 4 # acceptable placeholder azero for charged species
    
                    
                        dissrxn = self.logK_S_db["dissrxn"][i].split(" ")
                        n_metal = float(dissrxn[dissrxn.index(cation_name)-1])
                        if ligand_name in dissrxn:
                            n_ligand = float(dissrxn[dissrxn.index(ligand_name)-1])
                        elif ligand_basis in dissrxn:
                            n_ligand = float(dissrxn[dissrxn.index(ligand_basis)-1])
                        else:
                            # todo: throw error
                            pass
                        n_complex = -float(dissrxn[dissrxn.index(sp)-1])

                        complex_charge = n_metal*metal_charge + n_ligand*ligand_charge
                        complex_azero = self.logK_S_db["azero"][i]
                        
                        A=0.5114
                        B=0.3288
                        Bdot=0.041
                        If = 0 # what ionic strength to extrapolate to
                        
                        ari=[metal_azero, ligand_azero]
                        api=[complex_azero]
                        vri=[n_metal, n_ligand]
                        vpi=[n_complex]
                        zri=[metal_charge, ligand_charge]
                        zpi=[complex_charge]
                        
                        def loggamma(vparam, zparam, aparam, I):
                            x=[v*((-1*A*z**2*I**0.5)/(1+a*B*I**0.5)+Bdot*I) for v,z,a in zip(vparam, zparam, aparam)]
                            return x
                        
                        def f(vparam, zparam, aparam, I):
                            return sum(loggamma(vparam, zparam, aparam, I))
                        
                        logK_25C = -(-logK_25C+(f(vpi,zpi,api,IS_ref)-f(vri,zri,ari,IS_ref))-(f(vpi,zpi,api,If)-f(vri,zri,ari,If)))
                        
                    logK_list = self._est_logK_S(T_list, logK_25C, Delta_S)
                    
                    
                    if isinstance(self.logK_S_db["ligand_element"][i], str):
                        # modify element database with pseudoelements
                        pseudoelement = self.logK_S_db["ligand_element"][i]
                        if pseudoelement not in list(self.element_db["element"]):
                            e_df = pd.DataFrame(
                                {'element':[self.logK_S_db["ligand_element"][i]],
                                 'state':[self.logK_S_db["state"][i]],
                                 'source':[self.logK_S_db["ligand_name"][i]],
                                 'mass':[self.logK_S_db["ligand_mass"][i]],
                                 's':[self.logK_S_db["ligand_entropy"][i]],
                                 'n':[self.logK_S_db["ligand_n"][i]],
                                })

                            self.element_db = pd.concat([self.element_db, e_df], ignore_index=True)

                            # Update CHNOSZ's thermo object with the modified element database
                            pychnosz.thermo(element=self.element_db)

                    if isinstance(self.logK_S_db["ligand_basis"][i], str):
                        # add a basis species representing the pseudoelement
                        basis = self.logK_S_db["ligand_basis"][i]

                        if basis not in list(self.thermo_db["name"]):
                            b_df = pd.DataFrame(
                                {'name':[self.logK_S_db["ligand_basis"][i]],
                                 'abbrv':[""],
                                 'formula':[self.logK_S_db["ligand_formula"][i]],
                                 'state':[self.logK_S_db["state"][i]],
                                 'ref1':[self.logK_S_db["ref1"][i]],
                                 'ref2':[self.logK_S_db["ref2"][i]],
                                 'date':[self.logK_S_db["date"][i]],
                                 'E_units':["cal"],
                                 'G':[0], 'H':[0], 'S':[0],
                                 'Cp':[0], 'V':[0], 'a1.a':[0],
                                 'a2.b':[0], 'a3.c':[0], 'a4.d':[0],
                                 'c1.e':[0], 'c2.f':[0],
                                 'omega.lambda':[0],
                                 'z.T':[self.logK_S_db["ligand_charge"][i]],
                                 'azero':[self.logK_S_db["ligand_azero"][i]],
                                 'neutral_ion_type':[0],
                                 'dissrxn':[''],
                                 'tag':['basis'],
                                 'formula_ox':[self.logK_S_db["ligand_formula"][i]],
                                 'category_1':[self.logK_S_db["category_1"][i]],
                                })
                            
                            self.thermo_db = pd.concat([self.thermo_db, b_df], ignore_index=True)

                    if self.logK_S_db["name"][i] not in self.thermo_db["name"] and self.logK_S_db["name"][i] not in self.logK_db["name"]:

                        s_df = pd.DataFrame(
                                {'name':[self.logK_S_db["name"][i]],
                                 'abbrv':[""],
                                 'formula':[self.logK_S_db["formula"][i]],
                                 'state':[self.logK_S_db["state"][i]],
                                 'ref1':[self.logK_S_db["ref1"][i]],
                                 'ref2':[self.logK_S_db["ref2"][i]],
                                 'date':[self.logK_S_db["date"][i]],
                                 'logK1':[np.nan],'logK2':[np.nan],'logK3':[np.nan],'logK4':[np.nan],'logK5':[np.nan],'logK6':[np.nan],'logK7':[np.nan],'logK8':[np.nan],
                                 'T1':[np.nan],'T2':[np.nan],'T3':[np.nan],'T4':[np.nan],'T5':[np.nan],'T6':[np.nan],'T7':[np.nan],'T8':[np.nan],
                                 'P1':[np.nan],'P2':[np.nan],'P3':[np.nan],'P4':[np.nan],'P5':[np.nan],'P6':[np.nan],'P7':[np.nan],'P8':[np.nan],
                                 'azero':[self.logK_S_db["azero"][i]],
                                 'dissrxn':[self.logK_S_db["dissrxn"][i]],
                                 'tag':[''],
                                 'formula_ox':[self.logK_S_db["formula_ox"][i]],
                                 'category_1':[self.logK_S_db["category_1"][i]],
                                })
                        self.logK_db = pd.concat([self.logK_db, s_df], ignore_index=True)

                        for ti in range(0, len(T_list)):
                            if ti+1 > 8:
                                self.err_handler.raise_exception("Species ", sp, "in",
                                    self.logK_S_db_filename, "may only have up to",
                                    "eight temperature values in column T_vals")

                            self.logK_db.loc[self.logK_db.index[-1], "logK"+str(ti+1)] = logK_list[ti]
                            self.logK_db.loc[self.logK_db.index[-1], "T"+str(ti+1)] = T_list[ti]
                            self.logK_db.loc[self.logK_db.index[-1], "P"+str(ti+1)] = 'psat'

        def _est_logK_S(self, T_list, logK_25C, Delta_S):

            R = 8.31446261815324/4.184 # cal/mol/K
            T = 298.15

            logK_list = []
            # based on Eq3 of Prapaipong & Shock 2001:
            for T_C in T_list:
                #delta_logK = (Delta_S/(2.303*R*298.15))*(T_C-25)
                delta_logK = (Delta_S/(2.303*R*(T_C+273.15)))*(T_C-25)
                logK_list.append(logK_25C + delta_logK)
            
            return logK_list


        def _load_data0(self, db, source="URL"):
            """
            Load a data0 file from a file or URL.
            """

            if source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/data0.wrm"
                self.data0_db_filename, self.data0_db = self.__str_from_url(db)
                self.data0_db_type = "data0"
                self.data0_db_source = "URL"

            elif source == "file":
                # e.g., "data0.wrm"
                if os.path.exists(db) and os.path.isfile(db):
                    with open(db) as data0_content:
                        self.data0_db = data0_content.read()
                        self.data0_db_type = "data0"
                        self.data0_db_source = "file"
                        self.data0_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the data0 file '"+db+"'")


        def _load_csv(self, db, download_csv_files=False):
            """
            Load one or more WORM-styled thermodynamic database CSVs from files,
            URLs, or Pandas DataFrames. Multiple databases are concatenated.
            """

            self.csv_db_filename, self.csv_db, self.csv_db_source = \
                self._read_csv_dbs(db, description="thermodynamic database",
                                   download_csv_files=download_csv_files,
                                   key="name")
            self.csv_db_type = "CSV"

            # Check that thermodynamic database input files exist and are formatted correctly.
            self._check_csv_db()

            self.csv_db = self.csv_db.astype(WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT)

            # Special handling for missing 'model' column
            if 'model' not in self.csv_db.columns:
                # Create model column with proper values
                # - aqueous species (state == 'aq') get 'HKF'
                # - non-aqueous species get 'CGL'
                self.csv_db['model'] = self.csv_db['state'].apply(lambda x: 'HKF' if x == 'aq' else 'CGL')

                # Issue a warning to inform the user
                warnings.warn(
                    "The 'model' column was not found in the thermodynamic database CSV. "
                    "Auto-generating 'model' column: 'HKF' for aqueous species (state='aq'), "
                    "'CGL' for all other species.",
                    UserWarning
                )

            self._remove_missing_G_species()

            self.csv_db = self._exclude_category(df=self.csv_db, df_name=self.csv_db_filename)


        def _load_companion_csv(self, download_csv_files=False):
            """
            Load the WORM-styled thermodynamic database CSV that accompanies the
            active data0 or data1 file.

            A data0 or data1 file is all that EQ3/6 needs, but other functions,
            like plot_reaction_paths(), require the additional information found
            in a CSV database. The CSV becomes the `csv_db` attribute and does
            not replace the data0 or data1 file used for EQ3/6 calculations.
            """

            # only some data0 and data1 files have a companion CSV
            csv_names = WORM_COMPANION_CSVS.get(str(self.data0_lettercode).lower())
            if csv_names == None:
                return

            # a CSV in the working directory is used before one is downloaded
            local_csvs = [f for f in csv_names if os.path.isfile(f)]

            try:
                if len(local_csvs) > 0:
                    self._load_csv(local_csvs[0])
                else:
                    self._load_csv(WORM_DB_BASE_URL+"/"+csv_names[0],
                                   download_csv_files=download_csv_files)

            except Exception as e:
                # EQ3/6 calculations can go ahead without the companion CSV
                self.csv_db = None
                self.csv_db_type = None
                self.csv_db_source = None
                self.csv_db_filename = None

                if self.verbose > 0:
                    print("The CSV database that accompanies", self.thermo_db_filename,
                          "could not be loaded:", str(e))
                    print("Continuing anyway, but functions that require a "
                          "CSV database, like plot_reaction_paths(), will not "
                          "work. Consider placing", csv_names[0], "in your "
                          "working directory.")
                return

            if self.verbose > 0:
                print(self.csv_db_filename, "is loaded alongside",
                      str(self.thermo_db_filename), "to support functions that "
                      "require a CSV database.")


        def _apply_logK_overrides(self):
            """
            Species that appear in both the main CSV thermodynamic database and
            the logK database are taken from the logK database: the main
            database entries (including any polymorph rows) are removed so
            that the data0 file does not contain duplicate species blocks.
            Strict basis species cannot be overridden and are dropped from
            the logK database instead. If an overridden species was an
            auxiliary basis species, the logK entry inherits the 'aux' tag so
            that reactions written in terms of it remain valid.
            """
            if not isinstance(self.thermo_db, pd.DataFrame) or not isinstance(self.logK_db, pd.DataFrame):
                return
            if "name" not in self.logK_db.columns or "name" not in self.thermo_db.columns:
                return

            main_names = set(self.thermo_db["name"])
            overlap = [n for n in dict.fromkeys(self.logK_db["name"]) if n in main_names]
            if len(overlap) == 0:
                return

            basis_names = set(self.thermo_db.loc[self.thermo_db["tag"] == "basis", "name"])
            keep_basis = [n for n in overlap if n in basis_names]
            replace = [n for n in overlap if n not in basis_names]

            if len(keep_basis) > 0:
                if self.verbose > 0:
                    print("The logK database", self.logK_db_filename, "contains entries for the "
                          "strict basis species", keep_basis, "which cannot be overridden. "
                          "Ignoring these logK entries.")
                self.logK_db = self.logK_db[~self.logK_db["name"].isin(keep_basis)]

            if len(replace) > 0:
                # carry over the auxiliary basis tag where the main database used it
                aux_names = set(self.thermo_db.loc[self.thermo_db["tag"] == "aux", "name"])
                if "tag" in self.logK_db.columns:
                    self.logK_db["tag"] = self.logK_db["tag"].astype(object)
                    for n in replace:
                        if n in aux_names:
                            mask = (self.logK_db["name"] == n) & (self.logK_db["tag"].isnull() | (self.logK_db["tag"].astype(str).str.strip() == "") | (self.logK_db["tag"].astype(str) == "nan"))
                            self.logK_db.loc[mask, "tag"] = "aux"

                self.thermo_db = self.thermo_db[~self.thermo_db["name"].isin(replace)]
                if self.verbose > 0:
                    if len(replace) <= 20:
                        print("Entries in the logK database", self.logK_db_filename, "replace "
                              "the following species in", self.thermo_db_filename, ":", replace)
                    else:
                        print("Entries in the logK database", self.logK_db_filename, "replace",
                              len(replace), "species in", self.thermo_db_filename,
                              "(including", ", ".join(replace[:10]) + ", ...)")


        def _exclude_category(self, df, df_name):
            """
            Exclude entries from a df based on values in columns.
            e.g., {"category_1":["organic_aq", "organic_cr"]}
            Strict basis species are never excluded.
            """

            if isinstance(self.exclude_organics_except, list):
                organics_to_exclude = []
                for i,name in enumerate(df["name"]):
                    if df["category_1"].iloc[i] in ["organic_aq", "organic_cr"] and name not in self.exclude_organics_except:
                        organics_to_exclude.append(name)

                if not isinstance(self.exclude_category.get("name"), list) or "name" not in list(self.exclude_category.keys()):
                    self.exclude_category["name"] = organics_to_exclude
                else:
                    self.exclude_category["name"] += organics_to_exclude
                    self.exclude_category["name"] = list(set(self.exclude_category["name"]))

            
            exclude_keys = list(self.exclude_category.keys())
            if len(exclude_keys) > 0:
                for key in exclude_keys:
                    if self.verbose > 0:
                        if len(self.exclude_category[key]) <= 10:
                            print("Excluding", str(self.exclude_category[key]), "from column '" + str(key) + "'in", df_name)
                        else:
                            print("Excluding", len(self.exclude_category[key]), "different chemical species from column '" + str(key) + "' in", df_name)
                    
                    if isinstance(self.exclude_category[key], list):

                        to_exclude = df[key].isin(self.exclude_category[key])
                        # strict basis species are the foundation of every
                        # data0 file and are never excluded by category
                        if "tag" in df.columns and key != "name":
                            is_basis = df["tag"] == "basis"
                            if (to_exclude & is_basis).any() and self.verbose > 1:
                                print("Keeping strict basis species", list(df.loc[to_exclude & is_basis, "name"]),
                                      "even though they match an excluded category.")
                            to_exclude = to_exclude & ~is_basis

                        idx = list(df[to_exclude].index)
                        names = df["name"].loc[idx]

                        for name in names:
                            self._reject_species(name=name, reason="excluded by user")

                        df = df[~to_exclude]
                        

                    elif isinstance(self.exclude_category[key], str):
                        
                        idx = list(df[df[key] != self.exclude_category[key]].index)
                        names = df["name"].loc[idx]
                        
                        for name in names:
                            self._reject_species(name=name, reason="excluded by user")
                        
                        df = df[~df[key] != self.exclude_category[key]]

                    else:
                        self.err_handler.raise_exception("The parameter exclude_category must either be a string or a list.")
            return df


        def _check_csv_db(self):
            """
            Check for problems in the thermodynamic database CSV.
            """

            thermo_df = self.csv_db

            # does this file have the proper headers?
            required_headers = ["name", "abbrv", "formula", "state",
                                "ref1", "ref2", "date", "E_units",
                                "G", "H", "S", "Cp", "V",
                                "a1.a", "a2.b", "a3.c", "a4.d", "c1.e", "c2.f",
                                "omega.lambda", "z.T",
                                "azero", "neutral_ion_type",
                                "dissrxn", "tag", "formula_ox"]

            # columns that are needed to assign data types are required, too
            required_headers += [c for c in WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT
                                 if c not in required_headers]

            missing_headers = []
            for header in required_headers:
                if header not in thermo_df.columns:
                    missing_headers.append(header)
            if len(missing_headers) > 0:
                msg = ("The thermodynamic database "
                       "{}".format(self.csv_db_filename)+" is missing one or "
                       "more required columns: "
                       "{}".format(", ".join(missing_headers))+". "
                       "Are these headers spelled correctly in the file?")
                self.err_handler.raise_exception(msg)

            # does O2(g), and O2 exist in the file?
            required_species = ["O2", "O2(g)"]
            missing_species = []
            for species in required_species:
                if species not in list(thermo_df["name"]):
                    missing_species.append(species)
            if len(missing_species) > 0:
                msg = ("The thermodynamic database file "
                       "is missing required species:"
                       "{}".format(missing_species)+". Default thermodynamic values"
                       " will be used.")
                warnings.warn(msg)

            return


        def _suppress_redox_and_generate_dissrxns(self,
                                                  suppress_redox,
                                                  exceed_Ttr=True):
            
            thermo_df = self.thermo_db

            # if elements are being redox-suppressed, exclude all species with a
            # formula containing one or more of the redox-suppressed elements if the
            # species does not have a formula_ox.
            # e.g., if "methionine" does not have a formula_ox, ensure it is excluded
            #       if sulfur is redox-suppressed.
            if len(suppress_redox) > 0:
                thermo_db_no_formula_ox = thermo_df[thermo_df["formula_ox"].isnull()]
                if thermo_db_no_formula_ox.shape[0] > 0:
                    sp_names_to_exclude = []
                    for i,sp in enumerate(thermo_db_no_formula_ox["name"]):
                        f = thermo_db_no_formula_ox.iloc[i, thermo_db_no_formula_ox.columns.get_loc("formula")]
                        f_elems = list(parse_formula(f).keys())
                        for elem in suppress_redox:
                            if elem in f_elems:
                                sp_names_to_exclude.append(sp)
                    self.exclude_category["name"] = sp_names_to_exclude
                    if self.verbose > 0 and len(sp_names_to_exclude) > 0:
                        print("Excluding the following chemical species because "
                              "they contain redox-suppressed elements but do not "
                              "have element oxidation states given in the "
                              "'formula_ox' column of the thermodynamic database: "
                              ""+str(sp_names_to_exclude))

            # Use Python-based dissociation reaction processing instead of R
            thermo_df = assign_worm_db_col_dtypes(thermo_df)

            # Load the thermodynamic database into pyCHNOSZ before processing dissociation reactions
            # This ensures that CHNOSZ has access to all species in the database
            if self.verbose > 0:
                print("Loading thermodynamic database into pyCHNOSZ...")

            # Replace the default OBIGT database with user-supplied database
            # Keep only fixed species from the original OBIGT database
            # Filter out species that return NaN (not in default OBIGT)
            fixed_indices = pychnosz.info(FIXED_SPECIES, messages=False)
            fixed_indices_clean = [idx for idx in fixed_indices if not pd.isna(idx)]
            pychnosz.thermo(OBIGT=pychnosz.thermo().OBIGT.loc[fixed_indices_clean])

            # Add all species from the user-supplied thermodynamic database
            pychnosz.add_OBIGT(thermo_df, messages=False)

            # Call the Python function for processing dissociation reactions
            # Note: redox suppression will still be handled by R for now
            # This Python version handles dissociation reaction checking and generation
            # Get list of already-rejected species (e.g., from exclude_category)
            already_rejected = list(self.df_rejected_species['name']) if len(self.df_rejected_species) > 0 else []
            # A species excluded from the main database (e.g., by exclude_category)
            # but supplied again by the logK database is present in thermo_df and
            # must not cascade rejections onto the species that depend on it.
            present_names = set(thermo_df["name"])
            already_rejected = [n for n in already_rejected if n not in present_names]

            self.out_list = process_dissrxns(
                thermo_df=thermo_df,
                water_model=self.water_model,
                exceed_Ttr=exceed_Ttr,
                suppress_redox=suppress_redox,
                exclude_category=self.exclude_category,
                element_df=self.element_db,
                fixed_species=FIXED_SPECIES,
                verbose=self.verbose,
                already_rejected_species=already_rejected
            )

            thermo_df = self.out_list['thermo_df']

            # Handle rejected species
            rejected_species = self.out_list['dissrxns'].get('rejected_species', [])

            if rejected_species and len(rejected_species) > 0:
                for sp in rejected_species:
                    self._reject_species(sp, "A dissociation reaction could not be written with valid basis species.")

            thermo_df = assign_worm_db_col_dtypes(thermo_df)

            # convert E units and calculate missing GHS values
            self.thermo_db = OBIGT2eos(thermo_df, fixGHS=True, tocal=True)


    def append_pickup(self, input_file, pickup_file, output_file=None,
                       strip_minerals=False):
        """
        Update a .3i or .6i input file by replacing its bottom half with the
        bottom half of a .3p or .6p pickup file. The bottom half begins at the
        line containing "Start of the bottom half of the input file".

        Parameters
        ----------
        input_file : str
            Path to a .3i or .6i input file.

        pickup_file : str
            Path to a .3p or .6p pickup file whose bottom half will replace the
            bottom half of ``input_file``.

        output_file : str, optional
            Path to write the merged result. If not specified, ``input_file`` is
            overwritten in place.

        strip_minerals : bool, default False
            If True, remove all mineral (pure and solid solution) entries from
            the pickup file's bottom half before merging, and update the matrix
            index limits accordingly.
        """

        marker = "Start of the bottom half of the input file"

        with open(input_file, "r") as f:
            input_lines = f.readlines()

        with open(pickup_file, "r") as f:
            pickup_lines = f.readlines()

        # Keep everything before the marker in the input file.
        top_half = []
        found_marker_in_input = False
        for line in input_lines:
            if marker in line:
                found_marker_in_input = True
                break
            top_half.append(line)

        if not found_marker_in_input:
            self.err_handler.raise_exception(
                "Could not find '{}' in '{}'.".format(marker, input_file))

        # Keep everything from the marker onward in the pickup file.
        bottom_half = []
        found_marker_in_pickup = False
        for line in pickup_lines:
            if marker in line:
                found_marker_in_pickup = True
            if found_marker_in_pickup:
                bottom_half.append(line)

        if not found_marker_in_pickup:
            self.err_handler.raise_exception(
                "Could not find '{}' in '{}'.".format(marker, pickup_file))

        if strip_minerals:
            bottom_half = _strip_minerals_from_bottom_half(bottom_half)

        merged = top_half + bottom_half

        if output_file is None:
            output_file = input_file

        with open(output_file, "w") as f:
            f.writelines(merged)


    def prepare_6i_for_mixing(self, input_file, data0_file,
                               T_fluid1, T_fluid2,
                               mass_ratio_factor=0.0,
                               output_file=None,
                               eq_output_file=None):
        """
        Prepare a .6i input file for cross-database fluid mixing.

        After using ``append_pickup`` to combine a top half (from one database)
        with a bottom half (from another database), this method adjusts the .6i
        file so that it is compatible with the target thermodynamic database.
        Specifically, it:

        1. Sets the temperature option to fluid mixing tracking (jtemp=3).
        2. Adds any basis species that are referenced in the top half's reactant
           but missing from the bottom half.
        3. Adds basis species needed to cover every pseudo-element in the target
           data0 file, so that minerals using those pseudo-elements can be found
           by EQ6's data compression routine (cmpdat).
        4. Updates the matrix index limits to account for added species.
        5. When ``eq_output_file`` is given, uses the speciation results to
           compute physically meaningful mass balance totals for newly added
           pseudo-element species (e.g. splitting total Fe into Fe+2/Fe+3).

        Parameters
        ----------
        input_file : str
            Path to the .6i input file to modify.

        data0_file : str
            Path to the target data0 thermodynamic database file.

        T_fluid1 : float
            Temperature of fluid 1 in degrees Celsius.

        T_fluid2 : float
            Temperature of fluid 2 in degrees Celsius.

        mass_ratio_factor : float, default 0.0
            Starting mass ratio factor for fluid mixing.

        output_file : str, optional
            Path to write the modified file. If not specified, ``input_file``
            is overwritten in place.

        eq_output_file : str, optional
            Path to a .3o or .6o output file from the EQ3/6 run that produced
            the bottom half of the input. Used to compute the oxidation-state
            distribution when the target database uses pseudo-elements (e.g.
            splitting nrp's total Fe into vnt's Fejiip and Fejiiip based on
            the actual Fe+2/Fe+3 ratio from the speciation).
        """

        # --- read files ---------------------------------------------------- #
        with open(input_file, "r") as f:
            lines = f.readlines()
        with open(data0_file, "r") as f:
            data0_text = f.read()

        eq_output_text = None
        if eq_output_file is not None:
            with open(eq_output_file, "r") as f:
                eq_output_text = f.read()

        # --- split 6i into top and bottom ---------------------------------- #
        marker = "Start of the bottom half of the input file"
        split_idx = None
        for i, line in enumerate(lines):
            if marker in line:
                split_idx = i
                break
        if split_idx is None:
            self.err_handler.raise_exception(
                "Could not find '{}' in '{}'.".format(marker, input_file))

        top_lines = lines[:split_idx]
        bottom_lines = lines[split_idx:]

        # --- reconcile cross-database species ------------------------------ #
        top_lines, bottom_lines, _ = _reconcile_cross_database_mixing(
            top_lines, bottom_lines, data0_text, eq_output_text)

        # --- replace temperature block in top half ------------------------- #
        top_lines = self._replace_temperature_block(
            top_lines, T_fluid1, T_fluid2, mass_ratio_factor)

        # --- write output -------------------------------------------------- #
        merged = top_lines + bottom_lines
        if output_file is None:
            output_file = input_file
        with open(output_file, "w") as f:
            f.writelines(merged)


    @staticmethod
    def _parse_data0_species_elements(data0_text):
        """
        Parse a data0 file and return a dict mapping each basis species and
        auxiliary basis species to its list of element names.
        """
        lines = data0_text.split("\n")
        sep = "+--"  # section separator prefix

        species_elements = {}

        # Find section boundaries
        elem_start = None
        basis_start = None
        aux_start = None
        aq_start = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "elements":
                elem_start = i
            elif stripped == "basis species":
                basis_start = i
            elif stripped == "auxiliary basis species":
                aux_start = i
            elif stripped == "aqueous species":
                aq_start = i
                break

        # Parse basis species section
        if basis_start is not None:
            end = aux_start if aux_start is not None else (aq_start or len(lines))
            species_elements.update(
                _parse_species_section(lines, basis_start + 2, end))

        # Parse auxiliary basis species section
        if aux_start is not None:
            end = aq_start if aq_start is not None else len(lines)
            species_elements.update(
                _parse_species_section(lines, aux_start + 2, end))

        return species_elements


    @staticmethod
    def _parse_data0_all_species_stoichiometry(data0_text):
        """
        Parse basis, auxiliary basis, AND aqueous species sections from a
        data0 file, retaining stoichiometric coefficients.
        Returns dict: species_name -> {element_name: stoich_coefficient}.
        """
        lines = data0_text.split("\n")

        basis_start = None
        aux_start = None
        aq_start = None
        solids_start = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "basis species":
                basis_start = i
            elif stripped == "auxiliary basis species":
                aux_start = i
            elif stripped == "aqueous species":
                aq_start = i
            elif stripped in ("solids", "solid solutions"):
                solids_start = i
                break

        species_stoich = {}

        if basis_start is not None:
            end = aux_start if aux_start is not None else (
                aq_start or solids_start or len(lines))
            species_stoich.update(
                _parse_species_section_with_stoich(lines, basis_start + 2, end))

        if aux_start is not None:
            end = aq_start if aq_start is not None else (
                solids_start or len(lines))
            species_stoich.update(
                _parse_species_section_with_stoich(lines, aux_start + 2, end))

        if aq_start is not None:
            end = solids_start if solids_start is not None else len(lines)
            species_stoich.update(
                _parse_species_section_with_stoich(lines, aq_start + 2, end))

        return species_stoich


    @staticmethod
    def _parse_data0_dissociation_reactions(data0_text):
        """
        Parse dissociation reactions from auxiliary basis and aqueous species
        sections of a data0 file.

        Returns dict: species_name → {basis_species_name: coefficient}.
        Basis species themselves are included with {self: 1.0}.
        """
        import re
        lines = data0_text.split("\n")

        basis_start = None
        aux_start = None
        aq_start = None
        solids_start = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "basis species":
                basis_start = i
            elif stripped == "auxiliary basis species":
                aux_start = i
            elif stripped == "aqueous species":
                aq_start = i
            elif stripped in ("solids", "solid solutions"):
                solids_start = i
                break

        dissoc = {}

        # Basis species: trivial identity reactions
        if basis_start is not None:
            end = aux_start if aux_start is not None else (
                aq_start or solids_start or len(lines))
            i = basis_start + 2
            while i < end:
                line = lines[i].strip()
                if not line or line.startswith("+--"):
                    i += 1
                    continue
                sp_name = line
                dissoc[sp_name] = {sp_name: 1.0}
                i += 1
                while i < end and not lines[i].strip().startswith("+--"):
                    i += 1
                if i < end:
                    i += 1

        # Parse auxiliary basis and aqueous species sections
        for section_start, section_end in [
            (aux_start, aq_start or solids_start or len(lines)),
            (aq_start, solids_start or len(lines)),
        ]:
            if section_start is None:
                continue
            i = section_start + 2
            while i < section_end:
                line = lines[i].strip()
                if not line or line.startswith("+--"):
                    i += 1
                    continue
                sp_name = line
                i += 1
                rxn = {}
                while i < section_end:
                    mline = lines[i].strip()
                    if mline.startswith("+--"):
                        i += 1
                        break
                    m = re.match(
                        r"(\d+)\s+species in aqueous dissociation reaction:",
                        mline)
                    if m:
                        n_species = int(m.group(1))
                        i += 1
                        while i < section_end and len(rxn) < n_species:
                            rline = lines[i].strip()
                            if rline.startswith("+--") or rline.startswith(
                                    "****"):
                                break
                            tokens = rline.split()
                            j = 0
                            while j < len(tokens) - 1:
                                try:
                                    coeff = float(tokens[j])
                                    rxn[tokens[j + 1]] = coeff
                                    j += 2
                                except ValueError:
                                    j += 1
                            i += 1
                        while i < section_end:
                            if lines[i].strip().startswith("+--"):
                                i += 1
                                break
                            i += 1
                        break
                    i += 1
                dissoc[sp_name] = rxn

        return dissoc


    @staticmethod
    def _parse_data0_pseudo_elements(data0_text):
        """
        Parse the elements section of a data0 file and return a list of
        element names that are pseudo-elements (i.e. redox-suppression
        fake elements like Fejiiip, Sjiin, Cjivp, etc.), preserving the
        order they appear in the data0 file.
        """
        lines = data0_text.split("\n")
        pseudo = []

        # Find elements section boundaries
        elem_start = None
        elem_end = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped == "elements":
                elem_start = i
            elif stripped == "basis species":
                elem_end = i
                break

        if elem_start is None or elem_end is None:
            return pseudo

        for i in range(elem_start + 2, elem_end):
            line = lines[i].strip()
            if not line or line.startswith("+--"):
                continue
            elem_name = line.split()[0]
            _, ox = parse_pseudoelement_name(elem_name)
            if ox is not None:
                pseudo.append(elem_name)

        return pseudo


    @staticmethod
    def _parse_reactant_species(top_lines):
        """
        Extract species names from Reaction tables in the top half of a .6i
        file. Returns a list of unique species names (excluding reactant names
        like 'Fluid 2').
        """
        species = []
        in_reaction = False
        reactant_name = None

        # First pass: find all reactant names
        reactant_names = set()
        for line in top_lines:
            if "|Reactant" in line and "(ureac(n))" in line:
                parts = line.split("|")
                for p in parts:
                    p = p.strip()
                    if p and p != "Reactant" and "(ureac(n))" not in p:
                        reactant_names.add(p)

        for line in top_lines:
            text = line.rstrip("\n")

            if "|->|Reaction" in text:
                in_reaction = True
                continue

            if in_reaction:
                if text.startswith("|--->|") and "(ubsri(i,n)" in text:
                    # Extract species name from: |--->|Name  |coeff| ...
                    parts = text.split("|")
                    # parts[0] = '', parts[1] = '--->', parts[2] = species name
                    if len(parts) >= 3:
                        sp_name = parts[2].strip()
                        if sp_name and sp_name not in reactant_names:
                            if sp_name not in species:
                                species.append(sp_name)
                elif text.startswith("|---") and "table header" in text:
                    continue
                elif text.startswith("|---"):
                    # Separator line — could be end of reaction table
                    # Check if next relevant line is still reaction data
                    pass
                elif not text.startswith("|--->"):
                    in_reaction = False

        return species


    @staticmethod
    def _parse_bottom_basis_species(bottom_lines):
        """
        Extract current basis species names from the bottom half of a .6i file.
        Returns a list of species names.
        """
        species = []
        in_mass_balance = False

        for line in bottom_lines:
            text = line.rstrip("\n")

            if "Mass Balance Species (Matrix Row Variables)" in text:
                in_mass_balance = True
                continue

            if in_mass_balance:
                if text.startswith("* Valid jflag"):
                    break
                # Species lines look like:
                # |H2O                     Aqueous solution        |Moles   ...
                if (text.startswith("|") and "Aqueous solution" in text
                        and ("Moles" in text or "Make non-basis" in text)):
                    sp_name = text[1:25].strip()
                    if sp_name:
                        species.append(sp_name)

        return species


    @staticmethod
    def _parse_bottom_mass_balance_totals(bottom_lines, column="equilibrium"):
        """
        Parse the Mass Balance Totals table from the bottom half of a .6i
        file. Returns a dict mapping species name to its total in moles.

        Parameters
        ----------
        column : str
            "equilibrium" for the Equilibrium System column (mtbi(n)),
            "aqueous" for the Aqueous Solution column (mtbaqi(n)).
        """
        col_idx = 2 if column == "equilibrium" else 3
        totals = {}
        in_totals = False
        for line in bottom_lines:
            text = line.rstrip("\n")
            if "Mass Balance Totals (moles)" in text:
                in_totals = True
                continue
            if in_totals:
                if "Electrical imbalance" in text:
                    break
                if text.startswith("|") and "Aqueous" in text:
                    sp_name = text[1:25].strip()
                    parts = text.split("|")
                    if len(parts) >= (col_idx + 1):
                        val_str = parts[col_idx].strip()
                        try:
                            totals[sp_name] = float(val_str)
                        except ValueError:
                            pass
        return totals

    @staticmethod
    def _parse_output_numerical_composition(output_text):
        """
        Parse the "Numerical Composition of the Aqueous Solution" table
        from a .3o or .6o output file.  Returns a dict mapping species
        name to its molality (the last numeric column).

        For .6o files with multiple time steps, the *last* occurrence of
        the table is used (representing the final equilibrium state).
        """
        lines = output_text.split("\n")
        composition = {}
        i = 0
        while i < len(lines):
            if "Numerical Composition of the Aqueous Solution" in lines[i]:
                composition.clear()
                i += 1
                # Skip header lines until we hit species data
                while i < len(lines):
                    line = lines[i].strip()
                    if not line:
                        i += 1
                        continue
                    if line.startswith("Species") or line.startswith("---"):
                        i += 1
                        continue
                    if "physical significance" in line:
                        break
                    # Parse species line. Formats:
                    # .3o: " H2O    9.89E+05  9.66E+05  5.49E-02  5.55E+01"
                    # .6o: " H2O    9.66E+05  5.55E+01"
                    tokens = line.split()
                    if len(tokens) >= 2:
                        sp_name = tokens[0]
                        try:
                            molality = float(tokens[-1])
                            composition[sp_name] = molality
                        except ValueError:
                            pass
                    i += 1
                continue
            i += 1
        return composition

    @staticmethod
    def _parse_output_species_distribution(output_text):
        """
        Parse the "Distribution of Aqueous Solute Species" table from a
        .3o or .6o output file. Returns a dict mapping species name to its
        molality. For .6o files with multiple time steps, the last
        occurrence is used.
        """
        lines = output_text.split("\n")
        distribution = {}
        i = 0
        while i < len(lines):
            if "Distribution of Aqueous Solute Species" in lines[i]:
                distribution.clear()
                i += 1
                while i < len(lines):
                    line = lines[i].strip()
                    if not line:
                        if distribution:
                            break
                        i += 1
                        continue
                    if line.startswith("Species") or line.startswith("---"):
                        i += 1
                        continue
                    tokens = line.split()
                    if len(tokens) >= 2:
                        sp_name = tokens[0]
                        try:
                            molality = float(tokens[1])
                            distribution[sp_name] = molality
                        except ValueError:
                            pass
                    i += 1
                continue
            i += 1
        return distribution

    @staticmethod
    def _compute_missing_totals(missing, current_basis, all_species_elems,
                                existing_totals, output_composition,
                                species_distribution=None,
                                species_stoich=None,
                                h2o_aq_moles=None,
                                comp_moles=None):
        """
        Compute mass balance totals for species being added to the bottom
        half when switching databases (e.g. nrp → vnt).

        When ``species_distribution`` and ``species_stoich`` are provided,
        pseudo-element totals are computed from the full aqueous species
        distribution table cross-referenced against data0 element
        compositions.  Each missing species' pseudo-element total
        (molality summed over all aqueous species containing that
        pseudo-element) is converted to moles using the H2O mass from
        the pickup file.  The donor species (an existing basis species
        sharing the same real element) is reduced by the same amount to
        conserve total elemental mass.

        Falls back to ``output_composition`` ratio method, then to a
        small fraction (1E-10) of the donor total.

        Returns a tuple (missing_totals, donor_adjustments, allocation_log)
        where donor_adjustments maps existing species names to the amount
        that should be subtracted from their mass balance total to conserve
        the real-element mass.
        """
        missing_totals = {}
        donor_adjustments = {}
        allocation_log = {}
        fallback_fraction = 1e-10

        pe_totals = None
        pe_contributions = {}
        if species_distribution and species_stoich:
            pe_totals, pe_contributions = _compute_pseudo_element_totals(
                species_distribution, species_stoich)

        h2o_kg = None
        if h2o_aq_moles is not None and h2o_aq_moles > 0:
            h2o_kg = h2o_aq_moles * 18.01528 / 1000.0

        for sp in missing:
            sp_elems = all_species_elems.get(sp, [])
            computed_total = None
            donor_sp = None
            log_entry = {"method": None, "pseudo_element": None,
                         "molality_total": None, "h2o_kg": h2o_kg,
                         "moles": None, "contributors": []}

            for elem in sp_elems:
                real_elem, ox = parse_pseudoelement_name(elem)
                if ox is None:
                    continue

                if pe_totals is not None and h2o_kg is not None:
                    sp_pe_molality = pe_totals.get(elem, 0.0)
                    if sp_pe_molality > 0:
                        computed_total = sp_pe_molality * h2o_kg
                        log_entry["method"] = "species_distribution"
                    else:
                        pe_comp = (abs(comp_moles.get(elem, 0.0))
                                   if comp_moles else 0.0)
                        computed_total = pe_comp if pe_comp > 0 else 1e-30
                        log_entry["method"] = "pe_zero_floor"
                    log_entry["pseudo_element"] = elem
                    log_entry["molality_total"] = sp_pe_molality
                    log_entry["moles"] = computed_total
                    contribs = pe_contributions.get(elem, [])
                    contribs_sorted = sorted(
                        contribs, key=lambda x: x[3], reverse=True)
                    log_entry["contributors"] = [
                        {"species": c[0], "molality": c[1],
                         "stoich_coeff": c[2], "contribution": c[3]}
                        for c in contribs_sorted]
                    for existing_sp in current_basis:
                        for ex_elem in all_species_elems.get(
                                existing_sp, []):
                            ex_real, ex_ox = parse_pseudoelement_name(
                                ex_elem)
                            if (ex_ox is not None
                                    and ex_real == real_elem):
                                donor_sp = existing_sp
                                break
                        if donor_sp is not None:
                            break
                    break

                for existing_sp in current_basis:
                    for ex_elem in all_species_elems.get(existing_sp, []):
                        ex_real, ex_ox = parse_pseudoelement_name(ex_elem)
                        if ex_ox is not None and ex_real == real_elem:
                            donor_total = existing_totals.get(existing_sp)
                            if donor_total is None or abs(donor_total) <= 1e-30:
                                continue

                            sp_molal = output_composition.get(sp)
                            donor_molal = output_composition.get(existing_sp)
                            if (sp_molal is not None
                                    and donor_molal is not None
                                    and abs(donor_molal) > 0):
                                ratio = abs(sp_molal) / abs(donor_molal)
                                ratio = min(ratio, 0.5)
                                computed_total = abs(donor_total) * ratio
                                log_entry["method"] = "numerical_composition"
                            else:
                                computed_total = (abs(donor_total)
                                                  * fallback_fraction)
                                log_entry["method"] = "fallback_1e-10"
                            log_entry["moles"] = computed_total
                            donor_sp = existing_sp
                            break
                    if computed_total is not None:
                        break
                if computed_total is not None:
                    break

            if log_entry["method"] is None:
                log_entry["method"] = "no_pseudo_element"
            missing_totals[sp] = computed_total
            allocation_log[sp] = log_entry
            if computed_total is not None and donor_sp is not None:
                donor_adjustments[donor_sp] = (
                    donor_adjustments.get(donor_sp, 0.0) + computed_total)

        return missing_totals, donor_adjustments, allocation_log

    @staticmethod
    def _insert_missing_basis_species(bottom_lines, missing, flags,
                                      totals=None, donor_adjustments=None):
        """
        Insert missing basis species into the three relevant tables in the
        bottom half: Mass Balance Species, Mass Balance Totals, and Matrix
        Column Variables. Also updates the matrix index limits and adjusts
        donor species mass balance totals when donor_adjustments is given.
        """
        import math
        if donor_adjustments is None:
            donor_adjustments = {}
        n_add = len(missing)
        new_lines = []
        i = 0
        inserted_column_vars = False

        # Track which section we're in to avoid mismatched insertions
        SECTION_NONE = 0
        SECTION_MASS_BAL_SPECIES = 1
        SECTION_MASS_BAL_TOTALS = 2
        SECTION_COLUMN_VARS = 3
        current_section = SECTION_NONE

        while i < len(bottom_lines):
            text = bottom_lines[i].rstrip("\n")

            # --- Track current section ------------------------------------- #
            if "Mass Balance Species (Matrix Row Variables)" in text:
                current_section = SECTION_MASS_BAL_SPECIES
            elif "Mass Balance Totals (moles)" in text:
                current_section = SECTION_MASS_BAL_TOTALS
            elif "Matrix Column Variables and Values" in text:
                current_section = SECTION_COLUMN_VARS
            elif "Phases and Species in the PRS" in text:
                current_section = SECTION_NONE

            # --- Update matrix index limits -------------------------------- #
            if "(kbt)" in text:
                new_lines.append(_update_matrix_index(text, n_add))
                i += 1
                continue
            if "(kmt)" in text:
                new_lines.append(_update_matrix_index(text, n_add))
                i += 1
                continue
            if "(kxt)" in text:
                new_lines.append(_update_matrix_index(text, n_add))
                i += 1
                continue
            if "(kdim)" in text:
                new_lines.append(_update_matrix_index(text, n_add))
                i += 1
                continue

            # --- Insert into Mass Balance Species table -------------------- #
            if (current_section == SECTION_MASS_BAL_SPECIES
                    and text.startswith("* Valid jflag strings")):
                # Insert new rows just before the separator+comment block.
                if (new_lines and new_lines[-1].rstrip("\n").startswith(
                        "|----")):
                    sep_line = new_lines.pop()
                    for sp in missing:
                        flag = flags[sp]
                        row = "|{:<24s}Aqueous solution        |{:<16s}| --         |\n".format(
                            sp, flag)
                        new_lines.append(row)
                    new_lines.append(sep_line)
                new_lines.append(bottom_lines[i])
                i += 1
                continue

            # --- Adjust donor species in Mass Balance Totals table --------- #
            if (current_section == SECTION_MASS_BAL_TOTALS
                    and text.startswith("|") and "Aqueous" in text
                    and "Electrical imbalance" not in text):
                sp_name = text[1:25].strip()
                adj = donor_adjustments.get(sp_name)
                if adj is not None and adj > 0:
                    parts = text.split("|")
                    if len(parts) >= 4:
                        try:
                            eq_total = float(parts[2].strip())
                            aq_total = float(parts[3].strip())
                            new_eq = eq_total - adj
                            new_aq = aq_total - adj
                            row = "|{:<24s}Aqueous |{:22.15E}|{:22.15E}|\n".format(
                                sp_name, new_eq, new_aq)
                            new_lines.append(row)
                            i += 1
                            continue
                        except ValueError:
                            pass

            # --- Insert into Mass Balance Totals table --------------------- #
            if (current_section == SECTION_MASS_BAL_TOTALS
                    and "Electrical imbalance" in text):
                for sp in missing:
                    mt = (totals or {}).get(sp)
                    if mt is not None and mt > 0:
                        val = "{:22.15E}".format(mt)
                    else:
                        val = "{:22.15E}".format(1e-30)
                    row = "|{:<24s}Aqueous |{}|{}|\n".format(sp, val, val)
                    new_lines.append(row)
                new_lines.append(bottom_lines[i])
                i += 1
                continue

            # --- Insert into Matrix Column Variables table ----------------- #
            if (current_section == SECTION_COLUMN_VARS
                    and not inserted_column_vars
                    and "Aqueous solution" in text
                    and "(zvclgi(n))" not in text
                    and "table header" not in text
                    and text.startswith("|")):
                new_lines.append(bottom_lines[i])
                i += 1
                # Walk forward through remaining Aqueous solution entries
                while i < len(bottom_lines):
                    next_text = bottom_lines[i].rstrip("\n")
                    if ("Aqueous solution" in next_text
                            and next_text.startswith("|")
                            and "(zvclgi(n))" not in next_text):
                        new_lines.append(bottom_lines[i])
                        i += 1
                        continue
                    else:
                        # Reached end of Aqueous solution entries — insert here
                        for sp in missing:
                            mt = (totals or {}).get(sp)
                            if mt is not None and mt > 0:
                                log_val = math.log10(mt)
                            else:
                                log_val = -30.0
                            val = "{:22.15E}".format(log_val)
                            row = "|{:<24s}Aqueous solution        |{}| --   |\n".format(
                                sp, val)
                            new_lines.append(row)
                        inserted_column_vars = True
                        break
                continue

            new_lines.append(bottom_lines[i])
            i += 1

        return new_lines


    @staticmethod
    def _replace_temperature_block(top_lines, T_fluid1, T_fluid2,
                                    mass_ratio_factor):
        """
        Replace the temperature option block in the top half of a .6i file
        with the jtemp=3 fluid mixing tracking option.
        """
        temp_block = [
            "|Temperature option (jtemp):                                                   |\n",
            "|  [ ] ( 0) Constant temperature:                                              |\n",
            "|             Value (C)         | 0.00000E+00| (tempcb)                        |\n",
            "|  [ ] ( 1) Linear tracking in Xi:                                             |\n",
            "|             Base Value (C)    | 0.00000E+00| (tempcb)                        |\n",
            "|             Derivative        | 0.00000E+00| (ttk(1))                        |\n",
            "|  [ ] ( 2) Linear tracking in time:                                           |\n",
            "|             Base Value (C)    | 0.00000E+00| (tempcb)                        |\n",
            "|             Derivative        | 0.00000E+00| (ttk(1))                        |\n",
            "|  [x] ( 3) Fluid mixing tracking (fluid 2 = special reactant):                |\n",
            "|             T of fluid 1 (C)  |{:12.5E}| (tempcb)                        |\n".format(T_fluid1),
            "|             T of fluid 2 (C)  |{:12.5E}| (ttk(2))                        |\n".format(T_fluid2),
            "|             Mass ratio factor |{:12.5E}| (ttk(1))                        |\n".format(mass_ratio_factor),
            "|------------------------------------------------------------------------------|\n",
        ]

        new_lines = []
        i = 0
        skip = False
        while i < len(top_lines):
            text = top_lines[i].rstrip("\n")
            if "Temperature option (jtemp):" in text:
                skip = True
                new_lines.extend(temp_block)
                i += 1
                continue
            if skip:
                if "Pressure option (jpress):" in text:
                    skip = False
                    new_lines.append(top_lines[i])
                i += 1
                continue
            new_lines.append(top_lines[i])
            i += 1

        return new_lines


def join_input_files(input_files, filename="joined_input.csv",
                     fill_value=None, return_df=False):
    """
    Join two or more CSV input files into a single file.

    Each input file must contain a ``Temperature`` column and an ``H+``
    column with subheader ``pH``. If any file has a ``Pressure`` column
    then all files must have one. Species columns present in one file but
    not another are filled with *fill_value*.

    Parameters
    ----------
    input_files : list of str or pandas.DataFrame
        Inputs to join. Each element can be a file path (str) or a
        pandas DataFrame read from a CSV input file (where the first row
        contains unit subheaders and subsequent rows contain data).
        The list may mix strings and DataFrames.

    filename : str, default "joined_input.csv"
        Output CSV file path. Ignored when ``return_df`` is ``True``.

    fill_value : str or numeric, optional
        Value used to fill species columns that are missing from a file.
        If ``None`` (the default), missing values are left blank. A common
        alternative is ``"1E-18"``.

    return_df : bool, default False
        If ``True``, return a pandas DataFrame instead of writing a CSV
        file. The DataFrame has species columns as headers, with the
        first row containing unit subheaders. If ``False`` (the default),
        write a CSV file and return its path.

    Returns
    -------
    str or pandas.DataFrame
        The path to the created CSV file, or a DataFrame if
        ``return_df`` is ``True``.
    """
    import csv as csv_module

    if not isinstance(input_files, (list, tuple)) or len(input_files) < 2:
        raise Exception("input_files must be a list of at least two file "
                        "paths or DataFrames.")

    file_headers = []
    file_units = []
    file_rows = []

    for entry in input_files:
        if isinstance(entry, pd.DataFrame):
            df = entry
            if len(df) < 1:
                raise Exception("A DataFrame input must have at least "
                                "1 row (units) plus data rows.")
            if isinstance(df.columns, pd.MultiIndex):
                hdr = [str(c[0]) for c in df.columns]
                units = [str(c[1]) for c in df.columns]
                data = []
                for row_idx in range(len(df)):
                    data.append([str(v) for v in df.iloc[row_idx]])
            else:
                hdr = [str(c) for c in df.columns]
                units = [str(v) for v in df.iloc[0]]
                data = []
                for row_idx in range(1, len(df)):
                    data.append([str(v) for v in df.iloc[row_idx]])
            file_headers.append(hdr)
            file_units.append(units)
            file_rows.append(data)
        else:
            with open(entry, newline="") as f:
                reader = csv_module.reader(f)
                rows = list(reader)
            if len(rows) < 3:
                raise Exception("Input file '" + str(entry)
                                + "' must have at least 3 rows "
                                "(header, units, data).")
            file_headers.append(rows[0])
            file_units.append(rows[1])
            file_rows.append(rows[2:])

    input_labels = []
    for entry in input_files:
        if isinstance(entry, pd.DataFrame):
            input_labels.append("DataFrame")
        else:
            input_labels.append(str(entry))

    for i, (hdr, label) in enumerate(zip(file_headers, input_labels)):
        hdr_upper = [h.strip().upper() for h in hdr]
        if "TEMPERATURE" not in hdr_upper:
            raise Exception("Input '" + label
                            + "' is missing a Temperature column.")
        if "H+" not in [h.strip() for h in hdr]:
            raise Exception("Input '" + label
                            + "' is missing an H+ column.")
        units_for_h = file_units[i][
            [h.strip() for h in hdr].index("H+")].strip()
        if units_for_h.lower() != "ph":
            raise Exception("Input '" + label
                            + "' H+ column must have subheader 'pH', "
                            "got '" + units_for_h + "'.")

    has_pressure = []
    for hdr in file_headers:
        hdr_upper = [h.strip().upper() for h in hdr]
        has_pressure.append("PRESSURE" in hdr_upper)
    if any(has_pressure) and not all(has_pressure):
        raise Exception("Cannot join: some input files have a Pressure "
                        "column and others do not. Either all files must "
                        "have Pressure or none can.")

    meta_names_upper = {"SAMPLE", "H+", "PRESSURE", "TEMPERATURE",
                        "LOGFO2"}
    file_col_maps = []
    all_species = []
    for hdr, units in zip(file_headers, file_units):
        col_map = {}
        for j, (h, u) in enumerate(zip(hdr, units)):
            h_stripped = h.strip()
            col_map[h_stripped] = (j, u.strip())
            if h_stripped.upper() not in meta_names_upper:
                if h_stripped not in all_species:
                    all_species.append(h_stripped)
        file_col_maps.append(col_map)

    meta_order = ["Sample", "H+", "Pressure", "Temperature", "logfO2"]
    meta_units = {"Sample": "id", "H+": "pH", "Pressure": "bar",
                  "Temperature": "degC", "logfO2": "logfO2"}

    out_meta = []
    for m in meta_order:
        for cm in file_col_maps:
            if m in cm:
                out_meta.append((m, cm[m][1]))
                break

    out_header = [m for m, _ in out_meta] + all_species
    out_units = [u for _, u in out_meta]
    for sp in all_species:
        unit = "Molality"
        for cm in file_col_maps:
            if sp in cm:
                unit = cm[sp][1]
                break
        out_units.append(unit)

    fill_str = "" if fill_value is None else str(fill_value)

    out_rows = []
    for file_idx in range(len(input_files)):
        cm = file_col_maps[file_idx]
        for data_row in file_rows[file_idx]:
            row = []
            for m, _ in out_meta:
                if m in cm:
                    idx = cm[m][0]
                    row.append(data_row[idx] if idx < len(data_row)
                               else "")
                else:
                    row.append("")
            for sp in all_species:
                if sp in cm:
                    idx = cm[sp][0]
                    row.append(data_row[idx] if idx < len(data_row)
                               else fill_str)
                else:
                    row.append(fill_str)
            out_rows.append(row)

    if return_df:
        df_out = pd.DataFrame([out_units] + out_rows, columns=out_header)
        return df_out

    with open(filename, "w", newline="") as f:
        writer = csv_module.writer(f)
        writer.writerow(out_header)
        writer.writerow(out_units)
        for row in out_rows:
            writer.writerow(row)

    return filename


def drop_min_molal_bases(df, db, minimum_molality):
    """
    Drop species columns from a joined input DataFrame where all data
    values equal the minimum molality, provided no retained auxiliary
    basis species depends on them via its dissociation reaction.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame produced by ``join_input_files`` (first row is unit
        subheaders, subsequent rows are data).

    db : str
        Target database lettercode (e.g., ``"fcs"``) or path to a data0
        file.

    minimum_molality : float
        The minimum molality value. Columns where every data value equals
        this threshold are candidates for removal.

    Returns
    -------
    pandas.DataFrame
        A copy of *df* with eligible columns removed.
    """
    if os.path.exists(db):
        with open(db) as f:
            data0_text = f.read()
    else:
        candidate = db if db.startswith("data0.") else "data0." + db
        data0_text = None
        for search_dir in [os.getcwd()]:
            path = os.path.join(search_dir, candidate)
            if os.path.exists(path):
                with open(path) as f:
                    data0_text = f.read()
                break
        if data0_text is None:
            raise Exception(
                "Could not find target data0 file for '" + str(db)
                + "'. Place the file in the current working directory or "
                "specify the full path.")

    dissoc = AqEquil._parse_data0_dissociation_reactions(data0_text)

    basis_stoich, _ = _parse_data0_basis_species_full(data0_text)
    basis_set = set(basis_stoich.keys())

    _lines = data0_text.split("\n")
    _aux_start = None
    _aq_start = None
    for _i, _line in enumerate(_lines):
        _s = _line.strip()
        if _s == "auxiliary basis species":
            _aux_start = _i
        elif _s == "aqueous species":
            _aq_start = _i
            break
    aux_set = set()
    if _aux_start is not None:
        _end = _aq_start if _aq_start is not None else len(_lines)
        _j = _aux_start + 2
        while _j < _end:
            _sl = _lines[_j].strip()
            if not _sl or _sl.startswith("+--"):
                _j += 1
                continue
            aux_set.add(_sl)
            _j += 1
            while _j < _end and not _lines[_j].strip().startswith("+--"):
                _j += 1
            if _j < _end:
                _j += 1

    meta_upper = {"SAMPLE", "H+", "PRESSURE", "TEMPERATURE", "LOGFO2"}
    species_cols = [c for c in df.columns if c.upper() not in meta_upper]

    all_min_cols = set()
    for col in species_cols:
        vals = df[col].iloc[1:]
        try:
            if all(float(v) == minimum_molality for v in vals):
                all_min_cols.add(col)
        except (ValueError, TypeError):
            continue

    # Build dependency graph: for each aux species that is a column,
    # record which other columns its dissociation reaction depends on
    # (excluding system species).
    system_species = {"H2O", "H+", "O2(g)", "O2", "e-"}
    aux_deps = {}
    for sp in aux_set:
        if sp not in df.columns:
            continue
        rxn = dissoc.get(sp, {})
        deps = set()
        for product in rxn:
            if product == sp:
                continue
            if product in system_species:
                continue
            if product in df.columns:
                deps.add(product)
        if deps:
            aux_deps[sp] = deps

    # Iteratively determine which columns can be dropped:
    # A column can only be dropped if it is all-minimum AND no retained
    # aux species depends on it.
    to_drop = set()
    changed = True
    while changed:
        changed = False
        for col in list(all_min_cols - to_drop):
            depended_on = False
            for aux_sp, deps in aux_deps.items():
                if aux_sp in to_drop:
                    continue
                if col in deps:
                    depended_on = True
                    break
            if not depended_on:
                to_drop.add(col)
                changed = True

    return df.drop(columns=list(to_drop))


def _parse_species_section(lines, start, end):
    """
    Parse a basis species or auxiliary basis species section from a data0 file.
    Returns a dict mapping species name to a list of element names.
    """
    import re
    species_elements = {}
    i = start
    while i < end:
        line = lines[i].strip()
        # Skip separators and blank lines
        if not line or line.startswith("+--"):
            i += 1
            continue
        # This should be a species name
        sp_name = line
        i += 1
        # Skip metadata lines until we find the element count line
        elems = []
        while i < end:
            mline = lines[i].strip()
            if mline.startswith("+--"):
                # End of this species entry
                i += 1
                break
            elem_match = re.match(r"(\d+)\s+element\(s\):", mline)
            if elem_match:
                n_elems = int(elem_match.group(1))
                # Next line(s) have element data — read exactly n_elems elements
                i += 1
                while i < end and len(elems) < n_elems:
                    eline = lines[i].strip()
                    if eline.startswith("+--"):
                        break
                    if eline.startswith("*"):
                        i += 1
                        continue
                    # Parse element composition line:
                    # "1.0000 H   1.0000 O"
                    tokens = eline.split()
                    j = 0
                    while j < len(tokens) - 1 and len(elems) < n_elems:
                        try:
                            float(tokens[j])
                            elems.append(tokens[j + 1])
                            j += 2
                        except ValueError:
                            j += 1
                    i += 1
                # Skip remaining lines until separator
                while i < end:
                    if lines[i].strip().startswith("+--"):
                        i += 1
                        break
                    i += 1
                break
            i += 1
        species_elements[sp_name] = elems
    return species_elements


def _parse_species_section_with_stoich(lines, start, end):
    """
    Like _parse_species_section but retains stoichiometric coefficients.
    Returns a dict mapping species name to a dict of {element_name: coefficient}.
    """
    import re
    species_stoich = {}
    i = start
    while i < end:
        line = lines[i].strip()
        if not line or line.startswith("+--"):
            i += 1
            continue
        sp_name = line
        i += 1
        elems = {}
        while i < end:
            mline = lines[i].strip()
            if mline.startswith("+--"):
                i += 1
                break
            elem_match = re.match(r"(\d+)\s+element\(s\):", mline)
            if elem_match:
                n_elems = int(elem_match.group(1))
                i += 1
                while i < end and len(elems) < n_elems:
                    eline = lines[i].strip()
                    if eline.startswith("+--"):
                        break
                    if eline.startswith("*"):
                        i += 1
                        continue
                    tokens = eline.split()
                    j = 0
                    while j < len(tokens) - 1 and len(elems) < n_elems:
                        try:
                            coeff = float(tokens[j])
                            elems[tokens[j + 1]] = coeff
                            j += 2
                        except ValueError:
                            j += 1
                    i += 1
                while i < end:
                    if lines[i].strip().startswith("+--"):
                        i += 1
                        break
                    i += 1
                break
            i += 1
        species_stoich[sp_name] = elems
    return species_stoich


def _parse_data0_basis_species_full(data0_text):
    """
    Parse the basis species section of a data0 file and return dicts
    mapping each basis species name to:
      - its element stoichiometry  {element: coeff}
      - its charge (float)
    """
    import re
    lines = data0_text.split("\n")
    basis_start = None
    aux_start = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s == "basis species":
            basis_start = i
        elif s == "auxiliary basis species":
            aux_start = i
            break
    if basis_start is None:
        return {}, {}
    end = aux_start if aux_start is not None else len(lines)
    stoich = {}
    charges = {}
    i = basis_start + 2
    while i < end:
        line = lines[i].strip()
        if not line or line.startswith("+--"):
            i += 1
            continue
        sp_name = line
        i += 1
        elems = {}
        charge = 0.0
        while i < end:
            mline = lines[i].strip()
            if mline.startswith("+--"):
                i += 1
                break
            cm = re.match(r"charge\s*=\s*([-\d.]+)", mline)
            if cm:
                charge = float(cm.group(1))
            elem_match = re.match(r"(\d+)\s+element\(s\):", mline)
            if elem_match:
                n_elems = int(elem_match.group(1))
                i += 1
                while i < end and len(elems) < n_elems:
                    eline = lines[i].strip()
                    if eline.startswith("+--"):
                        break
                    if eline.startswith("*"):
                        i += 1
                        continue
                    tokens = eline.split()
                    j = 0
                    while j < len(tokens) - 1 and len(elems) < n_elems:
                        try:
                            coeff = float(tokens[j])
                            elems[tokens[j + 1]] = coeff
                            j += 2
                        except ValueError:
                            j += 1
                    i += 1
                while i < end:
                    if lines[i].strip().startswith("+--"):
                        i += 1
                        break
                    i += 1
                break
            i += 1
        stoich[sp_name] = elems
        charges[sp_name] = charge
    return stoich, charges


def _apply_pe_totals_to_composition(top_lines, pe_totals, h2o_kg,
                                    h2o_moles, target_elements,
                                    data0_text=None):
    """
    Replace pseudo-element entries in the Composition block with values
    computed from the species distribution.

    pe_totals maps pseudo-element name → total molality.  Values are
    converted to moles via ``molality * h2o_kg``.

    H and O entries are left unchanged because they represent total
    elemental amounts of the physical fluid; ``_rebuild_reaction_block``
    recomputes H2O, H+, and O2(g) in the Reaction block to close the
    hydrogen, oxygen, and charge balances against the updated
    pseudo-element values.
    """

    stoich_template = (
        "|--->|{elem}|{val}| (uesri(i,n), cesri(i,n))                |\n")

    # Parse the current Composition block
    comp = {}  # element -> moles
    comp_order = []  # preserve order
    in_comp = False
    comp_start_idx = None
    comp_end_idx = None
    for idx, line in enumerate(top_lines):
        text = line.rstrip("\n")
        if "|->|Composition" in text:
            in_comp = True
            comp_start_idx = idx
            continue
        if in_comp and "|->|Reaction" in text:
            comp_end_idx = idx
            break
        if in_comp and "|--->|" in text and "table header" not in text:
            parts = text.split("|")
            if len(parts) >= 4:
                elem = parts[2].strip()
                try:
                    val = float(parts[3])
                    comp[elem] = val
                    if elem not in [e for e, _ in comp_order]:
                        comp_order.append((elem, val))
                except ValueError:
                    pass

    if comp_start_idx is None or comp_end_idx is None:
        return top_lines

    # Update pseudo-element values from pe_totals, adding new ones if needed
    for pe, molality in pe_totals.items():
        moles = molality * h2o_kg
        if pe in comp:
            comp[pe] = moles
            comp_order = [(e, moles if e == pe else v)
                          for e, v in comp_order]
        else:
            comp[pe] = moles
            comp_order.append((pe, moles))

    sep = ("|----------------------------------------------"
           "--------------------------------|\n")

    # Rebuild Composition block: keep header lines, replace data lines
    new_comp_lines = []
    # Composition marker line
    new_comp_lines.append(top_lines[comp_start_idx])
    # separator + table header + separator
    new_comp_lines.append(sep)
    new_comp_lines.append(
        "|--->|Element |Stoich. Number        "
        "| (this is a table header)                |\n")
    new_comp_lines.append(sep)
    # Element entries
    for elem, val in comp_order:
        padded = elem.ljust(8)
        val_str = "{:22.15E}".format(val)
        new_comp_lines.append(stoich_template.format(
            elem=padded, val=val_str))
    # Closing separator
    new_comp_lines.append(sep)

    # Replace in top_lines
    new_top = (top_lines[:comp_start_idx] + new_comp_lines
               + top_lines[comp_end_idx:])
    return new_top


def _rebuild_reaction_block(top_lines, data0_text):
    """
    Rebuild the Reaction block in top_lines so that it uses the target
    database's basis species with molar amounts derived from the
    (already-remapped) Composition block.

    The Composition block encodes element → moles.  Each basis species in the
    target database carries exactly one "primary" element (i.e. the non-H,
    non-O element, or a pseudo-element).  The Reaction coefficient for that
    species is set to the moles of its primary element.  H2O, H+, and O2(g)
    are then calculated to close the hydrogen, oxygen, and charge balances.
    """
    from .redox_suppression import parse_pseudoelement_name

    basis_stoich, basis_charges = _parse_data0_basis_species_full(data0_text)

    # --- Parse the remapped Composition block from top_lines --- #
    comp = {}
    in_comp = False
    for line in top_lines:
        text = line.rstrip("\n")
        if "|->|Composition" in text:
            in_comp = True
            continue
        if in_comp and "|->|Reaction" in text:
            break
        if in_comp and "|--->|" in text and "table header" not in text:
            parts = text.split("|")
            if len(parts) >= 4:
                elem = parts[2].strip()
                try:
                    val = float(parts[3])
                    comp[elem] = val
                except ValueError:
                    pass

    if not comp:
        return top_lines

    # --- Build mapping: primary_element → basis_species --- #
    # "Primary element" is the non-H, non-O element of each species.
    # H2O, H+, and O2(g) are handled separately for balance closure.
    balance_species = {"H2O", "H+", "O2(g)"}
    elem_to_basis = {}
    for sp, sp_elems in basis_stoich.items():
        if sp in balance_species:
            continue
        for elem in sp_elems:
            if elem not in ("H", "O"):
                elem_to_basis[elem] = sp
                break

    # --- Assign reaction coefficients --- #
    rxn_coeffs = {}
    total_H_from_species = 0.0
    total_O_from_species = 0.0
    total_charge = 0.0

    for elem, moles in comp.items():
        if elem in ("H", "O"):
            continue
        sp = elem_to_basis.get(elem)
        if sp is None:
            continue
        sp_elems = basis_stoich[sp]
        elem_per_sp = sp_elems.get(elem, 1.0)
        if elem_per_sp == 0:
            elem_per_sp = 1.0
        sp_moles = moles / elem_per_sp
        rxn_coeffs[sp] = sp_moles
        h_per = sp_elems.get("H", 0.0)
        o_per = sp_elems.get("O", 0.0)
        total_H_from_species += sp_moles * h_per
        total_O_from_species += sp_moles * o_per
        total_charge += sp_moles * basis_charges.get(sp, 0.0)

    # --- Close balances for H2O, H+, O2(g) --- #
    # Total H from Composition:
    total_H = comp.get("H", 0.0)
    # Total O from Composition:
    total_O = comp.get("O", 0.0)

    # H+ carries 1 H and charge +1.  Its coefficient closes the charge
    # balance: sum of all species charges must be zero.
    # charge_from_H+ = c_Hplus * 1.0
    # total_charge + charge_from_H+ = 0
    c_Hplus = -total_charge
    rxn_coeffs["H+"] = c_Hplus
    total_H_from_species += c_Hplus * 1.0

    # H2O carries 2 H and 1 O.  It closes the hydrogen balance.
    # total_H_from_species + c_H2O * 2 = total_H
    c_H2O = (total_H - total_H_from_species) / 2.0
    rxn_coeffs["H2O"] = c_H2O
    total_O_from_species += c_H2O * 1.0

    # O2(g) carries 2 O.  It closes the oxygen balance.
    # total_O_from_species + c_O2g * 2 = total_O
    c_O2g = (total_O - total_O_from_species) / 2.0
    rxn_coeffs["O2(g)"] = c_O2g

    # --- Rebuild Reaction block in top_lines --- #
    rxn_template = (
        "|--->|{sp}|{val}| (ubsri(i,n), cbsri(i,n))|\n")

    new_top = []
    skip_old_rxn = False
    for line in top_lines:
        text = line.rstrip("\n")
        if "|->|Reaction" in text:
            new_top.append(line)
            # Add separator and header
            new_top.append(
                "|----------------------------------------------"
                "--------------------------------|\n")
            new_top.append(
                "|--->|Species                 |Reaction "
                "Coefficient  | (this is a table header)|\n")
            new_top.append(
                "|----------------------------------------------"
                "--------------------------------|\n")
            # Add Fluid 2 entry
            new_top.append(rxn_template.format(
                sp="Fluid 2".ljust(24),
                val="{:22.15E}".format(-1.0)))
            # Add reaction species in a stable order:
            # H2O first, then H+, then alphabetical, then O2(g) last
            ordered = []
            for sp in rxn_coeffs:
                if sp in ("H2O", "H+", "O2(g)"):
                    continue
                if rxn_coeffs[sp] == 0.0:
                    continue
                ordered.append(sp)
            ordered.sort()
            emit_order = ["H2O", "H+"] + ordered + ["O2(g)"]
            for sp in emit_order:
                if sp not in rxn_coeffs:
                    continue
                coeff = rxn_coeffs[sp]
                new_top.append(rxn_template.format(
                    sp=sp.ljust(24),
                    val="{:22.15E}".format(coeff)))
            new_top.append(
                "|----------------------------------------------"
                "--------------------------------|\n")
            skip_old_rxn = True
            continue
        if skip_old_rxn:
            if "|->|Surface area" in text:
                skip_old_rxn = False
                new_top.append(line)
            continue
        new_top.append(line)

    return new_top


def _rebuild_fluid2_block(lines_6i, speciation_pickup_top):
    """
    Replace the Composition and Reaction blocks in lines_6i (the top half)
    with those from a speciation's pickup top half.  This is used during
    cross-database mixing so that Fluid 2's Composition and Reaction blocks
    come from the target-database speciation rather than the cross-database
    Mixing_Fluid.

    The Prepare_Reaction options (temperature, pressure, Xi, mineral
    suppression, etc.) are preserved from lines_6i.
    """
    # Extract Composition and Reaction blocks from speciation pickup
    new_comp_lines = []
    new_rxn_lines = []
    in_comp = False
    in_rxn = False
    for line in speciation_pickup_top:
        if "|->|Composition" in line:
            in_comp = True
        if in_comp and "|->|Reaction" in line:
            in_comp = False
            in_rxn = True
        if in_rxn and "|->|Surface area" in line:
            in_rxn = False
            break
        if in_comp:
            new_comp_lines.append(line)
        elif in_rxn:
            new_rxn_lines.append(line)

    if not new_comp_lines or not new_rxn_lines:
        return lines_6i

    # Find and replace Composition+Reaction blocks in lines_6i
    result = []
    i = 0
    while i < len(lines_6i):
        text = lines_6i[i]
        if "|->|Composition" in text:
            # Replace from Composition through end of Reaction block
            result.extend(new_comp_lines)
            result.extend(new_rxn_lines)
            # Skip past original Composition+Reaction blocks
            while i < len(lines_6i):
                if "|->|Surface area" in lines_6i[i]:
                    break
                i += 1
            continue
        result.append(text)
        i += 1

    return result


def _parse_data0_elements(data0_text):
    """
    Parse the elements section of a data0 file and return the set of all
    element names (real and pseudo).
    """
    lines = data0_text.split("\n")
    elements = set()
    elem_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "elements":
            elem_start = i
        elif stripped == "basis species":
            break
    if elem_start is None:
        return elements
    for i in range(elem_start + 2, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith("+--"):
            continue
        if line == "basis species":
            break
        elements.add(line.split()[0])
    return elements


def _zero_composition_elements(top_lines, elements_to_zero):
    """
    Remove specified elements from the Composition block entirely.
    Elements absent from both fluids should not appear at all — even
    a zero-valued entry can cause EQ6 errors if the element isn't in
    the active data1 file.
    """
    new_top = []
    in_comp = False
    for line in top_lines:
        text = line.rstrip("\n")
        if "|->|Composition" in text:
            in_comp = True
        elif in_comp and "|->|Reaction" in text:
            in_comp = False

        if (in_comp and "|--->|" in text and "table header" not in text):
            parts = text.split("|")
            if len(parts) >= 4:
                elem = parts[2].strip()
                if elem in elements_to_zero:
                    continue

        new_top.append(line)
    return new_top


def _remap_composition_block(top_lines, target_elements, species_elements,
                             pseudo_elements):
    """
    Remap the Composition block in the top half so that its element names
    are compatible with the target data0 database.

    Elements that exist in the target database are kept as-is.  Elements
    that don't (e.g. plain ``C`` when the target uses ``Cjivp``/``Cjivn``)
    are replaced by the corresponding pseudo-elements.  The stoichiometric
    value is distributed among the replacement pseudo-elements by using
    the Reaction block coefficients of species that contain each
    pseudo-element.

    Parameters
    ----------
    top_lines : list of str
        Lines of the .6i top half.
    target_elements : set of str
        All element names present in the target data0 file.
    species_elements : dict
        Mapping of species name → list of element names in target data0.
    pseudo_elements : list of str
        Pseudo-element names in the target data0 (in data0 order).

    Returns
    -------
    list of str
        Modified top_lines with remapped Composition block.
    """
    from .redox_suppression import parse_pseudoelement_name

    pseudo_set = set(pseudo_elements)

    # --- Parse the Reaction block to get species→coefficient mapping --- #
    reaction_coeffs = {}
    in_reaction = False
    for line in top_lines:
        text = line.rstrip("\n")
        if "|->|Reaction" in text:
            in_reaction = True
            continue
        if in_reaction and "|->|Surface area" in text:
            break
        if in_reaction and "|--->|" in text and "table header" not in text:
            parts = text.split("|")
            if len(parts) >= 4:
                sp_name = parts[2].strip()
                try:
                    coeff = float(parts[3])
                    reaction_coeffs[sp_name] = coeff
                except ValueError:
                    pass

    # --- Build mapping: pseudo-element → sum of absolute reaction coeffs --- #
    pseudo_to_total = {}
    for sp, coeff in reaction_coeffs.items():
        for elem in species_elements.get(sp, []):
            if elem in pseudo_set:
                pseudo_to_total[elem] = (pseudo_to_total.get(elem, 0.0)
                                         + abs(coeff))

    # --- Build mapping: real element → list of replacement pseudo-elements --- #
    real_to_pseudos = {}
    for pe in pseudo_elements:
        real_elem, ox = parse_pseudoelement_name(pe)
        if ox is not None:
            if real_elem not in real_to_pseudos:
                real_to_pseudos[real_elem] = []
            real_to_pseudos[real_elem].append(pe)

    # --- Remap Composition block lines --- #
    new_top = []
    in_composition = False
    composition_done = False
    stoich_template = (
        "|--->|{elem}|{val}| (uesri(i,n), cesri(i,n))                |\n")

    for line in top_lines:
        text = line.rstrip("\n")

        if "|->|Composition" in text:
            in_composition = True
            new_top.append(line)
            continue

        if in_composition and "|->|Reaction" in text:
            in_composition = False
            composition_done = True
            new_top.append(line)
            continue

        if not in_composition or composition_done:
            new_top.append(line)
            continue

        # Inside Composition block
        if "|--->|" not in text or "table header" in text:
            new_top.append(line)
            continue

        # Parse element line
        parts = text.split("|")
        if len(parts) < 4:
            new_top.append(line)
            continue

        elem_name = parts[2].strip()
        try:
            elem_val = float(parts[3])
        except ValueError:
            new_top.append(line)
            continue

        if elem_name in target_elements:
            # Element exists in target — keep as-is
            new_top.append(line)
        elif elem_name in real_to_pseudos:
            # Real element needs splitting into pseudo-elements
            replacements = real_to_pseudos[elem_name]
            total_pseudo_moles = sum(
                pseudo_to_total.get(pe, 0.0) for pe in replacements)

            for pe in replacements:
                if total_pseudo_moles > 0:
                    fraction = pseudo_to_total.get(pe, 0.0) / total_pseudo_moles
                else:
                    fraction = 1.0 / len(replacements)
                pe_val = elem_val * fraction
                pe_padded = pe.ljust(8)
                val_str = "{:22.15E}".format(pe_val)
                new_top.append(stoich_template.format(
                    elem=pe_padded, val=val_str))
        else:
            # Check reverse: this might be a pseudo-element that needs
            # merging into a real element for the target database
            real_elem, ox = parse_pseudoelement_name(elem_name)
            if ox is not None and real_elem in target_elements:
                # This pseudo-element should be merged with others into
                # the real element.  Collect all pseudo-element lines for
                # the same real element, then emit one merged line.
                # We handle this by deferring — mark for merge.
                # For simplicity, just rename to the real element and
                # let duplicates accumulate; we'll merge below.
                re_padded = real_elem.ljust(8)
                val_str = "{:22.15E}".format(elem_val)
                new_top.append(stoich_template.format(
                    elem=re_padded, val=val_str))
            else:
                new_top.append(line)

    # --- Merge duplicate element entries (from pseudo→real collapse) --- #
    # Find the Composition block in new_top and merge duplicates
    merged_top = []
    in_comp = False
    comp_entries = []  # (elem_name, value, original_line)
    for line in new_top:
        text = line.rstrip("\n")
        if "|->|Composition" in text:
            in_comp = True
            merged_top.append(line)
            continue
        if in_comp and "|->|Reaction" in text:
            # Flush merged composition entries
            seen = {}
            order = []
            for en, ev, ol in comp_entries:
                if en is not None:
                    if en in seen:
                        seen[en] += ev
                    else:
                        seen[en] = ev
                        order.append(en)
                else:
                    # Non-element line (header/divider) — flush as-is
                    order.append(ol)
                    seen[ol] = None

            for key in order:
                if seen[key] is None:
                    merged_top.append(key)
                else:
                    padded = key.ljust(8)
                    val_str = "{:22.15E}".format(seen[key])
                    merged_top.append(stoich_template.format(
                        elem=padded, val=val_str))

            in_comp = False
            merged_top.append(line)
            continue

        if in_comp:
            if "|--->|" in text and "table header" not in text:
                parts = text.split("|")
                if len(parts) >= 4:
                    en = parts[2].strip()
                    try:
                        ev = float(parts[3])
                        comp_entries.append((en, ev, line))
                        continue
                    except ValueError:
                        pass
            comp_entries.append((None, 0, line))
            continue

        merged_top.append(line)

    return merged_top


def _compute_pseudo_element_totals(species_distribution, species_stoich):
    """
    Sum molality * stoichiometric_coefficient per pseudo-element across
    all aqueous species in the distribution table.

    Returns (totals, contributions) where:
      totals : dict  pseudo_element -> total_molality
      contributions : dict  pseudo_element -> list of
          (species_name, molality, stoich_coeff, contribution)
    """
    from .redox_suppression import parse_pseudoelement_name
    totals = {}
    contributions = {}
    for sp_name, molality in species_distribution.items():
        stoich = species_stoich.get(sp_name)
        if stoich is None:
            continue
        for elem, coeff in stoich.items():
            _, ox = parse_pseudoelement_name(elem)
            if ox is not None:
                contribution = molality * coeff
                totals[elem] = totals.get(elem, 0.0) + contribution
                if elem not in contributions:
                    contributions[elem] = []
                contributions[elem].append(
                    (sp_name, molality, coeff, contribution))
    return totals, contributions


def _reconcile_cross_database_mixing(top_lines, bottom_lines, data0_text,
                                     eq_output_text=None,
                                     fluid2_eq_output_text=None,
                                     bottom_half_matches_target=False,
                                     fluid2_basis=None):
    """
    Reconcile a .6i top and bottom half when they come from different
    thermodynamic databases.  Remaps the Composition block elements in
    the top half and adds missing basis species to the bottom half.

    Parameters
    ----------
    top_lines : list of str
        Lines of the .6i top half (will be modified with remapped elements).
    bottom_lines : list of str
        Lines of the .6i bottom half (will be modified with missing species).
    data0_text : str
        Full text of the target data0 database file (Fluid 1's database).
    eq_output_text : str, optional
        Full text of a .3o/.6o output file for oxidation-state distribution.
    fluid2_eq_output_text : str, optional
        Full text of the Mixing Fluid's .3o/.6o output file.
    bottom_half_matches_target : bool
        If True, the bottom half was produced by the same database as the
        target.  Existing mass balance totals are already correct and will
        not be overwritten; only missing species get new totals.
    fluid2_basis : list of str, optional
        Basis species names from the Mixing Fluid's pickup file.  When
        provided, species that appear in neither the bottom half's basis
        nor the mixing fluid's basis are excluded from the missing list
        (they were not present in either fluid and should not be created).

    Returns
    -------
    tuple of (list of str, list of str, dict)
        Modified (top_lines, bottom_lines, allocation_log).
    """
    from .redox_suppression import parse_pseudoelement_name

    species_elements = AqEquil._parse_data0_species_elements(data0_text)
    pseudo_elements = AqEquil._parse_data0_pseudo_elements(data0_text)
    target_elements = _parse_data0_elements(data0_text)

    # Determine which pseudo-elements are unwanted (their covering basis
    # species is absent from both fluids) BEFORE recomputing the
    # Composition block.  This prevents phantom H/O from being allocated
    # alongside species that will later be removed.
    current_basis = AqEquil._parse_bottom_basis_species(bottom_lines)
    unwanted_pes = set()
    if fluid2_basis is not None:
        either_fluid = set(current_basis) | set(fluid2_basis)
        for pe in pseudo_elements:
            covered_by = None
            for sp, elems in species_elements.items():
                if pe in elems:
                    covered_by = sp
                    break
            if covered_by is not None and covered_by not in either_fluid:
                unwanted_pes.add(pe)

    # --- Remap Composition block in top half --- #
    top_lines = _remap_composition_block(
        top_lines, target_elements, species_elements, pseudo_elements)

    # If we have the Mixing_Fluid's EQ3/6 output, recompute the
    # Composition block pseudo-element distribution from the species
    # distribution table (instead of relying on the source Reaction
    # block coefficients, which may use a different redox basis).
    if fluid2_eq_output_text is not None:
        target_stoich = AqEquil._parse_data0_all_species_stoichiometry(
            data0_text)
        fluid2_dist = AqEquil._parse_output_species_distribution(
            fluid2_eq_output_text)
        pe_totals, _ = _compute_pseudo_element_totals(
            fluid2_dist, target_stoich)
        # Remove unwanted pseudo-elements before applying to the
        # Composition block.  This ensures the O and H that would
        # accompany these species are never counted.
        for upe in unwanted_pes:
            pe_totals.pop(upe, None)
        if pe_totals:
            fluid2_num_comp = AqEquil._parse_output_numerical_composition(
                fluid2_eq_output_text)
            h2o_moles_f2 = fluid2_num_comp.get("H2O")
            if h2o_moles_f2 is None:
                h2o_moles_f2 = 55.51
            h2o_kg_f2 = h2o_moles_f2 * 18.01528 / 1000.0
            top_lines = _apply_pe_totals_to_composition(
                top_lines, pe_totals, h2o_kg_f2, h2o_moles_f2,
                target_elements, data0_text=data0_text)

    # Remove any unwanted pseudo-elements that may have been left by
    # _remap_composition_block (which runs before the pe_totals filter).
    if unwanted_pes:
        top_lines = _zero_composition_elements(
            top_lines, unwanted_pes)

    # --- Rebuild Reaction block to match remapped Composition --- #
    top_lines = _rebuild_reaction_block(top_lines, data0_text)

    output_composition = {}
    species_distribution = {}
    species_stoich = {}
    if eq_output_text is not None:
        output_composition = AqEquil._parse_output_numerical_composition(
            eq_output_text)
        species_distribution = AqEquil._parse_output_species_distribution(
            eq_output_text)
        species_stoich = AqEquil._parse_data0_all_species_stoichiometry(
            data0_text)

    reactant_species = AqEquil._parse_reactant_species(top_lines)

    # 1) Species referenced in the reactant but absent from the bottom half
    missing = [sp for sp in reactant_species if sp not in current_basis]

    # 2) Species needed to cover every pseudo-element in the target data0
    #    Only add a covering basis species when the pseudo-element actually
    #    carries nonzero mass in the Composition block.  If a pseudo-element
    #    has zero moles (e.g. because both fluids suppressed all species
    #    using that oxidation state), forcing in a basis species would
    #    fabricate mass that doesn't exist.
    comp_moles = {}
    in_comp = False
    for line in top_lines:
        text = line.rstrip("\n")
        if "|->|Composition" in text:
            in_comp = True
            continue
        if in_comp and "|->|Reaction" in text:
            break
        if in_comp and "|--->|" in text and "table header" not in text:
            parts = text.split("|")
            if len(parts) >= 4:
                try:
                    comp_moles[parts[2].strip()] = float(parts[3])
                except ValueError:
                    pass

    current_set = set(current_basis) | set(missing)
    pseudo_set = set(pseudo_elements)
    covered_pseudoelems = set()
    for sp in current_set:
        for elem in species_elements.get(sp, []):
            if elem in pseudo_set:
                covered_pseudoelems.add(elem)

    for pe in pseudo_elements:
        if pe in covered_pseudoelems:
            continue
        if abs(comp_moles.get(pe, 0.0)) == 0.0:
            continue
        for sp, elems in species_elements.items():
            if pe in elems and sp not in current_set:
                missing.append(sp)
                current_set.add(sp)
                for e in elems:
                    if e in pseudo_set:
                        covered_pseudoelems.add(e)
                break

    # Filter out species that don't exist in either fluid's basis.
    # current_basis is fluid 1's basis (the bottom half).  fluid2_basis
    # is the mixing fluid's basis.  If a species is in neither, it was
    # not present in either fluid and should not be fabricated.
    if fluid2_basis is not None:
        either_fluid = set(current_basis) | set(fluid2_basis)
        missing = [sp for sp in missing if sp in either_fluid]

    aq_totals = AqEquil._parse_bottom_mass_balance_totals(
        bottom_lines, column="aqueous")
    h2o_aq_moles = aq_totals.get("H2O")

    if not missing:
        # No missing species — recompute existing totals only if the
        # bottom half came from a different database than the target.
        existing_overrides = {}
        if species_distribution and not bottom_half_matches_target:
            existing_overrides = _recompute_mass_balance_totals(
                species_distribution, data0_text, current_basis,
                h2o_aq_moles)
        if existing_overrides:
            bottom_lines = _override_mass_balance_totals(
                bottom_lines, existing_overrides)
        override_log = {}
        for sp, new_val in existing_overrides.items():
            old_val = aq_totals.get(sp)
            override_log[sp] = {
                "method": "dissociation_recompute",
                "old_moles": old_val,
                "new_moles": new_val,
            }
        return top_lines, bottom_lines, override_log

    final_basis = list(current_basis) + missing

    # Recompute mass balance totals using the species distribution from
    # the source output and dissociation reactions from the target data0.
    # When the bottom half matches the target, only compute for missing
    # species (existing totals are already correct).  When the databases
    # differ, recompute everything.
    all_overrides = {}
    if species_distribution:
        recompute_basis = (final_basis if not bottom_half_matches_target
                           else missing)
        all_overrides = _recompute_mass_balance_totals(
            species_distribution, data0_text, recompute_basis, h2o_aq_moles)

    # Apply overrides for existing basis species to bottom half
    current_set = set(current_basis)
    existing_overrides = {sp: v for sp, v in all_overrides.items()
                         if sp in current_set}
    if existing_overrides:
        bottom_lines = _override_mass_balance_totals(
            bottom_lines, existing_overrides)

    override_log = {}
    for sp, new_val in existing_overrides.items():
        old_val = aq_totals.get(sp)
        override_log[sp] = {
            "method": "dissociation_recompute",
            "old_moles": old_val,
            "new_moles": new_val,
        }

    all_species_elems = {}
    for sp in final_basis:
        all_species_elems[sp] = species_elements.get(sp, [])

    missing_flags = {}
    for sp in missing:
        sp_elems = all_species_elems.get(sp, [])
        needs_moles = False
        for elem in sp_elems:
            _, ox = parse_pseudoelement_name(elem)
            if ox is None:
                continue
            covered = False
            for other_sp in final_basis:
                if other_sp == sp:
                    continue
                if elem in all_species_elems.get(other_sp, []):
                    covered = True
                    break
            if not covered:
                needs_moles = True
                break
        if needs_moles:
            missing_flags[sp] = "Moles"
        else:
            has_shared = False
            for elem in sp_elems:
                for other_sp in final_basis:
                    if other_sp == sp:
                        continue
                    if elem in all_species_elems.get(other_sp, []):
                        has_shared = True
                        break
                if has_shared:
                    break
            missing_flags[sp] = "Make non-basis" if has_shared else "Moles"

    missing_totals = {}
    donor_adjustments = {}
    allocation_log = {}

    # Use dissociation-recomputed totals when available and nonzero;
    # fall back to _compute_missing_totals for the rest.  This logic
    # is the same regardless of whether the bottom half matches the
    # target database — the difference is only in whether *existing*
    # basis species totals get overwritten (handled above).
    missing_with_recomputed = {
        sp for sp in missing if sp in all_overrides
        and all_overrides[sp] is not None
        and abs(all_overrides[sp]) > 0}
    missing_needing_fallback = [sp for sp in missing
                                if sp not in missing_with_recomputed]

    for sp in missing_with_recomputed:
        missing_totals[sp] = all_overrides[sp]
        allocation_log[sp] = {
            "method": "dissociation_recompute",
            "old_moles": None,
            "new_moles": all_overrides[sp],
        }

    if missing_needing_fallback:
        existing_totals = AqEquil._parse_bottom_mass_balance_totals(
            bottom_lines)
        fb_totals, fb_donors, fb_log = (
            AqEquil._compute_missing_totals(
                missing_needing_fallback, current_basis,
                all_species_elems, existing_totals, output_composition,
                species_distribution=species_distribution,
                species_stoich=species_stoich,
                h2o_aq_moles=h2o_aq_moles,
                comp_moles=comp_moles))
        missing_totals.update(fb_totals)
        allocation_log.update(fb_log)
        if not bottom_half_matches_target:
            for sp in all_overrides:
                fb_donors.pop(sp, None)
            donor_adjustments.update(fb_donors)

    for sp in allocation_log:
        allocation_log[sp]["flag"] = missing_flags.get(sp)

    bottom_lines = AqEquil._insert_missing_basis_species(
        bottom_lines, missing, missing_flags, missing_totals,
        donor_adjustments)

    allocation_log.update(override_log)

    return top_lines, bottom_lines, allocation_log


def _recompute_mass_balance_totals(species_distribution, data0_text,
                                   current_basis, h2o_aq_moles):
    """
    Recompute mass balance totals for all existing basis species using
    the species distribution from the source EQ3/6 output and the
    dissociation reactions from the target data0 file.

    For each aqueous species *i* with molality *m_i*, its dissociation
    reaction gives coefficients *a_ib* for each basis species *b*.
    The mass balance total for *b* is::

        MBT(b) = sum_i(a_ib * m_i) * kg_water

    H2O is special: the species distribution table only lists solute
    species, so H2O's own contribution (h2o_aq_moles) must be added
    directly to MBT(H2O).

    Returns a dict of {basis_species_name: new_total_moles}.
    """
    dissoc = AqEquil._parse_data0_dissociation_reactions(data0_text)

    if not dissoc or h2o_aq_moles is None or h2o_aq_moles <= 0:
        return {}

    kg_water = h2o_aq_moles * 18.01528 / 1000.0

    basis_set = set(current_basis)
    mbt = {b: 0.0 for b in current_basis}

    for sp_name, molality in species_distribution.items():
        rxn = dissoc.get(sp_name)
        if rxn is None:
            continue
        self_coeff = rxn.get(sp_name, None)
        if self_coeff is not None and self_coeff > 0:
            # Basis species: identity reaction {self: 1.0}
            if sp_name in basis_set:
                mbt[sp_name] += molality
            continue
        # Non-basis species: normalize by |self_coeff| (usually 1.0)
        norm = abs(self_coeff) if self_coeff is not None else 1.0
        if norm == 0:
            norm = 1.0
        for basis_sp, coeff in rxn.items():
            if basis_sp == sp_name:
                continue
            if basis_sp in basis_set:
                mbt[basis_sp] += (coeff / norm) * molality

    overrides = {}
    for b in current_basis:
        overrides[b] = mbt[b] * kg_water

    if "H2O" in overrides:
        overrides["H2O"] += h2o_aq_moles

    return overrides


def _override_mass_balance_totals(bottom_lines, overrides):
    """
    Replace both the equilibrium and aqueous mass balance totals for
    the species listed in *overrides* (a dict of species_name → new_moles).
    """
    new_lines = []
    in_totals = False
    for line in bottom_lines:
        text = line.rstrip("\n")
        if "Mass Balance Totals (moles)" in text:
            in_totals = True
        elif in_totals and "Electrical imbalance" in text:
            in_totals = False

        if (in_totals and text.startswith("|") and "Aqueous" in text
                and "Electrical imbalance" not in text
                and "(ubmtbi(n))" not in text and "(mtbi(n))" not in text):
            sp_name = text[1:25].strip()
            if sp_name in overrides:
                val = overrides[sp_name]
                row = "|{:<24s}Aqueous |{:22.15E}|{:22.15E}|\n".format(
                    sp_name, val, val)
                new_lines.append(row)
                continue

        new_lines.append(line)
    return new_lines


def _warn_if_collapsed_totals(bottom_lines):
    """
    Detect when a pickup file's Equilibrium System totals equal its Aqueous
    Solution totals for every species — a sign that clear_es_solids_end=True
    was used, which destroys mineral mass contributions.  This produces
    incorrect results in cross-database mixing because the total H+ and O2(g)
    budgets no longer account for mineral formation/dissolution.

    Skips the warning if kbt == kmt (no minerals in the equilibrium system),
    since Eq = Aq is expected when no minerals are present (e.g. a direct
    EQ3 speciation with no prior reaction).
    """
    import re, warnings

    kbt = kmt = None
    for line in bottom_lines:
        text = line.rstrip("\n")
        if "(kbt)" in text:
            m = re.search(r'\|\s*(\d+)\|', text)
            if m:
                kbt = int(m.group(1))
        elif "(kmt)" in text:
            m = re.search(r'\|\s*(\d+)\|', text)
            if m:
                kmt = int(m.group(1))
            break
    if kbt is not None and kmt is not None and kbt == kmt:
        return

    in_totals = False
    all_equal = True
    n_checked = 0
    for line in bottom_lines:
        if "Mass Balance Totals" in line:
            in_totals = True
            continue
        if in_totals and "Electrical imbalance" in line:
            break
        if not in_totals:
            continue
        parts = line.split("|")
        if len(parts) >= 4:
            eq_str = parts[2].strip()
            aq_str = parts[3].strip()
            try:
                eq_val = float(eq_str)
                aq_val = float(aq_str)
                n_checked += 1
                if eq_val != aq_val:
                    all_equal = False
                    break
            except ValueError:
                continue
    if all_equal and n_checked > 3:
        warnings.warn(
            "Cross-database mixing: the Mixing_Fluid's pickup file has "
            "identical Equilibrium System and Aqueous Solution totals for "
            "all species. This usually means clear_es_solids_end=True was "
            "used in the reaction that produced this fluid, which removes "
            "mineral mass contributions from the mass balance. This will "
            "produce incorrect pH and redox results during mixing. To fix "
            "this, set clear_es_solids_end=False in the Prepare_Reaction "
            "that creates the Mixing_Fluid's source speciation, and use "
            "strip_minerals=True on the Mixing_Fluid instead.",
            UserWarning,
            stacklevel=4,
        )


def _strip_minerals_from_bottom_half(bottom_lines):
    """
    Remove all mineral entries (pure minerals and solid solution end-members)
    from the bottom half of a .6p/.6i file, update the matrix index limits,
    and set equilibrium-system mass balance totals equal to aqueous-solution
    totals so that the removed mineral mass is not spuriously re-dissolved.
    """
    import re

    # First pass: find kbt so we know the target value for kmt, kxt, kdim
    # after all minerals are removed.
    kbt_value = None
    for line in bottom_lines:
        if "(kbt)" in line:
            m = re.search(r'\|\s*(\d+)\|', line)
            if m:
                kbt_value = int(m.group(1))
            break

    if kbt_value is None:
        return bottom_lines

    new_lines = []
    in_column_vars = False
    in_totals = False

    for line in bottom_lines:
        text = line.rstrip("\n")

        if "Matrix Column Variables and Values" in text:
            in_column_vars = True
        elif "Phases and Species in the PRS" in text:
            in_column_vars = False

        if "Mass Balance Totals (moles)" in text:
            in_totals = True
        elif in_totals and "Electrical imbalance" in text:
            in_totals = False

        # Set kmt, kxt, and kdim all equal to kbt (no minerals remain)
        if "(kmt)" in text or "(kxt)" in text or "(kdim)" in text:
            new_lines.append(_set_matrix_index(text, kbt_value))
            continue

        # In the Mass Balance Totals section, set equilibrium total = aqueous
        if (in_totals and text.startswith("|") and "Aqueous" in text
                and "Electrical imbalance" not in text
                and "(ubmtbi(n))" not in text and "(mtbi(n))" not in text):
            parts = text.split("|")
            if len(parts) >= 4:
                sp_name = parts[1][:24]
                aq_str = parts[3].strip()
                try:
                    aq_val = float(aq_str)
                    row = "|{:<24s}Aqueous |{:22.15E}|{:22.15E}|\n".format(
                        sp_name.strip(), aq_val, aq_val)
                    new_lines.append(row)
                    continue
                except ValueError:
                    pass

        # In the Matrix Column Variables section, drop non-aqueous entries
        # but keep the table header line containing "(zvclgi(n))"
        if (in_column_vars and text.startswith("|") and "| --   |" in text
                and "(zvclgi(n))" not in text):
            if "Aqueous solution" not in text:
                continue

        new_lines.append(line)

    return new_lines


def _set_matrix_index(line_text, new_val):
    """
    Set a matrix index value (e.g., kbt, kmt, kxt, kdim) to an absolute value.
    """
    import re
    match = re.search(r'\|\s*(\d+)\|', line_text)
    if match:
        inner = match.group(0)[1:-1]  # strip outer |
        width = len(inner)
        new_field = "|{:>{w}}|".format(new_val, w=width)
        new_line = line_text[:match.start()] + new_field + line_text[match.end():]
        return new_line + "\n"
    return line_text + "\n"


def _update_matrix_index(line_text, n_add):
    """
    Increment a matrix index value (e.g., kbt, kmt, kxt, kdim) by n_add.
    """
    import re
    match = re.search(r'\|\s*(\d+)\|', line_text)
    if match:
        old_val = int(match.group(1))
        new_val = old_val + n_add
        # Preserve the field width
        field = match.group(0)
        # The format is |   NN| with fixed width
        inner = match.group(0)[1:-1]  # strip outer |
        width = len(inner)
        new_field = "|{:>{w}}|".format(new_val, w=width)
        new_line = line_text[:match.start()] + new_field + line_text[match.end():]
        return new_line + "\n"
    return line_text + "\n"


def compare(*args):
    
    """
    Combine two or more speciations into a single speciation object for
    comparison. The speciation object returned by this function can produce
    scatterplots, barplots, and mass contribution plots, and contains a report
    that can be browsed with `lookup`. See documentation for the functions in
    the `Speciation` class for more detail.

    Parameters
    ----------
    *args : two or more objects of class `Speciation` to compare

    Returns
    ----------
    An object of class `Speciation`.
    """

    if all(["mass_contribution" in a.__dict__.keys() for a in args]):
        allow_mass_contribution = True
        mass_contribution_breaks = []
    else:
        allow_mass_contribution = False
    
    for i,sp in enumerate(args):
        if i == 0:
            sp_total = copy.deepcopy(sp)
            sp_total.sample_data_combined = {}
            sp_total.sample_data_combined[i] = sp_total.sample_data
            sp_total.sample_data = {}
            if allow_mass_contribution:
                mass_contribution_breaks.append(0)
        else:
            sp_i = copy.deepcopy(sp)
            sp_total.sample_data_combined[i] = sp_i.sample_data
            if allow_mass_contribution:
                mass_contribution_breaks.append(sp_total.report.shape[0])
            sp_total.report = pd.concat([sp_total.report, sp_i.report], axis=0, sort=False)
        

    sp_total.report.index = sp_total.report.index + ("_"+sp_total.report.groupby(level=0).cumcount().astype(str)).replace('_0','')

    # restore sample data and rename samples if necessary
    for i in sp_total.sample_data_combined.keys():
        for sample in sp_total.sample_data_combined[i].keys():
            if sample not in sp_total.sample_data.keys():
                sp_total.sample_data[sample] = sp_total.sample_data_combined[i][sample]
            else:
                prev_sample_name = copy.deepcopy(sample)
                ii = 0
                while sample in sp_total.sample_data.keys():
                    ii+=1
                    sample = prev_sample_name+"_"+str(ii)
                    if sample not in sp_total.sample_data.keys():
                        sp_total.sample_data[sample] = sp_total.sample_data_combined[i][prev_sample_name]
                        break
                    
    
    if allow_mass_contribution:
        mass_contribution_breaks.append(len(sp_total.report.index))
        mc_sample_names_with_suffixes = list(sp_total.report.index)
        for i,sp in enumerate(args):
            mc_i = copy.deepcopy(sp.mass_contribution)
            
            new_sample_names = copy.copy(mc_sample_names_with_suffixes[mass_contribution_breaks[i]:mass_contribution_breaks[i+1]])
            old_sample_names = list(args[i].sample_data.keys())
            
            old_new_sample_name_dict = {old:new for old,new in zip(old_sample_names, new_sample_names)}
                
            newsample = [old_new_sample_name_dict[old] for old in mc_i["sample"]]

            mc_i["sample"] = newsample
            
            if i == 0:
                mc_total = mc_i
            else:
                mc_total = pd.concat([mc_total, mc_i], axis=0, sort=False)
        
        sp_total.mass_contribution = mc_total
        
    else:
        def no_mass_contrib_message(*args, **kwargs):
            print("Mass contributions cannot be compared between these speciations "
                  "because one or more calculations lack mass contribution data.")
        sp_total.plot_mass_contribution = no_mass_contrib_message

    return sp_total


class Speciation(SpeciationGroups):
    
    """
    Stores the output of a speciation calculation.
    
    Parameters
    ----------
    args : dict
        Arguments inherited from class AqEquil.
    
    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this class?
        When True, error messages handled by this class will be short and to
        the point.
    
    Attributes
    ----------
    input : pd.DataFrame
        Pandas dataframe containing user-supplied sample chemistry data.
    
    mass_contribution : pd.DataFrame
        Pandas dataframe containing basis species contributions to mass balance
        of aqueous species.
    
    batch_3o : rpy2 ListVector
        An rpy2 ListVector (R object) containing speciation results, in case
        analysis in R is preferred.
    
    report : pd.DataFrame
        Pandas dataframe reporting major results of speciation calculation in
        across all samples.
    
    report_divs : dict
        An dictionary of column names within the different sections of the
        speciation report.
    
    sample_data : dict
        Dictionary with sample names as keys and speciation results as values.
    
    """
    
    # get functions from the AqEquil class
    _interpolate_logK = AqEquil._interpolate_logK 
    plot_logK_fit = AqEquil.plot_logK_fit
    
    def __init__(self, args, hide_traceback=True):
        self.err_handler = Error_Handler(clean=hide_traceback)
        
        self.reactions_for_plotting = None # stores formatted reactions for plotting results of affinity and energy supply calculations
        self.apply_redox_columns = None # stores a list of column names for output from apply_redox_reactions()
        self.redox_formatted_reactions = None # stores a dataframe of formatted affinity/energy supply reactions
        self.redox_reactions_table = None
        
        # stores a dictionary of speciation groups (e.g., "CO2", "HCO3-", "CO3-2"...)
        self.custom_grouping_filepath = None
        self.reactant_dict_scalar = None
        self.speciation_group_dict_unpacked = None
        
        for k in args:
            setattr(self, k, args[k])

        if 'verbose' not in list(self.__dict__.keys()):
            self.verbose = 1
        

    def __getitem__(self, item):
         return getattr(self, item)

    
    def add_new_half_reaction(self, oxidant, reductant, redox_couple=None):
        """
        Add a new half reaction to the bottom of the table of available half
        reactions that can be combined using `make_redox_reactions`.

        Parameters
        ----------
        oxidant, reductant : str
            Names of the oxidant and reductant as they appear in the
            thermodynamic database. It does not matter what speciated form of an
            oxidant or reductant is chosen. For example, choosing "acetate" as
            the oxidant is the same as choosing "acetic-acid" because they
            belong to the same speciation group.
        
        redox_couple : str, optional
            Desired name of the half reaction. This name is arbitrary and can be
            anything you'd like. If no name is provided, a name will be
            automatically generated based on the oxidant and reductant.

        Result
        ----------
        Appends a new row to the half_cell_reactions table.
        """

        if not isinstance(self.thermo.thermo_db, pd.DataFrame):
            # check to see if user is using a thermodynamic database CSV
            if self.verbose > 0:
                print("Did not add a new half reaction because a WORM-style thermodynamic database CSV is not available "
                      "to check the validity of the user-supplied oxidant and reductant or perform energy calculations. Perhaps you "
                      "are getting this message because you are using a calibrated database or data0 file (e.g. db='wrm' in AqEquil() ) "
                      "instead of a more complete CSV database (e.g. db='WORM' or db='wrm_data.csv' in AqEquil() ).")
            return
        else:
            # check to see if the oxidant and reductant are actually in the thermodynamic database
            if oxidant in list(self.thermo.thermo_db["name"]) and reductant in list(self.thermo.thermo_db["name"]):
                pass
            elif oxidant == "H2O" or reductant == "H2O":
                pass
            else:
                if oxidant not in list(self.thermo.thermo_db["name"]):
                    if self.verbose > 0:
                        print("Did not add a new half reaction because the oxidant '"+str(oxidant)+"' "
                              "was not found amongst the currently loaded chemical species.")
                else:
                    print("Did not add a new half reaction because the reductant '"+str(reductant)+"' "
                          "was not found amongst the currently loaded chemical species.")
                return
                
        # check to see if this half reaction already exists
        add_new_row = True
        dup_index = None
        for i in self.half_cell_reactions.index:
            if self.half_cell_reactions.iloc[i]["Oxidant"] == oxidant and self.half_cell_reactions.iloc[i]["Reductant"] == reductant:
                add_new_row = False
                dup_index = i
                break

        if add_new_row:

            if not isinstance(redox_couple, str):
                redox_couple = str(oxidant) + " to " + str(reductant)
            
            new_row = [redox_couple, oxidant, reductant]
            self.half_cell_reactions.loc[len(self.half_cell_reactions.index)] = new_row
            if self.verbose > 0:
                print("Added new half reaction to half_cell_reactions.")
        else:
            if self.verbose > 0:
                print("A row with this oxidant and reductant already exists with index", i)
    

    def delete_half_reaction(self, index):
        """
        Delete a half reaction from the half_cell_reactions table. Useful for
        removing unwanted or mistaken half reactions. This will reset the
        indices in the table to prevent a gap from the deleted row.

        Parameters
        ----------
        index : int or list of int
            Index or list of indices of the row in the half_cell_reactions table
            to be deleted.
        """
        if not isinstance(index, list):
            index = [index]
        self.half_cell_reactions.drop(index, inplace=True)
        self.half_cell_reactions.reset_index(inplace=True, drop=True)
        if self.verbose > 0:
                print("Rows", index, "have been deleted and indices have "
                      "been reset in the half_cell_reactions table.")

    
    def reset_half_reactions(self):
        """
        Reset the half_cell_reactions table back to its original default. Takes
        no parameters.
        """
        self.half_cell_reactions = copy.deepcopy(self._half_cell_reactions_original_copy)
        if self.verbose > 0:
                print("The table of half reactions half_cell_reactions has been reset.")

    
    def apply_redox_reactions(self, y_type="E", y_units="cal", limiting=None,
                                    grams_minerals=None,
                                    negative_energy_supplies=False,
                                    custom_grouping_filepath=None,
                                    append_report=True):

        """
        Calculate chemical affinities, energy supplies, and more in samples
        using the redox reactions generated by the function
        `make_redox_reactions`.

        Parameters
        ----------
        grams_minerals : float, int, dict, or pandas.DataFrame, optional
            Number of grams belonging to each mineral reactant when calculating
            the limiting reactant during an energy supply calculation. This
            parameter is only used when `y_type="E"`.
            For example, in the reaction
            4 goethite + 4 iron + 3 O2 = 4 magnetite + 2 H2O
            setting `grams_minerals = 0.001` would mean 0.001 grams of goethite
            is reacting with 0.001 grams of iron. If it is desirable to specify
            individual masses for each mineral reactant, then a dictionary can
            be provided. For example:
            `grams_minerals={"goethite": 0.001, "iron": 0.1},`
            To vary the mass of a mineral from sample to sample, provide a
            dataframe with a column per mineral and a row per sample.
            By default, mineral masses are left undefined and minerals are not
            considered when looking for a limiting reactant. Note that defining
            a mass of 0 grams is different: a mineral present at 0 grams is
            immediately the limiting reactant, giving an energy supply of 0.

        negative_energy_supplies : bool, default False
            Report negative energy supplies? If False, negative energy supplies
            are reported as 0. If True, negative energy supplies are
            reported. A 'negative energy supply' represents the energy cost of
            depleting the limiting reactant of a reaction. This metric is not
            always helpful when examing energy supply results, so this option is
            set to False by default.

        y_type : str, default "A"
            The variable to plot on the y-axis. Can be either 'A' (for chemical
            affinity), 'G' (for Gibbs free energy, ΔG), 'logK' (for the log
            of the equilibrium constant), 'logQ' (for the log of the reaction
            quotient), or 'E' for energy supply.

        y_units : str, default "kcal"
            The unit that energy will be reported in (per mol for G and A, or
            per kg fluid for energy supply, or unitless for logK and logQ).
            Can be 'kcal', 'cal', 'J', or 'kJ'.

        limiting : str, optional
            Name of the species to act as the limiting reactant when calculating
            energy supply. If this parameter is left undefined, then a
            limiting reactant will be chosen automatically based on
            concentration and stoichiometry. This parameter is ignored unless
            `y_type` is set to 'E' (energy supply).

        append_report : bool, default True
            Add or update calculated values to the speciation object's report
            dataframe?

        custom_grouping_filepath : str, optional
            Filepath for a TXT file containing customized speciation groups. Use
            to override the built-in speciation group file.

        raise_nonlimiting_exception : bool, default True
            This parameter can be ignored in almost all cases. Raise an
            exception when there are no available limiting reactants?
            The purpose of this parameter is toggle off error message
            interruptions when this function is called by
            `apply_redox_reactions`, which can test many different reactions at
            once, some of which do not have valid limiting reactants and would
            otherwise be interrupted by errors.
        
            
        Returns
        ----------
        Pandas dataframe
            Returns a multiindexed dataframe of samples and calculated results.
            If `append_report` is True, then the report attribute of the
            speciation object will be appended/updated.
        """
        
        self._set_speciation_groups(custom_grouping_filepath, reset=True)

        y_name_list = []
        val_list_list = []
        result_dict = {}
        result_lim_dict = {}

        if not isinstance(self.redox_reactions_table, pd.DataFrame):
            self.err_handler.raise_exception("Redox reactions cannot be applied "
                "because they have not yet been generated. Generate redox "
                "reactions with the function make_redox_reactions() and then "
                "try again.")
        
        coeff_colnames = [c for c in list(self.redox_reactions_table.columns) if "coeff_" in c]
        species_colnames = [c for c in list(self.redox_reactions_table.columns) if "species_" in c]

        # handle the progress bar
        max_count = len(list(self.redox_reactions_table.index))
        f = IntProgress(min=0, max=max_count) # instantiate the bar
        display(f) # display the bar

        df_rxn_list = []
        for i,rxn in enumerate(list(self.redox_reactions_table.index)):
        
            coeff_list = list(self.redox_reactions_table[coeff_colnames].loc[rxn])
            coeff_list = [v for v in coeff_list if not math.isnan(v)]
            coeff_list = [float(v) for v in coeff_list]
        
            species_list = list(self.redox_reactions_table[species_colnames].loc[rxn])
            species_list = [v for v in species_list if isinstance(v, str)]
            
            if len(coeff_list) != len(species_list):
                self.err_handler.raise_exception("There is a mismatch between "
                        "the number of coefficients and number of species in "
                        "the reaction. Coefficients: "+str(coeff_list)+" "
                        "Species: "+str(species_list))
        
            if y_type in ["A", "G"]:
                divisor = self.redox_reactions_table["mol_e-_transferred_per_mol_rxn"].loc[rxn]
                divisor = float(divisor)
            else:
                divisor = 1
        
            redox_pair = self.redox_reactions_table["redox_pairs"].loc[rxn]

            if y_type == "E":
                limiting_input = self._switch_limiting(limiting,
                                                  stoich=coeff_list,
                                                  species=species_list,
                                                  lenient=True)
            else:
                limiting_input = None
            
            df_rxn = self.calculate_energy(
                            species=species_list,
                            stoich=coeff_list,
                            divisor=divisor,
                            per_electron=True,
                            grams_minerals=grams_minerals,
                            y_type=y_type,
                            y_units=y_units,
                            limiting=limiting_input,
                            raise_nonlimiting_exception=False,
                            rxn_name=rxn,
                            append_report=append_report,
                            negative_energy_supplies=negative_energy_supplies,
                            )

            df_rxn_list.append(df_rxn)
            f.value += 1 # tick the counter

        df = pd.concat(df_rxn_list, axis=1)

        return df
    

    def show_redox_reactions(self, formatted=True,
                                   charge_sign_at_end=False,
                                   show=True):
        
        """
        Show a table of redox reactions generated with the function
        `make_redox_reactions`.
        
        Parameters
        ----------
        formatted : bool, default True
            Should reactions be formatted for html output?
            
        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-)? Ignored if
            `formatted` is False.
        
        show : bool, default False
            Show the table of reactions? Ignored if not run in a Jupyter
            notebook.
        
        Returns
        ----------
        A pandas dataframe containing balanced redox reactions written in full.
        """

        if isinstance(self.redox_reactions_table, pd.DataFrame):
            self.redox_formatted_reactions = copy.copy(self.redox_reactions_table.iloc[:, 0:1])
        else:
            self.err_handler.raise_exception("There are no redox reactions to display. "
                    "Try running make_redox_reactions() first.")
        
        df = copy.copy(self.redox_reactions_table)
        
        reactions = []
        for irow in range(0, df.shape[0]):
            redox_pair = df.loc[self.redox_reactions_table.index[irow], "redox_pairs"]

            oxidant = redox_pair[0]
            reductant = redox_pair[1]

            rxn_row = df.iloc[irow, 2:]
            rxn = rxn_row[rxn_row.notna()]
            coeffs = copy.copy(rxn[::2]).tolist()
            names = copy.copy(rxn[1::2]).tolist()

            reaction = self.format_reaction(coeffs=coeffs,
                                            names=names,
                                            formatted=formatted,
                                            charge_sign_at_end=charge_sign_at_end,
                                            show=False)
            
            reactions.append(reaction)
    
        self.redox_formatted_reactions["reaction"] = reactions
        

        df_out = copy.copy(self.redox_formatted_reactions)

        if isnotebook() and show:
            display(HTML(df_out.to_html(escape=False)))
        
        return df_out

    
    @staticmethod
    def format_reaction(coeffs, names, formatted=True,
                        charge_sign_at_end=True, show=True):
        
        react_grid = pd.DataFrame({"coeff":coeffs, "name":names})
        react_grid["coeff"] = pd.to_numeric(react_grid["coeff"])
        react_grid = react_grid.astype({'coeff': 'float'})

        reactants = " + ".join([(str(-int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else -react_grid["coeff"][i])+" " if -react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
        products = " + ".join([(str(int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else react_grid["coeff"][i])+" " if react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
        if formatted:
            reactants = " + ".join([format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
            products = " + ".join([format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
        reaction = reactants + " = " + products

        if isnotebook() and show:
            display(HTML(reaction))
        
        return reaction

    
    def make_redox_reactions(self,
                             idx_list="all",
                             show=True,
                             formatted=True,
                             charge_sign_at_end=False):

        """
        Generate an organized collection of redox reactions for calculating
        chemical affinity and energy supply values during speciation.
        
        Parameters
        ----------
        idx_list : list of int or "all", default "all"
            List of indices of half reactions in the half cell reaction table
            to be combined when generating full redox reactions.
            E.g. [0, 1, 4] will combine half reactions with indices 0, 1, and 4
            in the table stored in the `half_cell_reactions` attribute of the
            `Speciation` class.
            If "all", generate all possible redox reactions from available half
            cell reactions.

        show : bool, default True
            Show the table of reactions? Ignored if not run in a Jupyter
            notebook.

        formatted : bool, default True
            Should reactions be formatted for html output? Ignored if `show` is
            False.
            
        charge_sign_at_end : bool, default True
            Display charge with sign after the number (e.g. SO4 2-)? Ignored if
            `formatted` is False.
        
        Returns
        ----------
        Output is stored in the `redox_reactions_table` and
        `redox_formatted_reactions` attributes of the `Speciation` class.
        """
        
        # reset all redox variables stored in the AqEquil class
        self.redox_reactions_table = None
        self.redox_formatted_reactions = None

        if isinstance(idx_list, str):
            if idx_list == "all":
                pass
            else:
                self.err_handler.raise_exception("The parameter idx_list must "
                        "either be a list of indices corresponding to desired "
                        "half reactions in the half_cell_reactions table, or "
                        "'all' for all half reaction combinations.")
        elif isinstance(idx_list, list):
            if len(idx_list) == 0:
                self.err_handler.raise_exception("The list of desired half "
                        "reactions is empty. Please define the indices of at "
                        "least two half reactions to combine into a full redox "
                        "reaction. These indices correspond to rows in the "
                        "half_cell_reactions table.")
            if len(idx_list) == 1:
                self.err_handler.raise_exception("Only one half reaction index "
                        "was provided. At least two must be provided in order "
                        "to create a full redox reaction.")
        else:
            self.err_handler.raise_exception("The parameter idx_list must "
                    "either be a list of indices corresponding to desired "
                    "half reactions in the half_cell_reactions table, or "
                    "'all' for all half reaction combinations.")
        
        if self.verbose > 1:
            print("Generating redox reactions...")

        err_msg = ("redox_pairs can either be 'all' or a list of integers "
               "indicating the indices of half cell reactions in "
               "the half_cell_reactions table that should be combined into "
               "full redox reactions. For example, redox_pairs=[0, 1, 2, 6] "
               "will combine half cell reactions with indices 0, 1, 2, and 6 in "
               "the half_cell_reactions table. This table is an attribute in the "
               "class AqEquil.")
        if isinstance(idx_list, str):
            if idx_list == "all":
                idx_list = list(range(0, self.half_cell_reactions.shape[0]))
            else:
                self.err_handler.raise_exception(err_msg)
        elif isinstance(idx_list, list):
            if not all([isinstance(i, int) for i in idx_list]):
                self.err_handler.raise_exception(err_msg)
        else:
            self.err_handler.raise_exception(err_msg)
        
        self.idx_list = idx_list

        half_reaction_dict = {}
        bad_idx_list = []
        for idx in idx_list:
        
            #print(idx)

            if idx not in self.half_cell_reactions.index:
                self.err_handler.raise_exception("The index " + str(idx) + " was not "
                        "found in the half_cell_reactions table.")
            
            oxidant = self.half_cell_reactions["Oxidant"].iloc[idx]
            reductant = self.half_cell_reactions["Reductant"].iloc[idx]

            if oxidant == "O2" and reductant == "H2O":
                half_reaction_dict[idx] = {'O2': -1.0, 'e-': -4.0, 'H+': -4.0, 'H2O': 2.0}
                continue
            elif oxidant == "H2O" and reductant == "H2":
                half_reaction_dict[idx] = {'H2O': -2.0, 'e-': -2.0, 'H2': 1.0, 'OH-': 2.0}
                continue
            elif oxidant == "H2O2" and reductant == "H2O":
                half_reaction_dict[idx] = {'H2O2': -1.0, 'e-': -2.0, 'H2O': 2.0, "H+": -2.0}
                continue
            elif oxidant == "O2" and reductant == "H2O2":
                half_reaction_dict[idx] = {'O2': -1.0, 'e-': -2.0, 'H+': -2.0, 'H2O2': 1.0}
                continue
            elif oxidant == "H2O2" and reductant == "H2":
                half_reaction_dict[idx] = {'H2O2': -1.0, 'e-': -4.0, 'H+': -2.0, 'H2':1.0, 'OH-': 2.0}
                continue

            db_sp_names = list(self.thermo.thermo_db["name"])
            if oxidant == "H2O" or reductant == "H2O":
                pass
            elif oxidant not in db_sp_names or reductant not in db_sp_names:
                if oxidant not in db_sp_names:
                    problem_sp = oxidant
                else:
                    problem_sp = reductant
                if self.verbose > 0:
                    print("The species '"+str(problem_sp)+"' found in half cell reaction",
                          self.half_cell_reactions.iloc[idx][0], "( index", idx, ")",
                          "was not found in the thermodynamic database",
                          "used by the speciation. Check whether this species has",
                          "been excluded from the thermodynamic database prior to",
                          "the speciation step. Skipping this half reaction...")
                bad_idx_list.append(idx)
                continue
            
            # comparing common elements
            # requires "formula" and not "formula_modded" columns of thermo_db
            ox_formula_dict = parse_formula(self.thermo.thermo_db["formula"].loc[self.thermo.thermo_db["name"]==oxidant].values[0])
            ox_formula_ox = self.thermo.thermo_db["formula_ox"].loc[self.thermo.thermo_db["name"]==oxidant].values[0]
            ox_dissrxn = self.thermo.thermo_db["dissrxn"].loc[self.thermo.thermo_db["name"]==oxidant].values[0]
                
            
            red_formula_dict = parse_formula(self.thermo.thermo_db["formula"].loc[self.thermo.thermo_db["name"]==reductant].values[0])
            red_formula_ox = self.thermo.thermo_db["formula_ox"].loc[self.thermo.thermo_db["name"]==reductant].values[0]
            red_dissrxn = self.thermo.thermo_db["dissrxn"].loc[self.thermo.thermo_db["name"]==reductant].values[0]


            common_keys = []
            for key in list(ox_formula_dict.keys()):
                if key in list(red_formula_dict.keys()):
                    common_keys.append(key)
            
            ZOH_list = ["+", "-", "O", "H"]
            
            common_elems = [k for k in common_keys if k not in ZOH_list]
            elems_unique_to_oxidant = [k for k in list(ox_formula_dict.keys()) if k not in list(red_formula_dict.keys())+ZOH_list]
            elems_unique_to_reductant = [k for k in list(red_formula_dict.keys()) if k not in list(ox_formula_dict.keys())+ZOH_list]
            
            # print("common_elems")
            # print(common_elems)
            # print("unique to oxidant")
            # print(elems_unique_to_oxidant)
            # print("unique to reductant")
            # print(elems_unique_to_reductant)

            ox_red_formula_ox_dict = {}
            for i,ro in enumerate([ox_formula_ox, red_formula_ox]):
                ox_red_formula_ox_dict[[oxidant, reductant][i]] = self._formula_ox_to_dict(ro)

            common_elem_electron_dict = {}
            for ce in common_elems:
                ox_ce = ox_red_formula_ox_dict[list(ox_red_formula_ox_dict.keys())[0]][ce]
                red_ce = ox_red_formula_ox_dict[list(ox_red_formula_ox_dict.keys())[1]][ce]
            
                ox_ce_n = sum([v for v in list(ox_ce.values())])
                red_ce_n = sum([v for v in list(red_ce.values())])
            
                if ox_ce_n/red_ce_n < 1:
                    ox_coeff = red_ce_n/ox_ce_n
                    red_coeff = 1.0
                elif ox_ce_n/red_ce_n > 1:
                    ox_coeff = 1.0
                    red_coeff = ox_ce_n/red_ce_n
                else:
                    # oxidant and reductant have the same number of common element
                    ox_coeff = 1.0
                    red_coeff = 1.0
                
                total_ox_ce_oxstate = sum([k*v for k,v in zip(list(ox_ce.keys()), list(ox_ce.values()))])
                total_red_ce_oxstate = sum([k*v for k,v in zip(list(red_ce.keys()), list(red_ce.values()))])
                electron_coeff = red_coeff*total_red_ce_oxstate - ox_coeff*total_ox_ce_oxstate
                
                if electron_coeff > 0:
                    self.err_handler.raise_exception("The oxidant '"+oxidant+"' is more reduced than the "
                            "reductant '"+reductant+"' in the half reaction corresponding to index " + str(idx) + " "
                            "in the half_cell_reactions table. Please flip the species so that "
                            "the oxidant is in the Oxidant column.")
                    
                    # TODO: automatically flip it for the user...
                else:
                    # if e- transfer for this common elem is negative (e- is transferred)
                    # or zero (e.g., if K, Na, or Ca is a common element)
                    common_elem_electron_dict[ce] = {
                            "electron_coeff":electron_coeff,
                            "ox_coeff":-ox_coeff,
                            "total_ox_ce_oxstate":total_ox_ce_oxstate,
                            "red_coeff":red_coeff,
                            "total_red_ce_oxstate":total_red_ce_oxstate,
                            }

            # print("common_elem_electron_dict")
            # print(common_elem_electron_dict)
            
            if len(common_elem_electron_dict.keys()) > 1:
                # TODO: figure out what to do if there is more than one common element that
                # participates in electron transfer...
                pass
            
            if isinstance(ox_dissrxn, str):
                sp_diss_ox = ox_dissrxn.strip().split(" ")[1::2]
                sp_diss_ox = sp_diss_ox[1:]
                sp_diss_ox = [s if s != "O2(g)" else "O2" for s in sp_diss_ox]
            else:
                sp_diss_ox = [oxidant]
            
            sp_diss_ox = [s for s in sp_diss_ox if s not in ["H+", "H2O"]]
            
            basis_candidates = []
            basis_candidate_elems = []
            for f in sp_diss_ox:
                formula_dict = parse_formula(self.thermo.thermo_db["formula"].loc[self.thermo.thermo_db["name"]==f].values[0])
                elems = [k for k in list(formula_dict.keys()) if k not in ["+", "-", "H", "O"]]
                if len(elems) == 1:
                    if (elems[0] in common_elems or elems[0] in elems_unique_to_oxidant) and elems[0] not in basis_candidate_elems:
                        basis_candidates.append(f)
                        basis_candidate_elems.append(elems[0])
            
            if isinstance(red_dissrxn, str):
                sp_diss_red = red_dissrxn.strip().split(" ")[1::2]
                sp_diss_red = sp_diss_red[1:]
                sp_diss_res = [s if s != "O2(g)" else "O2" for s in sp_diss_red]
            else:
                sp_diss_red = [reductant]
            
            sp_diss_red = [s for s in sp_diss_red if s not in ["H+", "H2O"]]
            
            for f in sp_diss_red:
                formula_dict = parse_formula(self.thermo.thermo_db["formula"].loc[self.thermo.thermo_db["name"]==f].values[0])
                elems = [k for k in list(formula_dict.keys()) if k not in ["+", "-", "H", "O"]]
                if len(elems) == 1:
                    if elems[0] in elems_unique_to_reductant and elems[0] not in basis_candidate_elems:
                        basis_candidates.append(f)
                        basis_candidate_elems.append(elems[0])

            if len(common_elem_electron_dict) == 0:
                # no common elements between oxidant and reductant
                self.err_handler.raise_exception("A common element could not be found between "
                        "the oxidant " + str(oxidant) + " and the reductant " + str(reductant))
                
            
            unpacked_dict = common_elem_electron_dict[list(common_elem_electron_dict.keys())[0]]

            e_coeff = unpacked_dict["electron_coeff"]
            ox_coeff = unpacked_dict["ox_coeff"]
            red_coeff = unpacked_dict["red_coeff"]

            # print("BASIS CANDIDATES")
            # print(basis_candidates+["H2O", "e-", "H+"])
            
            b = pychnosz.basis(basis_candidates+["H2O", "e-", "H+"],
                           messages=False, global_state=False)
            
            sout = pychnosz.subcrt(species=[oxidant, reductant, "e-"],
                                   coeff=[ox_coeff, red_coeff, e_coeff],
                                   property="logK", T=25, show=False,
                                   messages=False, basis=b).reaction
            
            half_reaction_dict[idx] = self._create_sp_dict(sout)

        # print("HALF REACTION DICT")
        # print(half_reaction_dict)

        # find all possible combinations of idx pairs (but not reverse rxns to save time)
        good_idx_list = [idx for idx in idx_list if idx not in bad_idx_list]
        redox_pair_list = [list(p) for p in list(itertools.combinations(good_idx_list, 2))]


        # Weed out redox pairs that would have the same
        # electron-donating/accepting species as a product and reactant. E.g.
        # 'CO2': -1.0, 'CH4': 1.0, 'e-': -8.0, 'H2O': 2.0, 'H+': -8.0
        # 'CO': -1.0, 'CH4': 1.0, 'e-': -6.0, 'H2O': 1.0, 'H+': -6.0
        # In this example, CH4 would be both a product and reactant in the resulting
        # full redox reaction. Get rid of these pairs because they cause erroneous
        # calculation of electrons transferred when half reactions are summed.
        bad_redox_pair_index_list = []
        for i,p in enumerate(redox_pair_list):
            half_rxn_idx_1 = p[0]
            half_rxn_idx_2 = p[1]

            reactant_1 = [k for k,v in zip(half_reaction_dict[half_rxn_idx_1].keys(), half_reaction_dict[half_rxn_idx_1].values()) if v<0]
            reactant_2 = [k for k,v in zip(half_reaction_dict[half_rxn_idx_2].keys(), half_reaction_dict[half_rxn_idx_2].values()) if v<0]
            product_1 = [k for k,v in zip(half_reaction_dict[half_rxn_idx_1].keys(), half_reaction_dict[half_rxn_idx_1].values()) if v>0]
            product_2 = [k for k,v in zip(half_reaction_dict[half_rxn_idx_2].keys(), half_reaction_dict[half_rxn_idx_2].values()) if v>0]
            reactant_1 = [r for r in reactant_1 if r not in ["H2O", "H+", "e-"]]
            reactant_2 = [r for r in reactant_2 if r not in ["H2O", "H+", "e-"]]
            product_1 = [r for r in product_1 if r not in ["H2O", "H+", "e-"]]
            product_2 = [r for r in product_2 if r not in ["H2O", "H+", "e-"]]
            if len([ii for ii in reactant_1 if ii in reactant_2]) != 0:
                bad_redox_pair_index_list.append(i)
            if len([ii for ii in product_1 if ii in product_2]) != 0:
                bad_redox_pair_index_list.append(i)
        
        redox_pair_list = [i for j, i in enumerate(redox_pair_list) if j not in bad_redox_pair_index_list]
        
        if len(redox_pair_list) == 0:
            if self.verbose > 0:
                msg = ("No valid redox reactions could be made with these half reactions. "
                       "This is likely because a species appears as both an oxidant and a "
                       "reductant in the full redox reaction. For example, "
                       "combining the half reactions 'CO2 to CH4' and 'CO to CH4' "
                       "would produce a full redox reaction with CH4 as both a product "
                       "and a reactant. If you would like to model a "
                       "comproportionation or disproportionation reaction, be sure "
                       "to select half reactions that do not share reactants and products. "
                       "For example, 'CO2 to CO' and 'CO to CH4' is a valid combination that "
                       "would represent the comproportionation or disproportionation of CO in "
                       "the forward and backward redox reactions.")
                print(msg)
            return
        
        self.redox_pair_list = redox_pair_list

        # create 'reaction_dict': a dictionary of reactions keyed by their idx pairs
        reaction_dict = {}
        e_dict = {}
        for pair in redox_pair_list:

            # print("Pair:", str(pair))
            
            half_reaction_dict_1 = copy.deepcopy(half_reaction_dict[pair[0]])
            half_reaction_dict_2 = copy.deepcopy(half_reaction_dict[pair[1]])

            # print("half reactions")
            # print(half_reaction_dict_1)
            # print(half_reaction_dict_2)
            
            # find the lowest common multiple of e-
            # ensure e- is an integer value or else lcm() won't work
            assert int(half_reaction_dict_1["e-"]) == half_reaction_dict_1["e-"]
            assert int(half_reaction_dict_2["e-"]) == half_reaction_dict_2["e-"]

            # find lowest common multiple (lcm) of the electrons in the two half reactions
            e_lcm = math.lcm(int(half_reaction_dict_1["e-"]), int(half_reaction_dict_2["e-"]))
            
            # use the lcm to multiply half reaction coefficients to get the
            # same number of electrons transferred in each half reaction
            mult_1 = abs(e_lcm/half_reaction_dict_1["e-"])
            mult_2 = abs(e_lcm/half_reaction_dict_2["e-"])
            for k in list(half_reaction_dict_1.keys()):
                half_reaction_dict_1[k] = half_reaction_dict_1[k]*mult_1
            for k in list(half_reaction_dict_2.keys()):
                half_reaction_dict_2[k] = -half_reaction_dict_2[k]*mult_2

            # print("lcm multiplied")
            # print(half_reaction_dict_1)
            # print(half_reaction_dict_2)

            # sum the half reaction dicts to write the full balanced redox reaction
            full_rxn_dict = {k: half_reaction_dict_1.get(k, 0) + half_reaction_dict_2.get(k, 0) for k in set(half_reaction_dict_1) | set(half_reaction_dict_2)}

            # sum H+ and OH- to make H2O, and modify the full rxn dict appropriately
            if "OH-" in list(full_rxn_dict.keys()) and "H+" in list(full_rxn_dict.keys()):
                if ((full_rxn_dict["H+"]== full_rxn_dict["OH-"]) & (full_rxn_dict["H+"]==0)) or (full_rxn_dict["H+"]*full_rxn_dict["OH-"]>0):
                # check if the coefficients of OH- and H+ have the same sign

                    water_dict = {}
                    
                    if abs(full_rxn_dict["H+"]) == abs(full_rxn_dict["OH-"]):
                        water_dict = {"H2O": full_rxn_dict["H+"]}
            
                    elif abs(full_rxn_dict["H+"]) > abs(full_rxn_dict["OH-"]):
                        water_dict = {"H2O": full_rxn_dict["OH-"],
                                      "H+":full_rxn_dict["H+"] - full_rxn_dict["OH-"]}
                    else:
                        water_dict = {"H2O": full_rxn_dict["H+"],
                                      "OH-":full_rxn_dict["OH-"] - full_rxn_dict["H+"]}
                        
                    del full_rxn_dict["H+"]
                    del full_rxn_dict["OH-"]
        
                    # sum the full_rxn_dict and the water_dict to ensure OH- and H+ make H2O
                    full_rxn_dict = {k: full_rxn_dict.get(k, 0) + water_dict.get(k, 0) for k in set(full_rxn_dict) | set(water_dict)}

            
            for k in list(full_rxn_dict.keys()):
                if full_rxn_dict[k] == 0:
                    del full_rxn_dict[k]

            # print("full_rxn_dict before")
            # print(full_rxn_dict)
            
            # multiply all coefficients by their least common multiple
            denoms = [Fraction(x).limit_denominator().denominator for x in list(full_rxn_dict.values())]
            coeff_lcm = functools.reduce(lambda a,b: a*b//math.gcd(a,b), denoms)
            full_rxn_dict = {k:full_rxn_dict[k]*coeff_lcm for k in full_rxn_dict.keys()}

            # print("full_rxn_dict after")
            # print(full_rxn_dict)
            
            e_transferred = half_reaction_dict_1["e-"]*coeff_lcm
            e_dict[str(pair[0])+"_"+str(pair[1])] = abs(e_transferred)
            e_dict[str(pair[1])+"_"+str(pair[0])] = abs(e_transferred)

            # print("e_transferred")
            # print(e_transferred)

            reaction_dict[str(pair[0])+"_"+str(pair[1])] = full_rxn_dict
            reaction_dict[str(pair[1])+"_"+str(pair[0])] = {k:-v for k,v in zip(full_rxn_dict.keys(), full_rxn_dict.values())}

            # print("reaction_dict")
            # print(reaction_dict)

        for i,key in enumerate(list(reaction_dict.keys())):
            redox_pair = key.split("_")
            redox_pair = [int(v) for v in redox_pair]
    
            name = "rxn_"+str(key)
            
            species_dict_formatted = {"species_"+str(i+1):[k] for i,k in enumerate(reaction_dict[key].keys())}
            coeff_dict_formatted = {"coeff_"+str(i+1):[reaction_dict[key][k]] for i,k in enumerate(reaction_dict[key].keys())}
            
            reaction_dict[key] = {"reaction_name":[name],
                                  "redox_pairs":[redox_pair],
                                  "mol_e-_transferred_per_mol_rxn":[e_dict[key]]}
            reaction_dict[key].update(species_dict_formatted)
            reaction_dict[key].update(coeff_dict_formatted)

        self.reaction_dict = reaction_dict

        # make the affinity_energy_reactions_table
        for i,key in enumerate(list(self.reaction_dict.keys())):
            if i == 0:
                affinity_energy_reactions_table = pd.DataFrame(self.reaction_dict[key])
            else:
                affinity_energy_reactions_table = pd.concat([affinity_energy_reactions_table, pd.DataFrame(self.reaction_dict[key])])
        
        self.redox_reactions_table = affinity_energy_reactions_table.set_index("reaction_name")
        
        # rearrange column order
        non_coeff_non_sp_cols = [c for c in self.redox_reactions_table.columns if "coeff_" not in c and "species_" not in c]
        coeff_cols = [c for c in self.redox_reactions_table.columns if "coeff_" in c]
        species_cols = [c for c in self.redox_reactions_table.columns if "species_" in c]
        
        coeff_sp_cols = [None]*(len(coeff_cols)+len(species_cols))
        coeff_sp_cols[::2] = coeff_cols
        coeff_sp_cols[1::2] = species_cols
        new_col_order = non_coeff_non_sp_cols + coeff_sp_cols

        self.redox_reactions_table = self.redox_reactions_table[new_col_order]

        reverse_pair_list = [[v[1], v[0]] for v in self.redox_pair_list]
        pair_list_interleave = self.redox_pair_list + reverse_pair_list
        pair_list_interleave[::2] = self.redox_pair_list
        pair_list_interleave[1::2] = reverse_pair_list
        self.redox_pair_list = pair_list_interleave

        # sort by forward then backward reaction, e.g., [0, 1], [1, 0], [0, 2], [2, 0]...
        sort_order = ["rxn_"+str(v[0])+"_"+str(v[1]) for v in self.redox_pair_list]
        self.redox_reactions_table.reindex(sort_order)
        
        # # sort rows by ascending redox pairs
        # self.redox_reactions_table = self.redox_reactions_table.sort_values('redox_pairs', key=lambda col: col.map(lambda x: [x[0], x[1]]))

        if isinstance(self.reactions_for_plotting, pd.DataFrame):
            self.reactions_for_plotting = self.show_redox_reactions(
                    formatted=formatted,
                    charge_sign_at_end=charge_sign_at_end,
                    show=show).combine_first(self.reactions_for_plotting)
        else:
            self.reactions_for_plotting = self.show_redox_reactions(
                    formatted=formatted,
                    charge_sign_at_end=charge_sign_at_end,
                    show=show)

    
    
    @staticmethod
    def _create_sp_dict(sout):
        coeffs = list(sout["coeff"])
        species = list(sout["name"])
    
        species = [sp if sp != "water" else "H2O" for sp in species] # replace "water" with "H2O"
        
        sp_dict = {}
        for i,sp in enumerate(species):
            if sp in list(sp_dict.keys()):
                sp_dict[sp] = sp_dict[sp] + coeffs[i]
            else:
                sp_dict[sp] = coeffs[i]
        
        for sp in list(sp_dict.keys()):
            if sp_dict[sp] == 0:
                del sp_dict[sp]
    
        return sp_dict
    
    @staticmethod
    def _formula_ox_to_dict(f):
        f_split = f.strip().split(" ")
    
        out_dict = {}
        for f in f_split:
        
            num_elem_split = re.sub( r"([A-Z])", r" \1", f).split()
            if len(num_elem_split) == 1:
                # element has no coeff e.g., "S+6"
                n_elem = 1.0
                elem_ox = num_elem_split[0]
                split_pos = elem_ox.strip().split("+")
                split_neg = elem_ox.strip().split("-")
            else:
                # element has a coeff e.g., "4O-2"
                n_elem = float(num_elem_split[0])
                elem_ox = num_elem_split[1]
                split_pos = elem_ox.strip().split("+")
                split_neg = elem_ox.strip().split("-")
        
            if len(split_pos) != len([elem_ox]):
                # element has positive ox state
                elem = split_pos[0]
                if split_pos[1] == "":
                    split_pos[1] = 1.0
                ox_state = float(split_pos[1])
            elif len(split_neg) != len([elem_ox]):
                # element has negative ox state
                elem = split_neg[0]
                if split_neg[1] == "":
                    split_neg[1] = 1.0
                ox_state = -float(split_neg[1])
            else:
                # element has ox state of 0
                elem = elem_ox
                ox_state = 0.0
        
            # print("RESULTS")
            # print(f)
            # print(n_elem)
            # print(elem)
            # print(ox_state)
    
            if elem not in list(out_dict.keys()):
                out_dict[elem] = {ox_state:n_elem}
            else:
                out_dict[elem][ox_state] = n_elem
        

        return out_dict


    def calculate_energy(self, species, stoich,
                    divisor=1, per_electron=False,
                    grams_minerals=None,
                    rxn_name="custom reaction",
                    negative_energy_supplies=False,
                    y_type="A", y_units="kcal", 
                    limiting=None, charge_sign_at_end=False,
                    as_written=False,
                    simple_df_output=False,
                    append_report=True,
                    custom_grouping_filepath=None,
                    print_logK_messages=False,
                    raise_nonlimiting_exception=True):

        """
        Calculate Gibbs free energy, logK, logQ, chemical affinity, or energy
        supply for a user-defined reaction across all samples in a speciation.
        
        Parameters
        ----------
        species : list of str
            List of species in the reaction
            
        stoich : list of numeric
            List of stoichiometric reaction coefficients (reactants are negative)

        divisor : numeric, default 1, or list of numeric
            If a single numeric value is provided, divide all calculated values
            (Gibbs free energies, affinities, etc.) by that value. Synergizes
            with the parameter `per_electron` when normalizing calculated values
            to a per electron basis. See `per_electron` for more information.

        per_electron : bool, default False
            If False, values calculated by this function will be treated as
            'per mole of reaction'. If set to True, then the calculation will
            assume that `divisor` is normalizing calculated values to a per
            electron basis. For example, sulfide oxidation to sulfate has the
            following reaction: [2 O2 + H2S = SO4-2 + 2 H+]. The oxidation state
            of sulfur changes from S-2 in sulfide to S+6 in sulfate, a
            difference of 8 electrons. If you use `calculate_energy` to
            calculate the Gibbs free energy per mole of electrons transferred,
            you would set `divisor` to 8 and `per_electron` to True.

        grams_minerals : float, int, dict, or pandas.DataFrame, optional
            Number of grams belonging to each mineral reactant when calculating
            the limiting reactant during an energy supply calculation. This
            parameter is only used when `y_type="E"`.
            For example, in the reaction
            4 goethite + 4 iron + 3 O2 = 4 magnetite + 2 H2O
            setting `grams_minerals = 0.001` would mean 0.001 grams of goethite
            is reacting with 0.001 grams of iron. If it is desirable to specify
            individual masses for each mineral reactant, then a dictionary can
            be provided. For example:
            `grams_minerals={"goethite": 0.001, "iron": 0.1},`
            To vary the mass of a mineral from sample to sample, provide a
            dataframe with a column per mineral and a row per sample.
            By default, mineral masses are left undefined and minerals are not
            considered when looking for a limiting reactant. Note that defining
            a mass of 0 grams is different: a mineral present at 0 grams is
            immediately the limiting reactant, giving an energy supply of 0.

        rxn_name : str, default "custom reaction"
            Name for the reaction, e.g., "sulfide oxidation to sulfate".

        negative_energy_supplies : bool, default False
            Report negative energy supplies? If False, negative energy supplies
            are reported as 0. If True, negative energy supplies are
            reported. A 'negative energy supply' represents the energy cost of
            depleting the limiting reactant of a reaction. This metric is not
            always helpful when examing energy supply results, so this option is
            set to False by default.

        y_type : str, default "A"
            The variable to plot on the y-axis. Can be either 'A' (for chemical
            affinity), 'G' (for Gibbs free energy, ΔG), 'logK' (for the log
            of the equilibrium constant), 'logQ' (for the log of the reaction
            quotient), or 'E' for energy supply.

        y_units : str, default "kcal"
            The unit that energy will be reported in (per mol for G and A, or
            per kg fluid for energy supply, or unitless for logK and logQ).
            Can be 'kcal', 'cal', 'J', or 'kJ'.

        limiting : str, optional
            Name of the species to act as the limiting reactant when calculating
            energy supply. If this parameter is left undefined, then a
            limiting reactant will be chosen automatically based on
            concentration and stoichiometry. This parameter is ignored unless
            `y_type` is set to 'E' (energy supply).

        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-)?

        as_written : bool, default False
            If `as_written` is False, then built-in speciation groups
            will be used to calculate limiting reactants. For example, if CaCO3
            is a reactant, then the code tests whether the limiting reactant is
            the sum of CaCO3, CO2, HCO3-, CO3-2, NaCO3-... everything in the
            carbonate speciation group. If `as_written` is True, then speciation
            groups will not be used; only the species defined as reactants
            will be used to test whether a reactant is limiting
            (e.g., just CaCO3). This parameter is ignored unless `y_type` is set
            to 'E' (energy supply).

        simple_df_output : bool, default False
            By default, the dataframe returned by this function is multiindexed;
            i.e., there is a header column, a subheader defining the units,
            and then columns of values. If `simple_df_output` is set to True,
            then the dataframe returned will not be multiindexed; there will
            only be one column header that includes units.

        append_report : bool, default True
            Add or update calculated values to the speciation object's report
            dataframe?

        custom_grouping_filepath : str, optional
            Filepath for a TXT file containing customized speciation groups. Use
            to override the built-in speciation group file.

        print_logK_messages : bool, default False
            Print additional messages related to calculating logK values for
            each sample?

        raise_nonlimiting_exception : bool, default True
            This parameter can be ignored in almost all cases. Raise an
            exception when there are no available limiting reactants?
            The purpose of this parameter is toggle off error message
            interruptions when this function is called by
            `apply_redox_reactions`, which can test many different reactions at
            once, some of which do not have valid limiting reactants and would
            otherwise be interrupted by errors.
        
            
        Returns
        ----------
        Pandas dataframe
            Returns a multiindexed dataframe of samples and calculated results.
            If `append_report` is True, then the report attribute of the
            speciation object will be appended/updated. If `simple_df_output` is
            set to True, then the dataframe will not be multiindexed.
            
        """

        # check that y_type is recognized
        if y_type not in ["logK", "logQ", "G", "A", "E"]:
            self.err_handler.raise_exception("Valid options for y_type include "
                    "'logK', 'logQ', 'G' (Gibbs free energy), 'A' "
                    "(chemical affinity), and 'E' (energy suppy.")
        
        # check that a thermodynamic CSV is being used
        if not isinstance(self.thermo.csv_db, pd.DataFrame):
            self.err_handler.raise_exception("The calculate_energy() function requires "
                    "a thermodynamic database in a WORM-style CSV format, e.g., "
                    "'wrm_data.csv'. You may be getting this message because "
                    "a data0 or data1 file was used.")
        
        # check that the divisor is valid
        if isinstance(divisor, list) or isinstance(divisor, pd.Series):
            if len(divisor) != len(self.misc_params["Temp(C)"]):
                self.err_handler.raise_exception("The length of the divisor is "
                    "not equal to the number of samples.")

        # set the custom grouping filepath
        if isinstance(custom_grouping_filepath, str):
            self.custom_grouping_filepath = custom_grouping_filepath
        
        # check that the reaction is balanced
        formulas = []
        for s in species:
            if s == "H+":
                formulas.append("H+")
            elif s == "H2O":
                formulas.append("H2O")
            else:
                if s in list(self.thermo.csv_db["name"]):
                    formulas.append(list(self.thermo.csv_db[self.thermo.csv_db["name"]==s]["formula"])[0])
                else:
                    self.err_handler.raise_exception("Valid thermodynamic data "
                            "was not found for species "+str(s)+"")
                    
        missing_composition = check_balance(formulas, stoich)

        if y_type == "E":
            self._set_speciation_groups()

        # assign aq_distribution_logact table to Speciation
        sample_dict = {}
        for i,sample in enumerate(self.sample_data.keys()):
            sample_data = self.sample_data[sample]["aq_distribution"]["log_activity"]
            sample_dict[sample] = {name:logact for name, logact in zip(sample_data.index, sample_data)}
        df_aq_distribution_logact = pd.DataFrame(sample_dict)
        df_aq_distribution_logact = df_aq_distribution_logact.T
        df_aq_distribution_logact.insert(0, "Xi", 0)
        self.aq_distribution_logact = df_aq_distribution_logact

        # assign aq_distribution_molal table to Speciation
        sample_dict = {}
        for i,sample in enumerate(self.sample_data.keys()):
            sample_data = self.sample_data[sample]["aq_distribution"]["molality"]
            sample_dict[sample] = {name:molal for name, molal in zip(sample_data.index, sample_data)}
        df_aq_distribution_molal = pd.DataFrame(sample_dict)
        df_aq_distribution_molal = df_aq_distribution_molal.T
        df_aq_distribution_molal.insert(0, "Xi", 0)
        self.aq_distribution_molal = df_aq_distribution_molal

        # assign misc_params table to Speciation
        df_misc_param_row_list = []
        for sample in self.sample_data.keys():
            df_misc_param_row_list.append(pd.DataFrame({
                "Temp(C)" : [self.sample_data[sample]['temperature']],
                "Press(bars)" : [self.sample_data[sample]['pressure']],
                "pH" : [-self.sample_data[sample]['aq_distribution']["log_activity"]['H+']],
                "logfO2" : [self.sample_data[sample]["fugacity"]["log_fugacity"]["O2(g)"]],
                "ionic strength (molal)" : [self.sample_data[sample]["ionic_strength"]],
            }))
        df_misc_params = pd.concat(df_misc_param_row_list)
        df_misc_params.insert(0, "Xi", 0)
        df_misc_params.index = list(self.sample_data.keys())
        self.misc_params = df_misc_params
        
        # check that there are valid limiting reactants when calculating energy
        # e.g., prevent issue when the only reactant is a mineral, etc.
        reactant_idx = [1 if i<0 else 0 for i in stoich]
        reactants = [species[i] for i,idx in enumerate(reactant_idx) if idx == 1]
        invalid_limiting_reactants = []
        for r in reactants:
            if r not in ["H2O", "H+", "OH-"]:
                if list(self.thermo.csv_db[self.thermo.csv_db["name"]==r]["state"])[0] != "aq" and not self._mineral_amounts_known(grams_minerals):
                    invalid_limiting_reactants.append(r)
            else:
                invalid_limiting_reactants.append(r)

        no_limiting_reactants = False
        if reactants == invalid_limiting_reactants and y_type == "E" and raise_nonlimiting_exception:
            self.err_handler.raise_exception("Energy supply for this reaction "
                "cannot be calculated because none of the reactants are "
                "limiting. A limiting reactant must be aqueous and cannot be H+ "
                "or OH-.")
        elif reactants == invalid_limiting_reactants and y_type == "E" and not raise_nonlimiting_exception:
            no_limiting_reactants = True

        
        if limiting != None and y_type == "E":

            limiting = self._switch_limiting(limiting,
                                              stoich=stoich,
                                              species=species,
                                              lenient=False)
            
            # check that the limiting reactant is in the thermodynamic database
            if limiting not in list(self.thermo.csv_db["name"]):
                self.err_handler.raise_exception("Valid thermodynamic data was "
                        "not found for limiting reactant "+str(limiting)+"")
            
            # check that the limiting reactant is aqueous or gaseous
            if list(self.thermo.csv_db[self.thermo.csv_db["name"]==limiting]["state"])[0] != "aq":
                self.err_handler.raise_exception("The limiting reactant must "
                        "be an aqueous species.")
            
            # check that the limiting reactant is a reactant in the `species` parameter
            if limiting not in reactants:
                self.err_handler.raise_exception("The species specified as a "
                        "limiting reactant, '"+str(limiting)+"', is not a "
                        "reactant in this reaction.")
                
        # format reaction equation
        equation_to_display = format_equation(
                                      species,
                                      stoich,
                                      charge_sign_at_end=charge_sign_at_end,
                                      )
        
        # create a dictionary of species logacts across samples
        s_logact_dict = {}
        s_molal_dict = {}
        
        for s in species:
            if s == "H+":
                s_logact_dict[s] = [v for v in list(self.aq_distribution_logact["H+"])]
                s_molal_dict[s] = [float("NaN")]*len(list(self.aq_distribution_logact["H+"]))
            elif s == "H2O":
                s_logact_dict[s] = [v for v in list(self.aq_distribution_logact["H2O"])]
                s_molal_dict[s] = [float("NaN")]*len(self.aq_distribution_logact["H2O"])
            elif list(self.thermo.csv_db[self.thermo.csv_db["name"]==s]["state"])[0] not in ["cr", "liq"]:

                if s in self.aq_distribution_logact.columns:

                    # aqueous species
                    s_logact_dict[s] = list(self.aq_distribution_logact[s])
                    if isinstance(self.reactant_dict_scalar, dict):
                        s_molal_dict[s] = self._grouped_molality(
                                                s,
                                                self.aq_distribution_molal,
                                                as_written=as_written)
                    else:
                        s_molal_dict[s] = list(self.aq_distribution_molal[s])
                else:
                    s_logact_dict[s] = [float('NaN')]*self.aq_distribution_logact.shape[0]
                    s_molal_dict[s] = [float('NaN')]*self.aq_distribution_logact.shape[0]
                    # self.err_handler.raise_exception("The species "+str(s)+" is "
                    #         "not among the distribution of aqueous species in "
                    #         "this calculation.")
            else:
                # liq and cr species
                s_logact_dict[s] = [0]*len(self.misc_params["Temp(C)"])
                s_molal_dict[s] = self._mineral_molality(
                                        s,
                                        len(self.misc_params["Temp(C)"]),
                                        grams_minerals)

        if y_type in ["logK", "logQ"]:
            y_type_plain = copy.copy(y_type)
        elif y_type == "A":
            if per_electron:
                y_type_plain = "affinity per mole e-"
            else:
                y_type_plain = "affinity per mole rxn"
        elif y_type == "G":
            y_type_plain = "Gibbs free energy"
        elif y_type == "E":
            y_type_plain = "energy supply"
        else:
            # this shouldn't happen and will be caught by the y_type check above
            pass
        
        if y_type not in ["logK", "logQ"]:
            if y_units in ["cal", "kcal"]:
                r_div = 4.184
            elif y_units in ["J", "kJ"]:
                r_div = 1
            else:
                self.err_handler.raise_exception("The specified y_unit '"+y_units+"' "
                        "is not recognized. Try 'cal', 'kcal', 'J', or 'kJ'.")
            R = 8.314/r_div  # gas constant, unit = [cal/mol/K]

        if "k" in y_units:
            k_div = 1000
        else:
            k_div = 1
            
        y_list = []
        lr_name_list = []
        for i,T in enumerate(list(self.misc_params["Temp(C)"])):
            
            if isinstance(divisor, list):
                divisor_i = divisor[i]
            else:
                divisor_i = divisor
            
            if y_type != "logQ":
                logK = pychnosz.subcrt(
                              species,
                              stoich,
                              T=T,
                              P=list(self.misc_params["Press(bars)"])[i],
                              IS=0, # ionic strength of 0 ensures that the logK calculated is standard, as required
                              messages=print_logK_messages, show=False).out["logK"]

                logK = float(logK.iloc[0])

            if y_type == "logK":
                ylab_out = "log K"
                y_list.append(round(logK/divisor_i, 4))
                df_y_name = "logK"
                continue

            # print("")
            # print("stoich")
            # print([st for st in stoich])
            # print("species")
            # print([sp for sp in species])
            # print("logact")
            # print([s_logact_dict[sp][i] for sp in species])
            # print("result")
            # print(sum([st*s_logact_dict[sp][i] for st,sp in zip(stoich,species)]))
            

            
            logQ = sum([st*s_logact_dict[sp][i] for st,sp in zip(stoich,species)])

            if y_type == "logQ":
                ylab_out = "log Q"
                y_list.append(round(logQ/divisor_i, 4))
                df_y_name = "logQ"
                continue
            
            else:
                A = 2.303 * R * (273.15+T) * (logK - logQ)  # affinity, unit = [cal/mol]
                A = A/k_div
                
                if y_type=="G":
                    G = -A # gibbs free energy, unit = [cal/mol]
                    y_list.append(G/divisor_i)
                    ylab_out="ΔG, {}/mol".format(y_units)
                    y_units_out = y_units+"/mol"
                    if per_electron:
                        y_units_out = y_units_out+" e-"
                elif y_type=="A":
                    y_list.append(round(A/divisor_i, 4))
                    ylab_out="A, {}/mol".format(y_units)
                    y_units_out = y_units+"/mol"
                    if per_electron:
                        y_units_out = y_units_out+" e-"
                elif y_type=="E":

                    y_units_out = y_units+"/kg fluid"
                    ylab_out="Energy Supply, {}".format(y_units+"/kg fluid")

                    if no_limiting_reactants:
                        df_y_name = y_type_plain+", "+y_units_out
                        y_list.append(float('NaN'))
                        lr_name_list.append("None")
                        continue
                    
                    lr_name, lr_reported, lr_concentration, lr_stoich = \
                            self._select_limiting_reactant(
                                    species=species,
                                    stoich=stoich,
                                    s_molal_dict=s_molal_dict,
                                    step=i,
                                    limiting=limiting,
                                    allow_minerals=self._mineral_amounts_known(grams_minerals),
                                    charge_sign_at_end=charge_sign_at_end,
                                    )

                    if lr_name != None:
                        E = A * (lr_concentration/lr_stoich)
                        y_list.append(E/divisor_i)
                        lr_name_list.append(lr_reported)
                    else:
                        y_list.append(float('NaN'))
                        lr_name_list.append(float('NaN'))

                df_y_name = y_type_plain+", "+y_units_out
        
        if not negative_energy_supplies and y_type == "E":
            for i,v in enumerate(y_list):
                if v < 0:
                    y_list[i] = 0.0

        if y_type == "E":
            df_out = pd.DataFrame({df_y_name:y_list,
                                   "limiting reactant":lr_name_list},
                                   index=self.misc_params.index)
        else:
            df_out = pd.DataFrame({df_y_name:y_list},
                                   index=self.misc_params.index)

        df_out_simple = copy.deepcopy(df_out)

        if y_type in ["logK", "logQ"]:
            headers = [rxn_name+" "+y_type]
            subheaders = [y_type]
        
        elif y_type in ["A", "G"]:
            hs = df_out.columns[0].strip().split(", ")
            headers = [rxn_name+" "+hs[0]]
            subheaders = [hs[1]]
            
        else:
            # energy supplies have an extra 'limiting reactant' column to deal with
            hs = [h.strip().split(", ") for h in df_out.columns]
            hs = [[""]+h if len(h) == 1 else h for h in hs]
            hs = [[rxn_name+" "+h[0], h[1]] for h in hs]
            hs = [hs[0], [hs[1][0]+"limiting reactant", hs[1][1]]]
            headers = [hs[0][0], hs[1][0]]
            subheaders = [hs[0][1], hs[1][1]]
            
        multicolumns = pd.MultiIndex.from_arrays(
            [headers, subheaders], names=['Sample', ''])
        df_out.columns = multicolumns

        if append_report:
            col_order = [c[0] for c in self.report.columns] + headers # get order of columns in report and affinity/energy df
            col_order = list(collections.OrderedDict.fromkeys(col_order)) # remove duplicates
            # Suppress unorderable values warning when combining DataFrames with MultiIndex columns
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'The values in the array are unorderable')
                self.report = df_out.combine_first(self.report) # update the report with affinity/energy results
            self.report = self.report.reindex(level=0, columns=col_order) # restore column order

            # update the report category dictionary so y_type_plain appears as a category with relevant columns
            if y_type_plain not in list(self.report_divs.keys()):
                self.report_divs[y_type_plain] = headers
            else:
                self.report_divs[y_type_plain] = self.report_divs[y_type_plain]+headers
                self.report_divs[y_type_plain]= list(collections.OrderedDict.fromkeys(self.report_divs[y_type_plain]))

        # create a formatted reaction to add to self.reactions_for_plotting
        # so that it can be invoked in a plot
        formatted_rxn = reaction = self.format_reaction(
                                            coeffs=stoich,
                                            names=species,
                                            formatted=True,
                                            charge_sign_at_end=True,
                                            show=False)
        
        calc_energy_df = pd.DataFrame({"reaction_name":[rxn_name],
                                       "redox_pairs":[float('NaN')],
                                       "reaction":[formatted_rxn],
                                      }).set_index("reaction_name")
        
        if isinstance(self.reactions_for_plotting, pd.DataFrame):
            self.reactions_for_plotting = calc_energy_df.combine_first(self.reactions_for_plotting)
        else:
            self.reactions_for_plotting = calc_energy_df

        if simple_df_output:
            return df_out_simple # simple dataframe
        else:
            return df_out # multiindex

    
    @staticmethod
    def __unique(seq):
        """
        Provide a sequence, get a list of non-repeating elements in the same order.
        """
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    
    def save(self, filename, messages=True):
        """
        Save the speciation as a '.speciation' file to your current working
        directory. This file can be loaded with `AqEquil.load(filename)`.
        
        Parameters
        ----------
        filename : str
            The desired name of the file.
            
        messages : str
            Print a message confirming the save?
        """
        
        if filename[-11:] != '.speciation':
            filename = filename + '.speciation'

        with open(filename, 'wb') as handle:
            dill.dump(self, handle, protocol=dill.HIGHEST_PROTOCOL)
            if messages:
                print("Saved as '{}'".format(filename))

                
    @staticmethod
    def _save_figure(fig, save_as, save_format, save_scale, plot_width, plot_height, ppi):
        if isinstance(save_format, str) and save_format not in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', 'html']:
            self.err_handler.raise_exception("{}".format(save_format)+" is an unrecognized "
                            "save format. Supported formats include 'png', "
                            "'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', "
                            "'json', or 'html'")
            
        if isinstance(save_format, str):
            if not isinstance(save_as, str):
                save_as = "newplot"
            if save_format=="html":
                fig.write_html(save_as+".html")
                print("Saved figure as {}".format(save_as)+".html")
                save_format = 'png'
            elif save_format in ['pdf', 'eps', 'json']:
                pio.full_figure_for_development(fig, warn=False)
                pio.write_image(fig, save_as+"."+save_format, format=save_format, scale=save_scale,
                                width=plot_width*ppi, height=plot_height*ppi)
                print("Saved figure as {}".format(save_as)+"."+save_format)
                save_format = "png"
            else:
                pio.write_image(fig, save_as+"."+save_format, format=save_format, scale=save_scale,
                                width=plot_width*ppi, height=plot_height*ppi)
                print("Saved figure as {}".format(save_as)+"."+save_format)
        else:
            save_format = "png"
            
        return save_as, save_format
    
    
    @staticmethod
    def __get_unit_info(subheader):
        
        unit_name_dict = {
            "pH" : ("", "pH"),
            "ppm" : ("", "ppm"),
            "ppb" : ("", "ppb"),
            "mg/L" : ("", "mg/L"),
            "degC" : ("temperature", "°C"),
            "log_molality" : ("log molality", "log(mol/kg)"),
            "Molality" : ("molality", "mol/kg"),
            "molality" : ("molality", "mol/kg"),
            "molal" : ("molality", "mol/kg"),
            "log_activity" : ("log activity", ""),
            "Log activity" : ("log activity", ""),
            "mg/kg.sol" : ("", "mg solute per kg solution"),
            "Alk., eq/kg.H2O" : ("alkalinity", "eq/kg"),
            "Alk., eq/L" : ("alkalinity", "eq/L"),
            "Alk., eq/kg.sol" : ("alkalinity", "eq/kg solution"),
            "Alk., mg/L CaCO3" : ("alkalinity", "mg/L CaCO3"),
            "Alk., mg/L HCO3-" : ("alkalinity", "mg/L HCO3-"),
            "pX" : ("-(log activity)", "-log(mol/kg)"),
            "activity" : ("activity", ""),
            "log_gamma" : ("log gamma", ""),
            "gamma" : ("gamma", ""),
            "%" : ("", "%"),
            "Eh_volts" : ("Eh", "volts"),
            "eq/kg.H2O" : ("charge", "eq/kg"),
            "logfO2" : ("", ""),
            "cal/mol e-" : ("affinity", "cal/mol e-"),
            "kcal/mol e-" : ("affinity", "kcal/mol e-"),
            "J/mol e-" : ("affinity", "J/mol e-"),
            "kJ/mol e-" : ("affinity", "kJ/mol e-"),
            "cal/mol" : ("affinity", "cal/mol"),
            "kcal/mol" : ("affinity", "kcal/mol"),
            "J/mol" : ("affinity", "J/mol"),
            "kJ/mol" : ("affinity", "kJ/mol"),
            "cal/kg.H2O" : ("energy supply", "cal/kg fluid"), # deprecated
            "cal/kg fluid" : ("energy supply", "cal/kg fluid"),
            "kcal/kg fluid" : ("energy supply", "kcal/kg fluid"),
            "J/kg fluid" : ("energy supply", "J/kg fluid"),
            "kJ/kg fluid" : ("energy supply", "kJ/kg fluid"),
            "Log ion-H+ activity ratio" : ("Log ion-H+ activity ratio", ""),
            "log_fugacity" : ("log fugacity", "log(bar)"),
            "fugacity" : ("fugacity", "bar"),
            "bar" : ("", "bar"),
            "affinity_kcal" : ("saturation index", "kcal/mol"),
        }
        
        out = unit_name_dict.get(subheader)
        
        return out[0], out[1]

    
    def lookup(self, col=None):
        
        """
        Look up desired columns in the speciation report.
        
        Parameters
        ----------
        col : str or list of str
            Leave blank to get a list of section names in the report:
            ```speciation.lookup()```
            Provide the name of a section to look up the names of columns in
            that section of the report:
            ```speciation.lookup("aq_distribution")```
            Provide a column name (or a list of column names) to retrieve the
            column from the report:
            ```speciation.lookup(["Temperature", "O2"])```
            
        Returns
        ----------
        Pandas dataframe or list of str
            If a column name (or list of column names) is provided, returns the
            speciation report with only the desired column(s). Otherwise returns
            a list of section names (if no arguments are provided), or a list of
            columns in a section (if a section name is provided).
        """
        
        names_length = len(self.report_divs.keys())
        
        if col==None and names_length>0:
            return list(self.report_divs.keys())
        
        if names_length>0:
            if col in list(self.report_divs.keys()):
                return list(self.report_divs[col])
        
        if isinstance(col, str):
            col = [col]

        col = list(dict.fromkeys(col)) # remove duplicate search terms while preserving order
        
        df = self.report.iloc[:, self.report.columns.get_level_values(0).isin(set(col))]

        l = [c for c in col if c in df.columns]

        nonexistant_col = [c for c in col if c not in df.columns]
        
        if self.verbose > 0 and len(nonexistant_col) > 0:
            print("Column(s) not found:", nonexistant_col)
        
        return df[l]
    
    
    def __convert_aq_units_to_log_friendly(self, species, rows):

        col_data = self.lookup(species)
        
        col_data = col_data.loc[rows]
        
        if col_data.columns.get_level_values(1) == 'log_activity':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'activity'
        elif col_data.columns.get_level_values(1) == 'log_molality':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'molality'
        elif col_data.columns.get_level_values(1) == 'log_gamma':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'gamma'
        elif col_data.columns.get_level_values(1) == 'log_fugacity':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'fugacity'
        else:
            y = [float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = col_data.columns.get_level_values(1)[0]
        return y, out_unit
    
    
    def plot_mineral_saturation(self, sample_name, title=None,
                                mineral_sat_type="affinity",
                                plot_width=4, plot_height=3, ppi=122,
                                colors=["blue", "orange"],
                                save_as=None, save_format=None, save_scale=1,
                                interactive=True, plot_out=False):
        """
        Vizualize mineral saturation states in a sample as a bar plot.
        
        Parameters
        ----------
        sample_name : str
            Name of the sample to plot.
            
        title : str, optional
            Title of the plot.
        
        mineral_sat_type : str, default "affinity"
            Metric for mineral saturation state to plot. Can be "affinity" or
            "logQoverK".
        
        colors : list of two str, default ["blue", "orange"]
            Sets the color of the bars representing supersaturated
            and undersaturated states, respectively.
            
        save_as : str, optional
            Provide a filename to save this figure. Filetype of saved figure is
            determined by `save_format`.
            Note: interactive plots can be saved by clicking the 'Download plot'
            button in the plot's toolbar.

        save_format : str, default "png"
            Desired format of saved or downloaded figure. Can be 'png', 'jpg',
            'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', or 'html'. If 'html',
            an interactive plot will be saved. Only 'png', 'svg', 'jpeg',
            and 'webp' can be downloaded with the 'download as' button in the
            toolbar of an interactive plot.

        save_scale : numeric, default 1
            Multiply title/legend/axis/canvas sizes by this factor when saving
            the figure.
        
        interactive : bool, default True
            Return an interactive plot if True or a static plot if False.
            
        plot_out : bool, default False
            Return a plotly figure object? If True, a plot is not displayed as
            it is generated.
        """
        
        if sample_name not in self.report.index:
            msg = ("Could not find '{}'".format(sample_name)+" among sample "
                   "names in the speciation report. Sample names include "
                   "{}".format(list(self.report.index)))
            self.err_handler.raise_exception(msg)
        
        if isinstance(self.sample_data[sample_name].get('mineral_sat', None), pd.DataFrame):
            mineral_data = self.sample_data[sample_name]['mineral_sat'][mineral_sat_type].astype(float).sort_values(ascending=False)
            x = mineral_data.index
        else:
            msg = ("This sample does not have mineral saturation state data."
                   "To generate this data, ensure get_mineral_sat=True when "
                   "running speciate(), or ensure this sample has "
                   "mineral-forming basis species.")
            self.err_handler.raise_exception(msg)
        
        color_list = [colors[0] if m >= 0 else colors[1] for m in mineral_data]
            
        if mineral_sat_type == "affinity":
            ylabel = 'affinity, kcal/mol'
        if mineral_sat_type == "logQoverK":
            ylabel = 'logQ/K'
        
        if title==None:
            title = sample_name + " mineral saturation index"
        
        df = pd.DataFrame(mineral_data)

        fig = px.bar(df, x=df.index, y=mineral_sat_type,
            height=plot_height*ppi, width=plot_width*ppi,
            labels={mineral_sat_type: ylabel}, template="simple_white")
        
        fig.update_traces(hovertemplate = "%{x} <br>"+ylabel+": %{y}",
                          marker_color=color_list)
        
        fig.update_layout(xaxis_tickangle=-45, xaxis_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'},
                          margin={"t":40},
                          xaxis={'fixedrange':True},
                          yaxis={'fixedrange':True, 'exponentformat':'power'})
        
        save_as, save_format = self._save_figure(fig, save_as, save_format,
                                                  save_scale, plot_width,
                                                  plot_height, ppi)

        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d',
                                             'lasso2d', 'zoomIn2d', 'zoomOut2d',
                                             'autoScale2d', 'resetScale2d',
                                             'toggleSpikelines'],
                  
                  'toImageButtonOptions': {
                                             'format': save_format, # one of png, svg, jpeg, webp
                                             'filename': save_as,
                                             'height': plot_height*ppi,
                                             'width': plot_width*ppi,
                                             'scale': save_scale,
                                          },
                 }
        if not interactive:
            config['staticPlot'] = True

        if plot_out:
            return fig
        else:
            fig.show(config=config)

    

    def barplot(self, y="pH", title=None, convert_log=True, plot_zero=True,
                show_missing=True, plot_width=4, plot_height=3, ppi=122,
                colormap="WORM", save_as=None, save_format=None, save_scale=1,
                interactive=True, plot_out=False):
        
        """
        Show a bar plot to vizualize one or more variables across all samples.
        
        Parameters
        ----------
        y : str or list of str, default "pH"
            Name (or list of names) of the variables to plot. Valid variables
            are columns in the speciation report.

        title : str, optional
            Title of the plot.
            
        convert_log : bool, default True
            Convert units "log_activity", "log_molality", "log_gamma", and
            "log_fugacity" to "activity", "molality", "gamma", and "fugacity",
            respectively?

        plot_zero : bool, default True
            Plot zero values? Additionally, include series with all NaN (blank)
            values in the legend?
        
        show_missing : bool, default True
            Show samples that do not have bars?
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches.

        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
        
        colormap : str, default "WORM"
            Name of the colormap to color plotted data. Accepts "WORM",
            "colorblind", or matplotlib colormaps.
            See https://matplotlib.org/stable/tutorials/colors/colormaps.html
            The "colorblind" colormap is referenced from Wong, B. Points of view:
            Color blindness. Nat Methods 8, 441 (2011).
            https://doi.org/10.1038/nmeth.1618
            
        save_as : str, optional
            Provide a filename to save this figure. Filetype of saved figure is
            determined by `save_format`.
            Note: interactive plots can be saved by clicking the 'Download plot'
            button in the plot's toolbar.
        
        save_format : str, default "png"
            Desired format of saved or downloaded figure. Can be 'png', 'jpg',
            'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', or 'html'. If 'html',
            an interactive plot will be saved. Only 'png', 'svg', 'jpeg',
            and 'webp' can be downloaded with the 'download as' button in the
            toolbar of an interactive plot.

        save_scale : numeric, default 1
            Multiply title/legend/axis/canvas sizes by this factor when saving
            the figure.
        
        interactive : bool, default True
            Return an interactive plot if True or a static plot if False.
            
        plot_out : bool, default False
            Return a plotly figure object? If True, a plot is not displayed as
            it is generated.
            
        Returns
        -------
        fig : Plotly figure object
            A figure object is returned if `plot_out` is true. Otherwise, a
            figure is simply displayed.
        """
        
        if not isinstance(y, list):
            y = [y]

        colors = get_colors(colormap, len(y))

        # convert rgba to hex
        colors = [matplotlib.colors.rgb2hex(c) for c in colors]

        # map each species to its color, e.g.,
        # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
        dict_species_color = {sp:color for sp,color in zip(y, colors)}
        
        # html format color dict key names
        dict_species_color = {chemlabel(k):v for k,v in dict_species_color.items()}
            
        y_cols = self.lookup(y)

        if not show_missing:
            y_cols = y_cols.dropna(how='all') # this df will keep subheaders
        x = y_cols.index # names of samples

        y = [yi for yi in y if "limiting reactant" not in yi]

        if len(y) == 0:
            self.err_handler.raise_exception("There are no numeric variables to plot.")
        
        df = self.lookup(y).copy()
        if not show_missing:
            df = df.dropna(how='all') # this df will lose subheaders (flattened)
        df.loc[:, "name"] = df.index
        df.columns = df.columns.get_level_values(0)

        
        for i, yi in enumerate(y):

            y_col = y_cols.iloc[:, y_cols.columns.get_level_values(0)==yi]

            try:
                subheader = y_col.columns.get_level_values(1)[0]
            except:
                msg = ("Could not find '{}' ".format(yi)+"in the speciation "
                       "report. Available variables include "
                      "{}".format(list(set(self.report.columns.get_level_values(0)))))
                self.err_handler.raise_exception(msg)
            try:
                unit_type, unit = self.__get_unit_info(subheader)
            except:
                unit_type = ""
                unit = ""
                
            try:
                y_vals = [float(y0[0]) if y0[0] != 'NA' else float("nan") for y0 in y_col.values.tolist()]
            except:
                msg = ("One or more the values belonging to "
                       "'{}' are non-numeric and cannot be plotted.".format(y_col.columns.get_level_values(0)[0]))
                self.err_handler.raise_exception(msg)

            if convert_log and [abs(y0) for y0 in y_vals] != y_vals: # convert to bar-friendly units if possible
                if subheader in ["log_activity", "log_molality", "log_gamma", "log_fugacity"]:
                    y_plot, out_unit = self.__convert_aq_units_to_log_friendly(yi, rows=x)
                    unit_type, unit = self.__get_unit_info(out_unit)
                else:
                    y_plot = y_vals
            else:
                y_plot = y_vals

            if i == 0:
                subheader_previous = subheader
                unit_type_previous = unit_type
            if unit_type != unit_type_previous and i != 0:
                
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format(unit, yi_previous, unit_type_previous)+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)
            elif "activity" in subheader.lower() and "molality" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("activity", yi_previous, "molality")+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)
            elif "molality" in subheader.lower() and "activity" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("molality", yi_previous, "activity")+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)

            yi_previous = copy.deepcopy(yi)
            unit_type_previous = copy.deepcopy(unit_type)
            subheader_previous = copy.deepcopy(subheader)

            # Use pd.Series to avoid pandas 2.x string dtype assignment errors
            df[yi] = pd.Series(y_plot, index=df.index)


        if len(y) > 1:
            if unit != "":
                ylabel = "{} [{}]".format(unit_type, unit)
            else:
                ylabel = unit_type

        else:
            if 'pH' in y:
                ylabel = 'pH'
            elif 'Temperature' in y:
                ylabel = 'Temperature [°C]'
            else:
                if unit != "":
                    ylabel = "{} {} [{}]".format(chemlabel(y[0]), unit_type, unit)
                else:
                    ylabel = "{} {}".format(chemlabel(y[0]), unit_type)

        df = pd.melt(df, id_vars=["name"], value_vars=y)
        
        df = df.rename(columns={"Sample": "y_variable", "value": "y_value"})
        df = df.rename(columns={"variable": "y_variable"})
        
        df['y_variable'] = df['y_variable'].apply(chemlabel)

        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.reactions_for_plotting, pd.DataFrame):

            y_find = [yi.replace(" energy supply", "").replace(" affinity", "").replace(" Gibbs free energy", "") for yi in y]
            
            rxns = self.reactions_for_plotting.loc[y_find, :]["reaction"].tolist()
            
            # get the formatted reactions in the right order, then add as a
            # column in df
            formatted_rxn_list = []
            for rxn in rxns:
                for i in range(0,len(x)):
                    formatted_rxn_list.append(rxn)
            df["formatted_rxns"] = formatted_rxn_list

            if len(y) == 1:
                ylabel = "{}<br>{} [{}]".format(chemlabel(y_find[0]), unit_type, unit)

            if not plot_zero:
                df = df.dropna(subset=['y_value'])
                df = df[df.y_value != 0]
            
            # customdata for displaying reactions has to be here instead of in update_traces
            fig = px.bar(df, x="name", y="y_value",
                height=plot_height*ppi, width=plot_width*ppi,
                color='y_variable', barmode='group',
                labels={'y_value': ylabel}, template="simple_white",
                color_discrete_map=dict_species_color, custom_data=['formatted_rxns'])
            
            fig.update_traces(
                hovertemplate = "%{x} <br>"+ylabel+": %{y}<br>%{customdata}")

        else:

            if not plot_zero:
                df = df.dropna(subset=['y_value'])
                df = df[df.y_value != 0]
            
            fig = px.bar(df, x="name", y="y_value",
                height=plot_height*ppi, width=plot_width*ppi,
                color='y_variable', barmode='group',
                labels={'y_value': ylabel}, template="simple_white",
                color_discrete_map=dict_species_color)
            
            fig.update_traces(hovertemplate = "%{x} <br>"+ylabel+": %{y}")

        fig.update_layout(xaxis_tickangle=-45, xaxis_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'},
                          legend_title=None, margin={"t": 40},
                          xaxis={'fixedrange':True},
                          yaxis={'fixedrange':True, 'exponentformat':'power'})
        if len(y) == 1:
            fig.update_layout(showlegend=False)

        save_as, save_format = self._save_figure(fig, save_as, save_format,
                                                  save_scale, plot_width,
                                                  plot_height, ppi)
            
        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d',
                                             'lasso2d', 'zoomIn2d', 'zoomOut2d',
                                             'autoScale2d', 'resetScale2d',
                                             'toggleSpikelines'],
                  'toImageButtonOptions': {
                                             'format': save_format,
                                             'filename': save_as,
                                             'height': plot_height*ppi,
                                             'width': plot_width*ppi,
                                             'scale': save_scale,
                                           },
                  }
        if not interactive:
            config['staticPlot'] = True

        if plot_out:
            return fig
        else:
            fig.show(config=config)

        
    def scatterplot(self, x="pH", y="Temperature", samples=None, title=None,
                    log_x=False, log_y=False, plot_zero=True,
                    rxns_as_labels=True, charge_sign_at_end=False,
                    plot_width=4, plot_height=3, ppi=122,
                    fill_alpha=0.7, point_size=10,
                    ylab=None, lineplot=False, linemarkers=True,
                    colormap="WORM", save_as=None, save_format=None,
                    save_scale=1, interactive=True, plot_out=False):
        
        """
        Vizualize two or more sample variables with a scatterplot.
        
        Parameters
        ----------
        x, y : str, default for x is "pH", default for y is "Temperature"
            Names of the variables to plot against each other. Valid variables
            are columns in the speciation report. `y` can be a list of
            of variable names for a multi-series scatterplot.

        samples : list, optional
            List of samples to plot. By default, all samples in the speciation
            are plotted at once.

        title : str, optional
            Title of the plot.

        log_x, log_y : bool, default False
            Display the x_axis or y_axis in log scale?

        plot_zero : bool, default True
            Plot zero values? Additionally, include series with all NaN (blank)
            values in the legend?

        rxns_as_labels : bool, default True
            Display reactions as legend labels when plotting affinities and
            energy supplies?
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
        
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
        
        fill_alpha : numeric, default 0.7
            Transparency of scatterpoint area fill.
        
        point_size : numeric, default 10
            Size of scatterpoints.

        ylab : str, optional
            Custom label for the y-axis.
        
        lineplot : bool, default False
            Display a line plot instead of a scatterplot?

        linemarkers : bool, default True
            If `lineplot=True`, also plot markers?
        
        colormap : str, default "WORM"
            Name of the colormap to color the plotted data. Accepts "WORM",
            "colorblind", or matplotlib colormaps.
            See https://matplotlib.org/stable/tutorials/colors/colormaps.html
            The "colorblind" colormap is referenced from Wong, B. Points of view:
            Color blindness. Nat Methods 8, 441 (2011).
            https://doi.org/10.1038/nmeth.1618
            
        save_as : str, optional
            Provide a filename to save this figure. Filetype of saved figure is
            determined by `save_format`.
            Note: interactive plots can be saved by clicking the 'Download plot'
            button in the plot's toolbar.

        save_format : str, default "png"
            Desired format of saved or downloaded figure. Can be 'png', 'jpg',
            'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', or 'html'. If 'html',
            an interactive plot will be saved. Only 'png', 'svg', 'jpeg',
            and 'webp' can be downloaded with the 'download as' button in the
            toolbar of an interactive plot.
    
        save_scale : numeric, default 1
            Multiply title/legend/axis/canvas sizes by this factor when saving
            the figure.
        
        interactive : bool, default True
            Return an interactive plot if True or a static plot if False.
            
        plot_out : bool, default False
            Return a plotly figure object? If True, a plot is not displayed as
            it is generated.
            
        Returns
        -------
        fig : Plotly figure object
            A figure object is returned if `plot_out` is true. Otherwise, a
            figure is simply displayed.
        """

        if not isinstance(y, list) and not isinstance(y, str):
            self.err_handler.raise_exception("y must be a string or a list. "
                    "If you are applying the result of a lookup() to y, then "
                    "make sure the string inside of lookup() is valid.")
        
        if not isinstance(y, list):
            y = [y]
        
        if not isinstance(x, str):
            self.err_handler.raise_exception("x must be a string.")
        
        x_col = self.lookup(x)
        
        try:
            xsubheader = x_col.columns.get_level_values(1)[0]
        except:
            msg = ("Could not find '{}' ".format(x)+"in the speciation "
                   "report. Available variables include "
                   "{}".format(list(set(self.report.columns.get_level_values(0)))))
            self.err_handler.raise_exception(msg)
            
        try:
            x_plot = [float(x0[0]) if x0[0] != 'NA' else float("nan") for x0 in x_col.values.tolist()]
        except:
            msg = ("One or more the values belonging to "
                   "'{}' are non-numeric and cannot be plotted.".format(x_col.columns.get_level_values(0)[0]))
            self.err_handler.raise_exception(msg)
        
        try:
            xunit_type, xunit = self.__get_unit_info(xsubheader)
        except:
            xunit_type = ""
            xunit = ""

        for i, yi in enumerate(y):

            if "limiting reactant" in yi and len(y) > 1:
                continue
            
            y_col = self.lookup(yi)
            
            try:
                subheader = y_col.columns.get_level_values(1)[0]
            except:
                msg = ("Could not find '{}' ".format(yi)+"in the speciation "
                       "report. Available variables include "
                      "{}".format(list(set(self.report.columns.get_level_values(0)))))
                self.err_handler.raise_exception(msg)
            try:
                unit_type, unit = self.__get_unit_info(subheader)
            except:
                unit_type = ""
                unit = ""
            
            try:
                y_plot = [float(y0[0]) if y0[0] != 'NA' else float("nan") for y0 in y_col.values.tolist()]
            except:
                msg = ("One or more the values belonging to "
                       "'{}' are non-numeric and cannot be plotted.".format(y_col.columns.get_level_values(0)[0]))
                self.err_handler.raise_exception(msg)
                
            if i == 0:
                subheader_previous = subheader
                unit_type_previous = unit_type
            if unit_type != unit_type_previous and i != 0:
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format(unit_type, yi_previous, unit_type_previous)+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)
            elif "activity" in subheader.lower() and "molality" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("activity", yi_previous, "molality")+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)
            elif "molality" in subheader.lower() and "activity" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("molality", yi_previous, "activity")+""
                       "Plotted variables must share units.")
                self.err_handler.raise_exception(msg)
                
            yi_previous = copy.deepcopy(yi)
            unit_type_previous = copy.deepcopy(unit_type)
            subheader_previous = copy.deepcopy(subheader)

        if len(y) > 1:
            if unit != "":
                ylabel = "{} [{}]".format(unit_type, unit)
            else:
                ylabel = unit_type
        else:
            if 'pH' in y:
                ylabel = 'pH'
            elif 'Temperature' in y:
                ylabel = 'Temperature [°C]'
            else:
                y_formatted = chemlabel(y[0], charge_sign_at_end=charge_sign_at_end)
                if unit != "":
                    ylabel = "{} {} [{}]".format(y_formatted, unit_type, unit)
                else:
                    ylabel = "{} {}".format(y_formatted, unit_type)
        
        if x == 'pH':
            xlabel = 'pH'
        elif x == 'Temperature':
            xlabel = 'Temperature [°C]'
        else:
            x_formatted = chemlabel(x, charge_sign_at_end=charge_sign_at_end)
            if xunit != "":
                xlabel = "{} {} [{}]".format(x_formatted, xunit_type, xunit)
            else:
                xlabel = "{} {}".format(x_formatted, xunit_type)

        y = [yi for yi in y if "limiting reactant" not in yi]
        df = self.lookup([x]+y).copy() # TODO: is this where the "can't find name" message comes from?
        df.loc[:, "name"] = df.index
        df.columns = df.columns.get_level_values(0)
        
        if x in y:
            df[str(x)+"_copy"] = df[str(x)]
            df = df.melt(id_vars=["name", str(x)+"_copy"], value_vars=y).rename(columns={str(x)+'_copy': str(x)})
        else:
            df = pd.melt(df, id_vars=["name", x], value_vars=y)
            
        df = df.rename(columns={"Sample": "y_variable", "value": "y_value"})

        if not plot_zero:
            df = df.dropna(subset=['y_value'])
            df = df[df.y_value != 0]

        if isinstance(colormap, str):
            # get colors
            colors = get_colors(colormap, len(y), alpha=fill_alpha)
            
            # convert rgba to hex
            colors = [matplotlib.colors.rgb2hex(c) for c in colors]
    
            # map each species to its color, e.g.,
            # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
            dict_species_color = {sp:color for sp,color in zip(y, colors)}
            
            # html format color dict key names
            dict_species_color = {chemlabel(k, charge_sign_at_end=charge_sign_at_end):v for k,v in dict_species_color.items()}
        elif isinstance(colormap, list):
            colors = colormap
            dict_species_color = {sp:color for sp,color in zip(y, colors)}
        
        else:
            dict_species_color = {}

        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.reactions_for_plotting, pd.DataFrame):
            
            y_find = [yi.replace(" energy supply", "").replace(" affinity per mole rxn", "").replace(" affinity per mole e-", "").replace(" Gibbs free energy", "") for yi in y]
            y_find = [yi for yi in y_find if "limiting reactant" not in yi]

            rxns = self.reactions_for_plotting.loc[y_find, :]["reaction"].tolist()

            rxn_dict = {rxn_name:rxn for rxn_name,rxn in zip(y, rxns)}

            if len(y) == 1:
                ylabel = "{}<br>{} [{}]".format(chemlabel(y_find[0], charge_sign_at_end=charge_sign_at_end), unit_type, unit)
            
            df["formatted_rxn"] = df["y_variable"].map(rxn_dict)
        else:
            df["formatted_rxn"] = ""

        df['y_variable_original'] = df['y_variable']
        df['y_variable'] = df['y_variable'].apply(chemlabel, charge_sign_at_end=charge_sign_at_end)
        
        if ylab != None:
            ylabel=ylab

        if "log" in xlabel and log_x and xlabel != "pH":
            log_x = False
        if "log" in ylabel and log_y and ylabel != "pH":
            log_y = False

        if isinstance(samples, list):
            df = df.loc[df['name'].isin(samples)]
        
        if lineplot:
            df = df.sort_values(x).reset_index(drop=True)
            fig = px.line(df, x=x, y="y_value", color="y_variable",
                             log_x=log_x, log_y=log_y,
                             hover_data=[x, "y_value", "y_variable", "name", "formatted_rxn"],
                             width=plot_width*ppi, height=plot_height*ppi,
                             labels={x: xlabel,  "y_value": ylabel},
                             category_orders={"species": y},
                             color_discrete_map=dict_species_color,
                             custom_data=['name', 'formatted_rxn', 'y_variable_original'],
                             template="simple_white", markers=linemarkers)
        else:
            fig = px.scatter(df, x=x, y="y_value", color="y_variable",
                             log_x=log_x, log_y=log_y,
                             hover_data=[x, "y_value", "y_variable", "name", "formatted_rxn"],
                             width=plot_width*ppi, height=plot_height*ppi,
                             labels={x: xlabel,  "y_value": ylabel},
                             category_orders={"species": y},
                             color_discrete_map=dict_species_color,
                             opacity=fill_alpha,
                             custom_data=['name', 'formatted_rxn', 'y_variable_original'],
                             template="simple_white")
            
        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.reactions_for_plotting, pd.DataFrame):
            if rxns_as_labels:
                newnames = {y:r for y,r in zip(list(df["y_variable"]), list(df["formatted_rxn"]))}
                fig.for_each_trace(lambda t: t.update(name = newnames[t.name],
                                                      legendgroup = newnames[t.name],
                                                      hovertemplate = t.hovertemplate.replace(t.name, newnames[t.name])
                                                     ))
            fig.update_traces(marker=dict(size=point_size),
                              hovertemplate = "%{customdata[0]}<br>"+xlabel+": %{x} <br>"+ylabel+": %{y}<br>Reaction name: %{customdata[2]}<br>Reaction: %{customdata[1]}")
        
        else:
            fig.update_traces(marker=dict(size=point_size),
                              hovertemplate = "%{customdata[0]}<br>"+xlabel+": %{x} <br>"+ylabel+": %{y}<br>%{customdata[1]}")
        
        fig.update_layout(legend_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'},
                          margin={"t": 40},
                          yaxis={'exponentformat':'power'})
        if len(y) == 1:
            fig.update_layout(showlegend=False)

        fig.update_xaxes(exponentformat = 'E')
        fig.update_yaxes(exponentformat = 'E')
            
        save_as, save_format = self._save_figure(fig, save_as, save_format,
                                                  save_scale, plot_width,
                                                  plot_height, ppi)

        config = {'displaylogo': False, 'scrollZoom': True,
                  'modeBarButtonsToRemove': ['select2d', 'lasso2d', 'toggleSpikelines', 'resetScale2d'],
                  'toImageButtonOptions': {
                                           'format': save_format, # one of png, svg, jpeg, webp
                                           'filename': save_as,
                                           'height': plot_height*ppi,
                                           'width': plot_width*ppi,
                                           'scale': save_scale,
                                           },
                 }

        if not interactive:
            config['staticPlot'] = True
        
        if plot_out:
            return fig
        else:
            fig.show(config=config)

            
    def plot_mass_contribution(self, basis, title=None, sort_by=None,
                                     ascending=True, sort_y_by=None, width=0.9,
                                     colormap="WORM", sample_label="sample",
                                     colors=None,
                                     plot_width=4, plot_height=3, ppi=122,
                                     save_as=None, save_format=None,
                                     save_scale=1, interactive=True,
                                     plot_out=False):
        
        """
        Plot basis species contributions to mass balance of aqueous species
        across all samples.
        
        Parameters
        ----------
        basis : str
            Name of the basis species.

        title : str, optional
            Title of the plot.
            
        sort_by : str, optional
            Name of the variable used to sort samples. Variable names must be
            taken from the speciation report column names. No sorting is done by
            default.
        
        ascending : bool, default True
            Should sample sorting be in ascending order? Descending if False.
            Ignored unless `sort_by` is defined.
        
        sort_y_by : list of str or 'alphabetical', optional
            List of species names in the order that they should be stacked, from
            the bottom of the plot to the top. 'alphabetical' will sort species
            alphabetically.
        
        width : float, default 0.9
            Width of bars. No space between bars if width=1.0.
        
        colormap : str, default "WORM"
            Name of the colormap to color the scatterpoints. Accepts "WORM",
            "colorblind", or matplotlib colormaps.
            See https://matplotlib.org/stable/tutorials/colors/colormaps.html
            The "colorblind" colormap is referenced from Wong, B. Points of view:
            Color blindness. Nat Methods 8, 441 (2011).
            https://doi.org/10.1038/nmeth.1618
        
        sample_label : str, default "sample"
            Name of the label that appears when hovering over an element in the
            interactive mass contribution plot. By default, this is "sample".
            However, other words might be more appropriate to describe the
            calculations you are performing. For instance, if you are comparing
            reaction progress, `sample_label = "Xi"` might be more appropriate.
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
        
        save_as : str, optional
            Provide a filename to save this figure. Filetype of saved figure is
            determined by `save_format`.
            Note: interactive plots can be saved by clicking the 'Download plot'
            button in the plot's toolbar.

        save_format : str, default "png"
            Desired format of saved or downloaded figure. Can be 'png', 'jpg',
            'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', or 'html'. If 'html',
            an interactive plot will be saved. Only 'png', 'svg', 'jpeg',
            and 'webp' can be downloaded with the 'download as' button in the
            toolbar of an interactive plot.
    
        save_scale : numeric, default 1
            Multiply title/legend/axis/canvas sizes by this factor when saving
            the figure.
        
        interactive : bool, default True
            Return an interactive plot if True or a static plot if False.
            
        plot_out : bool, default False
            Return a plotly figure object? If True, a plot is not displayed as
            it is generated.
            
        Returns
        -------
        fig : Plotly figure object
            A figure object is returned if `plot_out` is true. Otherwise, a
            figure is simply displayed.
        """
        
        try:
            self.mass_contribution
        except:
            msg = ("Results for basis species contributions to aqueous mass "
                   "balance could not be found. Ensure that "
                   "get_mass_contribution = True when running speciate().")
            self.err_handler.raise_exception(msg)
            
        if basis not in set(self.mass_contribution['basis']):
            msg = ("The basis species {} ".format(basis)+"could not be found "
                   "among available basis species: "
                   "{}".format(str(list(set(self.mass_contribution['basis'])))))
            self.err_handler.raise_exception(msg)
            
        df_sp = copy.deepcopy(self.mass_contribution.loc[self.mass_contribution['basis'] == basis])
        
        if isinstance(sort_y_by, list):
            for species in sort_y_by:
                if species not in df_sp["species"]:
                    for sample in set(df_sp["sample"]):
                        df2 = pd.DataFrame({'sample':[sample], 'basis':[basis], 'species':[species], 'factor':[None], 'molality':[None], 'percent':[0]})
                        df_sp = pd.concat([df_sp, df2], ignore_index=True)
    
        if sort_by != None:
            if sort_by in self.report.columns.get_level_values(0):
                sort_col = self.lookup(sort_by)
                sort_by_unit = sort_col.columns.get_level_values(1)[0]
                sort_index = sort_col.sort_values([(sort_by, sort_by_unit)], ascending=ascending).index
                
                df_list = []
                for i in sort_index:
                    df_list.append(df_sp[df_sp['sample']==i])

                df_sp = pd.concat(df_list)
                
            else:
                msg = ("Could not find {}".format(sort_by)+" in the "
                       "speciation report. Available variables include "
                       "{}".format(list(self.report.columns.get_level_values(0))))
                self.err_handler.raise_exception(msg)
        
        df_sp['percent'] = df_sp['percent'].astype(float)
        
        unique_species = self.__unique(df_sp["species"])
        
        if "Other" in unique_species:

            unique_species.append(unique_species.pop(unique_species.index("Other")))
        
        labels = self.__unique(df_sp["sample"])

        bottom = np.array([0]*len(labels))

        if sort_y_by != None:
            if isinstance(sort_y_by, list):
                if len(unique_species) == len(sort_y_by):
                    if len([s for s in unique_species if s in sort_y_by]) == len(unique_species) and len([s for s in sort_y_by if s in unique_species]) == len(unique_species):
                        unique_species = sort_y_by
                    else:
                        valid_needed = [s for s in unique_species if s not in sort_y_by]
                        invalid = [s for s in sort_y_by if s not in unique_species]
                        msg = ("sort_y_by is missing the following species: "
                               "{}".format(valid_needed)+" and was provided "
                               "these invalid species: {}".format(invalid))
                        self.err_handler.raise_exception(msg)
                        
                elif len(sort_y_by) < len(unique_species):
                    msg = ("sort_y_by must have of all of the "
                           "following species: {}".format(unique_species)+". "
                           "You are missing {}".format([s for s in unique_species if s not in sort_y_by]))
                    self.err_handler.raise_exception(msg)
#                 else:
#                     msg = ("sort_y_by can only have the "
#                            "following species: {}".format(unique_species)+".")
#                     self.err_handler.raise_exception(msg)
            elif sort_y_by == "alphabetical":
                if "Other" in unique_species:
                    unique_species_no_other = [sp for sp in unique_species if sp != "Other"]
                    unique_species_no_other = sorted(unique_species_no_other)
                    unique_species = unique_species_no_other + ["Other"]
                else:
                    unique_species = sorted(unique_species)
            else:
                self.err_handler.raise_exception("sort_y_by must be either None, 'alphabetical', "
                                "or a list of species names.")
        
        if isinstance(colors, list):
            pass
        else:
            # get colormap
            colors = get_colors(colormap, len(unique_species))

            # convert rgba to hex
            colors = [matplotlib.colors.rgb2hex(c) for c in colors]

        df_sp["species"] = df_sp["species"].apply(chemlabel)
        unique_species = [chemlabel(sp) for sp in unique_species]
        
        if title == None:
            title = '<span style="font-size: 14px;">Species accounting for mass balance of {}</span>'.format(chemlabel(basis))
        
        
        # map each species to its color, e.g.,
        # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
        dict_species_color = {sp:color for sp,color in zip(unique_species, colors)}
        
        category_orders = {"species": unique_species, "sample": labels}


        fig = px.bar(df_sp, x="sample", y="percent", color="species",
                     width=plot_width*ppi, height=plot_height*ppi,
                     labels={"sample": sample_label,  "percent": "mole %", "species": "species"},
                     category_orders=category_orders,
                     color_discrete_map=dict_species_color,
                     template="simple_white",
                    )
        fig.update_layout(xaxis_tickangle=-45, xaxis_title=None, legend_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'},
                          margin={"t": 40}, bargap=0, xaxis={'fixedrange':True},
                          yaxis={'fixedrange':True})

        fig.update_traces(width=width, marker_line_width=0)
        
        save_as, save_format = self._save_figure(fig, save_as, save_format,
                                                  save_scale, plot_width,
                                                  plot_height, ppi)
            
        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d',
                                             'lasso2d', 'zoomIn2d', 'zoomOut2d',
                                             'autoScale2d', 'resetScale2d',
                                             'toggleSpikelines'],
                  'toImageButtonOptions': {
                                           'format': save_format, # one of png, svg, jpeg, webp
                                           'filename': save_as,
                                           'height': plot_height*ppi,
                                           'width': plot_width*ppi,
                                           'scale': save_scale,
                                           },
                 }
        
        if not interactive:
            config['staticPlot'] = True
        
        if plot_out:
            return fig
        else:
            fig.show(config=config)


    def plot_solid_solutions(self, sample, title=None,
                                   width=0.9, colormap="WORM",
                                   affinity_plot=True,
                                   affinity_plot_colors=["blue", "orange"],
                                   plot_width=4, plot_height=4, ppi=122,
                                   save_as=None, save_format=None,
                                   save_scale=1, interactive=True,
                                   plot_out=False):
        
        """
        Plot fractions of minerals of hypothetical solid solutions in a sample.
        
        Parameters
        ----------
        sample : str
            Name of the sample.

        title : str, optional
            Title of the plot.
        
        width : float, default 0.9
            Width of bars. No space between bars if width=1.0.
        
        colormap : str, default "WORM"
            Name of the colormap to color the scatterpoints. Accepts "WORM",
            "colorblind", or matplotlib colormaps.
            See https://matplotlib.org/stable/tutorials/colors/colormaps.html
            The "colorblind" colormap is referenced from Wong, B. Points of view:
            Color blindness. Nat Methods 8, 441 (2011).
            https://doi.org/10.1038/nmeth.1618
            
        affinity_plot : bool, default True
            Include the affinity subplot?
        
        affinity_plot_colors : list of two str, default ["blue", "orange"]
            Colors indicating positive and negative values in the affinity
            subplot, respectively.
            
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
        
        save_as : str, optional
            Provide a filename to save this figure. Filetype of saved figure is
            determined by `save_format`.
            Note: interactive plots can be saved by clicking the 'Download plot'
            button in the plot's toolbar.

        save_format : str, default "png"
            Desired format of saved or downloaded figure. Can be 'png', 'jpg',
            'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', or 'html'. If 'html',
            an interactive plot will be saved. Only 'png', 'svg', 'jpeg',
            and 'webp' can be downloaded with the 'download as' button in the
            toolbar of an interactive plot.
    
        save_scale : numeric, default 1
            Multiply title/legend/axis/canvas sizes by this factor when saving
            the figure.
        
        interactive : bool, default True
            Return an interactive plot if True or a static plot if False.
            
        plot_out : bool, default False
            Return a plotly figure object? If True, a plot is not displayed as
            it is generated.
            
        Returns
        -------
        fig : Plotly figure object
            A figure object is returned if `plot_out` is true. Otherwise, a
            figure is simply displayed.
        """

        if sample not in self.sample_data.keys():
            msg = ("The sample "+sample+" was not found in this speciation dataset."
                   " Samples with solid solutions in this dataset include:"+str([s for s in self.sample_data.keys() if "solid_solutions" in self.sample_data[s].keys()]))
            self.err_handler.raise_exception(msg)
        
        try:
            self.sample_data[sample]["solid_solutions"]
        except:
            msg = ("Results for solid solutions could not be found for this "
                   "sample. Samples with solid solutions in this speciation "
                   "dataset include:"+str([s for s in self.sample_data.keys() if "solid_solutions" in self.sample_data[s].keys()]))
            self.err_handler.raise_exception(msg)
        
        if title == None:
            title = "Hypothetical solid solutions in " + sample
        
        df_full = copy.deepcopy(self.sample_data[sample]["solid_solutions"])

        df = copy.deepcopy(df_full.dropna(subset=['x']))
        df = df[df['x'] != 0]

        unique_minerals = self.__unique(df["mineral"])
        
        # get colormap
        colors = get_colors(colormap, len(unique_minerals))
        
        # convert rgba to hex
        colors = [matplotlib.colors.rgb2hex(c) for c in colors]
        
        # map each species to its color, e.g.,
        # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
        dict_minerals_color = {sp:color for sp,color in zip(unique_minerals, colors)}

        solid_solutions = list(dict.fromkeys(df["solid solution"]))
        
        df_ss_only = df_full[df_full["x"].isnull()]
        
        mineral_dict = {m:[] for m in unique_minerals}
        for ss in solid_solutions:
            for m in unique_minerals:
                df_sub = df.loc[df["solid solution"] == ss,]
                frac = df_sub.loc[df_sub["mineral"] == m, "x"]
                if len(frac) > 0:
                    mineral_dict[m] = mineral_dict[m] + list(frac)
                else:
                    mineral_dict[m].append(0)

        if affinity_plot:
            rows = 2
            specs = [[{"type": "bar"}], [{"type": "bar"}]]
        else:
            rows = 1
            specs = [[{"type": "bar"}]]
                    
        fig = make_subplots(
            rows=rows, cols=1,
            specs=specs,
            vertical_spacing = 0.05
        )

        # subplot 1
        for m in unique_minerals[::-1]:
            fig.add_trace(go.Bar(name=m, x=solid_solutions, y=mineral_dict[m], marker_color=dict_minerals_color[m]), row=1, col=1)
        
        # subplot 2
        if affinity_plot:
            fig.add_trace(go.Bar(name="ss", x=solid_solutions, y=df_ss_only["Aff, kcal"],
                                 marker_color=[affinity_plot_colors[0] if val > 0 else affinity_plot_colors[1] for val in df_ss_only["Aff, kcal"]],
                                 showlegend=False),
                          row=2, col=1)

        fig.update_layout(barmode='stack', xaxis_tickangle=-45, xaxis_title=None, legend_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'}, autosize=False,
                          width=plot_width*ppi, height=plot_height*ppi,
                          margin={"t": 40}, bargap=0, xaxis={'fixedrange':True},
                          yaxis={'fixedrange':True}, template="simple_white")


        fig.update_xaxes(tickangle=-45)
        fig['layout']['yaxis']['title']='Mole Fraction'
        if affinity_plot:
            fig['layout']['yaxis2']['title']='Affinity, kcal/mol'
            fig.update_xaxes(showticklabels=False) # hide all the xticks
            fig.update_xaxes(showticklabels=True, row=2, col=1)
            
        
        save_as, save_format = self._save_figure(fig, save_as, save_format,
                                                  save_scale, plot_width,
                                                  plot_height, ppi)
            
        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'select2d',
                                             'lasso2d', 'zoomIn2d', 'zoomOut2d',
                                             'autoScale2d', 'resetScale2d',
                                             'toggleSpikelines'],
                  'toImageButtonOptions': {
                                           'format': save_format, # one of png, svg, jpeg, webp
                                           'filename': save_as,
                                           'height': plot_height*ppi,
                                           'width': plot_width*ppi,
                                           'scale': save_scale,
                                           },
                 }
        
        if not interactive:
            config['staticPlot'] = True

        if plot_out:
            return fig
        else:
            fig.show(config=config)
            
            
    def join_6i_p(self, filepath_6i, chain_mt, strip_minerals=False,
                  mixing_fluid_data0_lettercode=None,
                  mixing_fluid_pickup_bottom=None,
                  mixing_fluid_pickup_top=None,
                  mixing_fluid_eq_output_text=None,
                  data1_override=None):
        path='rxn_6i'
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            shutil.rmtree(path)
            os.makedirs(path)

        if chain_mt:
            raw_p_dict_bottom = self.raw_6_pickup_dict
        else:
            raw_p_dict_bottom = self.raw_3_pickup_dict_bottom

        # Detect whether cross-database reconciliation is needed
        fluid1_lettercode = getattr(self.thermo, 'data0_lettercode', None)
        needs_cross_db = (
            mixing_fluid_data0_lettercode is not None
            and fluid1_lettercode is not None
            and mixing_fluid_data0_lettercode != fluid1_lettercode
        )

        mix_with = getattr(filepath_6i, 'mix_with', None)
        mixing_sample_name = getattr(filepath_6i, 'mixing_sample_name', None)
        if mix_with is not None:
            samples_to_process = [s for s in raw_p_dict_bottom.keys()
                                  if s in mix_with]
        elif mixing_sample_name is not None:
            samples_to_process = [s for s in raw_p_dict_bottom.keys()
                                  if s != mixing_sample_name]
        else:
            samples_to_process = list(raw_p_dict_bottom.keys())

        for sample_name in samples_to_process:
            sample_filename = self.sample_data[sample_name]['filename'][:-3]

            if isinstance(filepath_6i, str):
                # if a string (filepath) is given
                with open(filepath_6i, "r") as f6i:
                    lines_6i = f6i.readlines()
            else:
                # if a Prepare_Reaction object is given
                all_lines = filepath_6i.formatted_reaction.split("\n")
                lines_6i = [e+"\n" for e in all_lines if e]

            # trim away any extra newlines at end of pre.6i, then add one.
            while lines_6i[-1] == "\n":
                lines_6i = lines_6i[:-1]

            if lines_6i[-1][-1:] != "\n": # \n counts as 1 character, not 2
                lines_6i[-1] = lines_6i[-1]+"\n"

            if needs_cross_db and mixing_fluid_pickup_bottom is not None:
                # Cross-database mixing.
                # Fluid 1 (bottom half) = self's pickup (the base speciation).
                # Fluid 2 (top half Composition/Reaction) = from the
                #   Prepare_Reaction template, which already contains the
                #   Mixing_Fluid's composition in the target database format.
                # The bottom half needs reconciliation to add basis species
                # that exist in the target database but not in self's.

                # Select the bottom half: prefer 6p (react result) over 3p.
                if (hasattr(self, 'raw_6_pickup_dict')
                        and sample_name in self.raw_6_pickup_dict):
                    lines_3p = list(self.raw_6_pickup_dict[sample_name])
                else:
                    lines_3p = list(raw_p_dict_bottom[sample_name])

                _warn_if_collapsed_totals(lines_3p)

                if strip_minerals:
                    lines_3p = _strip_minerals_from_bottom_half(lines_3p)

                # Determine the EQ3/6 output text for self (Fluid 1) to
                # use for species distribution when computing pseudo-element
                # mass balance totals.
                self_eq_output_text = None
                raw_6o = self.raw_6_output_dict.get(sample_name, [])
                if raw_6o:
                    self_eq_output_text = "".join(raw_6o)
                elif hasattr(self, 'raw_3_output_dict'):
                    raw_3o = self.raw_3_output_dict.get(sample_name, [])
                    if raw_3o:
                        self_eq_output_text = "\n".join(raw_3o)

                # Load the target database's data0 for reconciliation.
                target_lettercode = data1_override or fluid1_lettercode
                if (data1_override is not None
                        and data1_override != fluid1_lettercode):
                    data0_text = None
                else:
                    data0_text = getattr(self.thermo, 'data0_db', None)
                if data0_text is None and target_lettercode:
                    eq36da = getattr(self.thermo, 'eq36da', None)
                    if eq36da:
                        data0_path = os.path.join(eq36da,
                                                  "data0." + target_lettercode)
                        if os.path.exists(data0_path):
                            with open(data0_path) as f:
                                data0_text = f.read()
                    if data0_text is None:
                        data0_path = "data0." + target_lettercode
                        if os.path.exists(data0_path):
                            with open(data0_path) as f:
                                data0_text = f.read()
                if data0_text:
                    bottom_matches = (
                        fluid1_lettercode is not None
                        and fluid1_lettercode == target_lettercode)
                    f2_basis = None
                    if mixing_fluid_pickup_bottom is not None:
                        f2_basis = AqEquil._parse_bottom_basis_species(
                            mixing_fluid_pickup_bottom)
                    lines_6i, lines_3p, allocation_log = (
                        _reconcile_cross_database_mixing(
                            lines_6i, lines_3p, data0_text,
                            self_eq_output_text,
                            fluid2_eq_output_text=mixing_fluid_eq_output_text,
                            bottom_half_matches_target=bottom_matches,
                            fluid2_basis=f2_basis))
                    if allocation_log:
                        if not hasattr(self, 'cross_db_allocation_log'):
                            self.cross_db_allocation_log = {}
                        self.cross_db_allocation_log[sample_name] = (
                            allocation_log)
            else:
                lines_3p = raw_p_dict_bottom[sample_name]

                if strip_minerals:
                    lines_3p = _strip_minerals_from_bottom_half(lines_3p)

            lines_to_keep = []
            for line in lines_6i:
                if "Start of the bottom half of the input file" in line:
                    break
                else:
                    lines_to_keep.append(line)
            lines_to_keep += lines_3p
            
            if "{tval}" in "".join(lines_to_keep):
                # grab temperature
                for line in lines_3p:
                    if "Original temperature" in line:
                        o_t = line.split("|")[2]
                for i,line in enumerate(lines_to_keep):
                    if "{tval}" in line:
                        lines_to_keep[i] = line.format(tval=o_t)
                        
            if "{pval}" in "".join(lines_to_keep):
                # grab temperature
                for line in lines_3p:
                    if "Original pressure" in line:
                        o_p = line.split("|")[2]
                for i,line in enumerate(lines_to_keep):
                    if "{pval}" in line:
                        lines_to_keep[i] = line.format(tval=o_p)
            
            with open(path + "/" + sample_filename+".6i", "w") as f:
                f.writelines(lines_to_keep)


    def mt(self, sample):
        """
        Retrieve mass transfer results for a sample.

        Parameters
        ----------
        sample : str
            Name of the sample for which to retrieve mass transfer results.

        Returns
        -------
        An object of class `AqEquil.MassTransfer.Mass_Transfer`.
        """
        
        sample_data = getattr(self, "sample_data")

        if sample not in sample_data.keys():
            self.err_handler.raise_exception(("The sample '"+str(sample)+"'"
                    " was not found amongst samples with mass transfer"
                    " results: "+str(list(sample_data.keys()))))
        
        if "mass_transfer" in list(sample_data[sample].keys()):
            if sample_data[sample]["mass_transfer"] != None:
                return sample_data[sample]["mass_transfer"]
            
        msg = ("Mass transfer results are not stored for sample '"+sample+"'. "
              "This might be because the reaction calculation did not "
              "finish successfully or because the thermodynamic database "
              "is a data0 or data1 file without a supporting CSV file.")
        self.err_handler.raise_exception(msg)


    def create_input_file(self, db, filename="new_input.csv",
                          sample=None, xi=None, sample_name=None,
                          aux_basis=None, minimum_molality=None,
                          return_df=False):
        """
        Create a CSV input file from speciation or mass transfer results,
        compatible with a different thermodynamic database.

        Tallies the aqueous species distribution into the basis species of
        the target database by using its dissociation reactions.

        For mass transfer results, extracts the fluid at a chosen reaction
        progress (Xi). For basic speciation results (no mass transfer),
        extracts each sample directly. When ``sample`` is ``None``, all
        samples are included as separate rows.

        Parameters
        ----------
        db : str
            Target database lettercode (e.g., ``"fez"``) or path to a data0
            file (e.g., ``"data0.fez"``).

        filename : str, default "new_input.csv"
            Output CSV file path. Ignored when ``return_df`` is ``True``.

        sample : str, optional
            Name of the sample to extract. If ``None``, all samples are
            included. For mass transfer results, ``None`` auto-selects the
            sole sample or raises an error if multiple exist.

        xi : float, optional
            Reaction progress value at which to extract the fluid composition.
            Only used for mass transfer results. If ``None``, the last Xi
            value is used.

        sample_name : str, optional
            Name to write in the ``Sample`` column of the output CSV. Only
            used when a single sample is being extracted. If ``None``, the
            original sample name is used.

        aux_basis : list of str, optional
            Auxiliary basis species to include as separate columns in the
            output CSV. Their tallied molalities are subtracted from the
            corresponding donor basis species to conserve mass. For example,
            ``["NH4+"]`` would tally nitrogen species containing NH4+ into
            an NH4+ column and subtract that amount from the NO3- (donor)
            column. If ``None`` (the default), only strict basis species
            are included and no reallocation is performed.

        minimum_molality : float, optional
            Minimum molality for allocated basis or auxiliary basis species.
            Any tallied species value that is positive but below this
            threshold will be raised to this value. For example,
            ``minimum_molality=1E-18`` prevents extremely small values like
            ``4.25E-41`` from appearing in the output. If ``None`` (the
            default), no minimum is applied.

        return_df : bool, default False
            If ``True``, return a pandas DataFrame instead of writing a CSV
            file. The DataFrame has species columns as headers, with the
            first row containing unit subheaders (compatible with
            ``join_input_files``). If ``False`` (the default), write a CSV
            file and return its path.

        Returns
        -------
        str or pandas.DataFrame
            The path to the created CSV file, or a DataFrame if
            ``return_df`` is ``True``.
        """
        import csv as csv_module

        verbose = getattr(self, "verbose", 1)
        all_sample_data = getattr(self, "sample_data", {})

        # --- Locate and read the target data0 file ---
        target_data0_text = None
        if os.path.exists(db):
            with open(db) as f:
                target_data0_text = f.read()
        else:
            candidate = db if db.startswith("data0.") else "data0." + db
            for search_dir in [os.getcwd(),
                               getattr(self, "eq36da", None) or ""]:
                path = os.path.join(search_dir, candidate)
                if os.path.exists(path):
                    with open(path) as f:
                        target_data0_text = f.read()
                    break
        if target_data0_text is None:
            self.err_handler.raise_exception(
                "Could not find target data0 file for '" + str(db)
                + "'. Place the file in the current working directory or "
                "specify the full path.")

        dissoc = AqEquil._parse_data0_dissociation_reactions(target_data0_text)
        basis_stoich, _ = _parse_data0_basis_species_full(target_data0_text)
        basis_set = set(basis_stoich.keys())

        # Identify auxiliary basis species from the target data0 file
        all_aux_in_db = set()
        _lines = target_data0_text.split("\n")
        _aux_start = None
        _aq_start = None
        for _i, _line in enumerate(_lines):
            _s = _line.strip()
            if _s == "auxiliary basis species":
                _aux_start = _i
            elif _s == "aqueous species":
                _aq_start = _i
                break
        if _aux_start is not None:
            _end = _aq_start if _aq_start is not None else len(_lines)
            _j = _aux_start + 2
            while _j < _end:
                _sl = _lines[_j].strip()
                if not _sl or _sl.startswith("+--"):
                    _j += 1
                    continue
                all_aux_in_db.add(_sl)
                _j += 1
                while _j < _end and not _lines[_j].strip().startswith("+--"):
                    _j += 1
                if _j < _end:
                    _j += 1

        # Only include user-requested auxiliary basis species
        aux_set = set()
        aux_donor = {}
        if aux_basis is not None:
            for sp in aux_basis:
                if sp not in all_aux_in_db:
                    self.err_handler.raise_exception(
                        "'" + str(sp) + "' is not an auxiliary basis species "
                        "in the target database. Available: "
                        + str(sorted(all_aux_in_db)))
                aux_set.add(sp)

            # Map each requested aux species to its donor (the strict basis
            # species in its dissociation reaction, excluding system species).
            system_species = {"H2O", "H+", "O2(g)", "O2", "e-"}
            for aux_sp in aux_set:
                rxn = dissoc.get(aux_sp, {})
                for product, coeff in rxn.items():
                    if product == aux_sp:
                        continue
                    if product in basis_set and product not in system_species:
                        if coeff > 0:
                            aux_donor[aux_sp] = product
                            break

        tally_set = basis_set | aux_set

        # --- Determine which samples to process ---
        if sample is not None:
            if sample not in all_sample_data:
                self.err_handler.raise_exception(
                    "Sample '" + str(sample) + "' not found among: "
                    + str(list(all_sample_data.keys())))
            samples_to_process = [sample]
        else:
            samples_to_process = list(all_sample_data.keys())

        exclude = {"H2O", "H+", "O2(g)", "O2", "e-"}
        all_csv_species = set()
        rows = []

        for s in samples_to_process:
            sd = all_sample_data[s]
            has_mt = ("mass_transfer" in sd
                      and sd["mass_transfer"] is not None)

            if has_mt:
                if len(samples_to_process) > 1 and sample is None:
                    self.err_handler.raise_exception(
                        "Multiple samples with mass transfer results found. "
                        "Specify one with the 'sample' parameter: "
                        + str(samples_to_process))

                mt = sd["mass_transfer"]
                xi_values = mt.misc_params["Xi"].values
                if xi is None:
                    xi_idx = len(xi_values) - 1
                    xi_used = xi_values[xi_idx]
                else:
                    xi_idx = int(np.argmin(np.abs(xi_values - xi)))
                    xi_used = xi_values[xi_idx]

                aq_row = mt.aq_distribution_molal.iloc[xi_idx]
                species_distribution = {}
                for col in mt.aq_distribution_molal.columns:
                    if col == "Xi":
                        continue
                    val = aq_row[col]
                    if pd.notna(val) and val > 0:
                        species_distribution[col] = float(val)

                misc_row = mt.misc_params.iloc[xi_idx]
                temperature = misc_row["Temp(C)"]
                pressure = misc_row["Press(bars)"]
                pH = misc_row["pH"]
                logfO2 = misc_row["log fO2"]

            else:
                aq_dist = sd["aq_distribution"]
                molality_series = aq_dist["molality"]
                species_distribution = {
                    sp: float(mol) for sp, mol in molality_series.items()
                    if pd.notna(mol) and mol > 0}

                temperature = sd["temperature"]
                pressure = sd["pressure"]
                pH = -float(aq_dist["log_activity"]["H+"])
                fugacity = sd.get("fugacity")
                if fugacity is not None and "log_fugacity" in fugacity.columns:
                    logfO2 = float(fugacity["log_fugacity"]["O2(g)"])
                else:
                    logfO2 = ""

            # --- Tally molalities to target basis + aux basis species ---
            tally = {b: 0.0 for b in tally_set}
            unmatched = []

            for sp_name, molality in species_distribution.items():
                rxn = dissoc.get(sp_name)
                if rxn is None:
                    unmatched.append(sp_name)
                    continue
                # Basis and auxiliary basis species tally directly
                if sp_name in tally_set:
                    tally[sp_name] += molality
                    continue
                # Aqueous species: dissociate into basis/aux products
                self_coeff = rxn.get(sp_name, None)
                norm = abs(self_coeff) if self_coeff is not None else 1.0
                if norm == 0:
                    norm = 1.0
                for product_sp, coeff in rxn.items():
                    if product_sp == sp_name:
                        continue
                    if product_sp in tally_set:
                        tally[product_sp] += (coeff / norm) * molality

            # Subtract auxiliary basis tallies from their donors
            for aux_sp, donor_sp in aux_donor.items():
                if donor_sp in tally and tally[aux_sp] > 0:
                    tally[donor_sp] -= tally[aux_sp]

            if unmatched and verbose > 0:
                print("Warning (" + s + "): " + str(len(unmatched))
                      + " species not found in the target database"
                      " and skipped"
                      + (": " + str(unmatched[:10])
                         + ("..." if len(unmatched) > 10 else "")
                         if verbose > 1 else "")
                      + ".")

            csv_species = {}
            for sp, val in tally.items():
                if sp in exclude or val <= 0:
                    continue
                if minimum_molality is not None and val < minimum_molality:
                    val = minimum_molality
                csv_species[sp] = val
            all_csv_species.update(csv_species.keys())

            row_name = sample_name if sample_name is not None else s
            rows.append((row_name, pH, pressure, temperature, logfO2,
                         csv_species))

        # --- Build output ---
        sorted_species = sorted(all_csv_species)

        header = ["Sample", "H+", "Pressure", "Temperature", "logfO2"]
        units = ["id", "pH", "bar", "degC", "logfO2"]
        header += sorted_species
        units += ["Molality"] * len(sorted_species)

        out_data_rows = []
        for row_name, pH, pressure, temperature, logfO2, csv_sp in rows:
            data = [row_name, str(pH), str(pressure),
                    str(temperature), str(logfO2)]
            for sp in sorted_species:
                val = csv_sp.get(sp, 0.0)
                data.append("{:.4E}".format(val) if val > 0 else "0")
            out_data_rows.append(data)

        if return_df:
            df_out = pd.DataFrame([units] + out_data_rows, columns=header)
            return df_out

        with open(filename, "w", newline="") as f:
            writer = csv_module.writer(f)
            writer.writerow(header)
            writer.writerow(units)
            for data in out_data_rows:
                writer.writerow(data)

        if verbose > 0:
            has_mt = any("mass_transfer" in all_sample_data[s]
                         and all_sample_data[s]["mass_transfer"] is not None
                         for s in samples_to_process)
            xi_msg = (" at Xi=" + str(xi_used)) if has_mt else ""
            print("Created input file '" + filename + "' with "
                  + str(len(sorted_species)) + " species columns and "
                  + str(len(rows)) + " sample(s) for database '"
                  + str(db) + "'" + xi_msg + ".")

        return filename


    def rename_samples(self, names):
        """
        Rename samples throughout this Speciation object.

        Parameters
        ----------
        names : dict
            Mapping of ``{old_name: new_name}``. Only samples present in the
            dict are renamed; others are left unchanged.

        Raises
        ------
        Exception
            If an old name is not found, or if renaming would produce
            duplicate sample names.
        """
        if not isinstance(names, dict) or len(names) == 0:
            self.err_handler.raise_exception(
                "names must be a non-empty dict mapping old names to new names.")

        existing = set(self.sample_data.keys())
        for old in names:
            if old not in existing:
                self.err_handler.raise_exception(
                    "Sample '" + str(old) + "' not found among: "
                    + str(list(existing)))

        new_all = []
        for s in existing:
            new_all.append(names.get(s, s))
        if len(set(new_all)) != len(new_all):
            self.err_handler.raise_exception(
                "Renaming would produce duplicate sample names: "
                + str(new_all))

        def _rename_dict_keys(d):
            if not isinstance(d, dict):
                return
            for old, new in names.items():
                if old in d:
                    d[new] = d.pop(old)

        # 1. sample_data dict keys + nested "name" field
        for old, new in names.items():
            if old in self.sample_data:
                self.sample_data[old]["name"] = new
        _rename_dict_keys(self.sample_data)

        # 2. Raw EQ3/6 file dicts
        for attr in ["raw_3_input_dict", "raw_3_output_dict",
                      "raw_3_pickup_dict_top", "raw_3_pickup_dict_bottom",
                      "raw_6_input_dict", "raw_6_output_dict",
                      "raw_6_pickup_dict", "raw_6_pickup_dict_top"]:
            d = getattr(self, attr, None)
            if d is not None:
                _rename_dict_keys(d)

        # 3. DataFrames with sample names as row index
        for attr in ["report", "input", "aq_distribution_logact",
                      "aq_distribution_molal", "misc_params"]:
            df = getattr(self, attr, None)
            if df is not None and isinstance(df, pd.DataFrame):
                df.rename(index=names, inplace=True)

        # 4. mass_contribution "sample" column
        mc = getattr(self, "mass_contribution", None)
        if mc is not None and isinstance(mc, pd.DataFrame):
            if "sample" in mc.columns:
                mc["sample"] = mc["sample"].replace(names)
            for old, new in names.items():
                sd = self.sample_data.get(new)
                if sd is not None and "mass_contribution" in sd:
                    sd_mc = sd["mass_contribution"]
                    if (isinstance(sd_mc, pd.DataFrame)
                            and "sample" in sd_mc.columns):
                        sd_mc["sample"] = sd_mc["sample"].replace(names)

