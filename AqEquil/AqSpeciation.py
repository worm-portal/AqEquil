DEBUGGING_R = False
FIXED_SPECIES = ["H2O", "H+", "O2(g)", "water", "Cl-", "e-", "OH-", "O2", "H2O(g)"]
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

from ipywidgets import IntProgress
from IPython.display import display
import time

from urllib.request import urlopen
from io import StringIO

import warnings

warnings.simplefilter(action='ignore', category=FutureWarning) # TEMPORARY! Disable this once FutureWarning issues have been solved.

import subprocess
import pkg_resources
import pandas as pd
import numpy as np
from chemparse import parse_formula
from IPython.core.display import display, HTML
import periodictable
from collections import Counter

import pyCHNOSZ
from ._HKF_cgl import OBIGT2eos, calc_logK

# matplotlib for static plots
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotly.express as px
import plotly.io as pio

from plotly.subplots import make_subplots
import plotly.graph_objects as go

# rpy2 for Python and R integration
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # will display errors, but not warnings

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


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
    Return a list of duplicate elements in another list
    https://stackoverflow.com/questions/46554866/efficiently-finding-duplicates-in-a-list
    """
    c = Counter(array)
    return [k for k in c if c[k] > 1]


def _float_to_fraction (x, error=0.000001):
    
    """
    Convert a float into a fraction. Works with floats like 2.66666666.
    Solution from https://stackoverflow.com/a/5128558/8406195
    """
    n = int(math.floor(x))
    x -= n
    if x < error:
        return (n, 1)
    elif 1 - error < x:
        return (n+1, 1)

    # The lower fraction is 0/1
    lower_n = 0
    lower_d = 1
    # The upper fraction is 1/1
    upper_n = 1
    upper_d = 1
    while True:
        # The middle fraction is (lower_n + upper_n) / (lower_d + upper_d)
        middle_n = lower_n + upper_n
        middle_d = lower_d + upper_d
        # If x + error < middle
        if middle_d * (x + error) < middle_n:
            # middle is our new upper
            upper_n = middle_n
            upper_d = middle_d
        # Else If middle < x - error
        elif middle_n < (x - error) * middle_d:
            # middle is our new lower
            lower_n = middle_n
            lower_d = middle_d
        # Else middle is our best fraction
        else:
            return (n * middle_d + middle_n, middle_d)

    
def _float_to_formatted_fraction(x, error=0.000001):
    
    """
    Format a fraction for html.
    """
    f = _float_to_fraction(x, error=error)
    
    whole_number_float = int((f[0]-(f[0]%f[1]))/f[1])
    remainder_tuple = (f[0]%f[1], f[1])
    
    if remainder_tuple[0] == 0:
        return str(whole_number_float)
    else:
        if whole_number_float == 0:
            whole_number_float = ""
        return "{0}<sup>{1}</sup>&frasl;<sub>{2}</sub>".format(
                whole_number_float, remainder_tuple[0], remainder_tuple[1])


def _format_coeff(coeff):
    
    """
    Format a reaction coefficient for html.
    """
    if coeff == 1 or coeff == -1:
        coeff = ""
    elif coeff.is_integer() and coeff < 0:
        coeff = str(-int(coeff))
    elif coeff.is_integer() and coeff > 0:
        coeff = str(int(coeff))
    else:
        if coeff < 0:
            coeff = _float_to_formatted_fraction(-coeff)
        else:
            coeff = _float_to_formatted_fraction(coeff)

    if coeff != "":
        coeff = coeff + " "

    return coeff


def _convert_to_RVector(value, force_Rvec=True):
    
    """
    Convert a value or list into an R vector of the appropriate type.
    
    Parameters
    ----------
    value : numeric or str, or list of numeric or str
        Value to be converted.
    
    force_Rvec : bool, default True
        If `value` is not a list, force conversion into a R vector?
        False will return an int, float, or str if value is non-list.
        True will always return an R vector.
    
    Returns
    -------
    int, float, str, or an rpy2 R vector
        A value or R vector of an appropriate data type.
    """

    if not isinstance(value, list) and not force_Rvec:
        return value
    elif not isinstance(value, list) and force_Rvec:
        value = [value]
    else:
        pass

    if all(isinstance(x, bool) for x in value):
        return ro.BoolVector(value)
    elif all(isinstance(x, int) for x in value):
        return ro.IntVector(value)
    elif all(isinstance(x, float) or isinstance(x, int) for x in value):
        return ro.FloatVector(value)
    else:
        return ro.StrVector([str(v) for v in value])

    
def _clean_rpy2_pandas_conversion(
        df,
        float_cols=["G", "H", "S", "Cp",
                    "V", "a1.a", "a2.b",
                    "a3.c", "a4.d", "c1.e",
                    "c2.f", "omega.lambda", "z.T",
                    "azero", "neutral_ion_type", "regenerate_dissrxn",
                    "logK1", "logK2", "logK3", "logK4",
                    "logK5", "logK6", "logK7", "logK8",
                    "T1", "T2", "T3", "T4", "T5", "T6",
                    "T7", "T8"],
        str_cols=["name", "abbrv", "state", "formula",
                  "ref1", "ref2", "date",
                  "E_units", "tag", "dissrxn", "formula_ox",
                  "formula_modded", "formula_ox_modded", 
                  "P1", "P2", "P3", "P4", "P5", "P6",
                  "P7", "P8"],
        NA_string=""):

    df.replace(NA_string, np.nan, inplace=True)
    for col in float_cols:
        if col in df.columns:
            df[col] = df[col].astype(float)
    for col in str_cols:
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df
    

def _get_colors(colormap, ncol, alpha=1.0, hide_traceback=True):

    """
    Get a list of rgb values for a matplotlib colormap
    
    Parameters
    ----------
    colormap : str
        Name of the colormap to color the scatterpoints. Accepts "WORM",
        "colorblind", or matplotlib colormaps.
        See https://matplotlib.org/stable/tutorials/colors/colormaps.html
        The "colorblind" colormap is referenced from Wong, B. Points of view:
        Color blindness. Nat Methods 8, 441 (2011).
        https://doi.org/10.1038/nmeth.1618
    
    ncol : int
        Number of colors to return in the list.
    
    alpha : float, default 1.0
        An alpha value between 0.0 (transparent) and 1.0 (opaque).

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this function?
        When True, error messages handled by this class will be short and to
        the point.
    
    Returns
    -------
    colors : list
        A list of rgb color tuples
    """
    
    err_handler = Error_Handler(clean=hide_traceback)
    
    qualitative_cmaps = ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                         'Dark2', 'Set1', 'Set2', 'Set3',
                         'tab10', 'tab20', 'tab20b', 'tab20c']
    
    if colormap == "colorblind":
        # colors from Wong B. 2011, https://doi.org/10.1038/nmeth.1618
        colors = [(0, 0, 0, alpha), # black
                  (230/255, 159/255, 0, alpha), # orange
                  (86/255, 180/255, 233/255, alpha), # sky blue
                  (0, 158/255, 115/255, alpha), # bluish green
                  (240/255, 228/255, 66/255, alpha), # yellow
                  (0, 114/255, 178/255, alpha), # blue
                  (213/255, 94/255, 0, alpha), # vermillion
                  (204/255, 121/255, 167/255, alpha)] # reddish purple
        if ncol <= len(colors):
            return colors[:ncol]
        else:
            print("Switching from 'colorblind' colormap to 'viridis' because there are {} variables to plot.".format(ncol))
            colormap = "viridis"
    elif colormap == "WORM":
        colors = [(0, 0, 0, alpha), # black
                  (22/255, 153/255, 211/255, alpha), # blue
                  (232/255, 86/255, 66/255, alpha), # red
                  (245/255, 171/255, 80/255, alpha), # orange
                  (115/255, 108/255, 168/255, alpha), # purple
                  (151/255, 208/255, 119/255, alpha), # green
                  (47/255, 91/255, 124/255, alpha), # dark blue
                  (119/255, 119/255, 119/255, alpha)] # gray
        if ncol <= len(colors):
            return colors[:ncol]
        else:
            print("Switching from 'WORM' colormap to 'viridis' because there are {} variables to plot.".format(ncol))
            colormap = "viridis"
            
    if colormap in qualitative_cmaps:
        # handle qualitative (non-continuous) colormaps
        colors = [plt.cm.__getattribute__(colormap).colors[i] for i in range(ncol)]
        colors = [(c[0], c[1], c[2], alpha) for c in colors]
    else:
        # handle sequential (continuous) colormaps
        norm = matplotlib.colors.Normalize(vmin=0, vmax=ncol-1)
        try:
            cmap = cm.__getattribute__(colormap)
        except:
            valid_colormaps = [cmap for cmap in dir(cm) if "_" not in cmap and cmap not in ["LUTSIZE", "MutableMapping", "ScalarMappable", "functools", "datad", "revcmap"]]
            err_handler.raise_exception("'{}'".format(colormap)+" is not a recognized matplotlib colormap. "
                    "Try one of these: {}".format(valid_colormaps))
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = [m.to_rgba(i) for i in range(ncol)]
        colors = [(c[0], c[1], c[2], alpha) for c in colors]
    
    return colors


def _all_equal(iterable):
    # check that all elements of a list are equal
    g = itertools.groupby(iterable)
    return next(g, True) and not next(g, False)


def check_balance(formulas, stoich):
    """
    Check that a chemical reaction is balanced. If not, get missing composition.
    
    Parameters
    ----------
    formulas : list of str
        A list of species formulas that match the order of
        the stoichiometric reaction coefficients in the `stoich` parameter.
    
    stoich : list of numeric
        A list of stoichiometric reaction coefficients that match the order of
        the species formulas in the `formulas` parameter. Reactants are
        negative.
        
    Returns
    -------
    A printed warning and a dictionary of the missing composition if the
    reaction is unbalanced.
    """
    
    if len(formulas) != len(stoich):
        raise Exception("The number of species formulas does not match the "
              "number of stoichiometric coefficients in the reaction.")
    
    # sum all elements, +, and - by their reaction coefficient
    all_dict = {}
    for i,s in enumerate(formulas):
        s_dict = parse_formula(s)
        s_dict = {key: stoich[i]*s_dict[key] for key in s_dict.keys()}
        all_dict = {k: all_dict.get(k, 0) + s_dict.get(k, 0) for k in set(all_dict) | set(s_dict)}
    
    # sum + and - as Z (charge)
    if "+" not in list(all_dict.keys()):
        all_dict["+"] = 0
    if "-" not in list(all_dict.keys()):
        all_dict["-"] = 0
    all_dict["Z"] = all_dict["+"] - all_dict["-"]
    del all_dict["+"]
    del all_dict["-"]
    
    # delete all elements with a value of 0 (balanced)
    for key in list(all_dict.keys()):
        if all_dict[key] == 0:
            del all_dict[key]
    
    # print warnings, prepare missing composition dictionary
    if len(list(all_dict.keys())) > 0:
        missing_composition_dict = {k:[-all_dict[k]] for k in all_dict.keys()}
        print("Warning! The reaction is unbalanced. It is missing this composition:")
        print(pd.DataFrame(missing_composition_dict).to_string(index=False))
    else:
        missing_composition_dict = {}
    
    return missing_composition_dict
        

def chemlabel(name, charge_sign_at_end=False):
    
    """
    Format a chemical formula to display subscripts and superscripts in HTML
    (e.g., Plotly plots)
    Example, "CH3COO-" becomes "CH<sub>3</sub>COO<sup>-</sup>"
    
    Parameters
    ----------
    name : str
        A chemical formula.
    
    charge_sign_at_end : bool, default False
        Display charge with sign after the number (e.g. SO4 2-)?
        
    
    Returns
    -------
    A formatted chemical formula string.
    """
    
    # format only the first part of the name if it has "_(input)"
    if len(name.split("_(input)"))==2:
        if name.split("_(input)")[1] == '':
            name = name.split("_(input)")[0]
            input_flag=True
    else:
        input_flag = False
    
    name = _html_chemname_format(name, charge_sign_at_end=charge_sign_at_end)
    
    # add " (input)" to the end of the name
    if input_flag:
        name = name+" (input)"
    
    return(name)


def format_equation(species, stoich, charge_sign_at_end=False):
    """
    Format a chemical equation to display in HTML
    (e.g., Plotly plots)
    
    Parameters
    ----------
    species : list of str
        List of species in the reaction
        
    stoich : list of numeric
        List of stoichiometric reaction coefficients (reactants are negative)
    
    charge_sign_at_end : bool, default False
        Display charge with sign after the number (e.g. SO4 2-)?
        
    
    Returns
    -------
    A formatted chemical formula string.
    """
    reactants_list = []
    products_list = []
    for i,s in enumerate(species):
        s_f = chemlabel(s, charge_sign_at_end=charge_sign_at_end)
        if stoich[i] < 0:
            if stoich[i] != -1:
                entry = str(abs(stoich[i])) + " " + s_f
            else:
                entry = s_f
            reactants_list.append(entry)
        elif stoich[i] > 0:
            if stoich[i] != 1:
                entry = str(stoich[i]) + " " + s_f
            else:
                entry = s_f
            products_list.append(entry)
    
    reactants_together = " + ".join(reactants_list)
    products_together = " + ".join(products_list)
    
    equation_str = " → ".join([reactants_together, products_together])
    
    return equation_str


def _html_chemname_format(name, charge_sign_at_end=False):
    
    """
    Function duplicated from pyCHNOSZ
    """
    
    p = re.compile(r'(?P<sp>[-+]\d*?$)')
    name = p.sub(r'<sup>\g<sp></sup>', name)
    charge = re.search(r'<.*$', name)

    name_no_charge = re.match(r'(?:(?!<|$).)*', name).group(0)
    mapping = {"0": "<sub>0</sub>", "1": "<sub>1</sub>", "2": "<sub>2</sub>",
               "3": "<sub>3</sub>", "4": "<sub>4</sub>", "5": "<sub>5</sub>",
               "6": "<sub>6</sub>", "7": "<sub>7</sub>", "8": "<sub>8</sub>",
               "9": "<sub>9</sub>", ".":"<sub>.</sub>"}
    name_no_charge_formatted = "".join([mapping.get(x) or x
                                        for x in list(name_no_charge)])

    if charge != None:
        name = name_no_charge_formatted + charge.group(0)
    else:
        name = name_no_charge_formatted

    if charge_sign_at_end:
        if "<sup>-" in name:
            name = name.replace("<sup>-", "<sup>")
            name = name.replace("</sup>", "-</sup>")
        if "<sup>+" in name:
            name = name.replace("<sup>+", "<sup>")
            name = name.replace("</sup>", "+</sup>")

    return(name)
    

def _isnotebook():
    
    """
    Check if this code is running in a Jupyter notebook
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


class Error_Handler:
    
    """
    Handles how errors are printed in Jupyter notebooks. By default, errors that
    are handled by AqEquil are printed with an error message, but no traceback.
    Errors that are not handled by AqEquil, such as those thrown if the user
    encounters a bug, will display a full traceback.
    
    If the error handler prints an error message without traceback, all future
    errors regardless of origin will be shown without traceback until the
    notebook kernel is restarted.
    
    Parameters
    ----------
    clean : bool
        Report exceptions without traceback? If True, only the error message is
        shown. If False, the entire error message, including traceback, is
        shown. Ignored if AqEquil is not being run in a Jupyter notebook.
    
    """
    def __init__(self, clean=True):
        self.clean = clean # bool: hide traceback?
        pass
    
    
    @staticmethod
    def hide_traceback(exc_tuple=None, filename=None, tb_offset=None,
                       exception_only=False, running_compiled_code=False):
        
        """
        Return a modified ipython showtraceback function that does not display
        traceback when encountering an error.
        """
        
        ipython = get_ipython()
        etype, value, tb = sys.exc_info()
        value.__cause__ = None  # suppress chained exceptions
        return ipython._showtraceback(etype, value, ipython.InteractiveTB.get_exception_only(etype, value))
        

    def raise_exception(self, msg):
        
        """
        Raise an exception that displays the error message without traceback. This
        happens only when the exception is predicted by the AqEquil package
        (e.g., for common user errors).
        """
        if self.clean and _isnotebook():
            ipython = get_ipython()
            ipython.showtraceback = self.hide_traceback
            
        raise Exception(msg)
        
    
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
    
    db : str, default "WORM"
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
        "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
        
    solid_solutions : str
        Filepath of a CSV file containing parameters for solid solutions, e.g.,
        "my_solid_solutions.csv". If `db` is set to "WORM" and `solid_solutions`
        is not defined, then parameters for solid solutions will be retrieved
        from "Solid_solutions.csv" at https://github.com/worm-portal/WORM-db
    
    logK : str
        Filepath of a CSV file containing equilibrium constants for chemical
        species, e.g., "my_logK_entries.csv". If `db` is set to "WORM" and `logK`
        is not defined, then equilibrium constants will be retrieved from
        "wrm_data_logK.csv" at https://github.com/worm-portal/WORM-db
    
    logK_S : str
        Filepath of a CSV file containing equilibrium constants for chemical
        species, e.g., "my_logK_S_entries.csv". If `db` is set to "WORM" and `logK_S`
        is not defined, then equilibrium constants will be retrieved from
        "wrm_data_logK_S.csv" at https://github.com/worm-portal/WORM-db
    
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
        This is an experimental feature that is not yet fully supported.
        Desired water model. Can be either "SUPCRT92", "IAPWS95", or "DEW".
        These models are described here: http://chnosz.net/manual/water.html
        
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
    eq36da : str
        Path to directory where data1 files are stored.
        
    eq36co : str
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
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO'),
                 db="WORM",
                 elements=None,
                 solid_solutions=None,
                 logK=None,
                 logK_S=None,
                 logK_extrapolate="none",
                 download_csv_files=False,
                 exclude_organics=False,
                 exclude_category=None,
                 suppress_redox=[],
                 input_template="none",
                 water_model="SUPCRT92",
                 exceed_Ttr=True,
                 verbose=1,
                 load_thermo=True,
                 hide_traceback=True):
        
        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        self.water_model = water_model

        half_rxn_data = pkg_resources.resource_stream(__name__, "half_cell_reactions.csv")
        self.half_cell_reactions = pd.read_csv(half_rxn_data) #define the input file (dataframe of redox pairs)
        
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

        if load_thermo:
            
            # attributes to add to AqEquil class
            self.db = db
            self.elements = elements
            self.solid_solutions = solid_solutions

            if exclude_category is None:
                exclude_category = dict()
            
            if exclude_organics:
                if not isinstance(exclude_category.get("category_1"), list):
                    exclude_category["category_1"] = ["organic_aq", "organic_cr"]
                else:
                    print("category_1 exists in 'exclude_category' dict.")
                    if "organic_aq" not in exclude_category["category_1"]:
                        exclude_category["category_1"] = exclude_category["category_1"].append("organic_aq")
                    if "organic_cr" not in exclude_category["category_1"]:
                        exclude_category["category_1"] = exclude_category["category_1"].append("organic_cr")
            
            self.exclude_category = exclude_category
            self.logK = logK
            self.logK_S = logK_S
            self.logK_extrapolate = logK_extrapolate
            self.download_csv_files = download_csv_files
            self.suppress_redox = suppress_redox
            self.exceed_Ttr = exceed_Ttr
            self.input_template = input_template
            
            self.thermo = AqEquil.Thermodata(AqEquil_instance=self) # outer instance passed to inner instance
        
            self.data1 = self.thermo.data1

            
    def _capture_r_output(self):
        """
        Capture and create a list of R console messages
        """
        
        # Record output #
        self.stdout = []
        self.stderr = []
        
        # If DEBUGGING_R==False, uses python to print R lines after executing an R block 
        # If DEBUGGING_R==True, will ugly print from R directly. Allows printing from R to troubleshoot errors.
        if not DEBUGGING_R:
        
            # Dummy functions #
            def add_to_stdout(line): self.stdout.append(line)
            def add_to_stderr(line): self.stderr.append(line)

            # Keep the old functions #
            self.stdout_orig = rpy2.rinterface_lib.callbacks.consolewrite_print
            self.stderr_orig = rpy2.rinterface_lib.callbacks.consolewrite_warnerror

            # Set the call backs #
            rpy2.rinterface_lib.callbacks.consolewrite_print     = add_to_stdout
            rpy2.rinterface_lib.callbacks.consolewrite_warnerror = add_to_stderr

            
    def _print_captured_r_output(self):
        printable_lines = [line for line in self.stdout if line not in ['[1]', '\n']]
        printable_lines = [line for line in printable_lines if re.search("^\s*\[[0-9]+\]$", line) is None]
        printable_lines = [re.sub(r' \\n\"', "", line) for line in printable_lines]
        [print(line[2:-1]) for line in printable_lines]

        
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

    
    def _check_sample_input_file(self, input_filename, exclude, db,
                                       dynamic_db, charge_balance_on,
                                       suppress_missing,
                                       redox_suppression):
        """
        Check for problems in sample input file.
        """
        
        # does the input file exist? Is it a CSV?
        if self.__file_exists(input_filename):
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
            data_path = self.eq36da + "/data0." + db
        
        if self.thermo.thermo_db_type == "data0" and self.thermo.thermo_db_source == "URL":
            data_path = "data0." + self.thermo.data0_lettercode
        
        if not (os.path.exists(data_path) or os.path.isfile(data_path)) and self.thermo.thermo_db_source=="file":
            warn_no_data0 = ("Warning: Could not locate {}.".format(data_path) + " "
                "Unable to determine if column headers included in "
                "{} ".format(input_filename) + "match entries for species "
                "in the requested thermodynamic database '{}'.".format(db))
            if self.verbose > 0:
                print(warn_no_data0)
            
        if self.thermo.thermo_db_type == "data0":
            data0_lines = self.thermo.thermo_db.split("\n")
            recording_species = False
            for i,s in enumerate(data0_lines):
                if recording_species and "+---" in s:
                    end_index = i-1
                    recording_species=False
                    break
                if '*  species name' in s:
                    start_index = i+1
                    recording_species=True
            db_species = [i.split()[0] for i in data0_lines[start_index:end_index]]
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
        elif charge_balance_on != "none" and charge_balance_on not in list(set(df_in_headercheck.columns)):
            err_charge_balance_invalid_sp = ("The species chosen for charge balance"
                " '{}'".format(charge_balance_on)+""
                " was not found among the headers of the sample input file.")
            err_list.append(err_charge_balance_invalid_sp)
            
        if self.thermo.thermo_db_type in ["data0", "CSV"]:
            
            for species in list(dict.fromkeys(df_in_headercheck.columns)):
                if species not in db_species and species not in ['Temperature', 'logfO2', 'pH', 'Pressure', 'Eh', 'pe']+FIXED_SPECIES:
                    err_species_not_in_db = ("The species '{}'".format(species) + " "
                        "was not found in {}".format(data_path) + ". "
                        "If the column contains data that should not be "
                        "included in the speciation calculation, add the "
                        "column name to the 'exclude' argument. Try "
                        "help(AqEquil.AqEquil.speciate) "
                        "for more information about 'exclude'.")
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
                            "Log mean act", "pX", "pH", "pHCl", "pmH", "pmX",
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
        
        
        return sample_temps, sample_press
        

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
        
        args = ["cd", os.getcwd(), ";", self.eq36co+'/eqpt', "'"+os.getcwd()+"/data0."+db+"'"]
        args = " ".join(args)

        try:
            self.__run_script_and_wait(args) # run EQPT
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
            
        args = ["cd", "'" + cwdd+path_3i+"'", ";", # change directory to where 3i files are stored
                self.eq36co + '/eq3nr', # path to EQ3NR executable
                "'" + data1_path + "/data1." + db+"'", # path to data1 file
                "'"+cwdd + path_3i +"/"+ filename_3i+"'"] # path to 3i file
        
        args = " ".join(args)
        
        self.__run_script_and_wait(args) # run EQ3
        
        filename_3o = filename_3i[:-1] + 'o'
        filename_3p = filename_3i[:-1] + 'p'
        
        # The new eq36 build truncates names, e.g., MLS.Source.3i creates MLS.3o
        # Correct for this here:
        files_3o = [file for file in os.listdir(cwdd + path_3i) if ".3o" in file]
        files_3p = [file for file in os.listdir(cwdd + path_3i) if ".3p" in file]
        
        if len(files_3o) == 0:
            if self.verbose > 0:
                print('Error: EQ3 failed to produce output for ' + filename_3i)
        elif len(files_3o) == 1:
            file_3o = files_3o[0]
            try:
                # move output
                shutil.move(cwdd + path_3i+"/"+file_3o, cwdd + path_3o+"/"+filename_3o)
            except:
                self.err_handler.raise_exception("Error: could not move", path_3i+"/"+file_3o, "to", path_3o+"/"+filename_3o)
        else:
            # multiple 3o output files are present in the directory
            # this might happen when using runeq3() by itself in a directory with 3o files
            pass
            
        if len(files_3p) == 0:
            if self.verbose > 0:
                print('Error: EQ3 failed to produce output for ' + filename_3i)
        elif len(files_3p) == 1:
            file_3p = files_3p[0]
            try:
                # move output
                shutil.move(cwdd + path_3i+"/"+file_3p, cwdd + path_3p+"/"+filename_3p)
            except:
                self.err_handler.raise_exception("Error: could not move", path_3i+"/"+file_3p, "to", path_3p+"/"+filename_3p)
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

        args = ["cd", "'" + cwdd+path_6i+"'", ";", # change directory to 6i folder
                self.eq36co+'/eq6', # path of EQ6 executable
                "'" + data1_path + "/data1." + db+"'", # path to data1 file
                "'"+cwdd+path_6i + filename_6i+"'"] # path of 6i file
        
        args = " ".join(args)
        
        self.__run_script_and_wait(args) # run EQ6
        
                
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
        """
        
        # DEVNULL and STDOUT needed to suppress all warnings
        subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True).wait()

            
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
        
        
        f1_x = np.linspace(xvals[0], xvals[3], num=res)
        f2_x = np.linspace(xvals[3], xvals[7], num=res)
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
        
        # turns off poor polyfit warning
        # TODO: restore polyfit warning setting afterward
        warnings.simplefilter('ignore', np.RankWarning)
        
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
        
        
    def speciate(self,
                 input_filename,
                 db=None,
                 db_solid_solution=None,
                 db_logK=None,
                 logK_extrapolate=None,
                 activity_model="b-dot",
                 redox_flag="logfO2",
                 redox_aux="Fe+3",
                 default_logfO2=-6,
                 exclude=[],
                 suppress=[],
                 alter_options=[],
                 charge_balance_on="none",
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
                 get_affinity_energy=False, # deprecated
                 negative_energy_supplies=False, # deprecated
                 mineral_reactant_energy_supplies=False,
                 rxn_filename=None, # deprecated
                 not_limiting=["H+", "OH-", "H2O"], # deprecated
                 get_charge_balance=True,
                 custom_db=False, # deprecated
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
        
        db : str, default "wrm"
            Determines which thermodynamic database is used in the speciation
            calculation. There are several options available:
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
            "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
        
        db_solid_solution : str, optional
            Used only if `db` points to a thermodynamic data CSV file (or the
            URL of a CSV hosted online). Determines which thermodynamic database
            is used for idealized solid solutions in the speciation calculation.
            There are two options:
            - The name of a CSV file containing solid solution parameters
            located in the current working directory, e.g.,
            "wrm_solid_solutions.csv"
            - The URL of a CSV file containing solid solution parameters, e.g.,
            "https://raw.githubusercontent.com/worm-portal/WORM-db/master/solid_solutions.csv"
        
        db_logK : str, optional
            The name of the CSV file containing species with dissociation
            constants but no other properties or parameters. Used only if `db`
            points to a thermodynamic data CSV file (or the URL of a CSV hosted
            online).
        
        activity_model : str, default "b-dot"
            Activity model to use for speciation. Can be either "b-dot",
            or "davies". NOTE: the "pitzer" model is not yet implemented.
        
        redox_flag : str, default "O2(g)"
            Determines which column in the sample input file sets the overall
            redox state of the samples. Options for redox_flag include 'O2(g)',
            'pe', 'Eh', 'logfO2', and 'redox aux'. The code will search your
            sample spreadsheet file (see `filename`) for a column corresponding
            to the option you chose:
            
            * 'O2(g)' with a valid subheader for a gas
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
            
        charge_balance_on : str, default "none"
            If "none", will not balance electrical charge between cations and
            anions in the speciation calculation. If a name of a species is
            supplied instead, the activity of that species will be allowed to
            change until charge balance is obtained. For example,
            charge_balance_on = "H+" will calculate what pH a sample must have
            to have zero net charge.
        
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
        
        get_affinity_energy : bool, default False
            Deprecated; affinities and energy supplies are now calculated after
            speciation.
        
        negative_energy_supplies : bool, default False
            Deprecated.

        mineral_reactant_energy_supplies : bool, default False
            Report energy supplies for reactions with mineral reactants? This
            option is False by default because mineral reactants are considered
            to be unlimited. As a result, energy supplies from reactions with
            reactant minerals tend to be artificially high, especially in
            systems where the reactant minerals are unstable.
        
        rxn_filename : str, optional
            Name of .txt file containing reactions used to calculate affinities
            and energy supplies. Ignored if `get_affinity_energy` is False.
        
        not_limiting : list, default ["H+", "OH-", "H2O"]
            Deprecated.
        
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

        # deprecated!
        if get_affinity_energy:
            self.err_handler.raise_exception("Deprecation error: affinity and "
                    "energy supply calculations are now handled after speciation "
                    "using the speciation.apply_redox_reactions(...) or "
                    "speciation.calculate_energy(...) functions.")

        self.batch_T = []
        self.batch_P = []
        
        self.verbose = verbose
        
        if db != None:
            # load new thermodynamic database
            self.thermo._set_active_db(db)
        else:
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
        
        # check input sample file for errors
        if activity_model != 'pitzer': # TODO: allow check_sample_input_file() to handle pitzer
            sample_temps, sample_press = self._check_sample_input_file(
                                          input_filename, exclude, db,
                                          dynamic_db, charge_balance_on, suppress_missing,
                                          redox_suppression)
        
        if aq_dist_type not in ["molality", "log_molality", "log_gamma", "log_activity"]:
            self.err_handler.raise_exception("Unrecognized aq_dist_type. Valid "
                "options are 'molality', 'log_molality', 'log_gamma', 'log_activity'")
        if mineral_sat_type not in ["logQoverK", "affinity"]:
            self.err_handler.raise_exception("Unrecognized mineral_sat_type. Valid "
                "options are 'logQoverK' or 'affinity'")
        if redox_type not in ["Eh", "pe", "logfO2", "Ah"]:
            self.err_handler.raise_exception("Unrecognized redox_type. Valid "
                "options are 'Eh', 'pe', 'logfO2', or 'Ah'")
        
        if redox_flag == "O2(g)" or redox_flag == -3:
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
            self.err_handler.raise_exception("Unrecognized redox flag. Valid options are 'O2(g)'"
                                             ", 'pe', 'Eh', 'logfO2', 'redox aux'")
            
        # handle batch_3o naming
        if batch_3o_filename != None:
            if ".rds" in batch_3o_filename[-4:]:
                batch_3o_filename = batch_3o_filename
            else:
                batch_3o_filename = "batch_3o_{}.rds".format(data0_lettercode)
        else:
            batch_3o_filename = ro.r("NULL")
        
        # reset logK_models whenever speciate() is called
        # (prevents errors when speciations are run back-to-back)
        self.logK_models = {}
        
        # dynamic data0 creation per sample
        if dynamic_db:
            db_args["fill_data0"] = False
            db_args["dynamic_db"] = True
            db_args["verbose"] = self.verbose
            db_args["dynamic_db_sample_temps"] = sample_temps
            db_args["dynamic_db_sample_press"] = sample_press
            
            if db_logK != None:
                self.thermo._load_logK(db_logK, source="file")
            
            if logK_extrapolate != None:
                db_args["logK_extrapolate"] = logK_extrapolate
            elif self.thermo.logK_active:
                db_args["logK_extrapolate"] = self.thermo.logK_extrapolate
                logK_extrapolate = self.thermo.logK_extrapolate
            else:
                logK_extrapolate = "none"

            if db_solid_solution != None:
                if not (db_solid_solution[0:8].lower() == "https://" or db_solid_solution[0:7].lower() == "http://" or db_solid_solution[0:4].lower() == "www."):
                    if os.path.exists(db_solid_solution) and os.path.isfile(db_solid_solution):
                        db_args["filename_ss"] = db_solid_solution
                    else:
                        self.err_handler.raise_exception("Error: could not locate " + str(db_solid_solution))
                else:
                    db_solid_solution_csv_name = db_solid_solution.split("/")[-1].lower()

                    try:
                        # Download from URL and decode as UTF-8 text.
                        with urlopen(db_solid_solution) as webpage:
                            content = webpage.read().decode()
                    except:
                        self.err_handler.raise_exception("The webpage "+str(db_solid_solution)+" cannot"
                                " be reached at this time.")
                        
                    # Save to CSV file.
                    with open(db_solid_solution_csv_name, 'w') as output:
                        output.write(content)
                        
                    db_args["filename_ss"] = db_solid_solution_csv_name
                    
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
            
            data1_path = os.getcwd()+"/eqpt_files" # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA
            
            data0_path = "data0." + data0_lettercode
            
        elif dynamic_db:
            self.__mk_check_del_directory('eqpt_files')
            
        else:
            data0_path = self.eq36da + "/data0." + data0_lettercode
        
        # gather information from data0 file and perform checks
        if not dynamic_db:
            if os.path.exists(data0_path) and os.path.isfile(data0_path):
                with open(data0_path) as data0:
                    data0_lines = data0.readlines()
                    start_index = [i+1 for i, s in enumerate(data0_lines) if s == 'temperatures\n']
                    if activity_model == 'davies' or activity_model == 'b-dot':
                        end_index = [i for i, s in enumerate(data0_lines) if s == 'debye huckel a (adh)\n']
                    elif activity_model == 'pitzer':
                        end_index = [i for i, s in enumerate(data0_lines) if s == 'debye huckel aphi\n']
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
                              "15.5365", "39.7365", "85.8378", "165.2113"]
                
            grid_press_numeric = [float(n) for n in grid_press]
            if min(grid_press_numeric) == 1:
                P1=True
            else:
                P1=False
                
            self._capture_r_output()
        
            r_check_TP_grid = pkg_resources.resource_string(__name__, 'check_TP_grid.r').decode("utf-8")
        
            ro.r(r_check_TP_grid)
        
            list_tp = ro.r.check_TP_grid(grid_temps=_convert_to_RVector(grid_temps),
                                         grid_press=_convert_to_RVector(grid_press),
                                         P1=P1,
                                         water_model=water_model,
                                         check_for_errors=False,
                                         verbose=self.verbose)
        
            self._print_captured_r_output()
            
            grid_temps = list(list_tp.rx2("grid_temps"))
            grid_press = list(list_tp.rx2("grid_press"))
            poly_coeffs_1 = list_tp.rx2("poly_coeffs_1")
            poly_coeffs_2 = list_tp.rx2("poly_coeffs_2")
            
            
        else:
            grid_temps = ro.r("NULL")
            grid_press = ro.r("NULL")
            poly_coeffs_1 = ro.r("NULL")
            poly_coeffs_2 = ro.r("NULL")

        
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
                alter_options_dict[key] = _convert_to_RVector(list(ao[1:]))
        alter_options = ro.ListVector(alter_options_dict)
        
        input_dir = "rxn_3i"
        output_dir = "rxn_3o"
        pickup_dir = "rxn_3p"
            
        # preprocess for EQ3 using R scripts
        self._capture_r_output()
        
        r_prescript = pkg_resources.resource_string(
            __name__, 'preprocess_for_EQ3.r').decode("utf-8")
        ro.r(r_prescript)
        
        input_processed_list = ro.r.preprocess(input_filename=input_filename,
                                               exclude=_convert_to_RVector(exclude),
                                               grid_temps=_convert_to_RVector(grid_temps),
                                               grid_press=_convert_to_RVector(grid_press),
                                               strict_minimum_pressure=strict_minimum_pressure,
                                               dynamic_db=dynamic_db,
                                               poly_coeffs_1=poly_coeffs_1,
                                               poly_coeffs_2=poly_coeffs_2,
                                               water_model=water_model,
                                               verbose=self.verbose)
        
        
        self._print_captured_r_output()
        
        self.df_input_processed = ro.conversion.rpy2py(input_processed_list.rx2("df"))
        
        
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
        
        self.batch_T = list(input_processed_list.rx2("temp_degC"))
        self.batch_P = list(input_processed_list.rx2("pressure_bar"))
        
        # create and run a 3i file for each sample
        for sample_row_index in range(0, self.df_input_processed.shape[0]):
            
            temp_degC = list(input_processed_list.rx2("temp_degC"))[sample_row_index]
            pressure_bar = list(input_processed_list.rx2("pressure_bar"))[sample_row_index]

            df = self.df_input_processed.iloc[[sample_row_index]] # double brackets to keep as df row instead of series
            
            samplename = str(df.index[0])
            
            # handle dynamic data0 creation
            if dynamic_db:
                
                self.__fill_data0(thermo_df=ro.conversion.rpy2py(thermo_df),
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

                data1_path = os.getcwd()+"/eqpt_files" # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA

                data0_path = "data0." + data0_lettercode
                
            else:
                pressure_bar = list(input_processed_list.rx2("pressure_bar"))[sample_row_index]
                data1_path = self.thermo.eq36da
            
            # allowed aq block species are left after any category exclusion in db_args
            allowed_aq_block_species = ["all"]
            if dynamic_db:
                allowed_aq_block_species = list(thermo_df["name"]) + FIXED_SPECIES
            
            # write 3i files
            self._capture_r_output()

            warned_about_redox_column = ro.r.write_3i_file(df=ro.conversion.py2rpy(df),
                               temp_degC=temp_degC,
                               pressure_bar=pressure_bar,
                               minimum_pressure=input_processed_list.rx2("minimum_pressure"),
                               strict_minimum_pressure=strict_minimum_pressure,
                               pressure_override=dynamic_db,
                               suppress_missing=suppress_missing,
                               exclude=input_processed_list.rx2("exclude"),
                               allowed_aq_block_species=_convert_to_RVector(allowed_aq_block_species),
                               charge_balance_on=charge_balance_on,
                               suppress=_convert_to_RVector(suppress),
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

            self._print_captured_r_output()
        
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
                capture = False
                for line in lines:
                    if "Start of the bottom half of the input file" in line:
                        capture = True
                    if capture:
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
        
        df_input_processed_names = _convert_to_RVector(list(self.df_input_processed.columns))
        
        # mine output
        self._capture_r_output()
        
        r_3o_mine = pkg_resources.resource_string(
            __name__, '3o_mine.r').decode("utf-8")
        ro.r(r_3o_mine)
        
        batch_3o = ro.r.main_3o_mine(
            files_3o=_convert_to_RVector(files_3o),
            input_filename=input_filename,
            input_pressures=_convert_to_RVector(list(input_processed_list.rx2("pressure_bar"))),
            get_aq_dist=get_aq_dist,
            aq_dist_type=aq_dist_type,
            get_mass_contribution=get_mass_contribution,
            mass_contribution_other=mass_contribution_other,
            get_mineral_sat=get_mineral_sat,
            mineral_sat_type=mineral_sat_type,
            get_redox=get_redox,
            redox_type=redox_type,
            get_charge_balance=get_charge_balance,
            get_ion_activity_ratios=get_ion_activity_ratios,
            get_fugacity=get_fugacity,
            get_basis_totals=get_basis_totals,
            get_solid_solutions=get_solid_solutions,
            batch_3o_filename=batch_3o_filename,
            df_input_processed=ro.conversion.py2rpy(self.df_input_processed),
            # New rpy2 py2rpy2 conversion might not need the workaround below.
            # The old note regarding deprecated pandas2ri is shown below...
            # OLD NOTE:
            # Needed for keeping symbols in column names after porting
            #   df_input_processed in the line above. Some kind of check.names
            #   option for pandas2ri.py2ri would be nice. Workaround:
            df_input_processed_names=df_input_processed_names,
            verbose=self.verbose,
        )

        self._print_captured_r_output()
        
        if len(batch_3o) == 0:
            self.err_handler.raise_exception("Could not compile a speciation report. This is "
                            "likely because errors occurred during "
                            "the speciation calculation.")
            return
        
        if get_mass_contribution:
            mass_contribution = ro.conversion.rpy2py(batch_3o.rx2('mass_contribution'))
        df_report = ro.conversion.rpy2py(batch_3o.rx2('report'))
        
        #df_input = ro.conversion.rpy2py(batch_3o.rx2('input'))
        report_divs = batch_3o.rx2('report_divs')

        input_cols = list(report_divs.rx2('input'))
        df_input = df_report[input_cols].copy()
        
        # add a pressure column to df_input
        df_input["Pressure_bar"] = pd.Series(dtype='float')
        sample_data = batch_3o.rx2('sample_data')
        for sample in sample_data:
            df_input.loc[str(sample.rx2('name')[0]), "Pressure_bar"] = float(sample.rx2('pressure')[0])
        report_divs[0] = _convert_to_RVector(input_cols + ["Pressure_bar"])
            
        # handle headers and subheaders of input section
        headers = [col.split("_")[0] for col in list(df_input.columns)]
        headers = ["pH" if header == "H+" else header for header in headers]
        headers = [header+"_(input)" if header not in ["Temperature", "logfO2", "Pressure"]+exclude else header for header in headers]
        report_divs[0] = _convert_to_RVector(headers) # modify headers in the 'input' section, report_divs[0]
        subheaders = [subheader[1] if len(subheader) > 1 else "" for subheader in [
            col.split("_") for col in list(df_input.columns)]]
        multicolumns = pd.MultiIndex.from_arrays(
            [headers, subheaders], names=['Sample', ''])
        
        df_input.columns = multicolumns

        df_join = df_input

        if get_aq_dist:
            aq_distribution_cols = list(report_divs.rx2('aq_distribution'))
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
            aq_dist_indx = list(report_divs.names).index("aq_distribution")
            report_divs[aq_dist_indx] = _convert_to_RVector(list(headers))
            
            df_join = df_join.join(df_aq_distribution)

        if get_mineral_sat:
            mineral_sat_cols = list(report_divs.rx2('mineral_sat'))
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
            redox_cols = list(report_divs.rx2('redox'))
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
            charge_balance_cols = list(report_divs.rx2('charge_balance'))
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
            if type(report_divs.rx2('ion_activity_ratios')) != rpy2.rinterface_lib.sexp.NULLType:
                ion_activity_ratio_cols = list(report_divs.rx2('ion_activity_ratios'))

                df_ion_activity_ratios = df_report[ion_activity_ratio_cols]
                df_ion_activity_ratios = df_ion_activity_ratios.apply(pd.to_numeric, errors='coerce')

                # handle headers of df_ion_activity_ratios section
                headers = df_ion_activity_ratios.columns
                subheaders = ["Log ion-H+ activity ratio"]*len(headers)
                multicolumns = pd.MultiIndex.from_arrays(
                    [headers, subheaders], names=['Sample', ''])
                df_ion_activity_ratios.columns = multicolumns
                df_join = df_join.join(df_ion_activity_ratios)
            
        if get_fugacity:
            fugacity_cols = list(report_divs.rx2('fugacity'))
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
            sc_cols = list(report_divs.rx2('basis_totals'))
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

        sample_data = batch_3o.rx2('sample_data')

        # assemble sample data
        for i, sample in enumerate(sample_data):
            dict_sample_data = {
                "filename": str(sample.rx2('filename')[0]),
                "name": str(sample.rx2('name')[0]),
                "temperature": float(sample.rx2('temperature')[0]),
                "pressure": float(sample.rx2('pressure')[0]),
                "logact_H2O": float(sample.rx2('logact_H2O')[0]),
                "H2O_density": float(sample.rx2('H2O_density')[0]),
                "H2O_molality": float(sample.rx2('H2O_molality')[0]),
                "H2O_log_molality": float(sample.rx2('H2O_log_molality')[0]),
                }

            if get_aq_dist:
                sample_aq_dist = ro.conversion.rpy2py(sample.rx2('aq_distribution'))
                sample_aq_dist = sample_aq_dist.apply(pd.to_numeric, errors='coerce')
                
                sample_pH = -sample_aq_dist.loc["H+", "log_activity"]
                out_dict["report"].loc[str(sample.rx2('name')[0]), "pH"] = sample_pH
                
                dict_sample_data.update({"aq_distribution": sample_aq_dist})

            if get_mass_contribution:
                sample_mass_contribution = mass_contribution[mass_contribution["sample"] == sample.rx2('name')[0]]
                dict_sample_data.update(
                    {"mass_contribution": sample_mass_contribution})

            if get_mineral_sat:
                dict_sample_data.update(
                    {"mineral_sat": ro.conversion.rpy2py(sample.rx2('mineral_sat')).apply(pd.to_numeric, errors='coerce')})
                # replace sample mineral_sat entry with None if there is no mineral saturation data.
                if(len(dict_sample_data['mineral_sat'].index) == 1 and dict_sample_data['mineral_sat'].index[0] == 'None'):
                    dict_sample_data['mineral_sat'] = None

            if get_redox:
                dict_sample_data.update(
                    {"redox": ro.conversion.rpy2py(sample.rx2('redox')).apply(pd.to_numeric, errors='coerce')})

            if get_charge_balance:
                dict_sample_data.update({"charge_balance": df_charge_balance.loc[sample.rx2('name')[0], :]})
            
            if get_ion_activity_ratios:
                
                try:
                    dict_sample_data.update(
                        {"ion_activity_ratios": ro.conversion.rpy2py(sample.rx2('ion_activity_ratios'))})
                except:
                    dict_sample_data['ion_activity_ratios'] = None
            
            if get_fugacity:
                dict_sample_data.update(
                    {"fugacity": ro.conversion.rpy2py(sample.rx2('fugacity')).apply(pd.to_numeric, errors='coerce')})
                # replace sample fugacity entry with None if there is no fugacity data.
                if(len(dict_sample_data['fugacity'].index) == 1 and dict_sample_data['fugacity'].index[0] == 'None'):
                    dict_sample_data['fugacity'] = None
                else:
                    dict_sample_data["fugacity"]["fugacity"] = 10**dict_sample_data["fugacity"]["log_fugacity"]
                    
            if get_basis_totals:
                sc_dist = ro.conversion.rpy2py(sample.rx2('basis_totals'))
                sc_dist = sc_dist.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"basis_totals": sc_dist})

            if get_solid_solutions:
                sample_solid_solutions = batch_3o.rx2["sample_data"].rx2[str(sample.rx2('name')[0])].rx2["solid_solutions"]

                if not type(sample_solid_solutions.names) == rpy2.rinterface_lib.sexp.NULLType:

                    ss_df_list = []
                    for ss in list(sample_solid_solutions.names):
                        df_ss_ideal = ro.conversion.rpy2py(sample_solid_solutions.rx2[str(ss)].rx2["ideal solution"])
                        df_ss_mineral = ro.conversion.rpy2py(sample_solid_solutions.rx2[str(ss)].rx2["mineral"])
                        df_merged = pd.merge(df_ss_mineral, df_ss_ideal, left_on='mineral', right_on='component', how='left')
                        df_merged.insert(0, 'solid solution', ss)
                        del df_merged['component']
                        ss_df_list.append(df_merged)
                
                    dict_sample_data.update(
                        {"solid_solutions": pd.concat(ss_df_list)})

            out_dict["sample_data"].update(
                {sample_data.names[i]: dict_sample_data})

        out_dict.update({"batch_3o": batch_3o})
        
        out_dict.update({"water_model":water_model, "grid_temps":grid_temps, "grid_press":grid_press})
        
        speciation = Speciation(out_dict, hide_traceback=self.hide_traceback)

        speciation.half_cell_reactions = self.half_cell_reactions
        
        if report_filename != None:
            if ".csv" in report_filename[-4:]:
                out_dict["report"].to_csv(report_filename)
            else:
                out_dict["report"].to_csv(report_filename+".csv")

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
        
        
        self._capture_r_output()
        
        r_check_TP_grid = pkg_resources.resource_string(
            __name__, 'check_TP_grid.r').decode("utf-8")
        
        ro.r(r_check_TP_grid)
        
        list_tp = ro.r.check_TP_grid(grid_temps=_convert_to_RVector(grid_temps),
                                     grid_press=_convert_to_RVector(grid_press),
                                     P1=P1,
                                     water_model=water_model,
                                     check_for_errors=True,
                                     verbose=self.verbose)
        
        self._print_captured_r_output()
        
        grid_temps = list(list_tp.rx2("grid_temps"))
        grid_press = list(list_tp.rx2("grid_press"))
        
        if plot_poly_fit and len(grid_temps) == 8:
            self.__plot_TP_grid_polyfit(xvals=grid_temps,
                                        yvals=grid_press,
                                        poly_coeffs_1=list(list_tp.rx2("poly_coeffs_1")),
                                        poly_coeffs_2=list(list_tp.rx2("poly_coeffs_2")),
                                        res=500)

        self._print_captured_r_output()
        
        # calculate logK at each T and P for every species
        out_dfs = []
        for i,Tc in enumerate(grid_temps):
            out_dfs.append(calc_logK(thermo_df, Tc=Tc, P=grid_press[i], TP_i=i, water_model=water_model))
        
        dissrxn_logK_dict = {'name': out_dfs[0]["name"],
                             'logK_0': out_dfs[0]["dissrxn_logK_0"]}
        
        if len(grid_temps) == 8:
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
            logK_grid = list(dissrxn_logK.iloc[idx, 1:9])
            
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
                if (i+1) == 1 or (i+1) % 5 == 0: # first entry of a line
                    max_length = 11
                    end_char = ""
                elif (i+1) % 4 == 0 and (i+1) != len(logK_grid): # last entry of a line
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

        # handle data0 header section
        self._capture_r_output()
        
        r_fill_data0_header = pkg_resources.resource_string(
            __name__, 'fill_data0_header.r').decode("utf-8")
        
        ro.r(r_fill_data0_header)
        
        data0_file_lines = ro.r.fill_data0_head(data0_template=data0_file_lines,
                                       db=db,
                                       grid_temps=_convert_to_RVector(grid_temps),
                                       grid_press=_convert_to_RVector(grid_press),
                                       water_model=water_model,
                                       activity_model=activity_model)
        
        self._print_captured_r_output()
        
        with open("data0."+db, 'w') as f:
            for item in data0_file_lines:
                f.write("%s" % item)

                
    def plot_logK_fit(self, name, plot_out=False, res=200, internal=True, logK_extrapolate=None, T_vals=[]):
        """
        Plot the fit of logK values used in the speciation.

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
                    logK_label = "calculated LogK value(s)"
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

                if sp_press_grid == grid_press_list and sp_temps_grid == grid_temps:
                    # If pressures and temperature grid exactly matches that of the sp...
                    # need to test this!
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
                     activity_model="b-dot",
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
        
        self.batch_T = grid_temps
        self.batch_P = grid_press
        
        if not dynamic_db:
            if self.verbose >= 1:
                print("Creating data0.{}...".format(db), flush=True)
        
        if len(grid_temps) not in [1, 8]:
            self.err_handler.raise_exception("'grid_temps' must have either one or eight values.")
        if isinstance(grid_press, list):
            if len(grid_press) not in [1, 8]:
                self.err_handler.raise_exception("'grid_press' must have either one or eight values.")
        
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
            
        if water_model != "SUPCRT92":
            print("WARNING: water models other than SUPCRT92 are not yet fully supported.")
        
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
            
            free_logK_df = _clean_rpy2_pandas_conversion(self.thermo.logK_db)

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
            
            thermo_df = _clean_rpy2_pandas_conversion(thermo_df)
        
        if self.thermo.solid_solutions_active:
            solid_solution_df = ro.conversion.py2rpy(self.thermo.solid_solution_db)
        else:
            solid_solution_df = ro.r("NULL")
        
        template = pkg_resources.resource_string(
            __name__, 'data0.min').decode("utf-8")
        
        out_list = self.thermo.out_list
    
        self._capture_r_output()
    
        r_create_data0 = pkg_resources.resource_string(
            __name__, 'create_data0.r').decode("utf-8")
        
        ro.r(r_create_data0)
        
        # assemble data0 file
        data0_file_lines = ro.r.create_data0(
                          thermo_df=ro.conversion.py2rpy(thermo_df),
                          solid_solution_df=solid_solution_df,
                          db=db,
                          water_model=water_model,
                          template=template,
                          dissrxns=out_list.rx2("dissrxns"),
                          basis_pref=out_list.rx2("basis_pref"),
                          exceed_Ttr=exceed_Ttr,
                          fixed_species=_convert_to_RVector(FIXED_SPECIES),
                          verbose=self.verbose,
                          )
        
        self._print_captured_r_output()
        
        data0_file_lines = data0_file_lines[0].split("\n")
        
        if fill_data0:
            
            # begin TP-dependent processes
            self.__fill_data0(thermo_df=ro.conversion.rpy2py(thermo_df),
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
            

    def make_redox_reactions(self, *args, **kwargs):
        """
        Deprecated
        """
        self.err_handler.raise_exception("Deprecation error: make_redox_reactions "
                "now belongs to the Speciation class. Perform a speciation and then "
                "use: speciation.make_redox_reactions(...)")

    def show_redox_reactions(self, *args, **kwargs):
        """
        Deprecated
        """
        self.err_handler.raise_exception("Deprecation error: show_redox_reactions "
                "now belongs to the Speciation class. Perform a speciation and then "
                "use: speciation.show_redox_reactions(...)")


    
    class Thermodata(object):
        """
        Metaclass to store and load thermodynamic databases.
        Inherits attributes from its outer class, AqEquil.
        
        """

        def __init__(self, AqEquil_instance):

            self.AqEquil_instance = AqEquil_instance
            
            # attributes to add to AqEquil class
            self.db = self.AqEquil_instance.db
            self.elements = self.AqEquil_instance.elements
            solid_solutions = self.AqEquil_instance.solid_solutions
            self.exclude_category = self.AqEquil_instance.exclude_category
            self.water_model = self.AqEquil_instance.water_model
            self.elements = self.AqEquil_instance.elements
            logK = self.AqEquil_instance.logK
            logK_S = self.AqEquil_instance.logK_S
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

            self.verbose=verbose

            if isinstance(self.db, str):
                if self.db == "WORM":
                    if self.verbose > 0:
                        print("Loading Water-Organic-Rock-Microbe (WORM) thermodynamic databases...")
                    self.db = "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
                    self._set_active_db(db=self.db, download_csv_files=download_csv_files)
                    if self.elements == None:
                        self._load_elements("https://raw.githubusercontent.com/worm-portal/WORM-db/master/elements.csv", source="URL", download_csv_files=download_csv_files)
                    if solid_solutions == None:
                        self._load_solid_solutions("https://raw.githubusercontent.com/worm-portal/WORM-db/master/solid_solutions.csv", source="URL", download_csv_files=download_csv_files)
                    if logK == None:
                        self._load_logK("https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK.csv", source="URL", download_csv_files=download_csv_files)
                    if logK_S == None:
                        self._load_logK_S("https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK_S.csv", source="URL", download_csv_files=download_csv_files)
                else:
                    if self.verbose > 0:
                        print("Loading a user-supplied thermodynamic database...")
                    self._set_active_db(db=self.db, download_csv_files=download_csv_files)
            
            elif isinstance(self.db, pd.DataFrame):
                if self.verbose > 0:
                    print("Loading a user-supplied Pandas DataFrame thermodynamic database...")
                self._set_active_db(db=self.db, download_csv_files=download_csv_files)
                

            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                self._validate_thermodynamic_database()

            # elements must be loaded if thermo_db_type is a CSV
            if self.elements != None:
                self._load_elements(self.elements, source="file")
            if not self.element_active and self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                self._load_elements("https://raw.githubusercontent.com/worm-portal/WORM-db/master/elements.csv", source="URL", download_csv_files=download_csv_files)

            if solid_solutions != None:
                self._load_solid_solutions(solid_solutions, source="file")

            if logK != None:
                self._load_logK(logK, source="file")

            # must be loaded after the logK database
            if logK_S != None:
                self._load_logK_S(logK_S, source="file")
                
            if self.logK_active:
                self.thermo_db = pd.concat([self.thermo_db, self.logK_db], ignore_index=True)
                
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
            mineral_name_reject = list(set(self.csv_db[(self.csv_db["G"].isnull()) & (self.csv_db['state'].str.contains('cr'))]["name"]))
            
            idx = list(self.csv_db[self.csv_db["name"].isin(mineral_name_reject)].index)
            names = self.csv_db["name"].loc[idx]

            for name in names:
                self._reject_species(name=name, reason="missing Gibbs free energy for at least one polymorph")
            
            self.csv_db = self.csv_db[~self.csv_db["name"].isin(mineral_name_reject)]
    
            # TODO: other states besides minerals
                
        def _set_active_db(self, db=None, download_csv_files=False):
            """
            Set the main active thermodynamic database to a Pandas DataFrame,
            CSV file, data0 file, a data1 file on the server, a local file,
            or from a URL address.
            """

            if isinstance(db, pd.DataFrame):
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
                    if os.path.exists(self.eq36da + "/data1." + db) and os.path.isfile(self.eq36da + "/data1." + db):
                        self.thermo_db = None
                        self.thermo_db_type = "data1"
                        self.thermo_db_source = "file"
                        self.thermo_db_filename = "data1."+db
    
                        # store contents of data1 file in AqEquil object
                        with open(self.eq36da + "/data1." + db, mode='rb') as data1_file:
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
                        self.eq36da = os.getcwd()+"/eqpt_files"
    
                    elif os.path.exists("data1." + db) and os.path.isfile("data1." + db):
                        
                        if self.verbose > 0:
                            print("data1." + db + " was not found in the EQ36DA directory "
                                  "but a data1."+db+" was found in the current working "
                                  "directory. Using it...")
    
                        self.custom_data0 = True
                        self.thermo_db = None
                        self.eq36da = os.getcwd()+"/eqpt_files"
    
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
    
                    self._load_csv(db, source="file")
    
                    self.thermo_db = self.csv_db
                    self.thermo_db_filename = self.csv_db_filename
                    self.thermo_db_type = "CSV"
                    self.thermo_db_source = "file"
                    self.dynamic_db = True
                    self.custom_data0 = False
                    self.data0_lettercode = None
    
                elif db[-4:].lower() == ".csv" and (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
                    # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv
                    
                    
                    self._load_csv(db, source="URL", download_csv_files=download_csv_files)
    
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
                        "\n\t db='https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv'")

                self.db = db
            
            if self.verbose > 0:
                print(self.thermo_db_filename, "is now set as the active thermodynamic database.")

                if self.thermo_db_filename in ['data0.wrm', 'data1.wrm']:
                    print("This database is meant for rapid calculations between 0 and 350 °C at water saturation pressure.")
                elif self.thermo_db_filename == "wrm_data.csv":
                    print("This database is meant for calculations between 0 and 1000 °C and up to 5 kb pressure.")

            


        def __df_from_url(self, url, download_csv_files=False):
            """
            Get a filename and dataframe from a URL pointing to a CSV file.
            """

            filename = url.split("/")[-1].lower()

            try:
                # Download from URL and decode as UTF-8 text.
                with urlopen(url) as webpage:
                    content = webpage.read().decode()
            except:
                self.err_handler.raise_exception("The webpage "+str(url)+" cannot"
                        " be reached at this time.")

            if download_csv_files:
                if self.verbose > 0:
                    print("Downloading", filename, "from", url)
                with open(filename, 'w') as output:
                    output.write(content)

            return filename, pd.read_csv(StringIO(content), sep=",")


        def __str_from_url(self, url):
            """
            Get a filename and contents from a URL pointing to a txt file.
            """

            filename = url.split("/")[-1].lower()

            try:
                # Download from URL and decode as UTF-8 text.
                with urlopen(url) as webpage:
                    txt_content = webpage.read().decode()
            except:
                self.err_handler.raise_exception("The webpage "+str(url)+" cannot"
                        " be reached at this time.")

            if self.verbose > 0:
                print("Downloading", filename, "from", url)
            with open(filename, 'w') as output:
                output.write(txt_content)

            return filename, txt_content


        def _load_elements(self, db, source="url", download_csv_files=False):
            """
            Load an element database CSV file from a file or URL.
            """

            if source == "file":
                # e.g., "elements.csv"
                if os.path.exists(db) and os.path.isfile(db):
                    self.element_db = pd.read_csv(db)
                    self.element_db_source = "file"
                    self.element_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the CSV file '"+db+"'")

            elif source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/elements.csv"
                self.element_db_filename, self.element_db = self.__df_from_url(db, download_csv_files=download_csv_files)
                self.element_db_source = "URL"
            else:
                if self.verbose > 0:
                    print("No element database loaded.")

            if self.thermo_db_type in ["CSV", "Pandas DataFrame"]:
                if self.verbose > 0:
                    print("Element database", self.element_db_filename, "is active.")
                self.element_active = True
            else:
                if self.verbose > 0:
                    print("Element database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")


        def _load_solid_solutions(self, db, source="url", download_csv_files=False):
            """
            Load a solid solution database CSV file from a file or URL.
            """

            if source == "file":
                # e.g., "solid_solutions.csv"
                if os.path.exists(db) and os.path.isfile(db):
                    self.solid_solution_db = pd.read_csv(db)
                    self.solid_solution_db_source = "file"
                    self.solid_solution_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the CSV file '"+db+"'")

            elif source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/solid_solutions.csv"
                self.solid_solution_db_filename, self.solid_solution_db = self.__df_from_url(db, download_csv_files=download_csv_files)
                self.solid_solution_db_source = "URL"
            else:
                if self.verbose > 0:
                    print("No solid solution database loaded.")

            if self.thermo_db_type == "CSV":
                if self.verbose > 0:
                    print("Solid solution database", self.solid_solution_db_filename, "is active.")
                self.solid_solutions_active = True
            else:
                if self.verbose > 0:
                    print("Solid solution database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")


        def _load_logK(self, db, source="URL", download_csv_files=False):
            """
            Load a logK database CSV file from a file or URL.
            """

            if source == "file":
                # e.g., "logK.csv"
                if os.path.exists(db) and os.path.isfile(db):
                    self.logK_db = pd.read_csv(db)
                    self.logK_db_source = "file"
                    self.logK_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the CSV file '"+db+"'")

            elif source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK.csv"
                self.logK_db_filename, self.logK_db = self.__df_from_url(db, download_csv_files=download_csv_files)
                self.logK_db_source = "URL"

            else:
                if self.verbose > 0:
                    print("No logK database loaded.")

            if self.thermo_db_type == "CSV":
                if self.verbose > 0:
                    print("LogK database", self.logK_db_filename, "is active.")
                self.logK_active = True
            else:
                if self.verbose > 0:
                    print("LogK database is not active because the active thermodynamic database is a", self.thermo_db_type, "and not a CSV.")

            self.logK_db = self._exclude_category(df=self.logK_db, df_name=self.logK_db_filename)


        def _load_logK_S(self, db, source="URL", download_csv_files=False):
            """
            Load a logK_S database CSV file from a file or URL.
            """

            if source == "file":
                # e.g., "logK_S.csv"
                if os.path.exists(db) and os.path.isfile(db):
                    self.logK_S_db = pd.read_csv(db)
                    self.logK_S_db_source = "file"
                    self.logK_S_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the CSV file '"+db+"'")

            elif source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK_S.csv"
                self.logK_S_db_filename, self.logK_S_db = self.__df_from_url(db, download_csv_files=download_csv_files)
                self.logK_S_db_source = "URL"

            else:
                if self.verbose > 0:
                    print("No logK_S database loaded.")

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
                    
                    logK_25C = float(self.logK_S_db["logK_25"][i])
                    
                    IS_ref = float(self.logK_S_db["logK_25_IS"][i])
                    
                    T_list = self.logK_S_db["T_vals"][i].split(" ")
                    T_list = [float(T) for T in T_list]

                    Delta_S = float(self.logK_S_db["DeltaS"][i])
                    
                    if IS_ref > 0:
                        # extrapolate to ionic strength 0
                        
                        # collect azero values for metal, ligand, and complex
                        metal_name = self.logK_S_db["metal_name"][i]
                        
                        if self.thermo_db["name"].isin([metal_name]).any():
                            metal_azero = list(self.thermo_db[self.thermo_db["name"] == metal_name]["azero"])[0]
                            metal_charge = float(list(self.thermo_db[self.thermo_db["name"] == metal_name]["z.T"])[0])
#                         elif:
#                             # todo: get metal azero and charge from other databases, e.g., logK db
#                             pass
                        else:
                            # todo: throw error
                            pass
                        
                        ligand_name = self.logK_S_db["ligand_name"][i]
                        if isinstance(self.logK_S_db["ligand_azero"][i], float):
                            ligand_azero = float(self.logK_S_db["ligand_azero"][i])
                        elif self.thermo_db["name"].isin([ligand_name]).any():
                            ligand_azero = float(list(self.thermo_db[self.thermo_db["name"] == ligand_name]["azero"])[0])
                        else:
                            # todo: elif ligand_name in logK database names, get azero from there...
                            # or maybe this is not necessary if logK is merged with thermo_db at this point
                            pass
                        
                        if isinstance(self.logK_S_db["ligand_charge"][i], float):
                            ligand_charge = float(self.logK_S_db["ligand_charge"][i])
                        elif self.thermo_db["name"].isin([ligand_name]).any():
                            ligand_charge = float(list(self.thermo_db[self.thermo_db["name"] == ligand_name]["z.T"])[0])
                        else:
                            # todo: elif ligand_name in logK database names, get charge from there...
                            # or maybe this is not necessary if logK is merged with thermo_db at this point
                            pass
                    
                        dissrxn = self.logK_S_db["dissrxn"][i].split(" ")
                        n_metal = float(dissrxn[dissrxn.index(metal_name)-1])
                        n_ligand = float(dissrxn[dissrxn.index(ligand_name)-1])
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
                        if pseudoelement not in self.element_db["element"]:
                            e_df = pd.DataFrame(
                                {'element':[self.logK_S_db["ligand_element"][i]],
                                 'state':[self.logK_S_db["state"][i]],
                                 'source':[self.logK_S_db["ligand_name"][i]],
                                 'mass':[self.logK_S_db["ligand_mass"][i]],
                                 's':[self.logK_S_db["ligand_entropy"][i]],
                                 'n':[self.logK_S_db["ligand_n"][i]],
                                })

                            self.element_db = pd.concat([self.element_db, e_df], ignore_index=True)

                    if isinstance(self.logK_S_db["ligand_basis"][i], str):
                        # add a basis species representing the pseudoelement
                        basis = self.logK_S_db["ligand_basis"][i]
                        if basis not in self.thermo_db["name"]:
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

            R = 8.31446261815324/4.184 # cal/(mol K)

            # solve for G of reaction:
            # ∆_r G°= -2.303RT logK
            G_25 = -2.303*R*298.15*logK_25C # in cal/mol

            # solve for H of reaction:
            # ∆_r G°= ∆_r H°-T∆_r S°
            H = G_25 + 298.15*Delta_S # in cal/mol

            logK_list = []
            for T_C in T_list:

                T_K = T_C+273.15 # convert C to Kelvin

                # estimate G at temperature
                G_T = H - T_K*Delta_S

                # convert G to logK
                logK_T = G_T/(-2.303*R*T_K)
                logK_list.append(logK_T)

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


        def _load_csv(self, db, source="URL", download_csv_files=False):
            """
            Load a WORM-styled thermodynamic database CSV from a file or URL.
            """

            if source == "file":
                # e.g., "wrm_data.csv"
                if os.path.exists(db) and os.path.isfile(db):
                    self.csv_db = pd.read_csv(db)
                    self.csv_db_type = "CSV"
                    self.csv_db_source = "file"
                    self.csv_db_filename = db
                else:
                    self.err_handler.raise_exception("Could not locate the CSV file '"+db+"'")

            elif source == "URL":
                # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
                self.csv_db_filename, self.csv_db = self.__df_from_url(db, download_csv_files=download_csv_files)
                self.csv_db_type = "CSV"
                self.csv_db_source = "URL"

            self.csv_db = self.csv_db.astype(WORM_THERMODYNAMIC_DATABASE_COLUMN_TYPE_DICT)

            # Check that thermodynamic database input files exist and are formatted correctly.
            self._check_csv_db()
            self._remove_missing_G_species()

            self.csv_db = self._exclude_category(df=self.csv_db, df_name=self.csv_db_filename)


        def _exclude_category(self, df, df_name):
            """
            Exclude entries from a df based on values in columns.
            e.g., {"category_1":["organic_aq", "organic_cr"]}
            """
            
            exclude_keys = list(self.exclude_category.keys())
            if len(exclude_keys) > 0:
                for key in exclude_keys:
                    if self.verbose > 0:
                        print("Excluding", str(self.exclude_category[key]), "from column", str(key), "in", df_name)
                        
                    if isinstance(self.exclude_category[key], list):
                        
                        idx = list(df[df[key].isin(self.exclude_category[key])].index)
                        names = df["name"].loc[idx]
                        
                        for name in names:
                            self._reject_species(name=name, reason="excluded by user")
                        
                        df = df[~df[key].isin(self.exclude_category[key])]
                        

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

            missing_headers = []
            for header in required_headers:
                if header not in thermo_df.columns:
                    missing_headers.append(header)
            if len(missing_headers) > 0:
                msg = ("The thermodynamic database file "
                       "is missing one or more required columns: "
                       "{}".format(", ".join(missing_headers))+". "
                       "Are these headers spelled correctly in the file?")
                self.err_handler.raise_exception(msg)

            # does Cl-, O2(g), and O2 exist in the file?
            required_species = ["Cl-", "O2", "O2(g)"]
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
            
            suppress_redox = _convert_to_RVector(suppress_redox)

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

            if len(self.exclude_category) > 0:
                exclude_category_R =  {k:_convert_to_RVector(l) for k,l in zip(self.exclude_category.keys(), self.exclude_category.values())}
            else:
                exclude_category_R = {}
            exclude_category_R = ro.ListVector(exclude_category_R)

            self.AqEquil_instance._capture_r_output()

            r_redox_dissrxns = pkg_resources.resource_string(
                __name__, 'redox_and_dissrxns.r').decode("utf-8")

            ro.r(r_redox_dissrxns)
            
            thermo_df = _clean_rpy2_pandas_conversion(thermo_df)

            ro.conversion.py2rpy(thermo_df)

            self.out_list = ro.r.suppress_redox_and_generate_dissrxns(
                                   thermo_df=ro.conversion.py2rpy(thermo_df),
                                   water_model=self.water_model,
                                   exceed_Ttr=exceed_Ttr,
                                   suppress_redox=suppress_redox,
                                   exclude_category=exclude_category_R,
                                   element_df=ro.conversion.py2rpy(self.element_db),
                                   fixed_species=_convert_to_RVector(FIXED_SPECIES),
                                   verbose=self.verbose)
            
            self.AqEquil_instance._print_captured_r_output()

            thermo_df = self.out_list.rx2("thermo_df")
            thermo_df=ro.conversion.rpy2py(thermo_df)

            # Currently, species rejected by r.suppress_redox_and_generate_dissrxns()
            # are rejected because they cannot be written with valid basis species.
            # e.g., the mineral "iron" would be rejected when Fe is redox-isolated because
            # there is no aqueous basis species representing Fe with an oxidation state of 0.
            rejected_species = self.out_list.rx2("dissrxns").rx2("rejected_species")
            
            if type(rejected_species) != rpy2.rinterface_lib.sexp.NULLType:
                for i,sp in enumerate(rejected_species):
                    self._reject_species(sp, "A dissociation reaction could not be written with valid basis species.")

            thermo_df = _clean_rpy2_pandas_conversion(thermo_df)

            # convert E units and calculate missing GHS values
            self.thermo_db = OBIGT2eos(thermo_df, fixGHS=True, tocal=True)


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
            sp_total.sample_data = None
            if allow_mass_contribution:
                mass_contribution_breaks.append(0)
        else:
            sp_i = copy.deepcopy(sp)
            if allow_mass_contribution:
                mass_contribution_breaks.append(sp_total.report.shape[0])
            sp_total.report = pd.concat([sp_total.report, sp_i.report], axis=0, sort=False)
        

    sp_total.report.index = sp_total.report.index + ("_"+sp_total.report.groupby(level=0).cumcount().astype(str)).replace('_0','')
    
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


class Speciation(object):
    
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
    
    report_divs : rpy2 ListVector
        An rpy2 ListVector of column names within the different sections of the
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

        if 'report_divs' in list(self.__dict__.keys()):
            # create a dict of report categories and their child columns
            self.report_category_dict = {}
            for cat in [str(s) for s in list(self.report_divs.names)]:
                self.report_category_dict[cat] = list(self.report_divs.rx2(cat))
        

    def __getitem__(self, item):
         return getattr(self, item)

    
    def __make_speciation_group_dict(self):

        if isinstance(self.custom_grouping_filepath, str):
            with open(self.custom_grouping_filepath) as file:
                lines = [line.rstrip() for line in file]
        else:
            stream = pkg_resources.resource_stream(__name__, "speciation_groups_WORM.txt")
            with stream as s:
                content = s.read().decode()
            content = content.split("\n")
            lines = [line.rstrip() for line in content]

        lines = [l for l in lines if l != ""] # remove blank lines

        
        # Creates a dictionary with this format:
        #
        # {...
        #  "sulfides 1": ["H2S", "HS-"],
        #  "sulfides 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "sulfides 3": ["Pb(HS)3-"],
        #  "ferrous iron 1": [...],
        #  ...
        # }
        #
        # Used by calculate_energy() to calculate concentrations of limiting reactants.
        reactant_dict_scalar = {l.split(":")[0]:l.split(":")[1].strip() for l in lines}
        self.reactant_dict_scalar = {k:v.split(" ") for k,v in zip(reactant_dict_scalar.keys(), reactant_dict_scalar.values())}


        # Creates a dictionary with this format:
        # {...
        #  "H2S": ["H2S", "HS-"],
        #  "HS-": ["H2S", "HS-"],
        #  "Pb(HS)2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Ag(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Au(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  ...
        # }
        #
        # Used to switch a user-specified limiting reactant to one within the same scalar group.
        speciation_group_dict = {}
        for i,line in enumerate(lines):
            line = line.strip().split(":")
            assert len(line) == 2 # will fail if : in species names
            group_name = line[0].strip()
            group_species = line[1].strip().split(" ")
            group_species = list(collections.OrderedDict.fromkeys(group_species)) # remove duplicates
            speciation_group_dict[group_name] = group_species

        speciation_group_dict_unpacked = {}
        for key in list(speciation_group_dict.keys()):
            for sp in speciation_group_dict[key]:
                speciation_group_dict_unpacked[sp] = speciation_group_dict[key]

        self.speciation_group_dict_unpacked = speciation_group_dict_unpacked



    def __switch_limiting(self, limiting, stoich, species, lenient=False):
        if limiting != None:
            reactant_idx = [1 if i<0 else 0 for i in stoich]
            reactants = [species[i] for i,idx in enumerate(reactant_idx) if idx == 1]

            if limiting not in reactants and limiting in list(self.speciation_group_dict_unpacked.keys()):
                # if the user specifies a limiting reactant like "HCO3-" but the
                # reaction has "CO2" as a reactant, check the dict that
                # contains speciated groups and switch the limiting reactant to
                # the relevant limiting reactant to appear in reports,
                # e.g., "HCO3-" -> "CO2"
                lim_species_group_list = self.speciation_group_dict_unpacked[limiting]
                for s in lim_species_group_list:
                    if s in reactants:
                        if self.verbose > 0:
                            print("The specified limiting reactant", str(limiting),
                                  "has been switched to", str(s), "because the latter",
                                  "appears as a reactant in the reaction:")
                            if _isnotebook():
                                _ = self.format_reaction(coeffs=stoich,
                                   names=species,
                                   formatted=True,
                                   charge_sign_at_end=True,
                                   show=True)
                            else:
                                l = stoich + species
                                l[::2] = stoich
                                l[1::2] = species
                                l = [str(v) for v in l]
                                print(" ".join(l))
                        limiting = s
                        break
                            
            if limiting not in reactants and lenient:
                # if the user specifies "CO2" as the limiting reactant during a
                # batch calculation of many different reactions, some reactions
                # won't actually have "CO2" as a reactant. In this case, set
                # limiting to None so that CO2 will be limiting when applicable.
                limiting = None
                
        return limiting

    
    def apply_redox_reactions(self, y_type="E", y_units="cal", limiting=None,
                                    negative_energy_supplies=False,
                                    custom_grouping_filepath=None,
                                    append_report=True):
        
        self.custom_grouping_filepath = custom_grouping_filepath
        
        self.__make_speciation_group_dict()
        
        y_name_list = []
        val_list_list = []
        result_dict = {}
        result_lim_dict = {}

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
                limiting_input = self.__switch_limiting(limiting,
                                                  stoich=coeff_list,
                                                  species=species_list,
                                                  lenient=True)
            else:
                limiting_input = None
            
            df_rxn = self.calculate_energy(
                            species=species_list,
                            stoich=coeff_list,
                            divisor=divisor,
                            per_electron = True,
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

        if _isnotebook() and show:
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
            reactants = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
            products = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
        reaction = reactants + " = " + products

        if _isnotebook() and show:
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
            
            oxidant = self.half_cell_reactions["Oxidant"].iloc[idx]
            reductant = self.half_cell_reactions["Reductant"].iloc[idx]

            if oxidant == "O2" and reductant == "H2O":
                half_reaction_dict[idx] = {'O2': -1.0, 'e-': -4.0, 'H+': -4.0, 'H2O': 2.0}
                continue
            elif oxidant == "H2O" and reductant == "H2":
                half_reaction_dict[idx] = {'H2O': -2.0, 'e-': -2.0, 'H2': 1.0, 'OH-': 2.0}
                continue

            db_sp_names = list(self.thermo.thermo_db["name"])
            if oxidant not in db_sp_names or reductant not in db_sp_names:
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
                    print("Error: this half reaction's oxidant and reductant are switched. Electron coeff must be negative.")
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
            
            unpacked_dict = common_elem_electron_dict[list(common_elem_electron_dict.keys())[0]]

            e_coeff = unpacked_dict["electron_coeff"]
            ox_coeff = unpacked_dict["ox_coeff"]
            red_coeff = unpacked_dict["red_coeff"]

            # print("BASIS CANDIDATES")
            # print(basis_candidates+["H2O", "e-", "H+"])
            
            pyCHNOSZ.basis(basis_candidates+["H2O", "e-", "H+"],
                           messages=False)
            
            sout = pyCHNOSZ.subcrt(species=[oxidant, reductant, "e-"],
                                   coeff=[ox_coeff, red_coeff, e_coeff],
                                   property="logK", T=25,
                                   messages=False, show=False).reaction
            
            half_reaction_dict[idx] = self._create_sp_dict(sout)

        # print("HALF REACTION DICT")
        # print(half_reaction_dict)

        # find all possible combinations of idx pairs (but not reverse rxns to save time)
        good_idx_list = [idx for idx in idx_list if idx not in bad_idx_list]
        redox_pair_list = [list(p) for p in list(itertools.combinations(good_idx_list, 2))]

        self.redox_pair_list = redox_pair_list

        # create 'reaction_dict': a dictionary of reactions keyed by their idx pairs
        reaction_dict = {}
        e_dict = {}
        for pair in redox_pair_list:
            half_reaction_dict_1 = half_reaction_dict[pair[0]]
            half_reaction_dict_2 = half_reaction_dict[pair[1]]

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
            
            # divide all coefficients by their greatest common divisor
            argv = [int(c) for c in list(full_rxn_dict.values())]

            # all coefficients must be integers for math.gcd to work
            assert argv == list(full_rxn_dict.values())
            
            coeff_gcd = math.gcd(*argv)
            full_rxn_dict = {k:full_rxn_dict[k]/coeff_gcd for k in full_rxn_dict.keys()}
            
            # print("full_rxn_dict")
            # print(full_rxn_dict)

            e_transferred = half_reaction_dict_1["e-"]/coeff_gcd
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


    def __match_grouped_species(self, s):
        """
        Match whether a species is in a speciation group, and get a list of relevant groups and their scalars.
        
        e.g.,

        self.reactant_dict_scalar = {...
                                     "sulfides 1": ["H2S", "HS-"],
                                     "sulfides 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
                                     "sulfides 3": ["Pb(HS)3-"],
                                     "ferrous iron 1": [...],
                                     ...
                                     }
        If `s` = "Ag(HS)2-"
        then `scalars` = [1, 2, 3]
        and `groups` = ["sulfides 1", "sulfides 2", "sulfides 3"]
        
        """
        scalars = []
        groups = []
        for i,grp_list in enumerate(list(self.reactant_dict_scalar.values())):
            if s in grp_list:
                s_key = list(self.reactant_dict_scalar.keys())[i] #e.g., s_key can be "sulfates 1"
                s_key_split = s_key.split(" ")
                s_key_grp = " ".join(s_key_split[:-1])
                possible_scalars = [k.split(s_key_grp)[-1].strip() for k in list(self.reactant_dict_scalar.keys()) if s_key_grp in k]
                for k in possible_scalars:
                    try:
                        float(k) # test whether the scalar is a number. If so, append.
                        scalars.append(k)
                    except:
                        continue
                groups = [s_key_grp+" "+str(s) for s in scalars]
                break
        scalars = [float(s) for s in scalars]
        groups = [self.reactant_dict_scalar[g] for g in groups]
        return scalars, groups


    def calculate_energy(self, species, stoich,
                    divisor=1, custom_grouping_filepath=None,
                    per_electron=False, rxn_name="custom reaction",
                    negative_energy_supplies=False,
                    y_type="A", y_units="kcal", 
                    limiting=None, charge_sign_at_end=False,
                    simple_df_output=False,
                    append_report=True,
                    print_logK_messages=False,
                    raise_nonlimiting_exception=True):

        # check that y_type is recognized
        if y_type not in ["logK", "logQ", "G", "A", "E"]:
            self.err_handler.raise_exception("Valid options for y_type include "
                    "'logK', 'logQ', 'G' (Gibbs free energy), 'A' "
                    "(chemical affinity), and 'E' (energy suppy.")
        
        # check that a thermodynamic CSV is being used
        if not isinstance(self.thermo.csv_db, pd.DataFrame):
            self.err_handler.raise_exception("The plot_energy() function requires "
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
            if not isinstance(self.reactant_dict_scalar, dict):
                self.__make_speciation_group_dict()
        
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
                if list(self.thermo.csv_db[self.thermo.csv_db["name"]==r]["state"])[0] != "aq":
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

            limiting = self.__switch_limiting(limiting,
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

                        # check whether the species matches any of the groups in reactant_dict_scalar and retrieve scalars and groups
                        scalars, groups = self.__match_grouped_species(s)
                        
                        if len(scalars) > 0:
                            total_summed_scaled = [0]*self.aq_distribution_molal.shape[0]
                            
                            for i,scalar in enumerate(scalars):
                                col_subset = [col for col in groups[i] if col in self.aq_distribution_molal.columns]
                                scaled_df = self.aq_distribution_molal[col_subset].apply(lambda x: x*scalar)
                                summed_scaled = list(scaled_df.sum(axis=1, numeric_only=True))
                                total_summed_scaled = [ii+iii for ii,iii in zip(total_summed_scaled, summed_scaled)]
                                s_molal_dict[s] = total_summed_scaled

                        else:
                            s_molal_dict[s] = list(self.aq_distribution_molal[s])
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
                s_molal_dict[s] = [float("NaN")]*len(self.misc_params["Temp(C)"])

        if y_type in ["logK", "logQ"]:
            y_type_plain = copy.copy(y_type)
        elif y_type == "A":
            y_type_plain = "affinity"
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
                logK = pyCHNOSZ.subcrt(
                              species,
                              stoich,
                              T=T,
                              P=list(self.misc_params["Press(bars)"])[i],
                              show=False,
                              messages=print_logK_messages).out["logK"]

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
                    
                    if not isinstance(limiting, str):
                        lrc_dict = {}
                        for i_s,s in enumerate(species):
                            # identify valid limiting reactants and record concentrations
                            # 1. negative coefficient (reactant)
                            # 2. can't be OH-, H+, H2O
                            # 3. can't be cr or liq
                            if stoich[i_s] < 0 and s not in ["H2O", "H+", "OH-"] and list(self.thermo.csv_db[self.thermo.csv_db["name"]==s]["state"])[0] not in ["cr", "liq"]:
                                lrc_dict[s] = s_molal_dict[s][i]/abs(stoich[i_s])
                    
                    if not isinstance(limiting, str):
                        lr_name = min(lrc_dict, key=lrc_dict.get)
                        lr_val = lrc_dict[lr_name]
                        
                        # handle situations where there might be multiple limiting reactants
                        lr_list = []
                        for k,v in zip(lrc_dict.keys(), lrc_dict.values()):
                            if v == lr_val:
                                lr_list.append(str(k))
                    else:
                        lrc_dict = {}
                        lr_list = [limiting]

                    if len(lr_list) > 0 and sum([math.isnan(v) for v in list(lrc_dict.values())]) == 0:
                        # if there is a limiting reactant and no values of 'nan' for limiting reactant concentrations...
                        
                        lr_list_formatted = [chemlabel(lr_name, charge_sign_at_end=charge_sign_at_end) for lr_name in lr_list]
                        if len(lr_list_formatted) > 1:
                            lr_reported = ", ".join(lr_list_formatted)
                        else:
                            lr_reported = lr_list[0]
    
                        lr_name = lr_list[0] # doesn't matter which lr is used to calculate
                        lr_concentration = s_molal_dict[lr_name][i]
                        lr_name_list.append(lr_reported)
                        lr_stoich = -stoich[species.index(lr_name)]

                        E = A * (lr_concentration/lr_stoich)
    
                        y_list.append(E/divisor_i)
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
            self.report = df_out.combine_first(self.report) # update the report with affinity/energy results
            self.report = self.report.reindex(level=0, columns=col_order) # restore column order

            # update the report category dictionary so y_type_plain appears as a category with relevant columns
            if y_type_plain not in list(self.report_category_dict.keys()):
                self.report_category_dict[y_type_plain] = headers
            else:
                self.report_category_dict[y_type_plain] = self.report_category_dict[y_type_plain]+headers
                self.report_category_dict[y_type_plain]= list(collections.OrderedDict.fromkeys(self.report_category_dict[y_type_plain]))

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
        
        names_length = len(self.report_category_dict.keys())
        
        if col==None and names_length>0:
            return list(self.report_category_dict.keys())
        
        if names_length>0:
            if col in list(self.report_category_dict.keys()):
                return list(self.report_category_dict[col])
        
        if isinstance(col, str):
            col = [col]
        
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

        fig = px.bar(df, x=df.index, y="affinity",
            height=plot_height*ppi, width=plot_width*ppi,
            labels={'affinity': ylabel}, template="simple_white")
        
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

        colors = _get_colors(colormap, len(y))

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
            
            df.loc[:, yi] = y_plot


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

        
    def scatterplot(self, x="pH", y="Temperature", title=None, plot_zero=True,
                    rxns_as_labels=True, charge_sign_at_end=False,
                    plot_width=4, plot_height=3, ppi=122,
                    fill_alpha=0.7, point_size=10,
                    ylab=None, lineplot=False,
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

        title : str, optional
            Title of the plot.

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
        df = pd.melt(df, id_vars=["name", x], value_vars=y)
        df = df.rename(columns={"Sample": "y_variable", "value": "y_value"})

        if not plot_zero:
            df = df.dropna(subset=['y_value'])
            df = df[df.y_value != 0]


        if isinstance(colormap, str):
            # get colors
            colors = _get_colors(colormap, len(y), alpha=fill_alpha)
            
            # convert rgba to hex
            colors = [matplotlib.colors.rgb2hex(c) for c in colors]
    
            # map each species to its color, e.g.,
            # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
            dict_species_color = {sp:color for sp,color in zip(y, colors)}
            
            # html format color dict key names
            dict_species_color = {chemlabel(k, charge_sign_at_end=charge_sign_at_end):v for k,v in dict_species_color.items()}
        else:
            dict_species_color = {}

        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.reactions_for_plotting, pd.DataFrame):
            
            y_find = [yi.replace(" energy supply", "").replace(" affinity", "").replace(" Gibbs free energy", "") for yi in y]
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
        
        if lineplot:
            fig = px.line(df, x=x, y="y_value", color="y_variable",
                             hover_data=[x, "y_value", "y_variable", "name", "formatted_rxn"],
                             width=plot_width*ppi, height=plot_height*ppi,
                             labels={x: xlabel,  "y_value": ylabel},
                             category_orders={"species": y},
                             color_discrete_map=dict_species_color,
                             custom_data=['name', 'formatted_rxn', 'y_variable_original'],
                             template="simple_white")
        else:
            fig = px.scatter(df, x=x, y="y_value", color="y_variable",
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
            colors = _get_colors(colormap, len(unique_species))

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
        colors = _get_colors(colormap, len(unique_minerals))
        
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
            
            
    def join_6i_p(self, filepath_6i, chain_mt):
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
            
        for sample_name in raw_p_dict_bottom.keys():
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
                
            lines_3p = raw_p_dict_bottom[sample_name]
            
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
        
        if "mass_transfer" in list(sample_data[sample].keys()):
            if sample_data[sample]["mass_transfer"] != None:
                return sample_data[sample]["mass_transfer"]
            
        msg = ("Mass transfer results are not stored for sample '"+sample+"'. "
              "This might be because the reaction calculation did not "
              "finish successfully or because the thermodynamic database "
              "is a data0 or data1 file without a supporting CSV file.")
        self.err_handler.raise_exception(msg)

    
