DEBUGGING_R = True
FIXED_SPECIES = ["H2O", "H+", "O2(g)", "water", "Cl-", "e-", "OH-", "O2", "H2O(g)"]

import os
import re
import sys
import shutil
import copy
import collections
import dill
import math
from itertools import groupby

from urllib.request import urlopen
from io import StringIO

import warnings
import subprocess
import pkg_resources
import pandas as pd
import numpy as np
from chemparse import parse_formula
from IPython.core.display import display, HTML
import periodictable

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
        return "{0}<sup>{1}</sup>&frasl;<sub>{2}</sub>".format(whole_number_float, remainder_tuple[0], remainder_tuple[1])


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
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


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


def _html_chemname_format(name, charge_sign_at_end=False):
    
    """
    Function duplicated from pyCHNOSZ
    """
    
    p = re.compile(r'(?P<sp>[-+]\d*?$)')
    name = p.sub(r'<sup>\g<sp></sup>', name)
    charge = re.search(r'<.*$', name)

    name_no_charge = re.match(r'(?:(?!<|$).)*', name).group(0)
    mapping = {"0": "<sub>0</sub>", "1": "<sub>1</sub>", "2": "<sub>2</sub>", "3": "<sub>3</sub>", "4": "<sub>4</sub>", 
           "5": "<sub>5</sub>", "6": "<sub>6</sub>", "7": "<sub>7</sub>", "8": "<sub>8</sub>", "9": "<sub>9</sub>",
           ".":"<sub>.</sub>"}
    name_no_charge_formatted = "".join([mapping.get(x) or x for x in list(name_no_charge)])

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
    

class Thermodata:
    """
    Class to store and load thermodynamic databases.
    
    Parameters
    ----------
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
    
    logK_extrapolate : str, default "none"
        What method should be used to extrapolate equilibrium constants in the
        logK database (defined by parameter `logK`) as a function of
        temperature? Can be either "none", "flat", "poly", or "linear".
    
    download_csv_files : bool, default False
        Download copies of database CSV files to your current working directory?
    
    exclude_category : dict
        Exclude species from thermodynamic databases based on column values.
        For instance,
        `exclude_category={'category_1':["organic_aq", "organic_cr"]}`
        will exclude all species that have "organic_aq" or "organic_cr" in
        the column "category_1".
        Species are excluded from the main thermodynamic database CSV and the
        equilibrium constant (logK) CSV database. This parameter has no effect
        if the thermodynamic database is a data0 or data1 file.
    
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored. 
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.
        
    verbose : int, 0, 1, or 2, default 1
        Level determining how many messages are returned during a
        calculation. 2 for all messages, 1 for errors or warnings only,
        0 for silent.

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this class?
        When True, error messages handled by this class will be short and to
        the point.
    
    """
    
    def __init__(self,
                 db = "WORM",
                 solid_solutions=None,
                 logK=None,
                 logK_extrapolate="none",
                 download_csv_files=False,
                 exclude_category={},
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO'),
                 verbose=1,
                 hide_traceback=True):
    
        self.db = db
        self.exclude_category = exclude_category
    
        self.hide_traceback = hide_traceback
        self.err_handler = Error_Handler(clean=self.hide_traceback)
    
        self.eq36da = eq36da
        self.eq36co = eq36co
    
        # active thermo db attributes
        self.thermo_db = None
        self.thermo_db_type = None
        self.thermo_db_source = None
        self.thermo_db_filename = None
        self.custom_data0 = None
        self.data0_lettercode = None
        self.dynamic_db = None
        self.custom_obigt = None
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
        
        # solid solution attributes
        self.solid_solutions_active = False
        self.solid_solution_db = None
        self.solid_solution_db_source = None
        self.solid_solution_db_filename = None
        
        # logK attributes
        self.logK_active = False
        self.logK_extrapolate = logK_extrapolate
        self.logK_db = None
        self.logK_db_source = None
        self.logK_db_filename = None
        
        self.verbose=verbose
        
        if db == "WORM":
            if self.verbose > 0:
                print("Loading Water-Organic-Rock-Microbe (WORM) thermodynamic databases...")
            self.db = "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
            self._set_active_db(db=self.db, download_csv_files=download_csv_files)
            if solid_solutions == None:
                self._load_solid_solutions("https://raw.githubusercontent.com/worm-portal/WORM-db/master/solid_solutions.csv", source="URL", download_csv_files=download_csv_files)
            if logK == None:
                self._load_logK("https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data_logK.csv", source="URL", download_csv_files=download_csv_files)
            
        else:
            self._set_active_db(db=self.db, download_csv_files=download_csv_files)

        if solid_solutions != None:
            self._load_solid_solutions(solid_solutions, source="file")
        if logK != None:
            self._load_logK(logK, source="file")
        
    def _set_active_db(self, db=None, download_csv_files=False):
        """
        Set the main active thermodynamic database to a CSV file, a data0 file,
        or a data1 file on the server, a local file, or from a URL address.
        """
        
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
                self.custom_obigt = None

            elif os.path.exists("data1." + db) and os.path.isfile("data1." + db):
                if self.verbose > 0:
                    print("data1." + db + " was not found in the EQ36DA directory "
                          "but a data1."+db+" was found in the current working "
                          "directory. Using it...")

                self.custom_data0 = True
                self.thermo_db = None

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
            self.custom_obigt = None
            
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
            self.custom_obigt = None
        
        elif db[-4:].lower() == ".csv" and not (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
            # e.g., "wrm_data.csv"
            
            self._load_csv(db, source="file")
            
            self.thermo_db = self.csv_db
            self.thermo_db_filename = self.csv_db_filename
            self.thermo_db_type = "CSV"
            self.thermo_db_source = "file"
            self.dynamic_db = True
            self.custom_data0 = False
            self.custom_obigt = self.csv_db_filename
            self.data0_lettercode = None
            
        elif db[-4:].lower() == ".csv" and (db[0:8].lower() == "https://" or db[0:7].lower() == "http://" or db[0:4].lower() == "www."):
            # e.g., "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
            
            self._load_csv(db, source="URL", download_csv_files=download_csv_files)
            
            self.thermo_db = self.csv_db
            self.thermo_db_filename = self.csv_db_filename
            self.thermo_db_type = "CSV"
            self.thermo_db_source = "URL"
            self.dynamic_db = True
            self.custom_data0 = False
            self.custom_obigt = self.csv_db_filename
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
        
        if self.verbose > 0:
            print(self.thermo_db_filename, "is now set as the active thermodynamic database.")
        
            if self.thermo_db_filename in ['data0.wrm', 'data1.wrm']:
                print("This database is meant for rapid calculations between 0 and 350 °C at water saturation pressure.")
            elif self.thermo_db_filename == "wrm_data.csv":
                print("This database is meant for calculations between 0 and 1000 °C and up to 5 kb pressure.")
        
        self.db = db

        
    def __df_from_url(self, url, download_csv_files=False):
        """
        Get a filename and dataframe from a URL pointing to a CSV file.
        """
        
        filename = url.split("/")[-1].lower()
        
        # Download from URL and decode as UTF-8 text.
        with urlopen(url) as webpage:
            content = webpage.read().decode()
        
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
        
        # Download from URL and decode as UTF-8 text.
        with urlopen(url) as webpage:
            txt_content = webpage.read().decode()
        
        if self.verbose > 0:
            print("Downloading", filename, "from", url)
        with open(filename, 'w') as output:
            output.write(txt_content)
        
        return filename, txt_content
        
        
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

        self.csv_db = self.csv_db.astype({'name':'str', 'abbrv':'str', 'formula':'str',
                                      'state':'str', 'ref1':'str', 'ref2':'str',
                                      'date': 'str', 'E_units':'str',
                                      'G':'float', 'H':'float', 'S':'float',
                                      'Cp':'float', 'V':'float', 'a1.a':'float',
                                      'a2.b':'float', 'a3.c':'float', 'a4.d':'float',
                                      'c1.e':'float', 'c2.f':'float',
                                      'omega.lambda':'float', 'z.T':'float',
                                      'azero':'float', 'neutral_ion_type':'float',
                                      'dissrxn':'str', 'tag':'str',
                                      'formula_ox':'str', 'category_1':'str'})
        
        # Check that thermodynamic database input files exist and are formatted correctly.
        self._check_csv_db()
        
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
                    df = df[~df[key].isin(self.exclude_category[key])]
                elif isinstance(self.exclude_category[key], str):
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
            msg = ("The thermodynamic database file '{}'".format(filename)+" "
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
            msg = ("The thermodynamic database file '{}'".format(filename)+" "
                   "is missing required species:"
                   "{}".format(missing_species)+". Default thermodynamic values"
                   " will be used.")
            warnings.warn(msg)
        
        return
    
    
class AqEquil:

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
    
    logK_extrapolate : str, default "none"
        What method should be used to extrapolate equilibrium constants in the
        logK database (defined by parameter `logK`) as a function of
        temperature? Can be either "none", "flat", "poly", or "linear".
    
    download_csv_files : bool, default False
        Download copies of database CSV files to your current working directory?
    
    exclude_category : dict
        Exclude species from thermodynamic databases based on column values.
        For instance,
        `exclude_category={'category_1':["organic_aq", "organic_cr"]}`
        will exclude all species that have "organic_aq" or "organic_cr" in
        the column "category_1".
        Species are excluded from the main thermodynamic database CSV and the
        equilibrium constant (logK) CSV database. This parameter has no effect
        if the thermodynamic database is a data0 or data1 file.
        
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
        
    redox_pairs : list of int
        List of indices of half reactions in the `half_cell_reactions` table
        to be combined when generating full redox reactions.
            
    affinity_energy_reactions_raw : str
        A formatted TSV string of redox reactions for calculating chemical
        affinities and energy supplies during speciation.

    affinity_energy_reactions_table : pd.DataFrame
        A table of redox reactions for calculating chemical affinities and
        energy supplies during speciation.
    
    affinity_energy_formatted_reactions : pd.DataFrame
        A pandas dataframe containing balanced redox reactions written in full.
        
    """

    def __init__(self,
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO'),
                 db="WORM",
                 solid_solutions=None,
                 logK=None,
                 logK_extrapolate="none",
                 download_csv_files=False,
                 exclude_category={},
                 verbose=1,
                 load_thermo=True,
                 hide_traceback=True):

        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        
        half_rxn_data = pkg_resources.resource_stream(__name__, "half_cell_reactions.csv")
        self.half_cell_reactions = pd.read_csv(half_rxn_data) #define the input file (dataframe of redox pairs)
        self.redox_pairs = None
        self.affinity_energy_reactions_raw = None
        self.affinity_energy_reactions_table = None
        self.affinity_energy_formatted_reactions = None
        
        self.verbose = verbose
        self.hide_traceback = hide_traceback
        self.err_handler = Error_Handler(clean=self.hide_traceback)
        
        self.raw_input_dict = {}
        self.raw_output_dict = {}
        self.raw_pickup_dict = {}
        
        self.batch_T = []
        self.batch_P = []
        
        self.logK_models = {}
        self.df_rejected_species = pd.DataFrame({'database index':[], "name":[], "reason for rejection":[]})
        
        if load_thermo:
            self.thermo = Thermodata(
                     db=db,
                     solid_solutions=solid_solutions,
                     logK=logK,
                     logK_extrapolate=logK_extrapolate,
                     download_csv_files=download_csv_files,
                     exclude_category=exclude_category,
                     eq36da=self.eq36da,
                     eq36co=self.eq36co,
                     verbose=self.verbose,
                     hide_traceback=hide_traceback)
        
            self.data1 = self.thermo.data1

    def __capture_r_output(self):
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

    def __print_captured_r_output(self):
        printable_lines = [line for line in self.stdout if line not in ['[1]', '\n']]
        printable_lines = [line for line in printable_lines if re.search("^\s*\[[0-9]+\]$", line) is None]
        printable_lines = [re.sub(r' \\n\"', "", line) for line in printable_lines]
        [print(line[2:-1]) for line in printable_lines]

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
                                       dynamic_db, charge_balance_on, suppress_missing,
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
                err_bad_exclude = ("Could not exclude the header '{}'".format(exc)+". "
                                   "This header could not be found in {}".format(input_filename)+"")
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
        invalid_sample_names = [n for n in list(df_in.iloc[1:, 0]) if str(n[0])==" " or str(n[-1])==" "]
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
            start_index = [i+1 for i, s in enumerate(data0_lines) if '*  species name' in s]
            end_index = [i-1 for i, s in enumerate(data0_lines) if 'elements' in s]
            db_species = [i.split()[0] for i in data0_lines[start_index[0]:end_index[0]]]
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
                if species not in db_species and species not in ['Temperature', 'logfO2', 'pH', 'Pressure']+FIXED_SPECIES:
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
                            "logfO2", "Mineral", "bar"]
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
        
        args = ["cd", os.getcwd(), ";", self.eq36co+'/eqpt', os.getcwd()+"/data0."+db]
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
            
        args = ["cd", cwdd+path_3i, ";", # change directory to where 3i files are stored
                self.eq36co + '/eq3nr', # path to EQ3NR executable
                cwdd + data1_path + "/data1." + db, # path to data1 file
                cwdd + path_3i + "/" + filename_3i] # path to 3i file
        
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
            self.err_handler.raise_exception("Error: multiple output files detected for one speciation calculation.")
            
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
            self.err_handler.raise_exception("Error: multiple pickup files detected for one speciation calculation.")

                    
    def runeq6(self,
               filename_6i,
               db,
               samplename=None,
               path_6i="",
               path_6o="",
               path_6p="",
               path_extra_out="",
               data1_path="",
               dynamic_db_name=None):
        
        """
        Call EQ6 on a .6i input file.
        
        Parameters
        ----------
        filename_6i : str
            Name of 6i input file.
        
        db : str
            Three letter code of database.
        
        path_6i : path str, default current working directory
            Path of directory containing .6i input files.
            
        path_6o : path str, default current working directory
            Path of directory where .6o output files will be produced.
        
        path_6p : path str, default current working directory
            Path of directory where .6p pickup files will be produced.
            
        path_extra_out : path str, default current working directory
            Path of directory where additional output files, such as tab files,
            will be produced.
            
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

        # get current working dir
        cwd = os.getcwd()
        cwdd = cwd + "/"
        
        if samplename == None:
            samplename = filename_6i[:-3]
        
        if self.verbose > 0 and dynamic_db_name == None:
            print('Using ' + db + ' to react ' + samplename)
        elif self.verbose > 0 and isinstance(dynamic_db_name, str):
            print('Using ' + dynamic_db_name + ' to react ' + samplename)

        args = ["cd", cwdd+path_6i, ";", # change directory to 6i folder
                self.eq36co+'/eq6', # path of EQ6 executable
                cwdd+data1_path+"/data1."+db, # path of data1 file
                cwdd+path_6i + "/" + filename_6i] # path of 6i file
        
        args = " ".join(args)
        
        self.__run_script_and_wait(args) # run EQ6
        
        filename_6o = filename_6i[:-1] + 'o'
        filename_6p = filename_6i[:-1] + 'p'
        filename_6ba = filename_6i[:-1] + 'ba'
        filename_6bb = filename_6i[:-1] + 'bb'
        filename_6t = filename_6i[:-2] + 'csv'
        filename_6tx = filename_6i[:-1] + 'tx'

        # The new eq36 build truncates names, e.g., MLS.Source.3i creates MLS.3o
        # Correct for this here:
        files_6o = [file for file in os.listdir(cwdd+path_6i) if file[-3:] == ".6o"]
        files_6p = [file for file in os.listdir(cwdd+path_6i) if file[-3:] == ".6p"]
        files_6ba = [file for file in os.listdir(cwdd+path_6i) if file[-4:] == ".6ba"]
        files_6bb = [file for file in os.listdir(cwdd+path_6i) if file[-4:] == ".6bb"]
        files_6t = [file for file in os.listdir(cwdd+path_6i) if file[-3:] == ".6t"]
        files_6tx = [file for file in os.listdir(cwdd+path_6i) if file[-4:] == ".6tx"]
        
        if len(files_6o) == 0:
            if self.verbose > 0:
                print('Error: EQ6 failed to produce output for ' + filename_6i)
        elif len(files_6o) == 1:
            file_6o = files_6o[0]
            file_6ba = files_6ba[0]
            file_6bb = files_6bb[0]
            file_6t = files_6t[0]
            file_6tx = files_6tx[0]
            try:
                # move output
                shutil.move(cwdd+path_6i+"/"+file_6o, cwdd+path_6o+"/"+filename_6o)
                shutil.move(cwdd+path_6i+"/"+file_6ba, cwdd+path_extra_out+"/"+filename_6ba)
                shutil.move(cwdd+path_6i+"/"+file_6bb, cwdd+path_extra_out+"/"+filename_6bb)
                shutil.move(cwdd+path_6i+"/"+file_6t, cwdd+path_extra_out+"/"+filename_6t)
                shutil.move(cwdd+path_6i+"/"+file_6tx, cwdd+path_extra_out+"/"+filename_6tx)
            except:
                self.err_handler.raise_exception("Error: could not move", path_6i+"/"+file_6o, "to", path_6o+"/"+filename_6o)
        
        else:
            self.err_handler.raise_exception("Error: multiple output files detected for one mass transfer calculation.")
            
        if len(files_6p) == 0:
            if self.verbose > 0:
                print('Error: EQ6 failed to produce a pickup file for ' + filename_6i)
        elif len(files_6p) == 1:
            file_6p = files_6p[0]
            try:
                # move output
                shutil.move(cwdd+path_6i+"/"+file_6p, cwdd+path_6p+"/"+filename_6p)
            except:
                self.err_handler.raise_exception("Error: could not move", path_6i+"/"+file_6p, "to", path_6p+"/"+filename_6p)
        else:
            self.err_handler.raise_exception("Error: multiple pickup files detected for one mass transfer calculation.")
        
                
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

        print("R COEFFS")
        print(poly_coeffs_1)
        print(poly_coeffs_2)
        
        
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
                 get_affinity_energy=False,
                 negative_energy_supplies=False,
                 rxn_filename=None,
                 not_limiting=["H+", "OH-", "H2O"],
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
            Calculate affinities and energy supplies of reactions listed in a
            separate user-supplied file?
        
        negative_energy_supplies : bool, default False
            Report negative energy supplies? If False, negative energy supplies
            are reported as 0 cal/kg H2O. If True, negative energy supplies are
            reported. A 'negative energy supply' represents the energy cost of
            depleting the limiting reactant of a reaction. This metric is not
            always helpful when examing energy supply results, so this option is
            set to False by default.
        
        rxn_filename : str, optional
            Name of .txt file containing reactions used to calculate affinities
            and energy supplies. Ignored if `get_affinity_energy` is False.
        
        not_limiting : list, default ["H+", "OH-", "H2O"]
            List containing names of species that are not considered limiting
            when calculating energy supplies. Ignored if `get_affinity_energy`
            is False.
        
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
        
        if db != None:
            # load new thermodynamic database
            self.thermo._set_active_db(db, self.verbose)
        else:
            db = self.thermo.db
            
        if self.thermo.thermo_db_type == "CSV":
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
            db_args["generate_template"] = False
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
            
                    # Download from URL and decode as UTF-8 text.
                    with urlopen(db_solid_solution) as webpage:
                        content = webpage.read().decode()
                        
                    # Save to CSV file.
                    with open(db_solid_solution_csv_name, 'w') as output:
                        output.write(content)
                        
                    db_args["filename_ss"] = db_solid_solution_csv_name
                    
            if self.verbose > 0:
                print("Getting", self.thermo.thermo_db_filename, "ready. This will take a moment...")
    
            OBIGT_df, data0_file_lines, grid_temps, grid_press, data0_lettercode, water_model, P1, plot_poly_fit = self.create_data0(**db_args)
            
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
            
            data1_path = "eqpt_files" # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA
            
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
                
            self.__capture_r_output()
        
            r_check_TP_grid = pkg_resources.resource_string(__name__, 'check_TP_grid.r').decode("utf-8")
        
            ro.r(r_check_TP_grid)
        
            list_tp = ro.r.check_TP_grid(grid_temps=_convert_to_RVector(grid_temps),
                                         grid_press=_convert_to_RVector(grid_press),
                                         P1=P1,
                                         water_model=water_model,
                                         check_for_errors=False,
                                         verbose=self.verbose)
        
            self.__print_captured_r_output()
            
            grid_temps = list(list_tp.rx2("grid_temps"))
            grid_press = list(list_tp.rx2("grid_press"))
            poly_coeffs_1 = list_tp.rx2("poly_coeffs_1")
            poly_coeffs_2 = list_tp.rx2("poly_coeffs_2")
            
            
        else:
            grid_temps = ro.r("NULL")
            grid_press = ro.r("NULL")
            poly_coeffs_1 = ro.r("NULL")
            poly_coeffs_2 = ro.r("NULL")
            
            
        if get_affinity_energy:
            if rxn_filename == None and self.affinity_energy_reactions_raw==None:
                err = ("get_affinity_energy is set to True but a reaction TXT "
                       "file is not specified or redox reactions have not yet "
                       "been generated with make_redox_reactions()")
                self.err_handler.raise_exception(err)
            elif rxn_filename != None:
                self.__file_exists(rxn_filename, '.txt')
                
                self.affinity_energy_reactions_raw = pd.read_csv(rxn_filename, sep="\t", header=None, names=["col"+str(i) for i in range(1,50)])
                load_rxn_file = True
            else:
                if self.thermo.thermo_db_type != "CSV":
                    if self.verbose > 0:
                        warn_msg = ("Warning: get_affinity_energy is set to True but "
                            "the active thermodynamic database ("+self.thermo.db+") is not "
                            "in CSV format. This indicates a possible mismatch between "
                            "the thermodynamic database used to generate redox reactions "
                            "and the one used in this speciation calculation. Continuing anyway...")
                        print(warn_msg)
                rxn_filename = self.affinity_energy_reactions_raw
                load_rxn_file = False
            

            
        else:
            rxn_filename = ""
            load_rxn_file=False

        
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
        self.__capture_r_output()
        
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
        
        self.__print_captured_r_output()
        
        self.df_input_processed = ro.conversion.rpy2py(input_processed_list.rx2("df"))
        
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
                
                self.__fill_data0(OBIGT_df=ro.conversion.rpy2py(OBIGT_df),
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

                data1_path = "eqpt_files" # creating a folder name without spaces to store the data1 overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA

                data0_path = "data0." + data0_lettercode
                
            else:
                pressure_bar = list(input_processed_list.rx2("pressure_bar"))[sample_row_index]
            
            # allowed aq block species are left after any category exclusion in db_args
            allowed_aq_block_species = ["all"]
            if dynamic_db:
                allowed_aq_block_species = list(OBIGT_df["name"]) + FIXED_SPECIES
            
            # write 3i files
            self.__capture_r_output()

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

            self.__print_captured_r_output()
        
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
                self.raw_input_dict[samplename] = lines
            except:
                pass
            try:
                with open(output_dir + "/" + filename_3o, "r") as f:
                    lines=f.readlines()
                self.raw_output_dict[samplename] = lines
            except:
                pass
            try:
                with open(pickup_dir + "/" + filename_3p, "r") as f:
                    lines=f.readlines()
                self.raw_pickup_dict[samplename] = lines
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
        
        if self.thermo.thermo_db_type == "CSV":
            custom_obigt = self.thermo.thermo_db
        else:
            custom_obigt = ro.r("NULL")
        
        # mine output
        self.__capture_r_output()
        
        r_3o_mine = pkg_resources.resource_string(
            __name__, '3o_mine.r').decode("utf-8")
        ro.r(r_3o_mine)
        batch_3o = ro.r.main_3o_mine(
            files_3o=_convert_to_RVector(files_3o),
            input_filename=input_filename,
            input_pressures=_convert_to_RVector(list(input_processed_list.rx2("pressure_bar"))),
            rxn_filename=rxn_filename,
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
            get_affinity_energy=get_affinity_energy,
            negative_energy_supplies=negative_energy_supplies,
            load_rxn_file=load_rxn_file,
            not_limiting=_convert_to_RVector(not_limiting),
            batch_3o_filename=batch_3o_filename,
            df_input_processed=ro.conversion.py2rpy(self.df_input_processed),
            # New rpy2 py2rpy2 conversion might not need the workaround below.
            # The old note regarding deprecated pandas2ri is shown below...
            # OLD NOTE:
            # Needed for keeping symbols in column names after porting
            #   df_input_processed in the line above. Some kind of check.names
            #   option for pandas2ri.py2ri would be nice. Workaround:
            df_input_processed_names=df_input_processed_names,
            custom_obigt=custom_obigt,
            water_model=water_model,
            fixed_species=_convert_to_RVector(FIXED_SPECIES),
            verbose=self.verbose,
        )

        self.__print_captured_r_output()
        
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
            aq_dist_indx = report_divs.names.index("aq_distribution")
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
            
        if get_affinity_energy:
            affinity_cols = list(report_divs.rx2('affinity'))
            energy_cols = list(report_divs.rx2('energy'))
            
            df_affinity = df_report[affinity_cols]
            df_energy = df_report[energy_cols]
            
            df_affinity = df_affinity.apply(pd.to_numeric, errors='coerce')
            df_energy = df_energy.apply(pd.to_numeric, errors='coerce')
            
            # handle headers of df_affinity section
            headers = df_affinity.columns
            subheaders = ['cal/mol e-']*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_affinity.columns = multicolumns

            # handle headers of df_energy section
            headers = df_energy.columns
            subheaders = ['cal/kg.H2O']*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_energy.columns = multicolumns
            df_join = df_join.join(df_affinity)
            df_join = df_join.join(df_energy)

        out_dict = {'sample_data': {},
                    'report': df_join,
                    'input': df_input, 'report_divs': report_divs}
        
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
                sample_solid_solutions = batch_3o.rx2["sample_data"].rx2[sample.rx2('name')[0]].rx2["solid_solutions"]

                if not type(sample_solid_solutions.names) == rpy2.rinterface_lib.sexp.NULLType:

                    ss_df_list = []
                    for ss in list(sample_solid_solutions.names):
                        df_ss_ideal = ro.conversion.rpy2py(sample_solid_solutions.rx2[ss].rx2["ideal solution"])
                        df_ss_mineral = ro.conversion.rpy2py(sample_solid_solutions.rx2[ss].rx2["mineral"])
                        df_merged = pd.merge(df_ss_mineral, df_ss_ideal, left_on='mineral', right_on='component', how='left')
                        df_merged.insert(0, 'solid solution', ss)
                        del df_merged['component']
                        ss_df_list.append(df_merged)
                
                    dict_sample_data.update(
                        {"solid_solutions": pd.concat(ss_df_list)})
            
            if get_affinity_energy:
                dict_sample_data.update({"affinity_energy_raw": ro.conversion.rpy2py(
                    sample.rx2('affinity_energy_raw'))})
                dict_sample_data.update(
                    {"affinity_energy": ro.conversion.rpy2py(sample.rx2('affinity_energy'))})

            out_dict["sample_data"].update(
                {sample_data.names[i]: dict_sample_data})

        out_dict.update({"batch_3o": batch_3o})
        
        out_dict.update({"water_model":water_model, "grid_temps":grid_temps, "grid_press":grid_press})
        
        speciation = Speciation(out_dict, hide_traceback=self.hide_traceback)
        
        if get_affinity_energy:
            speciation.half_cell_reactions = self.half_cell_reactions
            speciation.affinity_energy_reactions_table = self.affinity_energy_reactions_table
            speciation.affinity_energy_formatted_reactions = self.affinity_energy_formatted_reactions
            speciation.show_redox_reactions = self.show_redox_reactions
        
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
        
        speciation.raw_input_dict = self.raw_input_dict
        speciation.raw_output_dict = self.raw_output_dict
        speciation.raw_pickup_dict = self.raw_pickup_dict
        speciation.thermo = self.thermo
        speciation.data1 = self.data1
        
        speciation.df_rejected_species = self.df_rejected_species
        speciation.logK_models = self.logK_models
        speciation.batch_T = self.batch_T
        speciation.batch_P = self.batch_P
        
        return speciation

    @staticmethod
    def __clean_rpy2_pandas_conversion(df,
                                       float_cols=["G", "H", "S", "Cp",
                                                    "V", "a1.a", "a2.b",
                                                    "a3.c", "a4.d", "c1.e",
                                                    "c2.f", "omega.lambda", "z.T",
                                                    "azero", "neutral_ion_type",
                                                    "logK1", "logK2", "logK3", "logK4",
                                                    "logK5", "logK6", "logK7", "logK8",
                                                    "T1", "T2", "T3", "T4", "T5", "T6",
                                                    "T7", "T8"],
                                        str_cols=["name", "abbrv", "state", "formula",
                                                  "ref1", "ref2", "date",
                                                  "E_units", "tag", "dissrxn", "formula_ox",
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
    

    @staticmethod
    def __s_d(x, k):
        # specify how many decimals are printed
        # e.g. 12.433 becomes "12.4330" if k=4
        kstr = '{:.'+str(k)+'f}'
        return kstr.format(round(x, k)).strip()
    
    
    def __fill_data0(self, OBIGT_df, data0_file_lines, grid_temps, grid_press, db,
                   water_model, activity_model, P1, plot_poly_fit, logK_extrapolate,
                   dynamic_db, verbose):
        
        
        self.__capture_r_output()
        
        r_check_TP_grid = pkg_resources.resource_string(
            __name__, 'check_TP_grid.r').decode("utf-8")
        
        ro.r(r_check_TP_grid)
        
        list_tp = ro.r.check_TP_grid(grid_temps=_convert_to_RVector(grid_temps),
                                     grid_press=_convert_to_RVector(grid_press),
                                     P1=P1,
                                     water_model=water_model,
                                     check_for_errors=True,
                                     verbose=self.verbose)
        
        self.__print_captured_r_output()
        
        grid_temps = list(list_tp.rx2("grid_temps"))
        grid_press = list(list_tp.rx2("grid_press"))
        
        if plot_poly_fit and len(grid_temps) == 8:
            self.__plot_TP_grid_polyfit(xvals=grid_temps,
                                        yvals=grid_press,
                                        poly_coeffs_1=list(list_tp.rx2("poly_coeffs_1")),
                                        poly_coeffs_2=list(list_tp.rx2("poly_coeffs_2")),
                                        res=500)

        self.__print_captured_r_output()
        
        # calculate logK at each T and P for every species
        out_dfs = []
        for i,Tc in enumerate(grid_temps):
            out_dfs.append(calc_logK(OBIGT_df, Tc=Tc, P=grid_press[i], TP_i=i, water_model=water_model))
        
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
        if "logK1" in OBIGT_df.columns:
            
            free_logK_df = OBIGT_df.dropna(subset=['logK1'])
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
        self.__capture_r_output()
        
        r_fill_data0_header = pkg_resources.resource_string(
            __name__, 'fill_data0_header.r').decode("utf-8")
        
        ro.r(r_fill_data0_header)
        
        data0_file_lines = ro.r.fill_data0_head(data0_template=data0_file_lines,
                                       db=db,
                                       grid_temps=_convert_to_RVector(grid_temps),
                                       grid_press=_convert_to_RVector(grid_press),
                                       water_model=water_model,
                                       activity_model=activity_model)
        
        self.__print_captured_r_output()
        
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
                if name not in list(self.df_rejected_species["name"]):
                    msg = "The chemical species " + str(name) + " is not recognized."
                    self.err_handler.raise_exception(msg)
                else:
                    reject_reason = list(self.df_rejected_species.loc[self.df_rejected_species['name'] == name, 'reason for rejection'])[0]
                    
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
                        print(sp)
                        print(logK_extrapolate)
                        
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
            reject_reasons =list(rejected_sp_i_dict.values())
            
            df_rejected_species = pd.DataFrame({'database index':reject_indices, "name":reject_names, "reason for rejection":reject_reasons})
                        
            return valid_sp_i, df_rejected_species
            
            
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
                     data0_formula_ox_name=None,
                     suppress_redox=[],
                     water_model="SUPCRT92",
                     activity_model="b-dot",
                     exceed_Ttr=True,
                     grid_temps=[0.0100, 50.0000, 100.0000, 150.0000,
                                 200.0000, 250.0000, 300.0000, 350.0000],
                     grid_press="Psat",
                     P1=True,
                     plot_poly_fit=False,
                     logK_extrapolate="none",
                     infer_formula_ox=False,
                     generate_template=True,
                     template_name=None,
                     template_type="strict",
                     exclude_category={},
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
        
        data0_formula_ox_name : str, optional
            Name of supplementary file containing data0 parameters and inferred
            formula oxidation states. Ignored if `infer_formula_ox` is False.
            See `infer_formula_ox` for more detail.
        
        suppress_redox : list of str, default []
            Suppress equilibrium between oxidation states of listed elements
            (Cl, H, and O cannot be included).

        water_model : str, default "SUPCRT92"
            This is an experimental feature that is not yet fully supported.
            Desired water model. Can be either "SUPCRT92", "IAPWS95", or "DEW".
            These models are described here: http://chnosz.net/manual/water.html

        exceed_Ttr : bool, default True
            Calculate Gibbs energies of mineral phases and other species
            beyond their transition temperatures?

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
        
        infer_formula_ox : bool, default False
            Create a supplementary file containing data0 parameters and
            inferred formula oxidation states? This option is useful for
            creating as many entries in the formula_ox column when creating a
            new supplementary file. Note that compounds like DySO4+ result in
            blank entries in formula_ox because the redox states of two
            elements, Dy and S, would have to be estimated together; S has many
            oxidation states and Dy's oxidation states are not hard-coded.
        
        generate_template : bool, default True
            Generate a CSV sample input template customized to this data0?
            Columns include 'Sample', 'Temperature', 'logfO2', and all strict
            basis species.
        
        template_name : str, optional
            Name of the sample input template file generated. If no name is
            supplied, defaults to 'sample_template_xyz.csv', where 'xyz' is
            the three letter code given to `db`. Ignored if `generate_template`
            is False.
        
        template_type : str, either 'strict', 'all basis', or 'all species'
            Determines which columns are written to the sample template.
            - 'strict' includes strict basis species
            - 'all basis' includes strict and auxiliary basis species
            - 'all species' includes all species in the thermodynamic database
            Ignored if `generate_template` is False.
        
        exclude_category : dict
            Exclude species from the custom thermodynamic dataset based on
            column values. For instance,
            `exclude_category={'category_1':["organic_aq", "organic_cr"]}`
            will exclude all species that have "organic_aq" or "organic_cr" in
            the column "category_1".
        
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
            if sum([P >= 10000 for P in grid_press]):
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
    
        if generate_template:
            if template_name == None:
                template_name = "sample_template_{}.csv".format(db)
            
            try:
                # check if template can be generated in specified location
                with open(template_name, 'w') as fp:
                    pass
            except:
                self.err_handler.raise_exception("The file {} could not be ".format(template_name)+""
                    "created. Is this a valid file path?")
        
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
            
            free_logK_df = self.__clean_rpy2_pandas_conversion(self.thermo.logK_db)

            valid_i, self.df_rejected_species = self.__get_i_of_valid_free_logK_sp(
                free_logK_df,
                grid_or_sample_temps,
                grid_or_sample_press,
                dynamic_db,
                logK_extrapolate,
                db_sp_names=thermo_df["name"],
                )
            free_logK_df_valid = copy.deepcopy(free_logK_df.iloc[valid_i])
            thermo_df = pd.concat([thermo_df, free_logK_df_valid], ignore_index=True)
            
            thermo_df = self.__clean_rpy2_pandas_conversion(thermo_df)
        
        template = pkg_resources.resource_string(
            __name__, 'data0.min').decode("utf-8")
        suppress_redox = _convert_to_RVector(suppress_redox)
        
        if self.thermo.solid_solutions_active:
            solid_solution_df = ro.conversion.py2rpy(self.thermo.solid_solution_db)
        else:
            solid_solution_df = ro.r("NULL")
            
        if data0_formula_ox_name == None:
            data0_formula_ox_name = ro.r("NULL")
        if template_name == None:
            template_name = "sample_template_{}.csv".format(db)
        if template_type not in ['strict', 'all basis', 'all species']:
            self.err_handler.raise_exception("template_type {} ".format(template_type)+"is not"
                            "recognized. Try 'strict', 'all basis', or 'all species'")

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
                exclude_category["name"] = sp_names_to_exclude
                if self.verbose > 0 and len(sp_names_to_exclude) > 0:
                    print("Excluding the following chemical species because "
                          "they contain redox-suppressed elements but do not "
                          "have element oxidation states given in the "
                          "'formula_ox' column of the thermodynamic database: "
                          ""+str(sp_names_to_exclude))
            
        if len(exclude_category) > 0:
            exclude_category_R =  {k:_convert_to_RVector(l) for k,l in zip(exclude_category.keys(), exclude_category.values())}
        else:
            exclude_category_R = {}
        exclude_category_R = ro.ListVector(exclude_category_R)
        
        self.__capture_r_output()

        r_redox_dissrxns = pkg_resources.resource_string(
            __name__, 'redox_and_dissrxns.r').decode("utf-8")
        
        ro.r(r_redox_dissrxns)
        
#         thermo_df["ref2"] = thermo_df["ref2"].astype(str)
#         thermo_df["dissrxn"] = thermo_df["dissrxn"].astype(str)
#         thermo_df["tag"] = thermo_df["tag"].astype(str)

        thermo_df = self.__clean_rpy2_pandas_conversion(thermo_df)

        ro.conversion.py2rpy(thermo_df)
        
        out_list = ro.r.suppress_redox_and_generate_dissrxns(thermo_df=ro.conversion.py2rpy(thermo_df),
                               db=db,
                               water_model=water_model,
                               template=template,
                               exceed_Ttr=exceed_Ttr,
                               data0_formula_ox_name=data0_formula_ox_name,
                               suppress_redox=suppress_redox,
                               infer_formula_ox=infer_formula_ox,
                               exclude_category=exclude_category_R,
                               fixed_species=_convert_to_RVector(FIXED_SPECIES),
                               verbose=self.verbose)
        
        self.__print_captured_r_output()
        
        OBIGT_df = out_list.rx2("OBIGT_df")
        OBIGT_df=ro.conversion.rpy2py(OBIGT_df)
        
#         regenerated_dissrxns = out_list.rx2("dissrxns")
        
#         regenerated_dissrxn_dict = {}
#         for name in regenerated_dissrxns.names:
#             if name != "basis_list":
#                 regenerated_dissrxn_dict[name] = regenerated_dissrxns.rx2(name)[0]


        OBIGT_df = self.__clean_rpy2_pandas_conversion(OBIGT_df)
        
        # convert E units and calculate missing GHS values
        OBIGT_df = OBIGT2eos(OBIGT_df, fixGHS=True, tocal=True)
        
        self.__capture_r_output()
        
        r_create_data0 = pkg_resources.resource_string(
            __name__, 'create_data0.r').decode("utf-8")
        
        ro.r(r_create_data0)
        
        ro.conversion.py2rpy(OBIGT_df)
        
        # assemble data0 file
        data0_file_lines = ro.r.create_data0(thermo_df=ro.conversion.py2rpy(OBIGT_df),
                          solid_solution_df=solid_solution_df,
                          db=db,
                          water_model=water_model,
                          template=template,
                          dissrxns=out_list.rx2("dissrxns"),
                          basis_pref=out_list.rx2("basis_pref"),
                          exceed_Ttr=exceed_Ttr,
                          fixed_species=_convert_to_RVector(FIXED_SPECIES),
                          verbose=self.verbose)
        
        self.__print_captured_r_output()
        
        data0_file_lines = data0_file_lines[0].split("\n")
#         print("lines:")
#         print(str(data0_file_lines))
#         with open('test.txt', 'w') as f:
#             for line in data0_file_lines:
#                 f.write(f"{line}\n")

        if generate_template:
            
            r_generate_template = pkg_resources.resource_string(
                __name__, 'generate_template.r').decode("utf-8")
        
            ro.r(r_generate_template)
            
            self.__capture_r_output()
            
            ro.r.generate_template(thermo_df=ro.conversion.py2rpy(OBIGT_df),
                                   template_name=template_name,
                                   template_type=template_type,
                                   fixed_species=_convert_to_RVector(FIXED_SPECIES))

            self.__print_captured_r_output()
        
        if fill_data0:
            
            # begin TP-dependent processes
            self.__fill_data0(OBIGT_df=ro.conversion.rpy2py(OBIGT_df),
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
            return OBIGT_df, data0_file_lines, grid_temps, grid_press, db, water_model, P1, plot_poly_fit

        if self.verbose > 0:
            print("Finished creating data0.{}.".format(db))
            

    def make_redox_reactions(self, db=None, redox_pairs="all", auto_load_db=True):
        
        """
        Generate an organized collection of redox reactions for calculating
        chemical affinity and energy supply values during speciation.
        
        Parameters
        ----------
        db : str
            Determines which thermodynamic database is used in the speciation
            calculation. The database must be a CSV file (not a data0file)
            because the code must look up properties of chemical species to
            calculate affinities and energy supplies of reactions.
            The `db` parameter can either be:
            - The name of a CSV file containing thermodynamic data located in
            the current working directory, e.g., "wrm_data.csv". The CSV file
            will be used to generate a data0 file for each sample (using
            additional arguments from `db_args` if desired).
            - The URL of a CSV file containing thermodynamic data, e.g.,
            "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
        
        redox_pairs : list of int or "all", default "all"
            List of indices of half reactions in the half cell reaction table
            to be combined when generating full redox reactions.
            E.g. [0, 1, 4] will combine half reactions with indices 0, 1, and 4
            in the table stored in the `half_cell_reactions` attribute of the
            `AqEquil` class.
            If "all", generate all possible redox reactions from available half
            cell reactions.
        
        auto_load_db : bool, default True
            Automatically download and use a WORM-styled CSV if the currently
            active thermodynamic database does not support affinity and energy
            supply calculations? If True, the most up-to-date copy of the
            wrm_data.csv will be downloaded from the URL
            https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv
            and set as the active thermodynamic database.
        
        Returns
        ----------
        Output is stored in the `affinity_energy_reactions_raw` and
        `affinity_energy_reactions_table` attributes of the `AqEquil` class.
        """
        
        if db != None:
            self.thermo._set_active_db(db, self.verbose)
            
        if self.thermo.thermo_db_type != "CSV":
            if self.verbose > 0:
                if auto_load_db:
                    print("Warning: Redox reactions require a WORM-styled thermodynamic database CSV file.")
                else:
                    self.err_handler.raise_exception("Error: Redox reactions require a WORM-styled CSV file as the active thermodynamic database.")
            
            if auto_load_db:
                if self.verbose > 0:
                    print("Warning: switching thermodynamic database from", str(self.thermo.thermo_db_filename), "to wrm_data.csv...")
                self.thermo._set_active_db(db="https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv",
                                    verbose=self.verbose)
            
        db = self.thermo.db
        
        # reset all redox variables stored in the AqEquil class
        self.affinity_energy_reactions_raw = None
        self.affinity_energy_reactions_table = None
        self.affinity_energy_formatted_reactions = None
        
        if self.verbose > 1:
            print("Generating redox reactions...")

        err_msg = ("redox_pairs can either be 'all' or a list of integers "
               "indicating the indices of half cell reactions in "
               "the half_cell_reactions table that should be combined into "
               "full redox reactions. For example, redox_pairs=[0, 1, 2, 6] "
               "will combine half cell reactions with indices 0, 1, 2, and 6 in "
               "the half_cell_reactions table. This table is an attribute in the "
               "class AqEquil.")
        if isinstance(redox_pairs, str):
            if redox_pairs == "all":
                redox_pairs = list(range(0, self.half_cell_reactions.shape[0]))
            else:
                self.err_handler.raise_exception(err_msg)
        elif isinstance(redox_pairs, list):
            if not all([isinstance(i, int) for i in redox_pairs]):
                self.err_handler.raise_exception(err_msg)
        else:
            self.err_handler.raise_exception(err_msg)
        
        self.redox_pairs = redox_pairs
        
        df = self.half_cell_reactions.iloc[redox_pairs].reset_index(drop=True)
        
        wrm_data = self.thermo.thermo_db
        basis_df = wrm_data.loc[wrm_data['tag'] == 'basis']
        
        db_names = []
        formulas = []
        for column in list(df.columns)[1:]:
            for item in list(df[column]):
                if item != 'nan':
                    if item in list(wrm_data.name) and isinstance(item, str):
                        index = wrm_data.loc[wrm_data['name'] == item].index[0]
                        formula = wrm_data.loc[wrm_data['name'] == item]['formula'][index]
                        df.replace(item, formula, inplace=True)
                        if item not in db_names:
                            db_names.append(item)
                            formulas.append(formula)
                    elif item in list(wrm_data.abbrv) and isinstance(item, str):
                        index = wrm_data.loc[wrm_data['abbrv'] == item].index[0]
                        formula = wrm_data.loc[wrm_data['abbrv'] == item]['formula'][index]
                        df.replace(item, formula, inplace=True)
                        if item not in db_names:
                            db_names.append(item)
                            formulas.append(formula)
                            
        df.replace('sulfur', 'S', inplace = True) ### remove this eventually
        db_names.append('sulfur') ### remove this eventually
        formulas.append('S') ### remove this eventually

        # append H+ and H2O to db for balancing
        if not 'H+' in db_names:
            db_names.append('H+')
            formulas.append('H+')
        if not 'H2O' in db_names:
            db_names.append('H2O')
            formulas.append('H2O')
        
        Oxidant_1 = df['Oxidant_1']
        Oxidant_2 = df['Oxidant_2']
        Oxidant_3 = df['Oxidant_3']
        Reductant_1 = df['Reductant_1']
        Reductant_2 = df['Reductant_2']
        
        #CREATING A LIST OF ALL SPECIES AND THEIR ELEMENT DICTIONARIES
        elements = [str(e) for e in list(periodictable.elements)[1:]]+['+','-']
        element_dictionary = dict()
        for i in formulas:
                parsed_formula = parse_formula(i)
                element_dictionary[i] = parsed_formula
                for e in elements:
                    if element_dictionary[i].get(e, 0) == 0:
                        element_dictionary[i][e] = 0
                        
        df_reax = pd.DataFrame() # empty df of reactions
        df_reax['rO_coeff'] = ''
        df_reax['rO'] = ''
        df_reax['rR_coeff'] = ''
        df_reax['rR'] = ''
        df_reax['pO_coeff'] = ''
        df_reax['pO'] = ''
        df_reax['pR_coeff'] = ''
        df_reax['pR'] = ''
        df_reax['Reaction'] = ''
        df_reax['redox_pair'] = ''
        index = 0

        reaction = [] 
        indices = np.arange(0, len(Oxidant_1), 1).tolist()*2 # how to loop back through the redox pairs to not run into out-of-range index
        rxn_num = 0 # counting unique reactions
        rxn_list = [] # list of unique reactions
        rxn_names = []
        rxn_pairs = [] # list of paired half reactions

        for i in range(0, len(Oxidant_1),1): # length of redox pairs - columns
            for n in range(0, len(Oxidant_1), 1):
                if Reductant_1[i] == Reductant_1[indices[i+n]] or Reductant_1[i] == Reductant_2[indices[i+n]]: # if both reductants are the same thing, skip
                    continue
                if Oxidant_1[i] == Oxidant_1[indices[i+n]] or Oxidant_1[i] == Oxidant_2[indices[i+n]] or Oxidant_1[i] == Oxidant_3[indices[i+n]]:
                    continue
#                 if Oxidant_1[i] == 'H2O' and Reductant_1[indices[i+n]] == 'H2O': #suppress the splitting of water
#                     continue
                else:

                    # GENERATING REACTIONS BETWEEN OXIDANT_1 AND REDUCTANT_1
                    reaction.append(Oxidant_1[i]+' \t ' + Reductant_1[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                    rxn_num+=1
                    rxn_list.append(rxn_num)
                    rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                    df_reax.loc[index, 'rO'] = Oxidant_1[i]
                    df_reax.loc[index, 'rR'] = Reductant_1[indices[i+n]]
                    df_reax.loc[index, 'pR'] = Reductant_1[i]
                    df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                    rxn_pairs.append([i, indices[i+n]])
                    index+=1

                # REACTIONS INVOLVING OTHER PH-DEPENDENT SPECIES
                if pd.isnull(Oxidant_2[i]) != True and Oxidant_2[i] != Reductant_2[indices[i+n]]:
                    reaction.append(Oxidant_2[i] + ' \t ' + Reductant_1[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                    rxn_list.append(rxn_num)
                    rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                    df_reax.loc[index, 'rO'] = Oxidant_2[i]
                    df_reax.loc[index, 'rR'] = Reductant_1[indices[i+n]]
                    df_reax.loc[index, 'pR'] = Reductant_1[i]
                    df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                    rxn_pairs.append([i, indices[i+n]])
                    index +=1

                    if pd.isnull(Reductant_2[indices[i+n]]) != True:
                        reaction.append(Oxidant_2[i] + ' \t ' + Reductant_2[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                        rxn_list.append(rxn_num)
                        rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                        df_reax.loc[index, 'rO'] = Oxidant_2[i]
                        df_reax.loc[index, 'rR'] = Reductant_2[indices[i+n]]
                        df_reax.loc[index, 'pR'] = Reductant_1[i]
                        df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                        rxn_pairs.append([i, indices[i+n]])
                        index +=1
                        
                if pd.isnull(Oxidant_2[i]) != True and Oxidant_2[i] == Reductant_2[indices[i+n]]:
                    reaction.append(Oxidant_2[i] + ' \t ' + Reductant_2[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                    rxn_list.append(rxn_num)
                    rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                    df_reax.loc[index, 'rO'] = Oxidant_2[i]
                    df_reax.loc[index, 'rR'] = Reductant_2[indices[i+n]]
                    df_reax.loc[index, 'pR'] = Reductant_1[i]
                    df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                    rxn_pairs.append([i, indices[i+n]])
                    index +=1

                if pd.isnull(Reductant_2[indices[i+n]]) != True and Oxidant_2[i] != Reductant_2[indices[i+n]]:
                    reaction.append(Oxidant_1[i] + ' \t ' + Reductant_2[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                    rxn_list.append(rxn_num)
                    rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                    df_reax.loc[index, 'rO'] = Oxidant_1[i]
                    df_reax.loc[index, 'rR'] = Reductant_2[indices[i+n]]
                    df_reax.loc[index, 'pR'] = Reductant_1[i]
                    df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                    rxn_pairs.append([i, indices[i+n]])
                    index +=1

                if pd.isnull(Oxidant_3[i]) != True:
                    reaction.append(Oxidant_3[i] + ' \t ' + Reductant_1[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                    rxn_list.append(rxn_num)
                    rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                    df_reax.loc[index, 'rO'] = Oxidant_3[i]
                    df_reax.loc[index, 'rR'] = Reductant_1[indices[i+n]]
                    df_reax.loc[index, 'pR'] = Reductant_1[i]
                    df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                    rxn_pairs.append([i, indices[i+n]])
                    index +=1

                    if pd.isnull(Reductant_2[indices[i+n]]) != True:
                        reaction.append(Oxidant_3[i] + ' \t ' + Reductant_2[indices[i+n]] + '=' + Reductant_1[i] + ' \t ' + Oxidant_1[indices[i+n]])
                        rxn_list.append(rxn_num)
                        rxn_names.append('red_'+Oxidant_1[i]+'_'+Reductant_1[i]+'_ox_'+Reductant_1[indices[i+n]]+'_'+Oxidant_1[indices[i+n]])
                        df_reax.loc[index, 'rO'] = Oxidant_3[i]
                        df_reax.loc[index, 'rR'] = Reductant_2[indices[i+n]]
                        df_reax.loc[index, 'pR'] = Reductant_1[i]
                        df_reax.loc[index, 'pO'] = Oxidant_1[indices[i+n]]
                        rxn_pairs.append([i, indices[i+n]])
                        index +=1
        
        df_reax['Reaction'] = rxn_list
        df_reax['Names'] = rxn_names
        df_reax['Temp_Pairs'] = rxn_pairs
        
        # if there are no reactions, return nothing
        if df_reax.shape[0] == 0:
            incompatible_half_reactions = pd.Series(self.half_cell_reactions["Redox Couple"], index=redox_pairs).tolist()
            redundant_reductant_or_oxidant = []
            for col in ["Oxidant_1", "Oxidant_2", "Oxidant_3", "Reductant_1", "Reductant_2"]:
                redox_col = pd.Series(self.half_cell_reactions[col], index=[0,1])
                if redox_col.eq(redox_col[0]).all():
                    redundant_reductant_or_oxidant.append(redox_col[0])
            err_no_rxns = ("Valid reactions could not be written between the half "
                "reactions {} ".format(incompatible_half_reactions)+"because "
                "{}".format(redundant_reductant_or_oxidant)+" is on both sides "
                "of all reactions.")
            print(err_no_rxns)
            return
        
        ### BALANCING NON-O, H ELEMENTS
        for r in range(0, len(df_reax['rO'])):
            count = 0 #to restart the loop through the elements
            temp_rO_coeff = [1] *(len(elements)-4) #loop through all elements except O, H, +, and -
            temp_rR_coeff = [1] *(len(elements)-4) #loop through all elements except O, H, +, and -
            temp_pO_coeff = [1] *(len(elements)-4) #loop through all elements except O, H, +, and -
            temp_pR_coeff = [1] *(len(elements)-4) #loop through all elements except O, H, +, and -

            
            for e in elements:
                if e in ['O','H','+','-']:
                    continue
                else:
            
                    temp1 = int(element_dictionary[df_reax['rO'][r]][e]) #count for the element in the list for rO at index r
                    temp2 = int(element_dictionary[df_reax['pR'][r]][e])
                    temp3 = int(element_dictionary[df_reax['rR'][r]][e])
                    temp4 = int(element_dictionary[df_reax['pO'][r]][e])
                    if temp1 == temp2:
                        temp_rO_coeff[count] = 1
                        temp_pR_coeff[count] = 1
                    if temp1 != temp2:
                        if temp1 ==0:
                            temp_rO_coeff[count] = 1
                        if temp1 != 0:
                            temp_rO_coeff[count] = np.lcm(temp1,temp2)/temp1
                        if temp2 == 0:
                            temp_pR_coeff[count] = 1
                        if temp2 != 0:
                            temp_pR_coeff[count] = np.lcm(temp1,temp2)/temp2
                    if temp3 == temp4:
                        temp_rR_coeff[count] = 1
                        temp_pO_coeff[count] = 1
                    if temp4 != temp3:
                        if temp3 == 0:
                            temp_rR_coeff[count] = 1
                        if temp3 != 0:
                            temp_rR_coeff[count] = np.lcm(temp3,temp4)/temp3
                        if temp4 == 0.0:
                            temp_pO_coeff[count] = 1
                        if temp4 !=0.0:
                            temp_pO_coeff[count] = np.lcm(temp3,temp4)/temp4
                    count +=1
            df_reax.loc[r, 'rO_coeff'] = -max(temp_rO_coeff)
            df_reax.loc[r, 'rR_coeff'] = -max(temp_rR_coeff)
            df_reax.loc[r, 'pR_coeff'] = max(temp_pR_coeff)
            df_reax.loc[r, 'pO_coeff'] = max(temp_pO_coeff)
        
        all_reax = df_reax.copy(deep=True)
        all_reax['rO_2_coeff'] = ''
        all_reax['rO_2'] = ''
        all_reax['rO_3_coeff'] = ''
        all_reax['rO_3'] = ''
        all_reax['rR_2_coeff'] = ''
        all_reax['rR_2'] = ''
        
        ### MAIN REACTION
        for r in range(1, max(all_reax['Reaction']+1)): # each reaction number once, 1 to 305
            if len(all_reax[all_reax['Reaction']==r].index.values) == 1: # if nothing to combine, skip
                continue
            else:
                temp = all_reax[all_reax['Reaction']==r].index.values[0] #index of first instance of this reaction which has multiple subreactions
                lst2 = []
                all_reax.loc[temp-0.5] = all_reax.loc[temp] # replicating the row to build on
                all_reax = all_reax.sort_index() # putting the replicated row above the first instance
                for i in all_reax[all_reax['Reaction']==r].index.values[2:]: #all but the first instance in the reactions (since that's copied already)
                    if all_reax.loc[i, 'rO'] != all_reax.loc[temp, 'rO'] and all_reax.loc[i, 'rO'] not in lst2: # if rO is new (and not the same as the first)
                        lst2.append(all_reax.loc[i, 'rO']) # list unique rO besides the first
                        for l in range(0, len(lst2)): # looping through unique rO
                            temp2 = 'rO_'+str(2+int(l)) # adding 0 or 1 to the rO number
                            all_reax.loc[temp-0.5,str(temp2)] = lst2[l] # add the unique rO to rO_2 or 3
                    if all_reax.loc[i, 'rR'] != all_reax.loc[temp, 'rR']: # if rR is new
                            all_reax.loc[temp-0.5,'rR_2'] = all_reax.loc[i, 'rR'] # add it to rR_2

                ##CHANGING COEFFICIENTS
                rO = all_reax.loc[temp-0.5,'rO'] #assigning easy variables
                rO_2 = all_reax.loc[temp-0.5,'rO_2']
                rO_3 = all_reax.loc[temp-0.5,'rO_3']
                rR = all_reax.loc[temp-0.5,'rR']
                rR_2 = all_reax.loc[temp-0.5,'rR_2']
                rO_coeff = all_reax.loc[temp-0.5,'rO_coeff'] #these are empty at the moment
                rO_2_coeff = all_reax.loc[temp-0.5,'rO_2_coeff']
                rO_3_coeff = all_reax.loc[temp-0.5,'rO_3_coeff']
                rR_2_coeff = all_reax.loc[temp-0.5,'rR_2_coeff']
                rR_coeff = all_reax.loc[temp-0.5,'rR_coeff']
                if rO_3 != '' and rO_2 != '': #if DIC is the oxidant
                    all_reax.loc[temp-0.5,'rO_coeff'] = rO_coeff/3 #this works fine
                    all_reax.loc[temp-0.5,'rO_2_coeff'] = rO_coeff/3
                    all_reax.loc[temp-0.5,'rO_3_coeff'] = rO_coeff/3

                    if rR_2 != '':
                        all_reax.loc[temp-0.5,'rR_coeff'] = rR_coeff/2 #this works fine
                        all_reax.loc[temp-0.5,'rR_2_coeff'] = rR_coeff/2

                        all_reax.loc[temp-0.4] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rR - this works fine : line 4
                        all_reax.loc[temp-0.4, 'rR_2_coeff'] = 0
                        all_reax.loc[temp-0.4, 'rR_coeff'] = rR_coeff

                        all_reax.loc[temp-0.3] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rR_2 - this works fine: line 5
                        all_reax.loc[temp-0.3, 'rR_2_coeff'] = rR_coeff
                        all_reax.loc[temp-0.3, 'rR_coeff'] = 0

                        #NEW ROWS WITH TWO DIC AND BOTH rR
                        all_reax.loc[temp-0.25] = all_reax.loc[temp-0.5] #eliminate CO2
                        all_reax.loc[temp-0.25, 'rO_coeff'] = 0
                        all_reax.loc[temp-0.25, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.25, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.24] = all_reax.loc[temp-0.5] #eliminate HCO3-
                        all_reax.loc[temp-0.24, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.24, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.24, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.23] = all_reax.loc[temp-0.5] #eliminate CO3-2
                        all_reax.loc[temp-0.23, 'rO_3_coeff'] = 0                
                        all_reax.loc[temp-0.23, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.23, 'rO_coeff'] = rO_coeff/2

                        all_reax.loc[temp-0.2] = all_reax.loc[temp-0.4] #NEW ROW WITH ONLY rR ELIMINATING CO2: line 6
                        all_reax.loc[temp-0.2, 'rO_coeff'] = 0
                        all_reax.loc[temp-0.2, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.2, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.1] = all_reax.loc[temp-0.4] #NEW ROW WITH ONLY rR ELIMINATING HCO3-: line 7
                        all_reax.loc[temp-0.1, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.1, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.1, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05] = all_reax.loc[temp-0.4] #NEW ROW WITH ONLY rR ELIMINATING CO3-2: line 8
                        all_reax.loc[temp-0.05, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05, 'rO_3_coeff'] = 0

                        all_reax.loc[temp-0.04] = all_reax.loc[temp-0.3] #NEW ROW WITH ONLY rR_2 ELIMINATING CO2: line 9
                        all_reax.loc[temp-0.04, 'rO_coeff'] = 0
                        all_reax.loc[temp-0.04, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.04, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.03] = all_reax.loc[temp-0.3] #NEW ROW WITH ONLY rR_2 ELIMINATING HCO3-: line 10
                        all_reax.loc[temp-0.03, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.03, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.03, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.02] = all_reax.loc[temp-0.3] #NEW ROW WITH ONLY rR_2 ELIMINATING CO3-2: line 11
                        all_reax.loc[temp-0.02, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.02, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.02, 'rO_3_coeff'] = 0

                        all_reax.loc[temp-0.01] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY CO2 and both RR
                        all_reax.loc[temp-0.01, 'rO_coeff'] = rO_coeff
                        all_reax.loc[temp-0.01, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.01, 'rO_3_coeff'] = 0

                        all_reax.loc[temp-0.009] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY HCO3- and both RR
                        all_reax.loc[temp-0.009, 'rO_coeff'] = 0 ###
                        all_reax.loc[temp-0.009, 'rO_2_coeff'] = rO_coeff
                        all_reax.loc[temp-0.009, 'rO_3_coeff'] = 0

                        all_reax.loc[temp-0.008] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY CO2 and both RR
                        all_reax.loc[temp-0.008, 'rO_coeff'] = 0
                        all_reax.loc[temp-0.008, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.008, 'rO_3_coeff'] = rO_coeff

                    if rR_2 == '':
                        all_reax.loc[temp-0.2] = all_reax.loc[temp-0.5] #NEW ROW ELIMINATING CO2
                        all_reax.loc[temp-0.2, 'rO_coeff'] = 0
                        all_reax.loc[temp-0.2, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.2, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.1] = all_reax.loc[temp-0.5] #NEW ROW ELIMINATING HCO3-
                        all_reax.loc[temp-0.1, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.1, 'rO_2_coeff'] = 0
                        all_reax.loc[temp-0.1, 'rO_3_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05] = all_reax.loc[temp-0.5] #NEW ROW ELIMINATING CO3-2
                        all_reax.loc[temp-0.05, 'rO_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05, 'rO_2_coeff'] = rO_coeff/2
                        all_reax.loc[temp-0.05, 'rO_3_coeff'] = 0

                if rO_2 != '' and rO_3 == '': # IF THERE ARE TWO OXIDANT OPTIONS
                    all_reax.loc[temp-0.5,'rO_coeff'] = rO_coeff/2
                    all_reax.loc[temp-0.5,'rO_2_coeff'] = rO_coeff/2

                    if rR_2 != '': #IF THERE ARE TWO REDUCTANT OPTIONS
                        all_reax.loc[temp-0.5,'rR_coeff'] = rR_coeff/2
                        all_reax.loc[temp-0.5,'rR_2_coeff'] = rR_coeff/2
                        
                        if rR_2 == rO_2:
                            continue
                        else:

                            all_reax.loc[temp-0.4] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rR
                            all_reax.loc[temp-0.4, 'rR_2_coeff'] = 0
                            all_reax.loc[temp-0.4, 'rR_coeff'] = rR_coeff

                            all_reax.loc[temp-0.3] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rR_2
                            all_reax.loc[temp-0.3, 'rR_2_coeff'] = rR_coeff
                            all_reax.loc[temp-0.3, 'rR_coeff'] = 0

                            all_reax.loc[temp-0.2] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rO_2
                            all_reax.loc[temp-0.2, 'rO_2_coeff'] = rO_coeff
                            all_reax.loc[temp-0.2, 'rO_coeff'] = 0

                            all_reax.loc[temp-0.1] = all_reax.loc[temp-0.5] #NEW ROW WITH ONLY rO
                            all_reax.loc[temp-0.1, 'rO_2_coeff'] = 0
                            all_reax.loc[temp-0.1, 'rO_coeff'] = rO_coeff

                if rO_2 == '' and rO_3 == '' and rR_2 != '': # IF THERE IS ONLY ONE OXIDANT BUT TWO REDUCTANTS
                    all_reax.loc[temp-0.5,'rR_coeff'] = rR_coeff/2
                    all_reax.loc[temp-0.5,'rR_2_coeff'] = rR_coeff/2

        all_reax = all_reax.sort_index().reset_index(drop=True)
        
        pair_list = []
        for i in range(0, len(all_reax['Temp_Pairs'])):
            pair_list.append([redox_pairs[all_reax.loc[i, 'Temp_Pairs'][0]], redox_pairs[all_reax.loc[i, 'Temp_Pairs'][1]]])
        all_reax['pairs'] = pair_list

        new_elements = []
        for r in range(0, len(all_reax['rO'])):
            for e in elements:
                if e in ['O','H','+','-']:
                    continue
                else:
                    temp1 = int(element_dictionary[all_reax['rO'][r]][e]) #count for the element in the list for rO at index r
                    temp2 = int(element_dictionary[all_reax['pR'][r]][e])
                    temp3 = int(element_dictionary[all_reax['rR'][r]][e])
                    temp4 = int(element_dictionary[all_reax['pO'][r]][e])
                    if temp1 != temp2:
                        if temp1 ==0:
                            all_reax.loc[r, 'red_'+e+'_coeff'] = -temp2
                            elmnt = basis_df.loc[basis_df['name'].str.contains(e)]['formula'].tolist()[0]
                            all_reax.loc[r, 'red_'+e] = elmnt
                            if elmnt not in new_elements:
                                new_elements.append(elmnt)
                        if temp2 == 0:
                            all_reax.loc[r, 'red_'+e+'_coeff'] = temp1
                            elmnt = basis_df.loc[basis_df['name'].str.contains(e)]['formula'].tolist()[0]
                            all_reax.loc[r, 'red_'+e] = elmnt
                            if elmnt not in new_elements:
                                new_elements.append(elmnt)
                    if temp4 != temp3:
                        if temp3 == 0:
                            all_reax.loc[r, 'ox_'+e+'_coeff'] = -temp4
                            elmnt = basis_df.loc[basis_df['name'].str.contains(e)]['formula'].tolist()[0]
                            all_reax.loc[r, 'ox_'+e] = elmnt
                            if elmnt not in new_elements:
                                new_elements.append(elmnt)
                        if temp4 == 0.0:
                            all_reax.loc[r, 'ox_'+e+'_coeff'] = temp3
                            elmnt = basis_df.loc[basis_df['name'].str.contains(e)]['formula'].tolist()[0]
                            all_reax.loc[r, 'ox_'+e] = elmnt
                            if elmnt not in new_elements:
                                new_elements.append(elmnt)


        for i in new_elements:
            if i not in db_names:
                db_names.append(i)
            if i not in formulas:
                formulas.append(i)
            parsed_formula = parse_formula(i)
            element_dictionary[i] = parsed_formula
            for e in elements:
                if element_dictionary[i].get(e, 0) == 0:
                    element_dictionary[i][e] = 0

        reax = all_reax.copy(deep=True)
        reax.drop('Temp_Pairs', axis=1, inplace=True)
        reax.drop('pairs', axis=1, inplace=True)
        reax.reset_index(drop=True, inplace=True)
        for s in ['O', 'H','-','+']:
            for i in range(0,len(reax['rO'])):
                red = 0
                ox=0
                for j in reax.columns.tolist():
                    if '_coeff' in j:
                        if 'rO_' in j or 'pR_' in j or 'red_' in j:
                            if str(reax[j][i]) != 'nan' and str(reax[j][i]) != '':
                                red_temp_coeff = reax[j][i]
                                red_temp = element_dictionary[reax[j.split('_coeff')[0]][i]][s]
                                red -= red_temp_coeff*red_temp

                        if 'rR_' in j or 'pO_' in j or 'ox_' in j:
                            if str(reax[j][i]) != 'nan' and str(reax[j][i]) != '':
                                ox_temp_coeff = reax[j][i]
                                ox_temp = element_dictionary[reax[j.split('_coeff')[0]][i]][s]
                                ox -= ox_temp_coeff*ox_temp      

                reax.loc[i, 'r_'+s] = red
                reax.loc[i, 'o_'+s] = ox

        reax['r_H'] = reax['r_H'] - 2*reax['r_O']
        reax['o_H'] = reax['o_H'] - 2*reax['o_O']
        reax['r_+'] = reax['r_+'] - reax['r_H']
        reax['o_+'] = reax['o_+'] - reax['o_H']
        reax['r_e-'] = reax['r_+'] - reax['r_-'] 
        reax['o_e-'] = reax['o_+'] - reax['o_-'] 
        reax.rename({'r_O': 'r_H2O', 'r_H': 'r_H+', 'o_O': 'o_H2O', 'o_H': 'o_H+'}, axis=1, inplace = True)
        
        ### MULTIPLYING SUB-REACTIONS
        lcm_charge = []
        electrons = []
        for i in range(0, len(reax['rO'])):
        # for i in range(0, 1):
            lcm_charge = np.lcm(round(reax['r_e-'][i]), round(reax['o_e-'][i]))
            electrons.append(str(lcm_charge)+'e')
            r_multiplier = abs(lcm_charge/int(reax['r_e-'][i]))
            o_multiplier = abs(lcm_charge/int(reax['o_e-'][i]))
            for s in list(reax.columns):
                if ('red_' in s and 'coeff' in s) or ('rO_' in s and 'coeff' in s) or ('pR_' in s and 'coeff' in s) or 'r_H2O' in s or 'r_H+' in s:
                    reax.loc[i, s] = reax.loc[i, s]*int(r_multiplier )
                if ('ox_' in s and 'coeff' in s) or ('rR_' in s and 'coeff' in s) or ('pO_' in s and 'coeff' in s) or 'o_H2O' in s or 'o_H+' in s:
                    reax.loc[i, s] = reax.loc[i, s]*int(o_multiplier)
        reax['H+'] = reax['r_H+'] + reax['o_H+']
        reax['protons'] = 'H+'
        reax['H2O'] = reax['r_H2O'] + reax['o_H2O']
        reax['water'] = 'H2O'
        reax.drop(columns = ['r_H2O', 'r_H+', 'r_-', 'r_+', 'r_e-', 'o_H2O', 'o_H+', 'o_-', 'o_+', 'o_e-'],axis = 1, inplace = True)

        count = 0
        for i in db_names:
            db_names[count] = 'start'+i+'end'
            count += 1

        real_reax = reax.replace(formulas,db_names)
        real_reax['rO_coeff'] = real_reax['rO_coeff'].astype('float') 
        real_reax['rR_coeff'] = real_reax['rR_coeff'].astype('float') 
        real_reax['pO_coeff'] = real_reax['pO_coeff'].astype('float') 
        real_reax['pR_coeff'] = real_reax['pR_coeff'].astype('float') 
        
        count = 0
        rxn_count = []
        rxn_number = []

        for i in real_reax['Reaction']:
            if i not in rxn_count:
                rxn_count.append(i)
                rxn_number.append(real_reax['Names'][count]+ '_'+str(count))
            else:
                rxn_number.append(real_reax['Names'][count]+ '_'+str(count)+'_sub')
            count += 1
        real_reax.insert(0, 'Reaction Number', rxn_number)
        real_reax.insert(1, 'electrons', electrons)

        lst3 = [] #list of reaction numbers
        lst4 = [] #list of reactions with issues
        for i in range(0, len(real_reax['Reaction'])):
            if real_reax['Reaction'][i] not in lst3:
                lst3.append(real_reax['Reaction'][i])
        for j in lst3: #looping through reaction numbers
            first_e = real_reax.loc[real_reax['Reaction'] == j]['electrons'].reset_index(drop=True)[0]
            for k in real_reax.loc[real_reax['Reaction'] == j]['electrons'].reset_index(drop=True):
                if k != first_e:
                    if j not in lst4:
                        lst4.append(j)
        for l in lst4:
            print(real_reax.loc[real_reax['Reaction'] ==  l])

        real_reax.drop(labels='Reaction', axis=1, inplace = True)
        real_reax.drop(columns = 'Names', inplace = True)
        pairs = real_reax['redox_pair']
        real_reax.drop(columns = 'redox_pair', inplace = True)
        
        # 2-16-2022 CHANGES START HERE
        for i in range(0, len(real_reax['Reaction Number'])):
            for j in range(2, len(real_reax.columns)):
                if str(real_reax.iloc[i, j]) == 'nan':
                    real_reax.iloc[i, j] = ''
        
        test_df = real_reax.copy(deep=True)
        count = 0
        for i in range(0, len(test_df['Reaction Number'])):
            coefficients = []
            species = []
            for j in range(2, len(test_df.columns)):
                if count % 2 == 0: 
                    if test_df.iloc[i, j] != '':
                        test_df.iloc[i, j] = round(test_df.iloc[i, j], 14)
                    if test_df.iloc[i, j] == 0 or test_df.iloc[i, j] == 0.0:
                        test_df.iloc[i, j] = ''
                        test_df.iloc[i, j+1] = ''
                    coefficient = test_df.iloc[i, j]
                    coefficients.append(test_df.iloc[i, j])
                if count % 2 != 0:
                    if test_df.iloc[i, j] != '':
                        test_df.iloc[i, j] = str(test_df.iloc[i, j]).split('start')[1].split('end')[0]
                    compound = test_df.iloc[i, j]
                    if compound in species:
                        og_location = species.index(compound) 
                        df_location = 3+og_location*2 
                        df_location_coeff = df_location - 1
                        old_coeff = coefficients[og_location] 
                        new_coeff = coefficient + old_coeff 
                        test_df.iloc[i, j] = ''
                        test_df.iloc[i, j-1] = ''
                        df_value = test_df.iloc[i, df_location] #values from first occurence remaining
                        test_df.iloc[i, df_location_coeff] = new_coeff
                    species.append(compound)
                count+=1
        for m in range(0, 7):
            for i in range(0, len(test_df['Reaction Number'])):
                line = []
                for j in range(2, len(test_df.columns)):
                    line.append(test_df.iloc[i, j])
                    if test_df.iloc[i, j] != '' and test_df.iloc[i, j-2] =='':
                        test_df.iloc[i, j-2] = test_df.iloc[i, j]
                        test_df.iloc[i, j] = ''

        file = test_df.to_csv(sep='\t', header=False, index=False, lineterminator='\n')

        file = file.split("\n") #not sure if I should keep this
        
        newlines = []
        for line in file:   
            line = line.strip()
            newlines.append(line)

        self.affinity_energy_reactions_raw = "\n".join(newlines)
        df_rxn = pd.DataFrame([x.split('\t') for x in self.affinity_energy_reactions_raw.split('\n')])
        df_rxn.columns = df_rxn.columns.map(str)
        df_rxn = df_rxn.rename(columns={"0": "reaction_name", "1": "mol_e-_transferred_per_mol_rxn"})
        df_rxn.insert(1, 'redox_pairs', all_reax['pairs'])
        df_rxn = df_rxn.set_index("reaction_name")
        df_rxn = df_rxn[df_rxn['mol_e-_transferred_per_mol_rxn'].notna()]
        self.affinity_energy_reactions_table = df_rxn
        
        prev_was_coeff = False
        n = 1
        for col in self.affinity_energy_reactions_table.iloc[:, 2:].columns:
            if not prev_was_coeff:
                new_col_name = "coeff_"+str(n)
                prev_was_coeff = True
            else:
                new_col_name = "species_"+str(n)
                prev_was_coeff = False
                n += 1
            self.affinity_energy_reactions_table = self.affinity_energy_reactions_table.rename(columns={col: new_col_name})
        
        nonsub_reaction_names = [name for name in self.affinity_energy_reactions_table.index if "_sub" not in name[-4:]]
        if self.verbose != 0:
            print("{} redox reactions have been generated.".format(len(nonsub_reaction_names)))

        
    def show_redox_reactions(self, formatted=True, charge_sign_at_end=False,
                                  hide_subreactions=True, simplify=True,
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
        
        hide_subreactions : bool, default True
            Hide subreactions?
        
        show : bool, default False
            Show the table of reactions? Ignored if not run in a Jupyter
            notebook.
        
        Returns
        ----------
        A pandas dataframe containing balanced redox reactions written in full.
        """
        
        self.affinity_energy_formatted_reactions = copy.copy(self.affinity_energy_reactions_table.iloc[:, 0:1])
        
        df = copy.copy(self.affinity_energy_reactions_table)
        
        if simplify:
            main_rxn_names = df.loc[[ind for ind in df.index if "_sub" not in ind[-4:]]].index
            df = df.iloc[[i-1 for i in range(0, len(df.index)) if "_sub" not in df.index[i][-4:]]]
            
            self.affinity_energy_formatted_reactions = copy.copy(df.iloc[:, 0:1])
            
            reactions = []
            for irow in range(0, df.shape[0]):
                redox_pair = df.loc[df.index[irow], "redox_pairs"]

                oxidant_1 = self.half_cell_reactions.loc[self.half_cell_reactions.index[redox_pair[0]], "Oxidant_1"]
                oxidant_2 = self.half_cell_reactions.loc[self.half_cell_reactions.index[redox_pair[0]], "Oxidant_2"]
                oxidant_3 = self.half_cell_reactions.loc[self.half_cell_reactions.index[redox_pair[0]], "Oxidant_3"]
                reductant_1 = self.half_cell_reactions.loc[self.half_cell_reactions.index[redox_pair[1]], "Reductant_1"]
                reductant_2 = self.half_cell_reactions.loc[self.half_cell_reactions.index[redox_pair[1]], "Reductant_2"]
                
                oxidants = [ox for ox in [oxidant_1, oxidant_2, oxidant_3] if str(ox) != 'nan']
                reductants = [rd for rd in [reductant_1, reductant_2] if str(rd) != 'nan']
                
                if len(oxidants) > 1:
                    oxidant_sigma_needed = True
                else:
                    oxidant_sigma_needed = False
                if len(reductants) > 1:
                    reductant_sigma_needed = True
                else:
                    reductant_sigma_needed = False
                    
                rxn_row = df.iloc[irow, 2:]
                rxn = rxn_row[rxn_row.notna()]
                coeffs = copy.copy(rxn[::2]).tolist()
                names = copy.copy(rxn[1::2]).tolist()
                
                if oxidant_sigma_needed or reductant_sigma_needed:

                    reactant_names = [names[i] for i in range(0, len(names)) if float(coeffs[i]) < 0]
                    for sp in reactant_names:
                        if sp in oxidants and oxidant_sigma_needed:
                            i = names.index(sp)
                            names[i] = u"\u03A3"+sp
                        if sp in reductants and reductant_sigma_needed:
                            if u"\u03A3"+sp not in names:
                                i = names.index(sp)
                                names[i] = u"\u03A3"+sp
                    
                react_grid = pd.DataFrame({"coeff":coeffs, "name":names})
                react_grid["coeff"] = pd.to_numeric(react_grid["coeff"])
                react_grid = react_grid.astype({'coeff': 'float'})

                reactants = " + ".join([(str(-int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else -react_grid["coeff"][i])+" " if -react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
                products = " + ".join([(str(int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else react_grid["coeff"][i])+" " if react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
                if formatted:
                    reactants = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
                    products = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
                reaction = reactants + " = " + products
                reactions.append(reaction)

            self.affinity_energy_formatted_reactions["reaction"] = reactions[1:] + reactions[:1] # because reactions got rotated with respect to reaction names, rotate the other way
            self.affinity_energy_formatted_reactions.index = main_rxn_names
            
        else:
            reactions = []
            for irow in range(0, df.shape[0]):
                redox_pair = df.loc[self.affinity_energy_reactions_table.index[irow], "redox_pairs"]

                oxidant = redox_pair[0]
                reductant = redox_pair[1]

                rxn_row = df.iloc[irow, 2:]
                rxn = rxn_row[rxn_row.notna()]
                coeffs = copy.copy(rxn[::2]).tolist()
                names = copy.copy(rxn[1::2]).tolist()
                react_grid = pd.DataFrame({"coeff":coeffs, "name":names})
                react_grid["coeff"] = pd.to_numeric(react_grid["coeff"])
                react_grid = react_grid.astype({'coeff': 'float'})

                reactants = " + ".join([(str(-int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else -react_grid["coeff"][i])+" " if -react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
                products = " + ".join([(str(int(react_grid["coeff"][i]) if react_grid["coeff"][i].is_integer() else react_grid["coeff"][i])+" " if react_grid["coeff"][i] != 1 else "") + react_grid["name"][i] for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
                if formatted:
                    reactants = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] < 0])
                    products = " + ".join([_format_coeff(react_grid["coeff"][i]) + chemlabel(react_grid["name"][i], charge_sign_at_end=charge_sign_at_end) for i in range(0, len(react_grid["name"])) if react_grid["coeff"][i] > 0])
                reaction = reactants + " = " + products
                reactions.append(reaction)
        
            self.affinity_energy_formatted_reactions["reaction"] = reactions
        

        df_out = copy.copy(self.affinity_energy_formatted_reactions)

        if hide_subreactions and not simplify:
            df_out = self.affinity_energy_formatted_reactions.loc[[ind for ind in self.affinity_energy_formatted_reactions.index if "_sub" not in ind[-4:]]]
        
        if _isnotebook() and show:
            display(HTML(df_out.to_html(escape=False)))
        
        return df_out
        

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
        
        self.df_rejected_species = pd.DataFrame({'database index':[], "name":[], "reason for rejection":[]})
        self.reactions_for_plotting = None # stores formatted reactions for plotting results of affinity and energy supply calculations
        for k in args:
            setattr(self, k, args[k])

    def __getitem__(self, item):
         return getattr(self, item)
    
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
    def __save_figure(fig, save_as, save_format, save_scale, plot_width, plot_height, ppi):
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
            "affinity_kcal" : ("affinity", "kcal/mol"),
            "%" : ("", "%"),
            "Eh_volts" : ("Eh", "volts"),
            "eq/kg.H2O" : ("charge", "eq/kg"),
            "logfO2" : ("", ""),
            "cal/mol e-" : ("affinity", "cal/mol e-"),
            "cal/kg.H2O" : ("energy supply", "cal/kg H2O"),
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
        
        names_length = len(self.report_divs.names)
        
        if col==None and names_length>0:
            return list(self.report_divs.names)
        
        if names_length>0:
            if col in list(self.report_divs.names):
                return list(self.report_divs.rx2(col))
        
        if isinstance(col, str):
            col = [col]
        
        return self.report.iloc[:, self.report.columns.get_level_values(0).isin(set(col))]
    
    
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
        
        save_as, save_format = self.__save_figure(fig, save_as, save_format,
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

    

    def barplot(self, y="pH", title=None, convert_log=True, show_missing=True,
                plot_width=4, plot_height=3, ppi=122, colormap="WORM",
                save_as=None, save_format=None, save_scale=1,
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

        df = self.lookup(["name"]+y).copy()
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

        df['y_variable'] = df['y_variable'].apply(chemlabel)
        
        
        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.affinity_energy_formatted_reactions, pd.DataFrame):
            
            # get formatted reactions to display
            if not isinstance(self.reactions_for_plotting, pd.DataFrame):

                self.reactions_for_plotting = self.show_redox_reactions(formatted=True,
                                                                       charge_sign_at_end=False,
                                                                       show=False, simplify=True)

            y_find = [yi.replace("_energy", "").replace("_affinity", "") for yi in y]
            
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
            
            # customdata for displaying reactions has to be here instead of in update_traces
            fig = px.bar(df, x="name", y="y_value",
                height=plot_height*ppi, width=plot_width*ppi,
                color='y_variable', barmode='group',
                labels={'y_value': ylabel}, template="simple_white",
                color_discrete_map=dict_species_color, custom_data=['formatted_rxns'])
            
            
            
            fig.update_traces(
                hovertemplate = "%{x} <br>"+ylabel+": %{y}<br>%{customdata}")

        else:
            
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

        save_as, save_format = self.__save_figure(fig, save_as, save_format,
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

        
    def scatterplot(self, x="pH", y="Temperature", title=None,
                          plot_width=4, plot_height=3, ppi=122,
                          fill_alpha=0.7, point_size=10,
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

        colors = _get_colors(colormap, len(y), alpha=fill_alpha)
        
        for i, yi in enumerate(y):
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
                y_formatted = chemlabel(y[0])
                if unit != "":
                    ylabel = "{} {} [{}]".format(y_formatted, unit_type, unit)
                else:
                    ylabel = "{} {}".format(y_formatted, unit_type)
        
        if x == 'pH':
            xlabel = 'pH'
        elif x == 'Temperature':
            xlabel = 'Temperature [°C]'
        else:
            x_formatted = chemlabel(x)
            if xunit != "":
                xlabel = "{} {} [{}]".format(x_formatted, xunit_type, xunit)
            else:
                xlabel = "{} {}".format(x_formatted, xunit_type)

        # convert rgba to hex
        colors = [matplotlib.colors.rgb2hex(c) for c in colors]

        # map each species to its color, e.g.,
        # {'CO2': '#000000', 'HCO3-': '#1699d3', 'Other': '#736ca8'}
        dict_species_color = {sp:color for sp,color in zip(y, colors)}
        
        # html format color dict key names
        dict_species_color = {chemlabel(k):v for k,v in dict_species_color.items()}
        
        df = self.lookup(["name", x]+y).copy()
        df.loc[:, "name"] = df.index
        df.columns = df.columns.get_level_values(0)
        df = pd.melt(df, id_vars=["name", x], value_vars=y)
        df = df.rename(columns={"Sample": "y_variable", "value": "y_value"})
        
        if (unit_type == "energy supply" or unit_type == "affinity") and isinstance(self.reactions_for_plotting, pd.DataFrame):
            
            # get formatted reactions to display
            if not isinstance(self.reactions_for_plotting, pd.DataFrame):
                self.reactions_for_plotting = self.show_redox_reactions(formatted=True,
                                                                       charge_sign_at_end=False,
                                                                       show=False, simplify=True)
            
            y_find = [yi.replace("_energy", "").replace("_affinity", "") for yi in y]
            
            
            rxns = self.reactions_for_plotting.loc[y_find, :]["reaction"].tolist()
            rxn_dict = {rxn_name:rxn for rxn_name,rxn in zip(y, rxns)}

            if len(y) == 1:
                ylabel = "{}<br>{} [{}]".format(chemlabel(y_find[0]), unit_type, unit)
            
            df["formatted_rxn"] = df["y_variable"].map(rxn_dict)
        else:
            df["formatted_rxn"] = ""
        
        df['y_variable'] = df['y_variable'].apply(chemlabel)
        
        fig = px.scatter(df, x=x, y="y_value", color="y_variable",
                         hover_data=[x, "y_value", "y_variable", "name", "formatted_rxn"],
                         width=plot_width*ppi, height=plot_height*ppi,
                         labels={x: xlabel,  "y_value": ylabel},
                         category_orders={"species": y},
                         color_discrete_map=dict_species_color,
                         opacity=fill_alpha,
                         custom_data=['name', 'formatted_rxn'],
                         template="simple_white")
        fig.update_traces(marker=dict(size=point_size),
                          hovertemplate = "%{customdata[0]}<br>"+xlabel+": %{x} <br>"+ylabel+": %{y}<br>%{customdata[1]}")
        fig.update_layout(legend_title=None,
                          title={'text':title, 'x':0.5, 'xanchor':'center'},
                          margin={"t": 40},
                          yaxis={'exponentformat':'power'})
        if len(y) == 1:
            fig.update_layout(showlegend=False)
            
        save_as, save_format = self.__save_figure(fig, save_as, save_format,
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
                                     colormap="WORM", sample_label = "sample",
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
        
        save_as, save_format = self.__save_figure(fig, save_as, save_format,
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
            
        
        save_as, save_format = self.__save_figure(fig, save_as, save_format,
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
            
            
    def join_6i_3p(self, filepath_6i):
        path='rxn_6i'
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            shutil.rmtree(path)
            os.makedirs(path)
            
        for sample_name in self.raw_pickup_dict.keys():
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
                
            lines_3p = self.raw_pickup_dict[sample_name]
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
        return sample_data[sample]["mass_transfer"]