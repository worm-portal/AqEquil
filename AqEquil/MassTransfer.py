import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pyCHNOSZ
from chemparse import parse_formula
import copy
import os
import shutil
from math import log10, isnan
import itertools
from operator import itemgetter
import re
import collections
from datetime import datetime
import numbers
import roman
from natsort import natsorted
from .AqSpeciation import Error_Handler, Speciation, AqEquil, chemlabel, check_balance, format_equation

FIXED_SPECIES = ["H2O", "H+", "O2(g)", "water", "Cl-", "e-", "OH-", "O2", "H2O(g)"]


def _format_pseudoelement_name(e):
    
    """
    Format a pseudoelement name
    E.g., "Fejiip" -> "Fe+2"
    """
    
    if len(e) > 2:
        e_split_list = e.split("j")
        
        if e_split_list[1][-1] == "p":
            charge_sign = "+"
        elif e_split_list[1][-1] == "n":
            charge_sign = "-"
        elif e_split_list[1][-1] == "z":
            charge_sign = ""
        else:
            print("ERROR in format_pseudoelement_name(): charge sign is not recognized")
        
        if charge_sign != "":
            charge_magnitude_roman = e_split_list[1][:-1].upper()
            charge_magnitude=roman.fromRoman(charge_magnitude_roman)
        else:
            charge_magnitude = 0
        
        formatted_elem_name = "".join([e_split_list[0], charge_sign,
                                       str(charge_magnitude)])
    else:
        formatted_elem_name = e
        
    return formatted_elem_name

        
def _get_ion_ratio_exponent(num, denom):

    num_formula = parse_formula(num)
    num_plus = num_formula.get("+", 0)
    num_minus = num_formula.get("-", 0)
    num_total_charge = num_plus - num_minus

    denom_formula = parse_formula(denom)
    denom_plus = denom_formula.get("+", 0)
    denom_minus = denom_formula.get("-", 0)
    denom_total_charge = denom_plus - denom_minus

    if num_total_charge == 0 or denom_total_charge == 0:
        return 0
    
    return num_total_charge/denom_total_charge


def __delete_file(file):
    if os.path.exists(file) and os.path.isfile(file):
        os.remove(file)

def __delete_dir(d):
    if os.path.exists(d) and os.path.isdir(d):
        shutil.rmtree(d)

def __move_file(file, destination_dir, silent=False):
    try:
        shutil.move(file, destination_dir+'/'+file)
    except:
        if not silent:
            print("Could not move", file, "to", destination_dir)

        
def react(speciation,
          reaction_setup,
          chain_mt=False,
          delete_generated_folders=False,
          hide_traceback=True,
          data1_override=None,
          format_element_names=True,
          eq36da=None,
          eq36co=None,
          verbose=1,
         ):
    
    """
    Calculate how speciated water reacts with minerals and/or gases.
    
    Parameters
    ----------
    speciation : Speciation object
        The output of a speciation calculation produced by the
        AqEquil.speciate function.
    
    reaction_setup : str or Prepare_Reaction object
        Defines how the reaction is to be set up. There are two ways to set up
        a reaction. The first way is to prepare the reaction with
        `Prepare_Reaction`. The second way is to prepare the first half of an
        EQ6 6i file (the part without the contents of the pickup file) and then
        pass the filename to `reaction_setup`.
    
    chain_mt : bool, default False
        Is the speciation the result of another mass transfer calculation?
        Choosing True will allow mass transfer calculations to be chained
        together.
    
    delete_generated_folders : bool, default False
        Delete the 'rxn_6i', 'rxn_6o', 'rxn_6p', and 'eq6_extra_out' folders
        containing raw EQ6 input, output, and pickup files once the
        reaction calculation is complete?

    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this class?
        When True, error messages handled by this class will be short and to
        the point.

    data1_override : str, optional
        The three letter code of a data1 file used to override the thermodynamic
        database used to speciate the sample(s). This is useful for chaining
        the results of one mass transfer calculation into another while
        simultaneously changing the temperature and pressure regime of the
        new system. See the description for the `chain_mt` parameter for more
        about chaining.

    format_element_names : bool, default True
        Format the chemical symbols of redox-isolated elements generated when
        elements are supplied to the `suppress_redox` parameter of
        `AqEquil.speciate`? Context: when elements are redox-isolated, new
        elements are designated representing different oxidations states in
        order to prevent redox reactions from occuring during calculations. The
        names of these elements appear like "Fejiip", where the first one or two
        letters represent the element (e.g. "Fe"), "j" denotes that the
        chemical symbol is finished, followed by roman numerals representing the
        magnitude of the charge (e.g., "ii" representing 2), finally followed by
        a letter that denotes the sign of the charge (e.g., "p" for positive).
        By setting this parameter to True, the element "Fejiip" will be
        converted into a more readable form, "Fe+2", in all figures and tables
        generated from calculation results. It may be useful to set this
        parameter to False if you are supplying your own data0 file containing
        custom pseudoelements that do not follow the convention used by AqEquil.
    
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored. 
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.

    Returns
    ----------
    An object of class `Speciation` modified with the results of the reaction
    (if a thermodynamic database formatted as a WORM-style CSV was used during
    speciation) or an unmodified speciation object (if a data0 or data1 file
    was used during speciation).
    """

    if not isinstance(eq36da, str):
        eq36da = os.environ.get('EQ36DA')
    if not isinstance(eq36co, str):
        eq36co = os.environ.get('EQ36CO')
    
    prev_wd = os.getcwd()
    ae = AqEquil(load_thermo=False, eq36co=eq36co, eq36da=eq36da, verbose=verbose)
    
    speciation.join_6i_p(reaction_setup, chain_mt)
    __delete_file("data1.dyn")
    
    paths=['rxn_6o', 'rxn_6p', 'eq6_extra_out']
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            shutil.rmtree(path)
            os.makedirs(path)
            
    for sample_name in list(speciation.sample_data.keys()):
        filename_6i = speciation.sample_data[sample_name]["filename"][:-3]+".6i"
        filename_6o = filename_6i[:-3]+".6o"
        filename_6p = filename_6i[:-3]+".6p"
        
        if data1_override != None:
            with open("data1."+data1_override, mode='rb') as data1:
                speciation.data1["all_samples"] = data1.read()
            speciation.thermo.thermo_db_filename = "data1."+data1_override
            
        if "all_samples" not in speciation.data1.keys():
            # each sample has a unique data1. e.g., with dynamic_db
            __delete_file("eq6_extra_out/data1.dyn")
            with open("eq6_extra_out/data1.dyn", 'wb') as f:
                f.write(speciation.data1[speciation.sample_data[sample_name]["filename"][:-3]])

        else:
            # all samples use the same data1.
            with open("eq6_extra_out/data1.dyn", 'wb') as f:
                f.write(speciation.data1["all_samples"])
        
        path_6i="rxn_6i/"
        path_6o="rxn_6o"
        path_6p="rxn_6p"
        path_extra_out="eq6_extra_out"
        
        ae.runeq6(filename_6i,
                  db="dyn",
                  path_6i=path_6i,
                  data1_path=os.getcwd()+"/eq6_extra_out", # ensuring data1 is read from a folder without spaces overcomes the problem where environment variables with spaces do not work properly when assigned to EQ36DA
                  dynamic_db_name=speciation.thermo.thermo_db_filename)

        # get current working dir
        cwd = os.getcwd()
        cwdd = cwd + "/"
        
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
            if ae.verbose > 0:
                print('Error: EQ6 failed to produce output for ' + filename_6i)
        elif len(files_6o) == 1:
            file_6o = files_6o[0]
            file_6ba = files_6ba[0]
            file_6bb = files_6bb[0]
            try:
                file_6t = files_6t[0]
                file_6tx = files_6tx[0]
            except:
                pass

            try:
                # report errors in output
                with open(cwdd+path_6i+"/"+file_6o) as file:
                    lines = [line.rstrip() for line in file]
                EQ6_errors_found = ae._report_3o_6o_errors(lines, sample_name)
            except:
                msg = ("Error: could not open "+path_6i+file_6o+" or there "
                      "is something wrong with EQ3/6 error reporting.")
                ae.err_handler.raise_exception(msg)
            
            try:
                # move output
                shutil.move(cwdd+path_6i+"/"+file_6o, cwdd+path_6o+"/"+filename_6o)
                shutil.move(cwdd+path_6i+"/"+file_6ba, cwdd+path_extra_out+"/"+filename_6ba)
                shutil.move(cwdd+path_6i+"/"+file_6bb, cwdd+path_extra_out+"/"+filename_6bb)
                try:
                    shutil.move(cwdd+path_6i+"/"+file_6t, cwdd+path_extra_out+"/"+filename_6t)
                    shutil.move(cwdd+path_6i+"/"+file_6tx, cwdd+path_extra_out+"/"+filename_6tx)
                except:
                    pass
            except:
                ae.err_handler.raise_exception(("Error: could not move", path_6i+"/"+file_6o, "to", path_6o+"/"+filename_6o))
        
        else:
            ae.err_handler.raise_exception("Error: multiple output files detected for one mass transfer calculation.")
            
        if len(files_6p) == 0:
            if ae.verbose > 0:
                print('Error: EQ6 failed to produce a pickup file for ' + filename_6i)
        elif len(files_6p) == 1:
            file_6p = files_6p[0]
            try:
                # move output
                shutil.move(cwdd+path_6i+"/"+file_6p, cwdd+path_6p+"/"+filename_6p)
            except:
                ae.err_handler.raise_exception(("Error: could not move", path_6i+"/"+file_6p, "to", path_6p+"/"+filename_6p))
        else:
            ae.err_handler.raise_exception("Error: multiple pickup files detected for one mass transfer calculation.")
        
        if not EQ6_errors_found:
            m = Mass_Transfer(thermo=speciation.thermo,
                              six_o_file='rxn_6o/'+filename_6o,
                              format_element_names=format_element_names,
                              tab_name="eq6_extra_out/"+filename_6i[:-2] + 'csv',
                              hide_traceback=hide_traceback)

            speciation.sample_data[sample_name]["mass_transfer"] = m
        else:
            speciation.sample_data[sample_name]["mass_transfer"] = None
            if ae.verbose > 0:
                print(("Mass transfer results for sample '"+sample_name+"' "
                       "could not be saved because the calculation did not "
                       "finish due to error(s).\n"))
        
        # store input, output, and pickup as dicts in speciation object
        try:
            with open(path_6i + "/" + filename_6i, "r") as f:
                lines=f.readlines()
            speciation.raw_6_input_dict[sample_name] = lines
        except:
            pass
        try:
            with open(path_6o + "/" + filename_6o, "r") as f:
                lines=f.readlines()
            speciation.raw_6_output_dict[sample_name] = lines
        except:
            pass
        try:
            with open(path_6p + "/" + filename_6p, "r") as f:
                lines=f.readlines()
            
            # Unlike 3p files, 6p files include headers that need to be removed
            bottom_half = []
            capture = False
            for line in lines:
                if "Start of the bottom half of the input file" in line:
                    capture = True
                if capture:
                    bottom_half.append(line)
            
            speciation.raw_6_pickup_dict[sample_name] = bottom_half
            
        except:
            pass
        
    if delete_generated_folders:
        __delete_dir("eq6_extra_out")
        __delete_dir("rxn_6i")
        __delete_dir("rxn_6p")
        __delete_dir("rxn_6o")
    
    return speciation


def join_mixes(m1, m2):
    """
    Join the results of two mixes, m1 and m2, so that results extend from
    1:0 to 0:1 m1:m2.
    
    Parameters
    ----------
    m1, m2 : objects of class Mass_Transfer
        The mass transfer results of two mixing calculations

    Returns
    ----------
    An object of class `Mass_Transfer` with results joined for plotting.
    """
    
    
    # tables belonging to the Mass Transfer class that are used by its plotting functions
    tabs = ["misc_params", "basis_logact", "dissolved_elements_molal", "dissolved_elements_ppm", 
            "aq_distribution_logact", "aq_distribution_molal", "aq_distribution_logmolal",
            "moles_minerals", "moles_product_minerals"]
    
    for tab in tabs:
        m2_reverse_rows = getattr(m2, tab).iloc[::-1]
        m2_reverse_rows["Xi"] = [1+(1-float(xi)) for xi in m2_reverse_rows["Xi"]]
        m2_reverse_rows = m2_reverse_rows[1:]
        m1_rows = getattr(m1, tab)
        setattr(m1, tab, pd.concat([m1_rows, m2_reverse_rows]).reset_index(drop=True))
        
    # do the same but for all EQ6 output tables
    if m1.tab != None and m2.tab != None:
        for tab in m1.tab.keys():
            if tab in m2.tab.keys():
                m2_reverse_rows = m2.tab[tab].iloc[::-1]
                m2_reverse_rows["Xi"] = [1+(1-float(xi)) for xi in m2_reverse_rows["Xi"]]
                m2_reverse_rows = m2_reverse_rows[1:]
                m1_rows = m1.tab[tab]
                m1.tab[tab] = pd.concat([m1_rows, m2_reverse_rows]).reset_index(drop=True)
    
    # do the same but for mass contribution tables
    for tab in m1.mass_contribution_dict.keys():
        if tab in m2.mass_contribution_dict.keys():
            m2_reverse_rows = m2.mass_contribution_dict[tab].iloc[::-1]
            m2_reverse_rows["Xi"] = [1+(1-float(xi)) for xi in m2_reverse_rows["Xi"]]
            m2_reverse_rows = m2_reverse_rows[1:]
            m1_rows = m1.mass_contribution_dict[tab]
            m1.mass_contribution_dict[tab] = pd.concat([m1_rows, m2_reverse_rows]).reset_index(drop=True)
    
    return m1


class Mass_Transfer:
    """
    Class containing functions to facilitate mass transfer and reaction path
    calculations and visualize results.
    
    Parameters
    ----------
    six_o_file : str
        Path name of the '6o' output file generated by EQ6.

    thermo : an object of class Thermodata
        The subclass containing thermodynamic data in a Speciation object.
    
    thermodata_csv : str
        Path name of the WORM-styled thermodynamic database CSV used in the EQ6
        calculation.
    
    tab_name : str
        Path name of the TAB file generated by EQ6.

    format_element_names : bool, default True
        Format the chemical symbols of redox-isolated elements generated when
        elements are supplied to the `suppress_redox` parameter of
        `AqEquil.speciate`? Context: when elements are redox-isolated, new
        elements are designated representing different oxidations states in
        order to prevent redox reactions from occuring during calculations. The
        names of these elements appear like "Fejiip", where the first one or two
        letters represent the element (e.g. "Fe"), "j" denotes that the
        chemical symbol is finished, followed by roman numerals representing the
        magnitude of the charge (e.g., "ii" representing 2), finally followed by
        a letter that denotes the sign of the charge (e.g., "p" for positive).
        By setting this parameter to True, the element "Fejiip" will be
        converted into a more readable form, "Fe+2", in all figures and tables
        generated from calculation results. It may be useful to set this
        parameter to False if you are supplying your own data0 file containing
        custom pseudoelements that do not follow the convention used by AqEquil.
    
    hide_traceback : bool, default True
        Hide traceback message when encountering errors handled by this class?
        When True, error messages handled by this class will be short and to
        the point.
    
    """
    def __init__(self, six_o_file, thermo=None, tab_name=None,
                 format_element_names=True, hide_traceback=True, verbose=1):
        
        self.err_handler = Error_Handler(clean=hide_traceback)
        
        self.six_o_file = six_o_file
        
        f=open(self.six_o_file, mode='r')
        self.six_o_file_lines=f.readlines()
        f.close()
        
        self.thermo=thermo
        self.tab_name = tab_name
        self.verbose = verbose
        
        self.inactive_species = self.__get_inactive_species()
        
        if isinstance(self.thermo.csv_db, pd.DataFrame):

            # these operations require a WORM-style thermodynamic database CSV
            obigt = pyCHNOSZ.thermo().OBIGT
            pyCHNOSZ.thermo(OBIGT = obigt.loc[ obigt.name.isin(FIXED_SPECIES), : ])
            _ = pyCHNOSZ.add_OBIGT(self.thermo.csv_db, force=True, messages=False)
            
            self.df = copy.deepcopy(self.thermo.csv_db)
            
            try:
                self.tab = self.process_tab(tab_name, self.thermo.csv_db)
            except:
                self.tab = None

            # remove species that do not have Gibbs free energy values
            self.df = self.df[~self.df["G"].isna()]
            
            # remove inactive species from the database to prevent them from
            # showing up as mineral fields or saturation lines
            if len(self.inactive_species) > 0:
                self.df = self.df[~self.df.name.isin(self.inactive_species)]
                
            basis_df = self.df[self.df["tag"] == "basis"]
            aux_df = self.df[self.df["tag"] == "aux"]
            
            # remove basis or aux species with no formula ox state col
            aux_df = aux_df[aux_df["formula_ox"] != ""]
            aux_df = aux_df[~aux_df['formula_ox'].isnull()]
            
            refstate_df = self.df[self.df["tag"] == "refstate"]
            self.basis_df = pd.concat([basis_df])
            self.basis_aux_df = pd.concat([basis_df, aux_df, refstate_df])
            self.df_cr = self.df[self.df["state"] == 'cr']
        else:
            self.df = None
            self.tab = None
            self.basis_df = None
            self.basis_aux_df = None
            self.df_cr = None
        
        self.misc_params = self.__get_misc_params()
        self.dissolved_elements_molal = self.__get_dissolved_elements(unit="molality")
        self.dissolved_elements_ppm = self.__get_dissolved_elements(unit="ppm")
 
        self.aq_distribution_logact = self.__get_aq_distribution(unit="log activity")
        self.aq_distribution_molal = self.__get_aq_distribution(unit="molality")
        self.aq_distribution_logmolal = self.__get_aq_distribution(unit="log molality")
        self.moles_minerals = self.__get_moles_minerals()

        self.saturation_states_pure_solids_log_Q_over_K = self.__get_saturation_states(unit="logQ/K")
        self.saturation_states_pure_solids_affinity = self.__get_saturation_states(unit="affinity")
        
        self.saturation_states_solid_solutions_log_Q_over_K = self.__get_ss_saturation_states(unit="logQ/K")
        self.saturation_states_solid_solutions_affinity = self.__get_ss_saturation_states(unit="affinity")
        self.solid_solution_names = [col for col in list(self.saturation_states_solid_solutions_affinity.columns) if col != "Xi"] # required by a few things below
        
        self.moles_product_minerals_and_solid_solutions = self.__get_moles_product_minerals(include_solid_solutions=True)
        self.moles_product_minerals = self.__get_moles_product_minerals(include_solid_solutions=False) # relies on self.solid_solution_names

        self.solid_solution_x_dict = self.__get_solid_solution_product_phases(unit="x") # relies on self.solid_solution_names
        self.moles_solid_solutions = self.moles_product_minerals_and_solid_solutions.loc[:, ["Xi"] + list(self.solid_solution_x_dict.keys())]
        self.solid_solution_log_x_dict = self.__get_solid_solution_product_phases(unit="log x") # relies on self.solid_solution_names
        self.solid_solution_log_lambda_dict = self.__get_solid_solution_product_phases(unit="log lambda") # relies on self.solid_solution_names
        self.solid_solution_log_lambda_dict = self.__get_solid_solution_product_phases(unit="log activity") # relies on self.solid_solution_names
        
        self.basis_molality = self.__get_basis_species(unit="molality")
        self.basis_ppm = self.__get_basis_species(unit="ppm")
        self.basis_logact = self.__get_basis_species(unit="logact") # this one needs to be after __get_aq_distribution()

        if format_element_names:
            try:
                # format element names in case there are redox-isolated elements
                self.dissolved_elements_molal.columns = ["Xi"]+[_format_pseudoelement_name(e) for e in self.dissolved_elements_molal.columns if e not in ["Xi", "t(days)"]]
                self.dissolved_elements_ppm.columns = ["Xi"]+[_format_pseudoelement_name(e) for e in self.dissolved_elements_ppm.columns if e not in ["Xi", "t(days)"]]
            except:
                self.err_handler.raise_exception("One of the following chemical symbols could not be formatted:" + str(self.dissolved_elements_molal.columns)+". Try setting the parameter format_element_names to False and try again.")
        
        if self.moles_minerals.shape[0] == 0 and self.moles_product_minerals.shape[0] > 0:
            # in the case of special reactants, there is no grand summary table in the 6o file
            # that combine reactant and product minerals (because there is no reactant mineral).
            # In this case, just assume these tables are equal.
            self.moles_minerals = self.moles_product_minerals

        self.mass_contribution_dict = self.__get_mass_contribution()
            
            
    def __get_misc_params(self):
        
        recording = False
        xi_vals = []
        t_vals = []
        p_vals = []
        pH_vals = []
        pmH_vals = []
        logfO2_vals = []
        Eh_vals = []
        pe_vals = []
        aw_vals = []
        for line in self.six_o_file_lines:
            if "                    Xi=" in line:
                recording = True
                splitstrings = line.strip().split(" ")
                xi = [float(v) for v in splitstrings if v not in ['', 'Xi=']][0]
                xi_vals.append(xi)
            elif " Temperature=" in line and recording:
                splitstrings = line.strip().split(" ")
                t = [float(v) for v in splitstrings if v not in ['', 'Temperature=', 'C']][0]
                t_vals.append(t)
            elif " Pressure=" in line and recording:
                splitstrings = line.strip().split(" ")
                p = [float(v) for v in splitstrings if v not in ['', 'Pressure=', 'bars']][0]
                p_vals.append(p)
            elif "NBS pH scale " in line and recording:
                splitstrings = line.strip().split(" ")
                multival = [v for v in splitstrings if v not in ['', 'NBS', 'pH', 'scale']]
                multival = [float(v) if v != "********" else float('nan') for v in multival]
                pH_vals.append(multival[0])
                Eh_vals.append(multival[1])
                pe_vals.append(multival[2])
            elif "Mesmer pH (pmH) scale " in line and recording:
                splitstrings = line.strip().split(" ")
                pmH = [v for v in splitstrings if v not in ['', 'Mesmer', 'pH', '(pmH)', 'scale']]
                pmH = [float(v) if v != "********" else float('nan') for v in pmH][0]
                pmH_vals.append(pmH)
            elif "  Log oxygen fugacity=" in line and recording:
                splitstrings = line.strip().split(" ")
                logfO2 = [float(v) for v in splitstrings if v not in ['', "Log", "oxygen", "fugacity="]][0]
                logfO2_vals.append(logfO2)
            elif "  Activity of water=" in line and recording:
                splitstrings = line.strip().split(" ")
                aw = [float(v) for v in splitstrings if v not in ['', "Activity", "of", "water="]][0]
                aw_vals.append(aw)
            if "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" in line:
                recording = False

        df = pd.DataFrame({
            "Xi":xi_vals,
            "t(days)":[0]*len(xi_vals), # TODO: make this actually reflect days, need an example output file
            "Temp(C)":t_vals,
            "Press(bars)":p_vals,
            "pH":pH_vals,
            "pmH":pmH_vals,
            "log fO2":logfO2_vals,
            "Eh(v)":Eh_vals,
            "pe":pe_vals,
            "aw":aw_vals,
        })
        
        return df
            
            
    def __get_saturation_states(self, unit="logQ/K"):
        
        if unit == "logQ/K":
            col_index = 1
        elif unit == "affinity":
            col_index = 2
        else:
            self.err_handler.raise_exception("Error in get_saturation_states()"
                    ". Unit not recognized.")
        
        return self.mine_6o_table(table_start="--- Saturation States of Pure Solids ---",
                table_stop="--- Saturation States of Pure Liquids ---",
                ignore = ["", "Phases", 'Phase', '---', '-'],
                col_index=col_index)
    
    
    def __get_ss_saturation_states(self, unit="logQ/K"):
        
        if unit == "logQ/K":
            col_index = 1
        elif unit == "affinity":
            col_index = 2
        else:
            self.err_handler.raise_exception("Error in get_ss_saturation_states()"
                    ". Unit not recognized.")
        
        return self.mine_6o_table(table_start="--- Saturation States of Solid Solutions ---",
                table_stop="--- Summary of Saturated and Supersaturated Phases ---",
                ignore = ["", "Phases", 'Phase', '---', '-'],
                col_index=col_index)
    
            
    def __get_dissolved_elements(self, unit="molality"):
        
        if unit == "molality":
            col_index = -1
        elif unit == "mg/kg.sol" or unit == "ppm":
            col_index = 1
        else:
            self.err_handler.raise_exception("Error in get_dissolved_elements()"
                    ". Unit not recognized.")
        
        return self.mine_6o_table(table_start="Elemental Composition of the Aqueous Solution",
                table_stop="Numerical Composition of the Aqueous Solution",
                ignore = ["", "Element", '---', '-'],
                col_index=col_index)
            
        
    def __get_basis_species(self, unit="logact"):
        
        if unit == "molality" or unit == "logact":
            col_index = -1
        elif unit == "mg/kg.sol" or unit == "ppm":
            col_index = 1
        else:
            self.err_handler.raise_exception("Error in get_basis_species()"
                    ". Unit not recognized.")
        
        basis_df = self.mine_6o_table(table_start="Numerical Composition of the Aqueous Solution",
                table_stop="Sensible Composition of the Aqueous Solution",
                ignore = ["", "Some", 'Species', '---', '-'],
                col_index=col_index)
        
        if unit in ["molality", "ppm", 'mg/kg.sol']:
            return basis_df
        elif unit == "logact":
            
            cols_to_harvest = [b for b in basis_df.columns if b not in ["O2(g)", "H2O"]]
            
            return self.aq_distribution_logact[cols_to_harvest]
        
        
        
    def __get_mass_contribution(self):
        
        bases = [b for b in self.basis_logact.columns if b not in ['Xi', 't(days)', 'H2O', "O2(g)", "H+"]]
        
        mass_contribution_dict = {}
        for basis in bases:

            df = self.mine_6o_table(
                            table_start="Species Accounting for 99% or More of Aqueous "+basis,
                            table_stop="Subtotal",
                            ignore = ["", "Species", '---', '-'],
                            col_index=-1)


            df['Xi'] = df['Xi'].astype(str)
            df['Other'] = 100 - df.sum(axis=1, numeric_only=True)
            df['Xi'] = df['Xi'].astype(float)

            df[df["Other"] < 0] = 0

            df["basis"] = basis
            df["factor"] = None
            df["molality"] = None
            
            mass_contribution_dict[basis] = df
        
        return mass_contribution_dict
            
        
            
    def mine_6o_table(self,
                      table_start="--- Distribution of Aqueous Solute Species ---",
                      table_stop="Species with molalities less than",
                      ignore = ["", "Species", '---'],
                      col_index=-1):
        
        """
        Mine a table in a '6o' EQ6 output file and consolidate results into a
        dataframe.
        
        Parameters
        ----------
        table_start : str
            A unique string that indicates the start of the table.

        table_stop : str, optional
            A unique string that indicates the end of the table.
        
        ignore : list of str
            A list of strings representing lines to ignore when mining a table.
            For example, it is prudent to ignore blank lines, or lines
            containing the table column headers.
            A line will be skipped if line.strip().split(' ')[0] matches any
            of the strings in the given list.
        
        col_index : int, default -1
            Integer representing the index of the table column to be mined.
            The default is -1, which is the last column in the table.
            
        Returns
        -------
        df : Pandas dataframe
            A dataframe with rows of the extent of reaction (Xi), and columns
            containing the values of chemical species mined from the file.
        """
        
        lines = self.six_o_file_lines
        
        species = []
        xi_vals=[]
        collect_values = False
        for i in lines:
            if len(i.strip().split(' ')) > 1 and i.strip().split(' ')[0] == "Xi=":
                this_xi_val = float(i.split(' ')[-1])
            if table_stop in i:
                collect_values = False
            if table_start in i:
                xi_vals.append(this_xi_val) # appending here prevents mismatch where there can be more Xi vals than tables to mine
                collect_values = True
            if collect_values:
                if i.strip().split(' ')[0] not in ignore:
                    species.append(i.strip().split(' ')[0])

        species = list(set(species))
        
        species_dict = {"Xi":xi_vals}
        
        if len(species) == 0:
            df = pd.DataFrame(species_dict)
            return df
        
        for s in species:

            vals=[]
            collect_values = False
            for i in lines:
                if collect_values and table_stop in i:
                    # stop collecting
                    collect_values = False
                    if not got_value:
                        vals.append(np.nan)
                if table_start in i:
                    # start collecting
                    collect_values = True
                    got_value = False
                if collect_values:
                    if len(i.strip().split(' ')) > 2 and i.strip().split(' ')[0] == s:
                        split_i = i.strip().split(' ')
                        split_i_clean = [v for v in split_i if v != '']
                        val = split_i_clean[col_index]
                        try:
                            val = float(val)
                        except:
                            # if a value is not a float, e.g., 3.1450-100
                            if "-" in val:
                                val_list = val.split("-")
                                val = "".join([val_list[0], "E-", val_list[1]])
                                val = float(val)
                            elif "+" in val:
                                val_list = val.split("+")
                                val = "".join([val_list[0], "E", val_list[1]])
                                val = float(val)
                            elif "SATD" in val:
                                pass
                            else:
                                self.err_handler.raise_exception(("Error: "
                                    "Encountered a non-numeric value when mining "
                                    "a .6o file: "+val))
                        vals.append(val)
                        got_value = True
            species_dict[s] = vals
        
        df = pd.DataFrame(species_dict)
        
        if 'None' in df.columns:
            df = df.drop(['None'], axis=1)
        
        return df

    
    def __get_inactive_species(self):
        
        lines=self.six_o_file_lines

        table_start = "--- Inactive Species ---"
        table_stop = "The activity coefficients of aqueous species"

        vals=[]
        collect_values=False
        for i in lines:
            if collect_values and table_stop in i:
                # stop collecting
                break
            if table_start in i:
                # start collecting
                collect_values = True
                got_value = False
            if collect_values:
                if table_start not in i:
                    split_i = i.strip().split(' ')
                    split_i_clean = [v for v in split_i if v != '' and v != "None"]
                    if len(split_i_clean) > 0:
                        vals.append(split_i_clean[0])
                        got_value = True
        
        return vals
        
        
    def __get_aq_distribution(self, unit="log activity"):
        
        if unit == "log activity":
            col_index = -1
        elif unit == "molality":
            col_index = 1
        elif unit == "log molality":
            col_index = 2
        
        return self.mine_6o_table(
                      table_start="--- Distribution of Aqueous Solute Species ---",
                      table_stop="Species with molalities less than",
                      ignore = ["", "Species", '---'],
                      col_index=col_index)

    def __get_moles_minerals(self):
        
        return self.mine_6o_table(
                      table_start="Grand Summary of Solid Phases",
                      table_stop="Mass, grams       Volume, cm3",
                      ignore = ["", "Phase/End-member", '---'],
                      col_index=2)

    
    def __get_moles_product_minerals(self, include_solid_solutions=True):
        
        if "Grand Summary of Solid Phases" in "\n".join(self.six_o_file_lines):
            table_stop = "Grand Summary of Solid Phases"
        else:
            table_stop = "Mass, grams       Volume, cm3"
        
        df = self.mine_6o_table(
                      table_start="--- Summary of Solid Phases (ES) ---",
                      table_stop=table_stop,
                      ignore = ["", "Phase/End-member", '---'],
                      col_index=2)
        
        if not include_solid_solutions:
            mineral_names_to_keep = []
            minerals_in_solid_solutions = []
            
            recording = True
            for col in df.columns:
                if col not in self.solid_solution_names:
                    mineral_names_to_keep.append(col)
            df = df[mineral_names_to_keep]
            
        return df
        
        
    def __get_solid_solution_product_phases(self, unit="x"):
        
        if unit == "x":
            col_index = 1
        elif unit == "log x":
            col_index = 2
        elif unit == "log lambda":
            col_index = 3
        elif unit == "log activity":
            col_index = 4
        
        ss_dict = {}
        for ss in self.solid_solution_names:
            df = self.mine_6o_table(
                         table_start="--- "+ss+" ---",
                         table_stop="Mineral",
                         ignore = ["", "Component", "Ideal", '---'],
                         col_index=col_index)
            
            if not df.empty:
                ss_dict[ss] = df
        
        
        return ss_dict
        
        
    def print_tabs(self):
        """
        Print the names of tables contained in a tab file processed by the
        the Mass_Transfer class.
        """
        if self.tab != None:
            [print(key) for key in self.tab.keys()]
        else:
            print("A processed TAB file is not associated with this sample.")

        
    @staticmethod
    def __is_all_same_value(s):
        a = s.to_numpy()
        return (a[0] == a).all()
        
        
    def plot_reaction_paths(self,
                            xyb=None,
                            path_margin=0.25,
                            flip_xy=False,
                            show_annotation=True,
                            annotation_coords=[0,0],
                            show_nonparticipating_mineral_lines=False,
                            minerals_to_show=[],
                            calculate_projected_points=True,
                            path_line_type = "markers+lines",
                            path_line_color = "red",
                            path_point_fill_color = "red",
                            path_point_line_color = "red",
                            projected_point_fill_color = "white",
                            projected_point_line_color = "red",
                            h_line_color="black",
                            v_line_color="black",
                            d_line_color="black",
                            res=300,
                            plot_width=4,
                            plot_height=3,
                            ppi=122,
                            borders=0,
                            save_as=None,
                            save_format=None,
                            save_scale=1,
                            colormap="bw"):
        
        """
        Create interactive plots of reaction paths in geochemical variable
        space.
        
        Parameters
        ----------
        xyb : list of three str, default None
            By default, this function will plot reaction paths in all possible
            dimensions.
            
            Optionally, if you want to produce only a specific plot,
            you can provide a list containing the basis species to be used for
            the x-axis and y-axis, followed by the basis species used for
            balance. For example, ["Fe+2", "Fe+3", "Mg+2"] will have the log
            activity of Fe+2 on the x-axis, the log activity of Fe+3 on the
            y-axis, and will be balanced on Mg+2.

        path_margin : float, default 0.25
            Controls the spacing between the reaction path and the plot axes.
            Increasing this value increases the spacing.
        
        flip_xy : bool, default False
            Transpose the plot so the x and y variables switch axes?

        show_annotation : bool, default True
            Show annotation in the bottom left of the figure? The annotation
            includes the species used for balance, the temperature, and
            the pressure.

        annotation_coords : list, default [0,0]
            List of two numeric values representing the X and Y coordinates of
            the annotation, where 0,0 is the bottom left, 0.5,0 is the bottom
            center, 1,0 is the bottom right, 1,1 is the top right, and so on.
            The annotation includes the species used for balance, the
            temperature and the pressure.

        show_nonparticipating_mineral_lines : bool, default False
            Depict lines for minerals even if those minerals do not participate
            in the reaction? This does not affect the mineral field of the
            diagram, only the mineral planes denoted as lines.
        
        path_line_type : str, default "markers+lines"
            Reaction path line type. Can be either "markers+lines", "lines", or
            "markers".
        
        path_line_color, str, default "black"
            Color of reaction path line.
        
        path_point_fill_color : str, default "black"
            Fill color of non-projected points along the reaction path. The
            fill color of projected points is handled by
            `projected_point_fill_color`.
        
        path_point_line_color : str, default "black"
            Color of the outlines of non-projected points along the reaction
            path. The outline color of projected points is handled by
            `projected_point_line_color`.
        
        projected_point_fill_color : str, default "white"
            Fill color of projected points along the reaction path. The
            fill color of non-projected points is handled by
            `path_point_fill_color`.
        
        projected_point_line_color : str, default "black"
            Color of the outlines of projected points along the reaction path.
            The outline color of non-projected points is handled by
            `path_point_line_color`.
            
        h_line_color : str, default "black"
            Color of horizontal lines representing minerals.
        
        v_line_color : str, default "black"
            Color of vertical lines representing minerals.
        
        d_line_color : str, default "black"
            Color of diagonal lines representing minerals.
        
        res : int, default 300
            Resolution, or number of calculations along each axis, for mineral
            stability fields. A lower number will be faster but will make appear
            boundaries blockier. A higher number takes longer to calculate, but
            will result in smoother boundaries.
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches.

        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
        
        borders : float, default 0
            Thickness of black lines forming boundaries between mineral
            stability regions. No lines appear if equal to 0.
        
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
        
        colormap : str, default "bw"
            Name of the colormap to color the scatterpoints. Accepts "bw"
            or names of matplotlib colormaps. If set to "bw", the plot will be
            set to black and white, except for the reaction path line itself.
            The colors of the reaction path line and its points are controlled
            by `path_line_color`, `path_point_fill_color`,
            `path_point_line_color`, `projected_point_fill_color`,
            and `projected_point_line_color`.
        
        Returns
        -------
        fig_list : a list of Plotly figure objects
            A list of interactive Plotly figures. If xyb equals None
            (the default), then the list will contain figures representing all
            combinations of geochemical variables. Optionally, if xyb is
            specified, fig_list will only contain the single figure of interest.
        """

        error_messages = []
        
        # check that there is only one temperature and pressure
        if not self.__is_all_same_value(self.misc_params["Temp(C)"]):
            error_messages.append("Reaction paths cannot be plotted when temperature changes with reaction progress.")
        if not self.__is_all_same_value(self.misc_params["Press(bars)"]):
            error_messages.append("Reaction paths cannot be plotted when pressure changes with reaction progress.")
        
        if isinstance(xyb, list):
            if len(xyb) != 3:
                error_messages.append(("Error in xyb={}".format(xyb)+". "
                        "The xyb parameter must either be None or a list of "
                        "three basis species to serve as x, y, and balance variables."))
        
        if not isinstance(self.df, pd.DataFrame):
            error_messages.append(("The plot_reaction_paths() function requires "
                    "a thermodynamic database in a WORM-style CSV format, e.g., "
                    "'wrm_data.csv'. You may be getting this message because "
                    "a data0 or data1 file was used."))
        
        if len(error_messages)>0:
            self.err_handler.raise_exception("\n".join(error_messages))
        
        self.T = float(self.misc_params["Temp(C)"][0])
        self.P = float(self.misc_params["Press(bars)"][0])
        self.path_margin = path_margin
        
        minerals_formed = [m for m in self.moles_minerals.columns if m != "Xi"]
        
        all_elements_of_interest = []
        for mineral in minerals_formed:
            
            if mineral not in list(self.df["name"]):
                if self.verbose > 0:
                    print("The mineral", mineral, "cannot be represented in a",
                           "reaction path diagram, likely because it is missing",
                           "a Gibbs free energy value in the thermodynamic",
                           "database. Continuing anyway, but be aware that",
                           "this mineral will not be represented in diagrams.")
                continue
            
            all_elements_of_interest += self.__get_elem_ox_of_interest_in_minerals(mineral)
        all_elements_of_interest_pre = list(set(all_elements_of_interest))
        
        # filter out elements like Fe+0, which has no aqueous species representative
        # for which to create an axis.
        bad_elem = []
        for elem in all_elements_of_interest_pre:
            if len(list(self.thermo.csv_db.loc[self.thermo.csv_db['name'] == self.__get_basis_from_elem(elem), 'state'])) == 0:
                bad_elem.append(elem)
            elif list(self.thermo.csv_db.loc[self.thermo.csv_db['name'] == self.__get_basis_from_elem(elem), 'state'])[0] != 'aq':
                bad_elem.append(elem)
                
        all_elements_of_interest = [elem for elem in all_elements_of_interest_pre if elem not in bad_elem]
        
        all_elements_of_interest = sorted(all_elements_of_interest)
        
        self.all_elements_of_interest = all_elements_of_interest
        
        fig_list = []
        
        # if there are only 2 elements of interest, these become the axes, and there is no
        # need to fuss with real vs projected points.
        if len(all_elements_of_interest) == 2:
            fig, _ , _ = self.__plot_reaction_path_main(
                                                triad = all_elements_of_interest,
                                                T=self.T, P=self.P,
                                                path_margin=self.path_margin,
                                                flip_xy=flip_xy,
                                                show_annotation=show_annotation,
                                                annotation_coords=annotation_coords,
                                                show_nonparticipating_mineral_lines=show_nonparticipating_mineral_lines,
                                                minerals_to_show=minerals_to_show,
                                                path_line_type=path_line_type,
                                                path_line_color=path_line_color,
                                                path_point_fill_color=path_point_fill_color,
                                                path_point_line_color=path_point_line_color,
                                                projected_point_fill_color=projected_point_fill_color,
                                                projected_point_line_color=projected_point_line_color,
                                                h_line_color=h_line_color,
                                                v_line_color=v_line_color,
                                                d_line_color=d_line_color,
                                                res=res,
                                                plot_width=plot_width,
                                                plot_height=plot_height,
                                                ppi=ppi,
                                                colormap=colormap,
                                                borders=borders,
                                                projected_points=["real"]*self.moles_product_minerals.shape[0],
                                                first_pass=False)
        
            fig_list = [fig]
        
        elif len(all_elements_of_interest) >= 3:

            # get a list of elem pairs for plotting
            alist = self.all_elements_of_interest
            element_plot_pairs = []
            for result in itertools.combinations(alist, 2):
                element_plot_pairs.append(list(result))

            element_plot_triad = []
            for pair in element_plot_pairs:
                elem_not_in_pair = [e for e in self.all_elements_of_interest if e not in pair]
                for e in elem_not_in_pair:
                    triad_to_append = pair + [e]
                    element_plot_triad.append(triad_to_append)

            if colormap == "bw":
                if borders == 0:
                    borders = 1
                colormap = "none"
                h_line_color = "black"
                v_line_color = "black"
                d_line_color = "black"

            fig_list = []
            pred_minerals_from_fields_list = []
            pred_minerals_from_lines_list = []
            

            if isinstance(xyb, list):
                try:
                    xyb_element_plot_triad = [[self.__get_elem_ox_of_interest_in_minerals(v)[0] for v in xyb]]
                except:
                    err = ("Plot axes cannot accomodate desired variables. "
                           "Available variables include {}".format([self.__get_basis_from_elem(elem) for elem in alist]))
                    self.err_handler.raise_exception(err)
                xyb_i = None
                # get index of triad that matches xyb:
                for i,triad in enumerate(element_plot_triad):
                    if set(xyb_element_plot_triad[0][0:2]) == set(triad[0:2]) and xyb_element_plot_triad[0][2] == triad[2]:
                        xyb_i = i
                if xyb_i == None:
                    err = ("Plot axes cannot accomodate desired variables. "
                           "Available variables include {}".format([self.__get_basis_from_elem(elem) for elem in alist]))
                    self.err_handler.raise_exception(err)
            
            if not calculate_projected_points or path_line_type=="lines":
                projected_points = ["real"]*self.moles_product_minerals.shape[0]
                fig_list_projected_points = [projected_points]*len(element_plot_triad)
                
            else:
                
                if len(element_plot_triad) > 20:
                    if self.verbose > 0:
                        print("Warning! There are {}".format(len(element_plot_triad)),
                              "different combinations of variables that must be considered",
                              "in order to plot markers.")
                        print("This might take a very long time or may not finish calculating at all.")
                        print("We recommend setting calculate_projected_points=False in",
                              "plot_reaction_paths() and then restarting the",
                              "calculation to avoid lengthy calculation times.")
                
                for triad in element_plot_triad:
                    
                    # do a quick first pass at making figures to see which points are projections.
                    fig, pred_minerals_from_fields, pred_minerals_from_lines = self.__plot_reaction_path_main(
                                                        triad, T=self.T, P=self.P,
                                                        show_nonparticipating_mineral_lines=False, # no need for this in first pass
                                                        minerals_to_show=[], # no need for this in first pass
                                                        path_margin=self.path_margin,
                                                        flip_xy=flip_xy,
                                                        first_pass=True, # flag for skipping certain calculations/plotting
                                                        res=1) # low res first pass

                    if pred_minerals_from_fields == None:
                        pred_minerals_from_fields=[None]*self.moles_product_minerals.shape[0]
                    
                    fig_list.append(fig)
                    pred_minerals_from_fields_list.append(pred_minerals_from_fields)
                    pred_minerals_from_lines_list.append(pred_minerals_from_lines)

                # determine which line segments in the reaction path are projections
                # and which are actually in the plane of the diagram.
                # This is the "first pass"
                fig_list_projected_points = []
                for i,triad in enumerate(element_plot_triad):

                    projected_points = ["projection"]*self.moles_product_minerals.shape[0]
                    for irow in range(0, self.moles_product_minerals.shape[0]):

                        # get names of minerals formed at this xi
                        xirow = list(self.moles_product_minerals.iloc[irow])
                        formed_minerals = [self.moles_product_minerals.columns[1:][ii] for
                                           ii,mineral in enumerate(xirow[1:]) if
                                           mineral>0]
                        
                        available_pred_minerals_from_fields = [l[irow] for l in pred_minerals_from_fields_list]

                        for mineral in formed_minerals:

                            if mineral == pred_minerals_from_fields_list[i][irow]:
                                # if this mineral is in pred_minerals_from_fields_list,
                                # then it is NOT a projection.
                                projected_points[irow] = "real"
                            if mineral in available_pred_minerals_from_fields and mineral in pred_minerals_from_lines_list[i]:
                                # if this mineral is in the irowth location of any of the
                                # lists in pred_minerals_from_lines_list, it is NOT a
                                # projection.
                                projected_points[irow] = "real"

                    fig_list_projected_points.append(projected_points)
                
            if isinstance(xyb, list):
                # if xyb is defined, make element_plot_triad have a length of 1
                element_plot_triad = xyb_element_plot_triad
                # give the list of projected points lists a length of 1
                fig_list_projected_points = [fig_list_projected_points[xyb_i]]
                
            fig_list = []
            for i,triad in enumerate(element_plot_triad):
                # re-run figure generation, passing in a list of which points are projected.
                fig, _ , _ = self.__plot_reaction_path_main(
                                    triad, T=self.T, P=self.P,
                                    path_margin=self.path_margin,
                                    flip_xy=flip_xy,
                                    show_annotation=show_annotation,
                                    annotation_coords=annotation_coords,
                                    show_nonparticipating_mineral_lines=show_nonparticipating_mineral_lines,
                                    minerals_to_show=minerals_to_show,
                                    path_line_type=path_line_type,
                                    path_line_color=path_line_color,
                                    path_point_fill_color=path_point_fill_color,
                                    path_point_line_color=path_point_line_color,
                                    projected_point_fill_color=projected_point_fill_color,
                                    projected_point_line_color=projected_point_line_color,
                                    h_line_color=h_line_color,
                                    v_line_color=v_line_color,
                                    d_line_color=d_line_color,
                                    res=res,
                                    plot_width=plot_width,
                                    plot_height=plot_height,
                                    ppi=ppi,
                                    colormap=colormap,
                                    borders=borders,
                                    projected_points=fig_list_projected_points[i],
                                    first_pass=False)

                fig_list.append(fig)
        
        if not fig_list and self.verbose > 0:
            print("Warning: a reaction path plot could not be generated for this system.")
        
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            for i,fig in enumerate(fig_list):
                
                if isinstance(xyb, list):
                    name_append = ""
                else:
                    name_append = "_{}".format(i+1)
                _, _ = dummy_sp._save_figure(fig,
                        save_as+name_append, save_format, save_scale,
                        plot_width, plot_height, ppi)
        
        return fig_list
        

    
    
    @staticmethod
    def process_tab(tab_name, thermodata_csv):
        
        """
        Process a TAB file (from EQ6) into a dictionary of Pandas dataframes.
        
        Parameters
        ----------
        tab_name : str
            Path name of the TAB file generated by EQ6.
        
        thermodata_csv : str or Pandas dataframe
            Path name of the WORM-styled thermodynamic database CSV used in the
            EQ6 calculation. Alternately, the thermodynamic database itself as a
            Pandas dataframe itself.

        Returns
        -------
        tab : a dict of Pandas dataframes
            A dictionary of dataframes representing tables mined from the TAB
            file.
        """
        
        if isinstance(thermodata_csv, str):
            thermo_db = pd.read_csv(thermodata_csv)
        else:
            thermo_db = thermodata_csv
        
        
        thermo_db_names = list(thermo_db["name"])

        with open(tab_name, "r") as tabfile:
            tab_lines = tabfile.readlines()

        tab = {}

        tables = ["B1", "B2", "C1", "C2", "C3", "C4",
                  "D1", "D2", "D3", "D4", "E1", "E2",
                  "E3", "J", "K", "P", "Q", "T", "W"]

        record_lines = False
        get_header = False
        recorded_lines = []
        for line in tab_lines:
            split_line = line.split(",")

            if True in [s in ["Table " + t for t in tables] for s in split_line]:
                get_header = True
                table_name = " ".join(split_line[1:-1])
                continue

            if get_header:
                header = split_line[1:-1]

                # handle instances where the tab file creates extra columns for things like "albite,low"
                if table_name in ["Table P Moles of product minerals",
                                  "Table Q Saturation indices of potential product phases",
                                  "Table J Moles of reactants destroyed/created",
                                  "Table K Affinities of reactants (kcal)"]:
                    
                    if table_name in ["Table P Moles of product minerals",
                                      "Table J Moles of reactants destroyed/created"]:
                        new_header = ["Xi", "t(days)"] # table P and J
                    elif table_name == "Table Q Saturation indices of potential product phases":
                        new_header = ["Xi", "t(days)", "H2O", "Gas"] # table Q
                    elif table_name == "Table K Affinities of reactants (kcal)":
                        new_header = ["Xi", "t(days)", "Total"]

                    for i,h in enumerate(header):

                        if h != "Xi" and h != "t(days)" and h != "H2O" and h != "Gas" and h != "Total":
                            if h not in thermo_db_names:
                                if header[i-1]+","+h in thermo_db_names:
                                    new_header = new_header[:-1]
                                    new_header.append(header[i-1]+","+h)
                                else:
                                    new_header.append(h)
                            else:
                                new_header.append(h)

                    header = new_header


                record_lines = True
                get_header = False
                continue

            if "EndTable:" in split_line:
                record_lines = False
                if len(recorded_lines) > 0:
                    df = pd.DataFrame(recorded_lines)
                    recorded_lines = []
                    df.columns = header
                    tab[table_name] = df

            if record_lines:
                recorded_lines.append(split_line[1:-1])

        return tab

    
    def __get_mineral_elem_ox(self, mineral):
        split_list = list(self.df[self.df["name"]==mineral]["formula_ox"])[0].split()
        split_list_clean = [s.replace(" ", "") for s in split_list]
        try:
            elem_ox_list = [re.findall(r"^(?:\d+|)([A-Z].*$)", s)[0] for s in split_list_clean]
        except:
            elem_ox_list = []
        return elem_ox_list

    
    def __get_elem_ox_of_interest_in_minerals(self, mineral_name):
        mineral_elements = self.__get_mineral_elem_ox(mineral_name)
        mineral_elements_of_interest = [e for e in mineral_elements if e not in ["H+", "O-2"]]
        return mineral_elements_of_interest


    def __get_mineral_elem_ox_dict(self, mineral):
        split_list = list(self.df[self.df["name"]==mineral]["formula_ox"])[0].split()
        split_list_clean = [s.replace(" ", "") for s in split_list]
        elem_ox_names = self.__get_mineral_elem_ox(mineral)

        elem_ox_list = []
        for s in split_list:
            coeff = re.findall(r"(\d+)[A-Z]", s)
            if len(coeff) == 0:
                coeff = 1
            else:
                coeff = float(coeff[0])
            elem_ox_list.append(coeff)

        return {key:val for key,val in zip(elem_ox_names, elem_ox_list)}

    
    def __get_mineral_elem_ox_dict_interest(self, mineral):
        mineral_dict = self.__get_mineral_elem_ox_dict(mineral)
        return {key:value for key,value in zip(mineral_dict.keys(), mineral_dict.values()) if key not in ["H+", "O-2"]}


    def __get_basis_from_elem(self, elem):

        basis_species_x = None

        for s in list(self.basis_df["name"]):
            if elem in self.__get_elem_ox_of_interest_in_minerals(s):
                basis_species_x = s

        if basis_species_x == None:
            for s in list(self.basis_aux_df["name"]):
                if elem in self.__get_elem_ox_of_interest_in_minerals(s):
                    basis_species_x = s

        return basis_species_x

    
    def __get_reaction_path(self, plot_basis_x, plot_basis_y, div_var_name):

        xi_vals=self.aq_distribution_logact["Xi"]
        proton_vals=self.aq_distribution_logact["H+"]
        x_vals=self.aq_distribution_logact[plot_basis_x]
        y_vals=self.aq_distribution_logact[plot_basis_y]

        assert len(proton_vals) == len(x_vals), f"number of proton values ({proton_vals}) should equal number of x values ({x_vals})"
        assert len(proton_vals) == len(y_vals), f"number of proton values ({proton_vals}) should equal number of y values ({y_vals})"

        x_vals = [log10((10**float(x))/(10**float(d))**_get_ion_ratio_exponent(plot_basis_x, "H+")) for x,d in zip(x_vals,proton_vals)]
        y_vals = [log10((10**float(y))/(10**float(d))**_get_ion_ratio_exponent(plot_basis_y, "H+")) for y,d in zip(y_vals,proton_vals)]
        
        return xi_vals, x_vals, y_vals

    
    def __get_plot_range(self, x_vals, y_vals):
        min_x_val = min(x_vals)
        min_y_val = min(y_vals)
        max_x_val = max(x_vals)
        max_y_val = max(y_vals)

        path_x_range = max_x_val - min_x_val
        path_y_range = max_y_val - min_y_val

        if len(list(set(x_vals))) == 1:
            # x values form a vertical line
            plot_x_range = [min_x_val-self.path_margin*(path_x_range+1), max_x_val+self.path_margin*(path_x_range+1)]
        else:
            plot_x_range = [min_x_val-self.path_margin*path_x_range, max_x_val+self.path_margin*path_x_range]
        
        if len(list(set(y_vals))) == 1:
            # y values form a horizontal line
            plot_y_range = [min_y_val-self.path_margin*(path_y_range+1), max_y_val+self.path_margin*(path_y_range+1)]
        else:
            plot_y_range = [min_y_val-self.path_margin*path_y_range, max_y_val+self.path_margin*path_y_range]

        return plot_x_range, plot_y_range

    
    @staticmethod
    def __get_xy_labs(plot_basis_x, plot_basis_y):
        try:
            xlab = pyCHNOSZ.ratlab(plot_basis_x)
        except:
            xlab = "log a"+chemlabel(plot_basis_x)
        try:
            ylab = pyCHNOSZ.ratlab(plot_basis_y)
        except:
            ylab = "log a"+chemlabel(plot_basis_y)

        return xlab, ylab


    def __plot_reaction_path_background(self, plot_basis_x, plot_basis_y,
                                        div_var_name, x_vals, y_vals,
                                        colormap="viridis", borders=0,
                                        field_minerals_exist=True, path_margin=0.25,
                                        plot_width=4, plot_height=3, ppi=122, res=300,
                                        annotation=None, annotation_coords=[0, 0],
                                        messages=False):
        
        plot_x_range, plot_y_range = self.__get_plot_range(x_vals, y_vals)

        xlab,ylab = self.__get_xy_labs(plot_basis_x, plot_basis_y)
        
        args = {plot_basis_x:plot_x_range+[res],
                plot_basis_y:plot_y_range+[res],
                "T":self.T, "P":self.P, "messages":messages}
        
        if field_minerals_exist:
            
            # check each value of Xi to see which mineral is most predominant
            pred_minerals_from_fields = []
            for i,val in enumerate(x_vals):
                
                args_temp = {plot_basis_x:[x_vals[i], x_vals[i], 1],
                             plot_basis_y:[y_vals[i], y_vals[i], 1],
                             "T":self.T, "P":self.P, "messages":messages}
                
                a = pyCHNOSZ.affinity(**args_temp)
                e = pyCHNOSZ.equilibrate(a, balance=self.__get_basis_from_elem(div_var_name), messages=messages)
                table = pyCHNOSZ.diagram(e, balance=self.__get_basis_from_elem(div_var_name), interactive=True, fig_out=False, plot_it=False, messages=messages)
                
                pred_minerals_from_fields.append(table["prednames"][0])
            
            a = pyCHNOSZ.affinity(**args)
            e = pyCHNOSZ.equilibrate(a, balance=self.__get_basis_from_elem(div_var_name), messages=messages)
            
            table,fig = pyCHNOSZ.diagram_interactive(e,
                                colormap=colormap, borders=borders,
                                balance=self.__get_basis_from_elem(div_var_name),
                                width=plot_width*ppi, height=plot_height*ppi,
                                xlab=xlab, ylab=ylab, annotation=annotation,
                                annotation_coords=annotation_coords,
                                plot_it=False, messages=messages)
            
            return table, fig, pred_minerals_from_fields
        
        else:
            # empty plot upon
            fig = go.Figure(go.Scatter(x=pd.Series(dtype=object),
                                       y=pd.Series(dtype=object),
                                       mode="markers",
                                       ),
                           layout_xaxis_range=plot_x_range,
                           layout_yaxis_range=plot_y_range,
                           )

            fig.add_annotation(x=annotation_coords[0],
                               y=annotation_coords[1],
                               xref="paper",
                               yref="paper",
                               align='left',
                               text=annotation,
                               bgcolor="rgba(255, 255, 255, 0.5)",
                               showarrow=False)

            fig.update_layout(
                width=plot_width*ppi, height=plot_height*ppi,
                xaxis={"title": xlab},
                yaxis={"title": ylab},
                template="simple_white",
            )

            fig.update_yaxes(autorange=True)

            return None,fig,None # a table, diagram without regions, and a list of predicted minerals at each Xi


    @staticmethod
    def __calc_dissrxn_logK(mineral, T, P):

        logK = pyCHNOSZ.subcrt([mineral], coeff=[-1], property='logK', T=T, P=P,
                      show=False, messages=False)["out"]["logK"].item()
        
        return logK

    
    def __add_reaction_path_to_plot(self, x_vals, y_vals, xi_vals, fig,
                                    basis_species_x, basis_species_y,
                                    path_margin=0.25, projected_points=[],
                                    path_line_type = "markers+lines",
                                    path_line_color = "black",
                                    path_point_fill_color = "black",
                                    path_point_line_color = "black",
                                    projected_point_fill_color = "white",
                                    projected_point_line_color = "black"):

        min_x_val = min(x_vals)
        min_y_val = min(y_vals)
        max_x_val = max(x_vals)
        max_y_val = max(y_vals)

        path_x_range = max_x_val - min_x_val
        path_y_range = max_y_val - min_y_val

        xlab,ylab = self.__get_xy_labs(basis_species_x, basis_species_y)
        
        with np.errstate(divide='ignore'):
            log_xi_vals = [round(np.log10(val), 4) for val in xi_vals]
            log_xi_vals = ["N/A" if np.isinf(val) else val for val in log_xi_vals]
        
        
        if path_line_type in ["markers+lines", "lines"]:
            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    line=go.scatter.Line(color=path_line_color),
                    mode='lines',
                    showlegend=True,
                    name='reaction path',
                    customdata = np.stack((xi_vals, log_xi_vals), axis=-1),
                    hovertemplate = 'Xi: %{customdata[0]}<br>log Xi: %{customdata[1]}<br>'+xlab+': %{x}<br>'+ylab+': %{y}<extra></extra>',
                    legendgroup='reaction path',
                )
            )
        
        
        if len(projected_points) > 0 and path_line_type in ["markers+lines", "markers"]:
            
            x_vals_real = []
            x_vals_projected = []
            y_vals_real = []
            y_vals_projected = []
            xi_vals_real = []
            xi_vals_projected = []
            for i,p in enumerate(projected_points):
                if p == "real":
                    x_vals_real.append(x_vals[i])
                    y_vals_real.append(y_vals[i])
                    xi_vals_real.append(xi_vals[i])
                else:
                    x_vals_projected.append(x_vals[i])
                    y_vals_projected.append(y_vals[i])
                    xi_vals_projected.append(xi_vals[i])
                    
            np.seterr(divide='ignore') # todo: reset np warnings
            log_xi_vals_real = [round(np.log10(val), 4) for val in xi_vals_real]
            log_xi_vals_real = ["N/A" if np.isinf(val) else val for val in log_xi_vals_real]
            log_xi_vals_projected = [round(np.log10(val), 4) for val in xi_vals_projected]
            log_xi_vals_projected = ["N/A" if np.isinf(val) else val for val in log_xi_vals_projected]
                    
            if len(x_vals_real) > 0:
                fig.add_trace(
                    go.Scatter(
                        x=x_vals_real,
                        y=y_vals_real,
                        marker=dict(
                            color=path_point_fill_color,
                            line=dict(
                                color=path_point_line_color,
                                width=2
                            )
                        ),
                        mode='markers',
                        showlegend=False,
                        name='reaction path',
                        customdata = np.stack((xi_vals_real, log_xi_vals_real), axis=-1),
                        hovertemplate = 'Xi: %{customdata[0]}<br>log Xi: %{customdata[1]}<br>'+xlab+': %{x}<br>'+ylab+': %{y}<extra></extra>',
                        legendgroup='reaction path',
                    )
                )
            if len(x_vals_projected) > 0:
                
                fig.add_trace(
                    go.Scatter(
                        x=x_vals_projected,
                        y=y_vals_projected,
                        marker=dict(
                            color=projected_point_fill_color,
                            line=dict(
                                color=projected_point_line_color,
                                width=2
                            )
                        ),
                        mode='markers',
                        showlegend=False,
                        name='reaction path',
                        customdata = np.stack((xi_vals_projected, log_xi_vals_projected), axis=-1),
                        hovertemplate = 'Xi: %{customdata[0]}<br>log Xi: %{customdata[1]}<br>'+xlab+': %{x}<br>'+ylab+': %{y}<extra></extra>',
                        legendgroup='reaction path',
                    )
                )
            

        fig.update_xaxes(range=self.__get_plot_range(x_vals, y_vals)[0], autorange=False)
        fig.update_yaxes(range=self.__get_plot_range(x_vals, y_vals)[1], autorange=False)

        return fig

    
    def __plot_reaction_path_main(self,
                                  triad, T=25, P=1, path_margin=0.25,
                                  flip_xy=False,
                                  show_annotation=False,
                                  annotation_coords=[0,0],
                                  show_nonparticipating_mineral_lines = False,
                                  minerals_to_show=[],
                                  path_line_type = "markers+lines",
                                  path_line_color = "black",
                                  path_point_fill_color = "black",
                                  path_point_line_color = "black",
                                  projected_point_color = "white",
                                  projected_point_fill_color= "white",
                                  projected_point_line_color = "black",
                                  h_line_color="red",
                                  v_line_color="blue",
                                  d_line_color="orange",
                                  res=300,
                                  plot_width=4,
                                  plot_height=3,
                                  ppi=122,
                                  colormap="viridis",
                                  borders=0,
                                  projected_points=[],
                                  first_pass=False):

        e_pair = triad[0:2]
        
        if len(triad) == 3:
            div_var_name = triad[2]
        elif len(triad) == 2:
            div_var_name = "None"

        if flip_xy:
            e_pair.reverse()
            
        basis_species_x = self.__get_basis_from_elem(e_pair[0])
        basis_species_y = self.__get_basis_from_elem(e_pair[1])

        basis_sp_list = list(set([self.__get_basis_from_elem(e) for e in triad] + ["H+","H2O"]))
    
        try:
            pyCHNOSZ.basis(basis_sp_list)
        except:
            pyCHNOSZ.basis(basis_sp_list + ["H2"])
        
        elems = []
        for elem in triad:
            elems.append(elem.split("+")[0].split("-")[0])
        
        mineral_names = []
        for elem in triad:
            elem = elem.split("+")[0].split("-")[0]
            m_idx = pyCHNOSZ.retrieve((elem), list(set(["O", "H"]+elems)), state=["cr"], messages=False)
            if len(m_idx) > 0:
                mineral_names += list(pyCHNOSZ.info(m_idx, messages=False)["name"])
        
        # exclude inactive minerals
        mineral_names_active = [m for m in mineral_names if m in list(self.df["name"])]
        mineral_names = list(set(mineral_names_active))
        
        mineral_formula_ox = [self.__get_elem_ox_of_interest_in_minerals(m) for m in mineral_names]
        mineral_formula_ox_singles = [e if len(e)==1 else [] for e in mineral_formula_ox]
        mineral_formula_ox_doubles = [e if len(e)==2 else [] for e in mineral_formula_ox]
        mineral_formula_ox_triples = [e if len(e)==3 else [] for e in mineral_formula_ox]
        
        retrieved_minerals = []
        for i,s in enumerate(mineral_formula_ox_singles):
            if [e_pair[0]] == s:
                retrieved_minerals.append(mineral_names[i])
        x_minerals_to_plot = retrieved_minerals

        retrieved_minerals = []
        for i,s in enumerate(mineral_formula_ox_singles):
            if [e_pair[1]] == s:
                retrieved_minerals.append(mineral_names[i])
        y_minerals_to_plot = retrieved_minerals

        xy_minerals_to_plot = []
        for i,s in enumerate(mineral_formula_ox_doubles):
            if e_pair[0] in s and e_pair[1] in s:
                xy_minerals_to_plot.append(mineral_names[i])

        if len(x_minerals_to_plot) > 1:
            pyCHNOSZ.species(x_minerals_to_plot, add=True)
        if len(y_minerals_to_plot) > 1:
            pyCHNOSZ.species(y_minerals_to_plot, add=True)
        if len(xy_minerals_to_plot) > 1:
            pyCHNOSZ.species(xy_minerals_to_plot, add=True)

#         print("x, y, xy")
#         print(x_minerals_to_plot)
#         print(y_minerals_to_plot)
#         print(xy_minerals_to_plot)
            
        field_minerals_to_plot = []
        for i,s in enumerate(mineral_formula_ox_triples):
            if e_pair[0] in s and e_pair[1] in s and div_var_name in s:
                field_minerals_to_plot.append(mineral_names[i])
        for i,s in enumerate(mineral_formula_ox_doubles):
            if e_pair[0] in s and div_var_name in s:
                field_minerals_to_plot.append(mineral_names[i])
            elif e_pair[1] in s and div_var_name in s:
                field_minerals_to_plot.append(mineral_names[i])
        for i,s in enumerate(mineral_formula_ox_singles):
            if div_var_name in s:
                field_minerals_to_plot.append(mineral_names[i])
                
        field_minerals_exist = True
        try:
            # get minerals relevant to plotted element pair
            pyCHNOSZ.species(field_minerals_to_plot)#, add=True)
        except:
            field_minerals_exist = False
        
#         print("field minerals")
#         print(field_minerals_to_plot)
    
        xi_vals, x_vals, y_vals = self.__get_reaction_path(basis_species_x, basis_species_y, div_var_name)
            
        
        if self.P <= 1:
            bar_bars = "bar"
        else:
            bar_bars = "bars"
        
        if show_annotation and div_var_name != "None":
            annotation = "Balanced on: "+chemlabel(self.__get_basis_from_elem(div_var_name))+"<br>"+'%g'%(self.T)+" °C, "+'%g'%(self.P)+" "+bar_bars
        elif show_annotation:
            annotation = '%g'%(self.T)+" °C, "+'%g'%(self.P)+" "+bar_bars
        else:
            annotation = None
        
        table,fig,pred_minerals_from_fields = self.__plot_reaction_path_background(
            basis_species_x, basis_species_y, div_var_name, x_vals, y_vals,
            plot_width=plot_width, plot_height=plot_height, ppi=ppi, res=res,
            colormap=colormap, borders=borders,
            annotation_coords=annotation_coords,
            field_minerals_exist=field_minerals_exist, path_margin=self.path_margin,
            annotation=annotation, messages=False)

        # plot minerals with a single element of interest as a line
        plot_x_range, plot_y_range = self.__get_plot_range(x_vals, y_vals)
        
        line_styles = ['dot', 'dash', 'dashdot', 'longdash', 'longdashdot']
        line_i = 0
        line_width = 1
        
        if first_pass:
            fig = None
        else:
            for mineral in x_minerals_to_plot + y_minerals_to_plot + xy_minerals_to_plot:
                
                if not show_nonparticipating_mineral_lines and mineral not in list(self.moles_product_minerals.columns) and mineral not in minerals_to_show:
                    continue
                
                # deal with mineral line style (dot, dashed, etc.)
                if line_i % 5 == 0:
                    if line_i != 0:
                        line_width += 0.5
                    line_i = 0

                    
                line_style = line_styles[line_i]
                    
                eoi = self.__get_elem_ox_of_interest_in_minerals(mineral)

                logK = self.__calc_dissrxn_logK(mineral, T, P)
                mineral_formula_dict = self.__get_mineral_elem_ox_dict_interest(mineral)

                xlab,ylab = self.__get_xy_labs(basis_species_x, basis_species_y)
                
                if len(eoi) == 1:

                    # if the element of interest is not in the current element pair, move on
                    if eoi[0] not in e_pair:
                        continue

                    if self.__get_elem_ox_of_interest_in_minerals(mineral)[0] == e_pair[0]:
                        # vertical line
                        x0, x1 = (1/mineral_formula_dict[e_pair[0]])*logK, (1/mineral_formula_dict[e_pair[0]])*logK
                        y0 = min(plot_y_range)
                        y1 = max(plot_y_range)
                        color = v_line_color
                        hovertemplate=mineral+'<br>'+xlab+' = '+str(round(x0, 3))
                    elif self.__get_elem_ox_of_interest_in_minerals(mineral)[0] == e_pair[1]:
                        # horizontal line
                        y0, y1 = (1/mineral_formula_dict[e_pair[1]])*logK, (1/mineral_formula_dict[e_pair[1]])*logK
                        x0 = min(plot_x_range)
                        x1 = max(plot_x_range)
                        color=h_line_color
                        hovertemplate=mineral+'<br>'+ylab+' = '+str(round(y0))

                if len(eoi) == 2:
                    
                    line_slope = mineral_formula_dict[e_pair[0]]/mineral_formula_dict[e_pair[1]]
                    
                    x0 = min(plot_x_range)
                    x1 = max(plot_x_range)
                    y0 = (1/mineral_formula_dict[e_pair[1]])*logK - line_slope*x0
                    y1 = (1/mineral_formula_dict[e_pair[1]])*logK - line_slope*x1
                    
                    intercept = (1/mineral_formula_dict[e_pair[1]])*logK
                    
                    color = d_line_color
                    
                    hovertemplate = mineral+'<br>slope = '+str(round(line_slope))+'<br>intercept = '+str(round(intercept))+'<extra></extra>'

                fig.add_trace(
                    go.Scatter(x=[x0, x1], y=[y0, y1], mode="lines",
                               name=mineral,
                               line=dict(color=color, width=line_width, dash=line_style),
                               hovertemplate=hovertemplate,
                              ),
                )
                
                line_i += 1
        
            fig = self.__add_reaction_path_to_plot(x_vals, y_vals, xi_vals, fig,
                                                   basis_species_x, basis_species_y,
                                                   path_margin=self.path_margin,
                                                   projected_points=projected_points,
                                                   path_line_type=path_line_type,
                                                   path_line_color=path_line_color,
                                                   path_point_fill_color=path_point_fill_color,
                                                   path_point_line_color=path_point_line_color,
                                                   projected_point_fill_color=projected_point_fill_color,
                                                   projected_point_line_color=projected_point_line_color)
        
        pred_minerals_from_lines = x_minerals_to_plot+y_minerals_to_plot+xy_minerals_to_plot
        return fig, pred_minerals_from_fields, pred_minerals_from_lines

    
    def __get_xlab_xvar(self, x_type):
    
        x_type_dict = {
            "logxi" : ["log Xi", "log Xi"],
            "xi" : ["Xi", "Xi"],
            "temperature" : ["Temp(C)", "Temperature, °C"],
            "pressure" : ["Press(bars)", "Pressure, bars"],
            "pH" : ["pH", "pH"],
            "pmH" : ["pmH", "pmH"],
            "logfO2" : ["log fO2", "log <i>f</i>O<sub>2</sub>"],
            "Eh" : ["Eh(v)", "Eh, volts"],
            "pe" : ["pe", "pe"],
            "aw" : ["aw", "aw"],
        }
        
        if x_type not in x_type_dict.keys():
            self.err_handler.raise_exception(("x_type must be set to either "
                "'logxi', 'xi', 'temperature', 'pressure', 'pH', 'pmH',"
                "'log fO2', 'Eh(v)', 'pe', or 'aw'."))
        
        xvar = x_type_dict[x_type][0]
        xlab = x_type_dict[x_type][1]
            
        return xlab, xvar
    
    
    def plot_elements(self, plot_elements=None, units="molality", log=True,
                      x_type="logxi", plot_width=4, plot_height=3, ppi=122,
                      ylim=None, show_legend=True, charge_sign_at_end=False,
                      save_as=None, save_format=None, save_scale=1):
        
        """
        Generate a line plot of the log activities of aqueous species as a
        function of the log of the extent of reaction (log Xi).
        
        Parameters
        ----------
        plot_elements : list of str, optional
            A list of elements to plot. If undefined, every element will be
            plotted at once.
        
        units : str, default "molality"
            Units of elemental abundance to plot. Can be "molality" or "ppm".
        
        log : bool, default True
            Display elemental abundances in log scale?
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        ylim : list of two numeric values, optional
            Minimum and maximum value of the y-axis.
            
        show_legend : bool, default True
            Show the legend?
            
        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-)?
            
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
            
        Returns
        -------
        fig : Plotly figure object
            A line plot.
        """
        
        title = "Concentrations of dissolved elements"
        if units == "molality":
            df = pd.concat([self.dissolved_elements_molal, self.misc_params[self.misc_params.columns[1:]]], axis=1)
            ylab = "{}molality"
        elif units == "ppm":
            df = pd.concat([self.dissolved_elements_ppm, self.misc_params[self.misc_params.columns[1:]]], axis=1)
            ylab = "{}ppm"
#             elif units == "molarity":

        plot_columns = [col for col in df.columns]
        if isinstance(plot_elements, list):
            plot_columns_temp = [col for col in plot_columns if col in plot_elements]
            plot_columns = plot_columns_temp
            
        df = pd.melt(df, id_vars=list(self.misc_params.columns), value_vars=plot_columns)
        df.columns = list(self.misc_params.columns)+["variable", "value"]
        df["variable"] = df["variable"].apply(chemlabel, charge_sign_at_end=charge_sign_at_end)

        df["Xi"] = pd.to_numeric(df["Xi"])

        df["value"] = pd.to_numeric(df["value"])
        df["value"] = df["value"].fillna(0)
        df["value"] = df["value"].replace(0, np.nan)

        with np.errstate(divide='ignore'):
            df['log Xi'] = np.log10(df['Xi'])
            if log:
                df['value'] = np.log10(df['value'])
                ylab = ylab.format("log ")
        if not log:
            ylab = ylab.format("")
        
        xlab, xvar = self.__get_xlab_xvar(x_type)

        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                              width=plot_width*ppi,  height=plot_height*ppi,
                              labels=dict(value=ylab, x=xlab), render_mode='svg',
                             )

        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          legend_title=None,
                          showlegend=show_legend,
                          yaxis={'showexponent': 'all',
                                 'exponentformat': 'power'},
                         )

        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})
            
        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)

        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
            
        return fig
    
    
    def plot_pH(self, x_type="logxi", show_neutrality=True, title=None,
                plot_width=4, plot_height=3, ppi=122, ylim=None, save_as=None,
                save_format=None, save_scale=1):
        
        """
        Generate a line plot of pH as a function of the log of the extent of
        reaction (log Xi) or some other variable.
        
        Parameters
        ----------
        x_type : str, default "logxi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".
        
        show_neutrality : bool, default True,
            Display a reference line representing neutral pH? Setting this
            option to True requires a thermodynamic database in a WORM-style CSV
            format, e.g., 'wrm_data.csv'.
        
        title : str
            Title of the plot to display.
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        ylim : list of two numeric values, optional
            Minimum and maximum value of the y-axis.
            
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
            
        Returns
        -------
        fig : Plotly figure object
            A line plot.
        """

        df = copy.deepcopy(self.misc_params)
        
        xlab, xvar = self.__get_xlab_xvar(x_type)
        
        if x_type == "logxi":
            with np.errstate(divide='ignore'):
                df['log Xi'] = np.log10(df['Xi'])

        fig = px.line(df, x=xvar, y="pH", template="simple_white",
                      width=plot_width*ppi,  height=plot_height*ppi,
                      labels=dict(value="pH", x=xlab), render_mode='svg',
                      )

        fig.update_layout(xaxis_title=xlab,
                          yaxis_title="pH",
                          showlegend=False)
        
        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})

        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)
        
        if show_neutrality and isinstance(self.thermo.csv_db, pd.DataFrame):
            _, df_pH = self.plot_energy(species=["H2O", "H+", "OH-"],
                                        stoich=[-1, 1, 1],
                                        divisor=-2,
                                        x_type=x_type,
                                        y_type="logK",
                                        df_out=True,
                                        )
            
            fig.add_trace(go.Scatter(x=df_pH[xvar],
                                     y=df_pH["logK"],
                                     mode='lines',
                                     name='neutral pH',
                                     showlegend=True,
                                     hovertemplate = xlab+': %{x}<br>pH: %{y}<extra></extra>',
                                     line=dict(color='silver',
                                               width=3,
                                               dash='dot')
                                    ))

            fig['data'][0]['showlegend']=True
            fig['data'][0]['name']='pH'
            fig['data'][1]['showlegend']=True
            fig['data'][1]['name']='neutral pH'
            fig.update_layout(showlegend=True)

            
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
            
        return fig
    
    
    def plot_product_minerals(self, show_reactant_minerals=False,
                              plot_minerals=None, x_type="logxi", y_type="mole",
                              log_y=True, df_out=False, markers=False,
                              plot_width=4, plot_height=3, ppi=122, ylim=None,
                              show_legend=True, save_as=None, save_format=None,
                              save_scale=1):
        
        """
        Generate a line plot of the log moles of product minerals as a
        function of the log of the extent of reaction (log Xi) or some other
        variable.
        
        Parameters
        ----------
        show_reactant_minerals : bool, default False
            Include log moles of reactant minerals?
            
        plot_minerals : list, optional
            List of minerals to plot. Useful for isolating one or more
            minerals.
            
        x_type : str, default "logxi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".
            
        y_type : str, default 'mole'
            The variable to plot on the y-axis. Can be either 'mole' (for moles
            of minerals), 'mass' (for masses of minerals), or 'volume' (for
            volumes of minerals).
        
        log_y : bool, default True
            Should the y-axis be logarithmic?
            
        df_out : bool, default False
            Should a dataframe of values also be returned? For example, if
            `y_type` is set to 'volume', should a table of mineral volumes be
            returned?
            
        markers : bool, default True
            Add circular markers to lines to indicate calculation steps?
            
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        ylim : list of two numeric values, optional
            Minimum and maximum value of the y-axis.
            
        show_legend : bool, default True
            Show the legend?
            
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
            
        Returns
        -------
        fig : Plotly figure object
            A line plot.
            
        df : a Pandas dataframe
            A dataframe is only returned if `df_out` is set to True (it is
            set to False by default).
        
        """
        
        xlab, xvar = self.__get_xlab_xvar(x_type)
        
        if show_reactant_minerals:
            df = pd.concat([self.moles_minerals, self.misc_params[self.misc_params.columns[1:]]], axis=1)
            title = "{}{} of reactant and product minerals"
        else:
            df = pd.concat([self.moles_product_minerals, self.misc_params[self.misc_params.columns[1:]]], axis=1)
            title = "{}{} of product minerals"
            
        if log_y:
            log_text = "log "
        else:
            log_text = ""

        # sort in order of appearance along Xi
        sort_order = list(self.moles_minerals.columns)
        
        sort_order = [v for v in sort_order if v not in self.misc_params.columns]
            
        if y_type == "mole":
            ylab = "{}moles".format(log_text)
            title = title.format(log_text, "moles")
        elif y_type == "mass": # not yet supported
            y_lab = "{}grams".format(log_text)
            title = title.format(log_text, "masses")
            self.err_handler.raise_exception("Plotting mineral masses is not yet "
                    "supported.")
        elif y_type == "volume":

            if not isinstance(self.df, pd.DataFrame):
                self.err_handler.raise_exception("Plotting mineral volume "
                        "requires a WORM-formatted CSV thermodynamic database. "
                        "You may be seeing this message because your speciation "
                        "used a data0-type thermodynamic database (e.g., 'wrm'), "
                        "which does not contain mineral volume data.")
            
            ylab = "{}cm<sup>3</sup>".format(log_text)
            title = title.format(log_text, "volumes")
            temps = df["Temp(C)"]
            
            minerals = [col for col in df.columns if col not in list(self.misc_params.columns)]

            for i,T in enumerate(temps):
                for ii,mineral in enumerate(minerals):
                    mineral_df = copy.deepcopy(self.df[self.df["name"]==mineral])
                    polymorph_idxs = []
                    for iii in range(0, mineral_df.shape[0]): # loop through mineral polymorphs

                        if float(T) < float(list(mineral_df["z.T"])[0]):
                            polymorph_idxs.append(iii)
                    if len(polymorph_idxs)==0:
                        polymorph_idx = iii
                    else:
                        polymorph_idx = polymorph_idxs[0]


                    partial_molal_volume = list(mineral_df["V"])[polymorph_idx]

                    df.at[i, mineral] = df[mineral][i]*partial_molal_volume
        else:
            self.err_handler.raise_exception("y_type must be either 'mole', "
                        "'mass', or 'volume'.")
            
        
        plot_columns = [col for col in df.columns if col not in list(self.misc_params.columns)]

        if isinstance(plot_minerals, list):
            plot_columns_temp = [col for col in plot_columns if col in plot_minerals]
            plot_columns = plot_columns_temp
            
        plot_columns = sorted(plot_columns, key=sort_order.index)
        
        df = pd.melt(df, id_vars=list(self.misc_params.columns), value_vars=plot_columns)
        df.columns = list(self.misc_params.columns)+["variable", "value"]
        
        df = df[df["variable"] != "None"]

        df["Xi"] = pd.to_numeric(df["Xi"])

        df["value"] = pd.to_numeric(df["value"])
        df["value"] = df["value"].fillna(0)

        with np.errstate(divide='ignore'):
            df['log Xi'] = np.log10(df['Xi'])
            if log_y:
                df["value"] = df["value"].replace(0, np.nan)
                df['value'] = np.log10(df['value'])

        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                      width=plot_width*ppi,  height=plot_height*ppi, markers=markers,
                      labels=dict(value=ylab, x=xlab), render_mode='svg',
                      )
        
        df_to_return = copy.deepcopy(df)
        
        # add lines that go to -9999 (representing -Inf) in log y plots
        if log_y:
            if not isinstance(ylim, list):
                # grab the automatic y-axis range from the plot, above
                full_fig = fig.full_figure_for_development(warn=False)
                ylim = list(full_fig.layout.yaxis.range)
                
            # add new rows
            for irow in range(0, df.shape[0]):
                if not isnan(df["value"].iloc[irow]) and irow != 0 and irow != df.shape[0]-1:
                    if isnan(df["value"].iloc[irow-1]):
                        
                        df_dict = {}
                        for col in list(self.misc_params.columns)+["log Xi"]:
                            df_dict[col] = df.loc[irow-1, col]
                        df_dict["variable"] = df.loc[irow-1, "variable"]
                        df_dict["value"] = df.loc[irow-1, "value"]
                            
                        new_row = pd.DataFrame(df_dict, index=[irow-1])
                        df = pd.concat([df.iloc[:irow], new_row, df.iloc[irow:]]).reset_index(drop=True)
                        irow=0
                        
                    if isnan(df["value"].iloc[irow+1]):
                        
                        df_dict = {}
                        for col in list(self.misc_params.columns)+["log Xi"]:
                            df_dict[col] = df.loc[irow+1, col]
                        df_dict["variable"] = df.loc[irow+1, "variable"]
                        df_dict["value"] = df.loc[irow+1, "value"]
                        
                        new_row = pd.DataFrame(df_dict, index=[irow+1])
                        df = pd.concat([df.iloc[:irow], new_row, df.iloc[irow:]]).reset_index(drop=True)
                        irow=0
            
            # fill new rows with -9999 value
            for irow in range(0, df.shape[0]):
                if not isnan(df["value"].iloc[irow]) and irow != 0 and irow != df.shape[0]-1:
                    if isnan(df["value"].iloc[irow-1]):
                        df.loc[irow-1, "value"] = -9999
                    if isnan(df["value"].iloc[irow+1]):
                        df.loc[irow+1, "value"] = -9999
                        
            # re-create the plot with log values down to -9999. This will screw up the
            # automatic y-axis range, but we grabbed it from the first generated plot.
            # We will reset the range in update_layout(yaxis_range) a little later.
            fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                          width=plot_width*ppi,  height=plot_height*ppi, markers=markers,
                          labels=dict(value=ylab, x=xlab), render_mode='svg',
                          )

        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          legend_title=None,
                          showlegend=show_legend)
        
        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})
        
        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)
        
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
        
        if df_out:
            df_to_return = pd.pivot_table(df_to_return, index='log Xi', columns='variable', values='value').reset_index()
            df_to_return.columns = [col for col in df_to_return.columns[:]] # make index column name blank
            return df_to_return, fig
        else:
            return fig


    def plot_mineral_saturation(self, solid_solutions=False, plot_minerals=None,
                                      x_type="logxi", y_type="affinity", log_y=True,
                                      df_out=False, markers=False, title=None,
                                      plot_width=4, plot_height=3, ppi=122, ylim=None,
                                      show_legend=True, save_as=None, save_format=None,
                                      save_scale=1):
        
        """
        Generate a line plot of the saturation indices of minerals or solid
        solutions as a function of the log of the extent of reaction (log Xi) or
        some other variable.
        
        Parameters
        ----------
        plot_minerals : list, optional
            List of minerals to plot. Useful for isolating one or more
            minerals.
            
        x_type : str, default "logxi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".
            
        y_type : str, default 'affinity'
            The variable to plot on the y-axis. Can be either 'affinity'
            or 'logQ/K'.
            
        df_out : bool, default False
            Should a dataframe of values also be returned? For example, if
            `y_type` is set to 'volume', should a table of mineral volumes be
            returned?
            
        markers : bool, default True
            Add circular markers to lines to indicate calculation steps?

        title : str, optional
            Used to customize the title of the plot.
            
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        ylim : list of two numeric values, optional
            Minimum and maximum value of the y-axis.
            
        show_legend : bool, default True
            Show the legend?
            
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
            
        Returns
        -------
        fig : Plotly figure object
            A line plot.
            
        df : a Pandas dataframe
            A dataframe is only returned if `df_out` is set to True (it is
            set to False by default).
        
        """
        
        xlab, xvar = self.__get_xlab_xvar(x_type)

        if not isinstance(title, str):
            title = "Saturation indices of "
        
        if not solid_solutions:
            if y_type == "affinity":
                df = pd.concat([self.saturation_states_pure_solids_affinity,
                                self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "affinity, kcal/mol"
            else:
                df = pd.concat([self.saturation_states_pure_solids_log_Q_over_K,
                                self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "logQ/K"
            
            title = title + "minerals"
        else:
            if y_type == "affinity":
                df = pd.concat([self.saturation_states_solid_solutions_affinity,
                                self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "affinity, kcal/mol"
            else:
                df = pd.concat([self.saturation_states_solid_solutions_log_Q_over_K,
                                self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "logQ/K"
            title = title + "solid solutions"
            
        plot_columns = [col for col in df.columns if col not in list(self.misc_params.columns)]

        if isinstance(plot_minerals, list):
            plot_columns_temp = [col for col in plot_columns if col in plot_minerals]
            plot_columns = plot_columns_temp
        
        df = pd.melt(df, id_vars=list(self.misc_params.columns), value_vars=plot_columns)
        df.columns = list(self.misc_params.columns)+["variable", "value"]
        
        df = df[df["variable"] != "None"]

        df["Xi"] = pd.to_numeric(df["Xi"])

        df["value"] = pd.to_numeric(df["value"])
        df["value"] = df["value"].fillna(0)

        with np.errstate(divide='ignore'):
            df['log Xi'] = np.log10(df['Xi'])

        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                      width=plot_width*ppi,  height=plot_height*ppi, markers=markers,
                      labels=dict(value=ylab, x=xlab), render_mode='svg',
                      )
        

        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          legend_title=None,
                          showlegend=show_legend)
        
        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})
        
        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)
        
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
        
        return fig







    

    
    def plot_aqueous_species(self, plot_basis=False, plot_species=None,
                             x_type="logxi", y_type="log activity",
                             initially_visible=None, show_legend=True,
                             charge_sign_at_end=False,
                             plot_width=4, plot_height=3, ppi=122, xlim=None,
                             ylim=None, save_as=None, save_format=None,
                             save_scale=1):
        
        """
        Generate a line plot of the log activities of aqueous species as a
        function of the log of the extent of reaction (log Xi) or some other
        variable.
        
        Parameters
        ----------
        plot_basis : bool, default False
            Plot basis species only?
        
        plot_species : list of str, optional
            A list of aqueous species to plot. If undefined, every species at
            will be plotted at once.
        
        x_type : str, default "logxi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".
        
        y_type : str, default 'log activity'
            The variable to plot on the y-axis. Can be either 'log activity',
            'molality', or 'log molality'.
        
        initially_visible : list of str, optional
            A list of aqueous species that will be visible on the plot
            initially. All other species will be hidden, but can still be
            toggled back on in the legend.
        
        show_legend : bool, default True
            Show the legend?
        
        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-)?
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        xlim, ylim : list of two numeric values, optional
            Minimum and maximum value of the x-axis and y-axis, respectively.
            
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
            
        Returns
        -------
        fig : Plotly figure object
            A line plot.
        """
        
        if plot_basis:

            df = pd.concat([self.basis_logact, self.misc_params[self.misc_params.columns[1:]]], axis=1)
            title = "Solute basis species"
            ylab = "log activity"

        else:
            if y_type == "log activity":
                df = pd.concat([self.aq_distribution_logact, self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "log activity"
            elif y_type == "molality":
                df = pd.concat([self.aq_distribution_molal, self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "molality"
            elif y_type == "log molality":
                df = pd.concat([self.aq_distribution_logmolal, self.misc_params[self.misc_params.columns[1:]]], axis=1)
                ylab = "log molality"
            else:
                self.err_handler.raise_exception("The chosen 'y_type' parameter "
                    "is not recognized. 'y_type' can be 'log activity', "
                    "'molality', or 'log molality'")
            title = "Solute species"
            
        plot_columns = [col for col in df.columns]
        if isinstance(plot_species, list):
            plot_columns_temp = [col for col in plot_columns if col in plot_species]
            plot_columns = plot_columns_temp
            
        df = pd.melt(df, id_vars=list(self.misc_params.columns), value_vars=plot_columns)
        df.columns = list(self.misc_params.columns)+["variable", "value"]
        
        df["variable"] = df["variable"].apply(chemlabel, charge_sign_at_end=charge_sign_at_end)

        df["Xi"] = pd.to_numeric(df["Xi"])

        df["value"] = pd.to_numeric(df["value"])
        df["value"] = df["value"].fillna(0)
        df["value"] = df["value"].replace(0, np.nan)

        with np.errstate(divide='ignore'):
            df['log Xi'] = np.log10(df['Xi'])
        
        xlab, xvar = self.__get_xlab_xvar(x_type)
        
        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                              width=plot_width*ppi,  height=plot_height*ppi,
                              labels=dict(value=ylab, x=xlab), render_mode='svg',
                             )

        if isinstance(initially_visible, list):
            initially_visible_html = [chemlabel(sp, charge_sign_at_end=charge_sign_at_end) for sp in initially_visible]
            for trace in fig['data']: 
                if (not trace['name'] in initially_visible_html):
                    trace['showlegend'] = True
                    trace['visible'] = 'legendonly'
        
        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          legend_title=None,
                          showlegend=show_legend)

        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})

        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)
            
        if isinstance(xlim, list):
            fig.update_layout(xaxis_range=xlim)
            
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
            
        return fig

    
    def plot_energy(self, species, stoich,
                    divisor=1, x_type="logxi", y_type="A", y_units="kcal", 
                    show_zero_line=False, limiting=None, xlab=None, ylab=None,
                    title=None, charge_sign_at_end=False, log_y=False,
                    plot_width=4, plot_height=3, ppi=122,
                    xlim=None, ylim=None, df_out=False,
                    save_as=None, save_format=None,
                    save_scale=1, print_logK_messages=False):
        
        """
        Generate a line plot of the energy profile for a reaction as a
        function of the log of the extent of reaction (log Xi) or some other
        variable.
        
        Parameters
        ----------
        species : list of str
            A list of species that match the order of the stoichiometric
            reaction coefficients in the `stoich` parameter.

        stoich : list of numeric
            A list of stoichiometric reaction coefficients that match the order
            of the species in the `species` parameter.
        
        divisor : float, default 1
            Divide all values in the energy profile by this number. Useful for
            calculating energy per electron transferred or similar.
        
        x_type : str, default "logxi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".
        
        y_type : str, default 'A'
            The variable to plot on the y-axis. Can be either 'A' (for chemical
            affinity), 'G' (for Gibbs free energy, ΔG), 'logK' (for the log
            of the equilibrium constant), 'logQ' (for the log of the reaction
            quotient), or 'E' for energy supply.
        
        y_units : str, default 'kcal'
            The unit that energy will be reported in (per mol for G and A, or
            per kg fluid for energy supply, or unitless for logK and logQ).
            Can be 'kcal', 'cal', 'J', or 'kJ'.
        
        show_zero_line : bool, default False
            If True, displays a dotted line where affinity or ΔG equals 0 (at
            equilibrium).
        
        limiting : str, optional
            Name of the species to act as the limiting reactant when calculating
            energy supply. If this parameter is left undefined, then a
            limiting reactant will be chosen automatically based on
            concentration and stoichiometry. This parameter is ignored unless
            `y_type` is set to 'E' (energy supply).
        
        xlab, ylab : str, optional
            Custom x and y axis labels.
        
        title : str, optional
            Title of the plot to display.
        
        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2- instead of
            SO4-2) in species names when the reaction is displayed in the plot
            title?
        
        plot_width, plot_height : numeric, default 4 by 3
            Width and height of the plot, in inches. Size of interactive plots
            is also determined by pixels per inch, set by the parameter `ppi`.
            
        ppi : numeric, default 122
            Pixels per inch. Along with `plot_width` and `plot_height`,
            determines the size of interactive plots.
            
        xlim, ylim : list of two numeric values, optional
            Minimum and maximum value of the x-axis and y-axis, respectively.
            
        df_out = bool, default False
            Return a pandas dataframe in addition to a figure?
            
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
            
        print_logK_messages : bool, default False
            Print pyCHNOSZ messages while the logK of the reaction is
            calculated?
            
        Returns
        -------
        fig : Plotly figure object, optionally a Pandas Dataframe
            A line plot. If `df_out` is True, also returns a dataframe.
        """
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
                    "not equal to the number of reported xi steps.")

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
        if reactants == invalid_limiting_reactants and y_type == "E":
            self.err_handler.raise_exception("Energy supply for this reaction "
                "cannot be calculated because none of the reactants are "
                "limiting. A limiting reactant must be aqueous and cannot be H+ "
                "or OH-.")
        
        if limiting != None:
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
        
        # create a dictionary of species logacts across xi
        s_logact_dict = {}
        s_molal_dict = {}
        for s in species:
            if s == "H+":
                s_logact_dict[s] = [-pH for pH in list(self.misc_params["pH"])]
                s_molal_dict[s] = [float("NaN")]*len(self.misc_params["pH"])
            elif s == "H2O":
                s_logact_dict[s] = [0]*len(self.misc_params["Temp(C)"])
                s_molal_dict[s] = [float("NaN")]*len(self.misc_params["Temp(C)"])
            elif list(self.thermo.csv_db[self.thermo.csv_db["name"]==s]["state"])[0] not in ["cr", "liq"]:
                if s in self.aq_distribution_logact.columns:
                    # aqueous species
                    s_logact_dict[s] = list(self.aq_distribution_logact[s])
                    s_molal_dict[s] = list(self.aq_distribution_molal[s])
                else:
                    self.err_handler.raise_exception("The species "+str(s)+" is "
                            "not among the distribution of aqueous species in "
                            "this calculation.")
            else:
                # liq and cr species
                s_logact_dict[s] = [0]*len(self.misc_params["Temp(C)"])
                s_molal_dict[s] = [float("NaN")]*len(self.misc_params["Temp(C)"])
                
            
        xlab, xvar = self.__get_xlab_xvar(x_type)
            
        if y_type not in ["logK", "logQ"]:
            if y_units in ["cal", "kcal"]:
                r_div = 4.184
            elif y_units in ["J", "kJ"]:
                r_div = 1
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
                if title == None:
                    title = "Equilibrium constant for the reaction<br>"+equation_to_display
                y_list.append(round(logK/divisor_i, 4))
                df_y_name = "logK"
                continue
                

            logQ = sum([st*s_logact_dict[sp][i] for st,sp in zip(stoich,species)])

            if y_type == "logQ":
                ylab_out = "log Q"
                if title == None:
                    title = "Reaction quotient for the reaction<br>"+equation_to_display
                y_list.append(round(logQ/divisor_i, 4))
                df_y_name = "logQ"
                continue
            
            else:
                A = 2.303 * R * (273.15+T) * (logK - logQ)  # affinity, unit = [cal/mol]
                A = A/k_div
                
                if title == None:
                    title = "Energy profile for the reaction<br>"+equation_to_display
                
                if y_type=="G":
                    G = -A # gibbs free energy, unit = [cal/mol]
                    y_list.append(G/divisor_i)
                    ylab_out="ΔG, {}/mol".format(y_units)
                    y_units_out = y_units+"/mol"
                elif y_type=="A":
                    y_list.append(round(A/divisor_i, 4))
                    ylab_out="A, {}/mol".format(y_units)
                    y_units_out = y_units+"/mol"
                elif y_type=="E":
                    
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
                        lr_list = [limiting]
                        
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

                    y_list.append(round(E/divisor_i, 4))
                    y_units_out = y_units+"/kg fluid"
                    ylab_out="Energy Supply, {}".format(y_units+"/kg fluid")
#                     else:
#                         lr_name_list.append("NA")
#                         y_list.append(float('NaN'))
#                         y_units_out = y_units+"/kg fluid"
#                         ylab_out="Energy Supply, {}".format(y_units+"/kg fluid")

                df_y_name = y_type+", "+y_units_out

            
        if xlab != None:
            xlab_out = xlab
        if ylab != None:
            ylab_out = ylab

        df = copy.deepcopy(self.misc_params)
        with np.errstate(divide='ignore'):
            df['log Xi'] = np.log10(df['Xi'])
            
        df[df_y_name] = y_list

        
        fig = px.line(df, x=xvar, y=df_y_name, log_y=log_y,
                      width=plot_width*ppi, height=plot_height*ppi,
                      template="simple_white")
        
        if y_type=="E":
            fig.add_trace(
                go.Scatter(
                    x=df[xvar],
                    y=df[df_y_name],
                    mode='lines',
                    customdata = lr_name_list,
                    hovertemplate = xlab+': %{x}<br>'+ylab_out+': %{y}<br>Limiting : %{customdata}<extra></extra>',
                )
            )
            fig.data[0].visible=False
        
        fig.update_layout(xaxis_title=xlab_out,
                          yaxis_title=ylab_out)
        
        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})

        if isinstance(ylim, list):
            fig.update_layout(yaxis_range=ylim)
            
        if isinstance(xlim, list):
            fig.update_layout(xaxis_range=xlim)
            
        if isinstance(save_as, str):
            dummy_sp = Speciation({})
            save_as, save_format = dummy_sp._save_figure(fig,
                    save_as, save_format, save_scale,
                    plot_width, plot_height, ppi)
            
        if show_zero_line:
            fig.add_hline(y=0, line_width=3, line_dash="dash", line_color="black")
            
        if df_out:
            return fig, df
        else:
            return fig
    
    
    def plot_mass_contribution(self, *args, x_type="xi", x_decimals=3,
                                     track_steps=True, keep_xi_order=False,
                                     **kwargs):
        
        """
        Generate a bar plot of mass contributions (in mole percent) of aqueous
        species formed as a function of reaction progress Xi or some other
        variable.
        
        Parameters
        ----------
        *args : iterable
            Arguments to be passed to `Speciation.plot_mass_contribution`.

        x_type : str, default "xi"
            Variable to appear on the x-axis. Can be "logxi", "xi",
            "temperature", "pressure", "pH", "pmH", "logfO2", "Eh", "pe", or
            "aw".

        x_decimals : int
            Number of decimals to display in the numeric values of the x-axis
            variable defined by `x_type`.

        track_steps : bool, default True
            Show reported xi steps on x-axis ticks? Useful for plotting as a
            function of an x-axis variable that can be out-of-order, like
            Eh or temperature. This parameter will become True automatically
            if it will prevent column stacking (which happens if there are
            duplicate x-axis values).

        **kargs : dict
            Keyword arguments to be passed to `Speciation.plot_mass_contribution`.
            
        Returns
        -------
        fig : Plotly figure object
            A mass contribution bar plot.
        """
        
        basis = args[0]
        
        xlab, kwargs["sample_label"] = self.__get_xlab_xvar(x_type)
        
        df_sp = pd.concat([self.mass_contribution_dict.get(basis, "error"),
                           self.misc_params[self.misc_params.columns[1:]]], axis=1)
        
        with np.errstate(divide='ignore'):
            df_sp['log Xi'] = np.log10(df_sp['Xi'])
        
        # handle display of the x axis variable
        # if number of decimals to display is too low, columns will stack
        # check to see if this happens, then increment x_decimals until the stacking problem is solved
        if x_decimals < 0:
            msg = "The parameter x_decimals must be greater than or equal to 0."
            self.err_handler.raise_exception(msg)
        original_x_decimals = copy.copy(x_decimals)
        
        len_unique_labels_rounded = len(set(df_sp[kwargs["sample_label"]].apply(lambda x: ('%.'+str(x_decimals)+'e') % x)))
        len_unique_labels_unrounded = len(set(df_sp[kwargs["sample_label"]]))
        
        if len_unique_labels_rounded < len_unique_labels_unrounded:
            x_decimals += 1
            solved_decimals = False
            for i in range(x_decimals, x_decimals+10):
                len_unique_labels_rounded = len(set(df_sp[kwargs["sample_label"]].apply(lambda x: ('%.'+str(x_decimals)+'e') % x)))
                if len_unique_labels_rounded == len_unique_labels_unrounded:
                    if self.verbose > 0:
                        print("Number of decimals to display for x-axis variable",
                              "increased to", x_decimals, "to prevent column stacking.")
                    solved_decimals = True
                    break
                else:
                    x_decimals += 1
            if not solved_decimals:
                msg = ("X-axis value decimal formatting is resulting in column "
                       "stacking even after attempting 10 increments of x_decimals.")
                self.err_handler.raise_exception(msg)

        if isinstance(df_sp, str):
            msg = ("The basis species {} ".format(basis)+"could not be found "
                   "among available basis species: "
                   "{}".format(str(list(self.mass_contribution_dict.keys()))))
            self.err_handler.raise_exception(msg)
        
        df_sp["position"] = list(range(0, df_sp.shape[0]))
        
        df_sp_melt = df_sp.melt(
                id_vars=list(self.misc_params.columns)+["log Xi", "basis", "factor", "molality", "position"])
        
        df_sp_melt.rename(columns={kwargs["sample_label"] : "sample",
                                   df_sp_melt.columns[-2] : "species",
                                   df_sp_melt.columns[-1] : "percent",
                                   }, inplace=True)
                                      
        df_sp_melt = df_sp_melt[df_sp_melt['percent'].notna()]
        
        if keep_xi_order:
            sample_order = ["position", "species", "percent"]
        else:
            sample_order = ["sample", "position", "species", "percent"]
            
        df_sp_melt.sort_values(sample_order, axis=0, ascending=True, inplace=True)

        if df_sp_melt.dtypes["sample"] != "O": # if the column isn't formatted as a string
            if any(["e" in v2 for v2 in [str(v1) for v1 in df_sp_melt['sample']]]): # if value is in scientific notation
                df_sp_melt['sample'] = df_sp_melt['sample'].apply(lambda x: ('%.'+str(x_decimals)+'e') % x) # converts numeric to string
            else:
                df_sp_melt['sample'] = [('{0:.'+str(x_decimals)+'f}').format(v) for v in df_sp_melt['sample']]
            
            
        if len(list(set(df_sp_melt["sample"]))) != len(list(set(df_sp_melt["position"]))) or track_steps:
            # handle duplicate x-axis values to prevent stacking
            temp_col = []
            for i,v in enumerate(df_sp_melt["sample"]):
                temp_col.append(" (step "+str(list(df_sp_melt["position"])[i])+")")
            df_sp_melt["sample"] = [str(v)+a for v,a in zip(df_sp_melt["sample"], temp_col)]
            if not track_steps and self.verbose > 0:
                print("Reported Xi step tracking has been added to x-axis ticks to prevent column stacking.")
            
        sp = Speciation(args={})
        sp.mass_contribution = df_sp_melt

        if not kwargs.get("plot_out", False):
            plot_out = False
        else:
            plot_out = True
        kwargs["plot_out"] = True
        
        fig = sp.plot_mass_contribution(*args, **kwargs)
        
        fig.update_layout(
            xaxis_title=xlab, # add an x axis title
        )
        
        if plot_out:
            return fig
        else:
            fig.show()


template = """|------------------------------------------------------------------------------|
| Main Title             | (utitl1(n))                                         |
|------------------------------------------------------------------------------|
|Sample: {{sample_name}}                                                         |
|Date created:  {date_created}                                       |
|EQ6 input file generated by AqEquil                                           |
|                                                                              |
|------------------------------------------------------------------------------|
|Temperature option (jtemp):                                                   |
|  [{t_checkbox_1}] ( 0) Constant temperature:                                              |
|             Value (C)         |{tval1}| (tempcb)                        |
|  [{t_checkbox_2}] ( 1) Linear tracking in Xi:                                             |
|             Base Value (C)    |{tval2}| (tempcb)                        |
|             Derivative        |{tval3}| (ttk(1))                        |
|  [{t_checkbox_3}] ( 2) Linear tracking in time:                                           |
|             Base Value (C)    |{tval4}| (tempcb)                        |
|             Derivative        |{tval5}| (ttk(1))                        |
|  [{t_checkbox_4}] ( 3) Fluid mixing tracking (fluid 2 = special reactant):                |
|             T of fluid 1 (C)  |{tval6}| (tempcb)                        |
|             T of fluid 2 (C)  |{tval7}| (ttk(2))                        |
|             Mass ratio factor |{tval8}| (ttk(1))                        |
|------------------------------------------------------------------------------|
|Pressure option (jpress):                                                     |
|  [{p_checkbox_1}] ( 0) Follow the data file reference pressure curve                      |
|  [{p_checkbox_2}] ( 1) Follow the 1.013-bar/steam-saturation curve                        |
|  [{p_checkbox_3}] ( 2) Constant pressure:                                                 |
|             Value (bars)      |{pval1}| (pressb)                        |
|  [{p_checkbox_4}] ( 3) Linear tracking in Xi:                                             |
|             Base Value (bars) |{pval2}| (pressb)                        |
|             Derivative        |{pval3}| (ptk(1))                        |
|  [{p_checkbox_5}] ( 4) Linear tracking in time:                                           |
|             Base Value (bars) |{pval4}| (pressb)                        |
|             Derivative        |{pval5}| (ptk(1))                        |
|------------------------------------------------------------------------------|
|Reactants (Irreversible Reactions) | (nrct)                                   |
|------------------------------------------------------------------------------|{reactant_blocks}
* Valid reactant type strings (urcjco(jcode(n))) are:                          *
*    Pure mineral                Solid solution                                *
*    Special reactant            Aqueous species                               *
*    Gas species                 Generic ion exchanger                         *
*------------------------------------------------------------------------------*
* Valid reactant status strings (urcjre(jreac(n))) are:                        *
*    Saturated, reacting         Reacting                                      *
*    Exhausted                   Saturated, not reacting                       *
*------------------------------------------------------------------------------*
* Valid forward rate law strings (urcnrk(nrk(1,n))) are:                       *
*    Use backward rate law       Relative rate equation                        *
*    TST rate equation           Linear rate equation                          *
*------------------------------------------------------------------------------*
* Valid backward rate law strings (urcnrk(nrk(2,n))) are:                      *
*    Use forward rate law        Partial equilibrium                           *
*    Relative rate equation      TST rate equation                             *
*    Linear rate equation                                                      *
*------------------------------------------------------------------------------*
|Starting, minimum, and maximum values of key run parameters.                  |
|------------------------------------------------------------------------------|
|Starting Xi value        |{start_xi}| (xistti)                              |
|------------------------------------------------------------------------------|
|Maximum Xi value         |{max_xi}| (ximaxi)                              |
|------------------------------------------------------------------------------|
|Starting time (seconds)  |{start_time}| (tistti)                              |
|------------------------------------------------------------------------------|
|Maximum time (seconds)   |{max_time}| (timmxi)                              |
|------------------------------------------------------------------------------|
|Minimum value of pH      |{min_pH}| (phmini)                              |
|------------------------------------------------------------------------------|
|Maximum value of pH      |{max_pH}| (phmaxi)                              |
|------------------------------------------------------------------------------|
|Minimum value of Eh (v)  |{min_Eh}| (ehmini)                              |
|------------------------------------------------------------------------------|
|Maximum value of Eh (v)  |{max_Eh}| (ehmaxi)                              |
|------------------------------------------------------------------------------|
|Minimum value of log fO2 |{min_fO2}| (o2mini)                              |
|------------------------------------------------------------------------------|
|Maximum value of log fO2 |{max_fO2}| (o2maxi)                              |
|------------------------------------------------------------------------------|
|Minimum value of aw      |{min_aw}| (awmini)                              |
|------------------------------------------------------------------------------|
|Maximum value of aw      |{max_aw}| (awmaxi)                              |
|------------------------------------------------------------------------------|
|Maximum number of steps  |{max_n_steps}| (kstpmx)                              |
|------------------------------------------------------------------------------|
|Print interval parameters.                                                    |
|------------------------------------------------------------------------------|
|Xi print interval        |{xi_print_int}| (dlxprn)                              |
|------------------------------------------------------------------------------|
|Log Xi print interval    |{log_xi_print_int}| (dlxprl)                              |
|------------------------------------------------------------------------------|
|Time print interval      |{time_print_int}| (dltprn)                              |
|------------------------------------------------------------------------------|
|Log time print interval  |{log_time_print_int}| (dltprl)                              |
|------------------------------------------------------------------------------|
|pH print interval        |{pH_print_interval}| (dlhprn)                              |
|------------------------------------------------------------------------------|
|Eh (v) print interval    |{Eh_print_interval}| (dleprn)                              |
|------------------------------------------------------------------------------|
|Log fO2 print interval   |{logfO2_print_interval}| (dloprn)                              |
|------------------------------------------------------------------------------|
|aw print interval        |{aw_print_interval}| (dlaprn)                              |
|------------------------------------------------------------------------------|
|Steps print interval     |{n_steps_print_interval}| (ksppmx)                              |
|------------------------------------------------------------------------------|
|Plot interval parameters.                                                     |
|------------------------------------------------------------------------------|
|Xi plot interval         | 1.00000E+38| (dlxplo)                              |
|------------------------------------------------------------------------------|
|Log Xi plot interval     | 1.00000E+38| (dlxpll)                              |
|------------------------------------------------------------------------------|
|Time plot interval       | 1.00000E+38| (dltplo)                              |
|------------------------------------------------------------------------------|
|Log time plot interval   | 1.00000E+38| (dltpll)                              |
|------------------------------------------------------------------------------|
|pH plot interval         | 1.00000E+38| (dlhplo)                              |
|------------------------------------------------------------------------------|
|Eh (v) plot interval     | 1.00000E+38| (dleplo)                              |
|------------------------------------------------------------------------------|
|Log fO2 plot interval    | 1.00000E+38| (dloplo)                              |
|------------------------------------------------------------------------------|
|aw plot interval         | 1.00000E+38| (dlaplo)                              |
|------------------------------------------------------------------------------|
|Steps plot interval      |       10000| (ksplmx)                              |
|------------------------------------------------------------------------------|
|Iopt Model Option Switches ("( 0)" marks default choices)                     |
|------------------------------------------------------------------------------|
|iopt(1) - Physical System Model Selection:                                    |
|  [{i1_checkbox_1}] ( 0) Closed system                                                      |
|  [{i1_checkbox_2}] ( 1) Titration system                                                   |
|  [{i1_checkbox_3}] ( 2) Fluid-centered flow-through open system                            |
|------------------------------------------------------------------------------|
|iopt(2) - Kinetic Mode Selection:                                             |
|  [{i2_checkbox_1}] ( 0) Reaction progress mode (arbitrary kinetics)                        |
|  [{i2_checkbox_2}] ( 1) Reaction progress/time mode (true kinetics)                        |
|------------------------------------------------------------------------------|
|iopt(3) - Phase Boundary Searches:                                            |
|  [{i3_checkbox_1}] ( 0) Search for phase boundaries and constrain the step size to match   |
|  [{i3_checkbox_2}] ( 1) Search for phase boundaries and print their locations              |
|  [{i3_checkbox_3}] ( 2) Don't search for phase boundaries                                  |
|------------------------------------------------------------------------------|
|iopt(4) - Solid Solutions:                                                    |
|  [{i4_checkbox_1}] ( 0) Ignore                                                             |
|  [{i4_checkbox_2}] ( 1) Permit                                                             |
|------------------------------------------------------------------------------|
|iopt(5) - Clear the ES Solids Read from the INPUT File:                       |
|  [{i5_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i5_checkbox_2}] ( 1) Do it                                                              |
|------------------------------------------------------------------------------|
|iopt(6) - Clear the ES Solids at the Initial Value of Reaction Progress:      |
|  [{i6_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i6_checkbox_2}] ( 1) Do it                                                              |
|------------------------------------------------------------------------------|
|iopt(7) - Clear the ES Solids at the End of the Run:                          |
|  [{i7_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i7_checkbox_2}] ( 1) Do it                                                              |
|------------------------------------------------------------------------------|
|iopt(9) - Clear the PRS Solids Read from the INPUT file:                      |
|  [{i9_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i9_checkbox_2}] ( 1) Do it                                                              |
|------------------------------------------------------------------------------|
|iopt(10) - Clear the PRS Solids at the End of the Run:                        |
|  [{i10_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i10_checkbox_2}] ( 1) Do it, unless numerical problems cause early termination           |
|------------------------------------------------------------------------------|
|iopt(11) - Auto Basis Switching in pre-N-R Optimization:                      |
|  [{i11_checkbox_1}] ( 0) Turn off                                                           |
|  [{i11_checkbox_2}] ( 1) Turn on                                                            |
|------------------------------------------------------------------------------|
|iopt(12) - Auto Basis Switching after Newton-Raphson Iteration:               |
|  [{i12_checkbox_1}] ( 0) Turn off                                                           |
|  [{i12_checkbox_2}] ( 1) Turn on                                                            |
|------------------------------------------------------------------------------|
|iopt(13) - Calculational Mode Selection:                                      |
|  [{i13_checkbox_1}] ( 0) Normal path tracing                                                |
|  [{i13_checkbox_2}] ( 1) Economy mode (if permissible)                                      |
|  [{i13_checkbox_3}] ( 2) Super economy mode (if permissible)                                |
|------------------------------------------------------------------------------|
|iopt(14) - ODE Integrator Corrector Mode Selection:                           |
|  [{i14_checkbox_1}] ( 0) Allow Stiff and Simple Correctors                                  |
|  [{i14_checkbox_2}] ( 1) Allow Only the Simple Corrector                                    |
|  [{i14_checkbox_3}] ( 2) Allow Only the Stiff Corrector                                     |
|  [{i14_checkbox_4}] ( 3) Allow No Correctors                                                |
|------------------------------------------------------------------------------|
|iopt(15) - Force the Suppression of All Redox Reactions:                      |
|  [{i15_checkbox_1}] ( 0) Don't do it                                                        |
|  [{i15_checkbox_2}] ( 1) Do it                                                              |
|------------------------------------------------------------------------------|
|iopt(16) - BACKUP File Options:                                               |
|  [ ] (-1) Don't write a BACKUP file                                          |
|  [x] ( 0) Write BACKUP files                                                 |
|  [ ] ( 1) Write a sequential BACKUP file                                     |
|------------------------------------------------------------------------------|
|iopt(17) - PICKUP File Options:                                               |
|  [ ] (-1) Don't write a PICKUP file                                          |
|  [x] ( 0) Write a PICKUP file                                                |
|------------------------------------------------------------------------------|
|iopt(18) - TAB File Options:                                                  |
|  [{i18_checkbox_1}] (-1) Don't write a TAB file                                             |
|  [{i18_checkbox_2}] ( 0) Write a TAB file                                                   |
|  [{i18_checkbox_3}] ( 1) Write a TAB file, prepending TABX file data from a previous run    |
|------------------------------------------------------------------------------|
|iopt(20) - Advanced EQ6 PICKUP File Options:                                  |
|  [{i20_checkbox_1}] ( 0) Write a normal EQ6 PICKUP file                                     |
|  [{i20_checkbox_2}] ( 1) Write an EQ6 INPUT file with Fluid 1 set up for fluid mixing       |
|------------------------------------------------------------------------------|
|Iopr Print Option Switches ("( 0)" marks default choices)                     |
|------------------------------------------------------------------------------|
|iopr(1) - Print All Species Read from the Data File:                          |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print                                                              |
|------------------------------------------------------------------------------|
|iopr(2) - Print All Reactions:                                                |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print the reactions                                                |
|  [ ] ( 2) Print the reactions and log K values                               |
|  [ ] ( 3) Print the reactions, log K values, and associated data             |
|------------------------------------------------------------------------------|
|iopr(3) - Print the Aqueous Species Hard Core Diameters:                      |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print                                                              |
|------------------------------------------------------------------------------|
|iopr(4) - Print a Table of Aqueous Species Concentrations, Activities, etc.:  |
|  [ ] (-3) Omit species with molalities < 1.e-8                               |
|  [ ] (-2) Omit species with molalities < 1.e-12                              |
|  [ ] (-1) Omit species with molalities < 1.e-20                              |
|  [x] ( 0) Omit species with molalities < 1.e-100                             |
|  [ ] ( 1) Include all species                                                |
|------------------------------------------------------------------------------|
|iopr(5) - Print a Table of Aqueous Species/H+ Activity Ratios:                |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print cation/H+ activity ratios only                               |
|  [ ] ( 2) Print cation/H+ and anion/H+ activity ratios                       |
|  [ ] ( 3) Print ion/H+ activity ratios and neutral species activities        |
|------------------------------------------------------------------------------|
|iopr(6) - Print a Table of Aqueous Mass Balance Percentages:                  |
|  [ ] (-1) Don't print                                                        |
|  [x] ( 0) Print those species comprising at least 99% of each mass balance   |
|  [ ] ( 1) Print all contributing species                                     |
|------------------------------------------------------------------------------|
|iopr(7) - Print Tables of Saturation Indices and Affinities:                  |
|  [ ] (-1) Don't print                                                        |
|  [x] ( 0) Print, omitting those phases undersaturated by more than 10 kcal   |
|  [ ] ( 1) Print for all phases                                               |
|------------------------------------------------------------------------------|
|iopr(8) - Print a Table of Fugacities:                                        |
|  [x] (-1) Don't print                                                        |
|  [ ] ( 0) Print                                                              |
|------------------------------------------------------------------------------|
|iopr(9) - Print a Table of Mean Molal Activity Coefficients:                  |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print                                                              |
|------------------------------------------------------------------------------|
|iopr(10) - Print a Tabulation of the Pitzer Interaction Coefficients:         |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print a summary tabulation                                         |
|  [ ] ( 2) Print a more detailed tabulation                                   |
|------------------------------------------------------------------------------|
|iopr(17) - PICKUP file format ("W" or "D"):                                   |
|  [x] ( 0) Use the format of the INPUT file                                   |
|  [ ] ( 1) Use "W" format                                                     |
|  [ ] ( 2) Use "D" format                                                     |
|------------------------------------------------------------------------------|
|Iodb Debugging Print Option Switches ("( 0)" marks default choices)           |
|------------------------------------------------------------------------------|
|iodb(1) - Print General Diagnostic Messages:                                  |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print Level 1 diagnostic messages                                  |
|  [ ] ( 2) Print Level 1 and Level 2 diagnostic messages                      |
|------------------------------------------------------------------------------|
|iodb(2) - Kinetics Related Diagnostic Messages:                               |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print Level 1 kinetics diagnostic messages                         |
|  [ ] ( 2) Print Level 1 and Level 2 kinetics diagnostic messages             |
|------------------------------------------------------------------------------|
|iodb(3) - Print Pre-Newton-Raphson Optimization Information:                  |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |
|  [ ] ( 3) Print more detailed information (including matrix equations)       |
|  [ ] ( 4) Print most detailed information (including activity coefficients)  |
|------------------------------------------------------------------------------|
|iodb(4) - Print Newton-Raphson Iteration Information:                         |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |
|  [ ] ( 3) Print more detailed information (including the Jacobian)           |
|  [ ] ( 4) Print most detailed information (including activity coefficients)  |
|------------------------------------------------------------------------------|
|iodb(5) - Print Step-Size and Order Selection:                                |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|  [ ] ( 2) Print detailed information                                         |
|------------------------------------------------------------------------------|
|iodb(6) - Print Details of Hypothetical Affinity Calculations:                |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|  [ ] ( 2) Print detailed information                                         |
|------------------------------------------------------------------------------|
|iodb(7) - Print General Search (e.g., for a phase boundary) Information:      |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|------------------------------------------------------------------------------|
|iodb(8) - Print ODE Corrector Iteration Information:                          |
|  [x] ( 0) Don't print                                                        |
|  [ ] ( 1) Print summary information                                          |
|  [ ] ( 2) Print detailed information (including the betar and delvcr vectors)|
|------------------------------------------------------------------------------|
|Mineral Sub-Set Selection Suppression Options | (nxopt)                       |
|------------------------------------------------------------------------------|
|Option  |Sub-Set Defining Species| (this is a table header)                   |
|------------------------------------------------------------------------------|{mineral_suppress_lines}
|------------------------------------------------------------------------------|
* Valid mineral sub-set selection suppression option strings (uxopt(n)) are:   *
*    None        All         Alwith      Allwith                               *
*------------------------------------------------------------------------------*
|Exceptions to the Mineral Sub-Set Selection Suppression Options | (nxopex)    |
|------------------------------------------------------------------------------|
|Mineral                 | (this is a table header)                            |
|------------------------------------------------------------------------------|
|None                    | (uxopex(n))                                         |
|------------------------------------------------------------------------------|
|Fixed Fugacity Options | (nffg)                                               |
|------------------------------------------------------------------------------|
|Gas                     |Moles to Add |Log Fugacity | --                      |
| (uffg(n))              | (moffg(n))  | (xlkffg(n)) | --                      |
|------------------------------------------------------------------------------|{gas_lines}
|------------------------------------------------------------------------------|
|Numerical Parameters                                                          |
|------------------------------------------------------------------------------|
|Max. finite-difference order               |{max_finite_difference_order}| (nordmx)            |
|Beta convergence tolerance                 |{beta_convergence_tolerance}| (tolbt)             |
|Del convergence tolerance                  |{del_convergence_tolerance}| (toldl)             |
|Max. No. of N-R iterations                 |{max_n_NR_iter}| (itermx)            |
|Search/find convergence tolerance          |{search_find_convergeance_tolerance}| (tolxsf)            |
|Saturation tolerance                       |{saturation_tolerance}| (tolsat)            |
|Max. No. of Phase Assemblage Tries         |{max_n_phase_assemblage_tries}| (ntrymx)            |
|Zero order step size (in Xi)               |{zero_order_step_size}| (dlxmx0)            |
|Max. interval in Xi between PRS transfers  |{max_interval_in_xi_between_PRS_transfers}| (dlxdmp)            |
|------------------------------------------------------------------------------|"""


rb_template = """
|Reactant        |{reactant_name}| (ureac(n))                         |
|------------------------------------------------------------------------------|
|->|Type         |{reactant_type}| (urcjco(jcode(n)))                 |
|------------------------------------------------------------------------------|
|->|Status       |{reactant_status}| (urcjre(jreac(n)))                 |
|------------------------------------------------------------------------------|
|->|Amount remaining (moles) |{amount_remaining}| (morr(n))                          |
|------------------------------------------------------------------------------|
|->|Amount destroyed (moles) |{amount_destroyed}| (modr(n))                          |
|------------------------------------------------------------------------------|{sr_block}
|->|Surface area option (nsk(n)):                                              |
|->|  [{sa_checkbox_1}] ( 0) Constant surface area:                                          |
|->|             Value (cm2)       |{sa_val_1}| (sfcar(n))                   |
|->|  [{sa_checkbox_2}] ( 1) Constant specific surface area:                                 |
|->|             Value (cm2/g)     |{sa_val_2}| (ssfcar(n))                  |
|->|  [{sa_checkbox_3}] ( 2) n**2/3 growth law- current surface area:                        |
|->|             Value (cm2)       |{sa_val_3}| (sfcar(n))                   |
|------------------------------------------------------------------------------|
|->|Surface area factor      |{sa_factor}| (fkrc(n))                          |
|------------------------------------------------------------------------------|
|->|Forward rate law          |{f_rate_law}| (urcnrk(nrk(1,n)))    |
|------------------------------------------------------------------------------|{f_rate_block}
|->|Backward rate law         |{b_rate_law}| (urcnrk(nrk(2,n)))    |
|------------------------------------------------------------------------------|{b_rate_block}"""


relative_rate_equation_template = """
|--->|dXi(n)/dXi (mol/mol)      |{eq1}| (rkb(1,1,n))                    |
|------------------------------------------------------------------------------|
|--->|d2Xi(n)/dXi2 (mol/mol2)   |{eq2}| (rkb(2,1,n))                    |
|------------------------------------------------------------------------------|
|--->|d3Xi(n)/dXi3 (mol/mol3)   |{eq3}| (rkb(3,1,n))                    |
|------------------------------------------------------------------------------|"""

# not yet supported
TST_equation_template = """
|--->|Mechanism  1                                                             |
|------------------------------------------------------------------------------|
|----->|sigma(i,+,n)            |{eq1}| csigma(i,1,n)                   |
|------------------------------------------------------------------------------|
|----->|k(i,+,n) (mol/cm2/sec)  |{eq2}| rkb(i,1,n)                      |
|------------------------------------------------------------------------------|
|----->|Ref. Temperature (C)    |{eq3}| trkb(i,1,n)                     |
|------------------------------------------------------------------------------|
|----->|Temperature dependence option (iact(i,1,n)):                           |
|----->|  [{rate_checkbox_1}] ( 0) No temperature dependence                                   |
|----->|  [{rate_checkbox_2}] ( 1) Constant activation energy:                                 |
|----->|             Value (kcal/mol)  |{eq4}| (eact(i,1,n))            |
|----->|  [{rate_checkbox_3}] ( 2) Constant activation enthalpy:                               |
|----->|             Value (kcal/mol)  |{eq5}| (hact(i,1,n))            |
|------------------------------------------------------------------------------|
|----->|Kinetic activity product species (ndact(i,1,n))                        |
|------------------------------------------------------------------------------|
|------->|Species                                         |-N(j,i,+,n)         |
|------->| (udac(j,i,1,n))                                | (cdac(j,i,1,n))    |
|------------------------------------------------------------------------------|{k_act_prod_species_block}
|------------------------------------------------------------------------------|"""
# part of TST rate equation; not yet supported
k_act_prod_species_template = """
|------->|{sp_name}}|{sp_N}|"""

gl_template = """
|{gas_name}|{gas_moles}|{gas_log_fugacity}| --                      |"""

mineral_suppress_template = """
|{mineral_suppress_option}|                        | (uxopt(n), uxcat(n))                       |"""

srb_template = """
|->|Molar volume (cm3/mol)   |{molar_volume}| (vreac(n))                         |
|------------------------------------------------------------------------------|
|->|Composition                                                                |
|------------------------------------------------------------------------------|
|--->|Element |Stoich. Number        | (this is a table header)                |
|------------------------------------------------------------------------------|{srb_stoich}
|------------------------------------------------------------------------------|
|->|Reaction                                                                   |
|------------------------------------------------------------------------------|
|--->|Species                 |Reaction Coefficient  | (this is a table header)|
|------------------------------------------------------------------------------|
|------------------------------------------------------------------------------|"""

srb_stoich_template = """
|--->|{elem}|{elem_val}| (uesri(i,n), cesri(i,n))                |"""


class Mixing_Fluid:
    def __init__(self,
                 speciation,
                 sample_name,
                 amount_remaining=1,
                 amount_destroyed=0,
                 molar_volume=1,
                 mass_ratio=1,
                 hide_traceback=True,
                ):

        """
        Class used to define the fluid to be mixed with other fluids in
        `Prepare_Reaction`.

        Parameters
        ----------
        speciation : object of class Speciation
            The speciation object containing the fluid to be mixed.

        sample_name : str
            The name of the fluid sample that will be mixed with all other
            speciated fluids.

        amount_remaining : float, default 1
            Number of moles of the fluid to be mixed with all others.

        amount_destroyed : float, default 0
            Number of moles of the mixing fluid that has been destroyed.

        molar_volume : float, default 1
            Molar volume of the mixing fluid, in moles/cm3.

        mass_ratio : float, default 1
            Ratio of mass of the mixing fluid to all other fluids.

        hide_traceback : bool, default True
            Hide traceback message when encountering errors handled by this
            class? When True, error messages handled by this class will be short
            and to the point.
            
        """
        
        self.err_handler = Error_Handler(clean=hide_traceback)

        # Prepare a special reactant to be used in a mixing calculation.
        if isinstance(speciation, Speciation):
            
            self.sample_name = sample_name

            if not sample_name in speciation.sample_data.keys():
                
                self.err_handler.raise_exception(("The sample '"+str(sample_name)+"'"
                        " was not found amongst the samples in this speciation"
                        " calculation: "+str(list(speciation.sample_data.keys()))))
                
            self.speciation_sample_data = speciation.sample_data[sample_name]
            self.T = self.speciation_sample_data["temperature"]
            self.mass_ratio = mass_ratio
            
            elemental_composition_lines = []
            capture = False
            for i,line in enumerate(speciation.raw_3_pickup_dict_top[sample_name]):
                if "|->|Composition" in line:
                    capture = True
                    i_start = i
                if "|->|Reaction" in line:
                    capture = False
                if capture and i > i_start + 3:
                    elemental_composition_lines.append(line)
            
            # ignore last line, which is a divider "|----------..."
            elemental_composition_lines = elemental_composition_lines[:-1]

            fluid_2_dict = {}
            for line in elemental_composition_lines:
                split_line = line.split("|")
                element = split_line[2].strip()
                value = float(split_line[3])
                fluid_2_dict[element] = value
            
            self.reactant = Reactant(reactant_name="Fluid 2",
                                     reactant_type="Special reactant",
                                     special_reactant_dict=fluid_2_dict,
                                     amount_remaining=amount_remaining,
                                     amount_destroyed=amount_destroyed,
                                     molar_volume=molar_volume,
                                     hide_traceback=hide_traceback)
            
            self.formatted_block = self.reactant.formatted_block
            self.reactant_type = "Special reactant"
            
        else:
            self.err_handler.raise_exception(("The speciation parameter was"
                    " not given a Speciation object. A Speciation object is"
                    " produced by the AqEquil.speciate() function."))

        # handle fluid mixing reaction block

        lines_to_keep = self.formatted_block.split("\n")
        raw_p_dict_top = speciation.raw_3_pickup_dict_top

        # grab this fluid's reaction block
        reaction_block_lines = []
        capture = False
        for line in raw_p_dict_top[sample_name]:
            if "|->|Reaction" in line:
                capture = True
            if "|->|Surface area" in line:
                capture = False
            if capture:
                reaction_block_lines.append(line)

        # insert this fluid's reaction block
        before_lines = []
        after_lines = []
        is_before = True
        for i,line in enumerate(lines_to_keep):
            if "|->|Reaction" in line:
                is_before=False
            if is_before:
                before_lines.append(line)
            else:
                after_lines.append(line)

        # trim redundant lines
        after_lines = after_lines[5:]

        lines_to_keep = before_lines + reaction_block_lines + after_lines

        self.formatted_block = "\n".join(lines_to_keep)

class Reactant:
    def __init__(self,
                 reactant_name,
                 reactant_type="Pure mineral",
                 reactant_status="Reacting",
                 amount_remaining=1,
                 amount_destroyed=0,
                 molar_volume=0,
                 surface_area_option=0,
                 surface_area_value=0,
                 surface_area_factor=0,
                 f_rate_law="Relative rate equation",
                 f_eq1=1,
                 f_eq2=0,
                 f_eq3=0,
                 b_rate_law="Partial equilibrium",
                 b_eq1=1,
                 b_eq2=0,
                 b_eq3=0,
                 special_reactant_dict={},
                 hide_traceback=True,
                 ):
        
        """
        Class used to define reactants for `Prepare_Reaction`.

        Parameters
        ----------
        reactant_name : str
            Name of the reactant.

        reactant_type : str, default "Pure mineral"
            Reactant type. Valid types include:
            
            - Pure mineral
            - Solid solution
            - Special reactant
            - Aqueous species
            - Gas species
            - Generic ion exchanger

        reactant_status : str, default "Reacting"
            Status of the reactant. Valid statuses include:
            
            - Reacting
            - Saturated, reacting
            - Exhausted
            - Saturated, not reacting

        amount_remaining : float, default 1
            Moles of reactant remaining.

        amount_destroyed : float, default 0
            Moles of reactant destroyed.

        surface_area_option : int, default 0
            Option for reactant surface area. Valid options include:
            
            - 0, for Constant surface area (in cm2)
            - 1, for Constant specific surface area (in cm2/g). Surface area
            changes in proportion to the reactant mass.
            - 2, for n**2/3 growth law- current surface area (in cm2)

        surface_area_value : float, default 0
            The value assigned to choice of surface area option according to:
            
            - If `surface_area_option` is 0, value is in cm2.
            - If `surface_area_option` is 1, value is in cm2/g.
            - If `surface_area_option` is 2, value is in cm2.

        surface_area_factor : float, default 0
            Value assigned to surface area factor.

        f_rate_law : str, default "Relative rate equation"
            Type of forward rate law. Valid types include:
            
            - Use backward rate law
            - Relative rate equation
            The TST rate equation is not yet supported.

        f_eq1, f_eq2, f_eq3 : float, default 1, 0, 0, respectively
            Coefficients of the forward rate law defined for `f_rate_law`.
            If `f_rate_law` is "Relative rate equation", then:
            
            - f_eq1 is dXi(n)/dXi (mol/mol)
            - f_eq2 is d2Xi(n)/dXi2 (mol/mol2)
            - f_eq3 is d3Xi(n)/dXi3 (mol/mol3)

        b_rate_law : str, default "Partial equilibrium"
            Type of backward rate law. Valid types include:
            
            - Use forward rate law
            - Partial equilibrium
            - Relative rate equation
            The TST rate equation is not yet supported.
            
        b_eq1, b_eq2, b_eq3 : float, default 1, 0, 0, respectively
            Coefficients of the backward rate law defined for `b_rate_law`.
            If `b_rate_law` is "Relative rate equation", then:
            
            - the value of b_eq1 represents dXi(n)/dXi (mol/mol)
            - the value of b_eq2 represents d2Xi(n)/dXi2 (mol/mol2)
            - the value of b_eq3 represents d3Xi(n)/dXi3 (mol/mol3)
        
        hide_traceback : bool, default True
            Hide traceback message when encountering errors handled by this class?
            When True, error messages handled by this class will be short and to
            the point.

        """
    
        self.err_handler = Error_Handler(clean=hide_traceback)
        
        if f_rate_law not in ["Use backward rate law", "Relative rate equation"]:
            self.err_handler.raise_exception(("f_rate_law must be either "
                    "'Use backward rate law' or 'Relative rate equation'."))
        if b_rate_law not in ["Use forward rate law", "Partial equilibrium", "Relative rate equation"]:
            self.err_handler.raise_exception(("b_rate_law must be either "
                    "'Use forward rate law', 'Partial equilibrium', or 'Relative rate equation'."))
        
        self.reactant_name=reactant_name
        self.reactant_type=reactant_type
        self.reactant_status=reactant_status
        self.amount_remaining=amount_remaining
        self.amount_destroyed=amount_destroyed
        self.molar_volume=molar_volume
        self.sa_val_1=0
        self.sa_val_2=0
        self.sa_val_3=0
        self.sa_checkbox_1= " "
        self.sa_checkbox_2= " "
        self.sa_checkbox_3= " "
        
        self.special_reactant_dict=special_reactant_dict
        
        if surface_area_option == 0:
            self.sa_checkbox_1 = "x"
            self.sa_val_1 = surface_area_value
        elif surface_area_option == 1:
            self.sa_checkbox_2 = "x"
            self.sa_val_2 = surface_area_value
        elif surface_area_option == 2:
            self.sa_checkbox_3 = "x"
            self.sa_val_3 = surface_area_value
        self.sa_factor=surface_area_factor
        
        self.f_rate_law=f_rate_law
        self.f_eq1=f_eq1
        self.f_eq2=f_eq2
        self.f_eq3=f_eq3
        self.b_rate_law=b_rate_law
        self.b_eq1=b_eq1
        self.b_eq2=b_eq2
        self.b_eq3=b_eq3
        
        self.__format_rate_block("forward")
        self.__format_rate_block("backward")

        self.__format_block()


    def __format_block(self):
        
        rb_dict_formatted = dict(
                reactant_name = f"{self.reactant_name:<24}",
                reactant_type = f"{self.reactant_type:<24}",
                reactant_status = f"{self.reactant_status:<24}",
                amount_remaining = f"{'{:.5E}'.format(self.amount_remaining):>12}",
                amount_destroyed = f"{'{:.5E}'.format(self.amount_destroyed):>12}",
                sa_checkbox_1 = self.sa_checkbox_1,
                sa_checkbox_2 = self.sa_checkbox_2,
                sa_checkbox_3 = self.sa_checkbox_3,
                sa_val_1 = f"{'{:.5E}'.format(self.sa_val_1):>12}",
                sa_val_2 = f"{'{:.5E}'.format(self.sa_val_2):>12}",
                sa_val_3 = f"{'{:.5E}'.format(self.sa_val_3):>12}",
                sa_factor = f"{'{:.5E}'.format(self.sa_factor):>12}",
                f_rate_law = f"{self.f_rate_law:<24}",
                b_rate_law = f"{self.b_rate_law:<24}",
                f_rate_block = self.formatted_f_rate_block,
                b_rate_block = self.formatted_b_rate_block,
                )

        if self.reactant_type.lower() == "special reactant":
            elem_lines = []
            for key in list(self.special_reactant_dict.keys()):
                sp_elem_line = copy.copy(srb_stoich_template)
                elem_name = f"{str(key):<8}"
                elem_value = f"{'{:.15E}'.format(float(self.special_reactant_dict[key])):>22}"
                sp_elem_line = sp_elem_line.format(**{"elem":elem_name, "elem_val":elem_value})
                elem_lines.append(sp_elem_line)
            elem_lines = "".join(elem_lines)
            
            srb_dict_formatted = dict(
                    molar_volume = f"{'{:.5E}'.format(self.molar_volume):>12}",
                    srb_stoich = elem_lines,
                    )
            
            sp_reactant_block_template = copy.copy(srb_template)
            sp_reactant_block_formatted = sp_reactant_block_template.format(**srb_dict_formatted)
            rb_dict_formatted["sr_block"] = sp_reactant_block_formatted
        else:
            rb_dict_formatted["sr_block"] = ""
   
        reactant_block_template = copy.copy(rb_template)
        self.formatted_block = reactant_block_template.format(**rb_dict_formatted)
    
    def __format_rate_block(self, direction):
        if direction == "forward":
            d="f"
            od="b"
            odirection="backward"
        else:
            d="b"
            od="f"
            odirection="forward"

        eq1 = getattr(self, d+"_eq1")
        eq2 = getattr(self, d+"_eq2")
        eq3 = getattr(self, d+"_eq3")

        if getattr(self, d+"_rate_law") == "Relative rate equation":
            rate_block = copy.deepcopy(relative_rate_equation_template)
            
            rate_options_formatted = dict(
                eq1=f"{'{:.5E}'.format(eq1):>12}",
                eq2=f"{'{:.5E}'.format(eq2):>12}",
                eq3=f"{'{:.5E}'.format(eq3):>12}",
            )
            
            rate_block = rate_block.format(**rate_options_formatted)
            setattr(self, "formatted_"+d+"_rate_block", rate_block)
        elif getattr(self, d+"_rate_law") == "Linear rate equation":
            # Note: "Linear rate equation" is valid in EQ6, but I cannot find any
            # EQ6 input files formatted for this option. My guess is that linear
            # rate equations can be imposed by only specifying the first parameter
            # of a relative rate equation.
            self.err_handler.raise_exception("The 'Linear rate equation' option is not supported at this time.")
        elif getattr(self, d+"_rate_law") == "TST rate equation":
            self.err_handler.raise_exception("TST rate equations are not supported at this time.")
        elif getattr(self, d+"_rate_law") == "Partial equilibrium":
            setattr(self, "formatted_"+d+"_rate_block", "")
        elif getattr(self, d+"_rate_law") == "Use backward rate law":
            setattr(self, "formatted_"+d+"_rate_block", "")
        elif getattr(self, d+"_rate_law") == "Use forward rate law":
            setattr(self, "formatted_"+d+"_rate_block", "")
        else:
            msg = ("Valid rate laws for the "+direction+" rate law includes "
                "'Relative rate equation', 'TST rate equation', "
                "'Linear rate equation', or 'Use "+odirection+" rate law'")
            self.err_handler.raise_exception(msg)

class Gas:
    def __init__(self,
                 gas_name="None",
                 gas_moles=0,
                 gas_log_fugacity=0,
                ):
        
        """
        Class used to define gases for `Prepare_Reaction`.

        Parameters
        ----------
        gas_name : str, default "None"
            Name of the gas.
        
        gas_moles : float, default 0
            Moles of gas.
        
        gas_log_fugacity : float, default 0
            The log10 fugacity of the gas.

        """
        
        self.gas_name=gas_name
        self.gas_moles=gas_moles
        self.gas_log_fugacity=gas_log_fugacity
        
        self.__format_line()
        
    def __format_line(self):
        
        gas_options_formatted = dict(
            gas_name = f"{self.gas_name:<24}",
            gas_moles = f"{'{:.5E}'.format(self.gas_moles):>13}",
            gas_log_fugacity = f"{'{:.5E}'.format(self.gas_log_fugacity):>13}",
        )
    
        gas_line_template = copy.copy(gl_template)
        self.formatted_line = gas_line_template.format(**gas_options_formatted)
        

class Prepare_Reaction:
    def __init__(self,
                 reactants,
                 gases=[],
                 t_option=None,
                 t_value_1=None, # temp
                 t_value_2=None,  # temp or deriv
                 t_value_3=None,  # mass ratio factor
                 p_option=0,
                 p_value_1=None,  # pressure
                 p_value_2=0,  # deriv
                 xi_range=[0, 1],
                 time_range=[0, 1e38],
                 pH_range=[-1e38, 1e38],
                 Eh_range=[-1e38, 1e38],
                 fO2_range=[-1e38, 1e38],
                 aw_range=[-1e38, 1e38],
                 max_n_steps=900,
                 xi_print_int=1,
                 log_xi_print_int=1,
                 time_print_int=1e38,
                 log_time_print_int=1e38,
                 pH_print_interval=1e38,
                 Eh_print_interval=1e38,
                 logfO2_print_interval=1e38,
                 aw_print_interval=1e38,
                 n_steps_print_interval=100,
                 physical_system_model="closed",
                 kinetic_mode="arbitrary",
                 phase_boundary_search=0,
                 permit_solid_solutions=False,
                 clear_es_solids_read=False,
                 clear_es_solids_initial=False,
                 clear_es_solids_end=False,
                 clear_prs_solids_read=False,
                 clear_prs_solids_end=False,
                 auto_basis_switching_pre_NR=False,
                 auto_basis_switching_post_NR=False,
                 calc_mode_selection=0,
                 ODE_corrector_mode=0,
                 mineral_suppression_option="None",
                 write_tab=-1, # do not write a TAB file by default because EQ6 can encounter an access violation when writing a TAB file
                 fluid_mixing_setup=False,
                 max_finite_difference_order=6,
                 beta_convergence_tolerance=0,
                 del_convergence_tolerance=0,
                 max_n_NR_iter=500,
                 search_find_convergeance_tolerance=0,
                 saturation_tolerance=0,
                 max_n_phase_assemblage_tries=0,
                 zero_order_step_size=0,
                 max_interval_in_xi_between_PRS_transfers=0,
                 filename=None,
                 hide_traceback=True,
                ):
        
        """
        Class used to set the parameters of a reaction between the results of a
        speciation calculation and minerals, gases, etc.

        Parameters
        ----------
        reactants : list
            List of reactants defined by the 'Reactant' class. Can be an empty
            list if no reactants are desired.
             
        gases : list, default []
            List of gases defined by the 'Gas' class. Can be an empty
            list if no gases are desired.
         
        t_option : int, default 0
            Desired option for handling temperature of the reaction. Valid
            choices include:
            
            - 0 for constant temperature
            - 1 for linear tracking in Xi
            - 2 for linear tracking in time
            
            Fluid mixing tracking is not yet supported.
         
        t_value_1, t_value_2, t_value_3 : default None, 0, 0, respectively
            By default, the temperature of samples in the speciation
            calculation will be used, so the user does not need to specify
            temperature values here. However, if a user wishes, temperature
            values can be defined here that will be applied to all samples in
            the speciation. Note that doing so may result in incongruous
            results. That said, the values specified here depend on which
            option is selected for `t_option`:
            
            - If 't_option' is 0, then t_value_1 is the value of the constant
            temperature (in degrees C), and t_value_2 and t_value_3 are ignored.
            - If 't_option' is 1, then t_value_1 is the base value
            temperature (in degrees C), t_value_2 is the derivative, and
            t_value_3 is ignored.
            - If 't_option' is 2, then t_value_1 is the base value
            temperature (in degrees C), t_value_2 is the derivative, and
            t_value_3 is ignored.
         
        p_option : int, default 0
            Desired option for handling pressure of the reaction. Valid
            choices include:
            
            - 0 to follow the data file reference pressure curve
            - 1 to follow the 1.013-bar/steam-saturation curve
            - 2 for constant pressure
            - 3 for linear tracking in Xi
            - 4 for linear tracking in time
        
        p_value_1, p_value_2 : float default None, 0
            Values assigned to desired `p_option`.
            
            - If `p_option` is 0 or 1, p_value_1 and p_value_2 are ignored.
            - If `p_option` is 2, p_value_1 represents a constant pressure,
            in bars, and p_value_2 is ignored.
            - If `p_option` is 3 or 4, p_value_1 represents the base pressure
            value, in bars, and p_value_2 is the derivative.
         
        xi_range : list of two float, default [0, 1]
            A list containing the starting and maximum value of Xi,
            respectively.
         
        time_range : list of two float, default [0, 1e38]
            A list containing the starting and maximum time, respectively.
         
        pH_range : list of two float, default [-1e38, 1e38]
            A list containing the minimum and maximum values of pH.
         
        Eh_range : list of two float, default [-1e38, 1e38]
            A list containing the minimum and maximum values of Eh.
         
        fO2_range : list of two float, default [-1e38, 1e38]
            A list containing the minimum and maximum values of the fugacity
            of oxygen, fO2.
         
        aw_range : list of two float, default [-1e38, 1e38]
            A list containing the minimum and maximum values of water activity.
         
        max_n_steps : int, default 900
            Maximum number of steps of Xi allowed.
         
        xi_print_int : int, default 1
            Xi print interval.
         
        log_xi_print_int : int, default 1
            Log Xi print interval.
         
        time_print_int : int, default 1e38
            Time print interval.
         
        log_time_print_int : int, default 1e38
            Log time print interval.
         
        pH_print_interval : int, default 1e38
            pH print interval.
         
        Eh_print_interval : int, default 1e38
            Eh (v) print interval.
         
        logfO2_print_interval : int, default 1e38
            Log fO2 print interval.
         
        aw_print_interval : int, default 1e38
            Activity of water (aw) print interval.
         
        n_steps_print_interval : int, default 100
            Steps print interval.

        physical_system_model : str, default "closed"
            Selection for the physical system model. Valid options include:
            
            - "closed"
            - "titration"
            - "fluid-centered flow-through open"
        
        kinetic_mode : str, default "arbitrary"
            Selection for kinetic mode. Valid options include:
            - "arbitrary" for arbitrary kinetics, reaction progress mode
            - "true" for true kinetics, reaction progress/time mode 
         
        phase_boundary_search : int, default 0
            Selection for phase boundary searches. Valid options include:
            
            - 0 to search for phase boundaries and constrain the step size to
            match.
            - 1 to search for phase boundaries and print their locations.
            - 2 to not search for phase boundaries.
        
        permit_solid_solutions : bool, default False
            Permit solid solutions? If False, solid solutions are ignored.
        
        clear_es_solids_read : bool, default False
            Clear the ES solids read from the input file?
        
        clear_es_solids_initial : bool, default False
            Clear the ES solids at the initial value of reaction progress?
        
        clear_es_solids_end : bool, default False
            Clear the ES solids at the end of the run?
        
        clear_prs_solids_read : bool, default False
            Clear the PRS solids read from the input file?
        
        clear_prs_solids_end : bool, default False
            Clear the PRS solids at the end of the run? If True, PRS solids will
            be cleared unless numerical problems cause early termination.
        
        auto_basis_switching_pre_NR : bool, default False
            Turn on auto basis switching in pre-Newton-Raphson optimization?
        
        auto_basis_switching_post_NR : bool, default False
            Turn on auto basis switching after Newton-Raphson iteration?
        
        calc_mode_selection : int, default 0
            Calculational mode selection. Valid options include:
            
            - 0 for normal path tracing
            - 1 for economy mode (if permissible)
            - 2 for super economy mode (if permissible)
        
        ODE_corrector_mode : int, default 0
            ODE integrator corrector mode selection. Valid options include:
            
            - 0 to allow stiff and simple correctors
            - 1 to allow only simple corrector
            - 2 to allow only stiff corrector
            - 3 to allow no correctors
        
        mineral_suppression_option : str, default "None"
            Option to suppress formation of minerals. Can be either "None" (no
            minerals are suppressed) or "All" (all minerals are suppressed).
        
        write_tab : int, default -1
            Option to write a TAB file. Valid options include:
            
            - -1 do not write a TAB file (default)
            - 0 write a TAB file
            - 1 write a TAB file prepending TABX file data from a previous run
        
        fluid_mixing_setup : bool, default False
            If True, will write an EQ6 input file with Fluid 1 set up for
            fluid mixing. If False, a normal EQ6 pickup file will be written.
         
        max_finite_difference_order : int, default 6
            Maximum finite-difference order (numerical parameter).
         
        beta_convergence_tolerance : float, default 0
            Beta convergence tolerance (numerical parameter).
         
        del_convergence_tolerance : float, default 0
            Delta convergence tolerance (numerical parameter).
         
        max_n_NR_iter : int, default 500
            Maximum number of N-R iterations (numerical parameter).
         
        search_find_convergeance_tolerance : float, default 0
            Search/find convergence tolerance (numerical parameter).
         
        saturation_tolerance : float, default 0
            Saturation tolerance (numerical parameter).
         
        max_n_phase_assemblage_tries : int, default 0
            Maximum number of phase assemblage tries (numerical parameter).
         
        zero_order_step_size : int, default 0
            Zero order step size in Xi (numerical parameter).
         
        max_interval_in_xi_between_PRS_transfers : int, default 0
            Maximum interval in Xi between PRS transfers (numerical parameter).
         
        filename : str, default None
            Filename where the results of `Prepare_Reaction` will be written.
            This is equivalent to the top half of an EQ3/6 6i file. If None,
            no file will be written.
             
        hide_traceback : bool, default True
            Hide traceback message when encountering errors handled by this
            class? When True, error messages handled by this class will be short
            and to the point.

        """
        
        self.err_handler = Error_Handler(clean=hide_traceback)
        
        if len(reactants) == 1 and not isinstance(reactants, list):
            reactants = list(reactants)
            
        if len(gases) == 0:
            gases = [Gas()]
        elif len(gases) == 1 and not isinstance(gases, list):
            gases = list(gases)
        
        self.reactants=reactants
        self.gases=gases
        self.t_option=t_option
        self.t_value_1=t_value_1
        self.t_value_2=t_value_2
        self.t_value_3=t_value_3
        self.p_option=p_option
        self.p_value_1=p_value_1
        self.p_value_2=p_value_2
        self.start_xi=xi_range[0]
        self.max_xi=xi_range[1]
        self.start_time=time_range[0]
        self.max_time=time_range[1]
        self.min_pH=pH_range[0]
        self.max_pH=pH_range[1]
        self.min_Eh=Eh_range[0]
        self.max_Eh=Eh_range[1]
        self.min_fO2=fO2_range[0]
        self.max_fO2=fO2_range[1]
        self.min_aw=aw_range[0]
        self.max_aw=aw_range[1]
        self.max_n_steps=max_n_steps
        self.xi_print_int=xi_print_int
        self.log_xi_print_int=log_xi_print_int
        self.time_print_int=time_print_int
        self.log_time_print_int=log_time_print_int
        self.pH_print_interval=pH_print_interval
        self.Eh_print_interval=Eh_print_interval
        self.logfO2_print_interval=logfO2_print_interval
        self.aw_print_interval=aw_print_interval
        self.n_steps_print_interval=n_steps_print_interval
        self.max_finite_difference_order=max_finite_difference_order
        self.beta_convergence_tolerance=beta_convergence_tolerance
        self.del_convergence_tolerance=del_convergence_tolerance
        self.max_n_NR_iter=max_n_NR_iter
        self.search_find_convergeance_tolerance=search_find_convergeance_tolerance
        self.saturation_tolerance=saturation_tolerance
        self.max_n_phase_assemblage_tries=max_n_phase_assemblage_tries
        self.zero_order_step_size=zero_order_step_size
        self.max_interval_in_xi_between_PRS_transfers=max_interval_in_xi_between_PRS_transfers
        self.mineral_suppression_option=mineral_suppression_option
        
        self.t_checkbox_1=" "
        self.t_checkbox_2=" "
        self.t_checkbox_3=" "
        self.t_checkbox_4=" "
        self.tval1=0
        self.tval2=0
        self.tval3=0
        self.tval4=0
        self.tval5=0
        self.tval6=0
        self.tval7=0
        self.tval8=0
        self.p_checkbox_1=" "
        self.p_checkbox_2=" "
        self.p_checkbox_3=" "
        self.p_checkbox_4=" "
        self.p_checkbox_5=" "
        self.pval1=0
        self.pval2=0
        self.pval3=0
        self.pval4=0
        self.pval5=0
        self.i1_checkbox_1 = " "
        self.i1_checkbox_2 = " "
        self.i1_checkbox_3 = " "
        self.i2_checkbox_1 = " "
        self.i2_checkbox_2 = " "
        self.i3_checkbox_1 = " "
        self.i3_checkbox_2 = " "
        self.i3_checkbox_3 = " "
        self.i4_checkbox_1 = " "
        self.i4_checkbox_2 = " "
        self.i5_checkbox_1 = " "
        self.i5_checkbox_2 = " "
        self.i6_checkbox_1 = " "
        self.i6_checkbox_2 = " "
        self.i7_checkbox_1 = " "
        self.i7_checkbox_2 = " "
        self.i9_checkbox_1 = " "
        self.i9_checkbox_2 = " "
        self.i10_checkbox_1 = " "
        self.i10_checkbox_2 = " "
        self.i11_checkbox_1 = " "
        self.i11_checkbox_2 = " "
        self.i12_checkbox_1 = " "
        self.i12_checkbox_2 = " "
        self.i13_checkbox_1 = " "
        self.i13_checkbox_2 = " "
        self.i13_checkbox_3 = " "
        self.i14_checkbox_1 = " "
        self.i14_checkbox_2 = " "
        self.i14_checkbox_3 = " "
        self.i14_checkbox_4 = " "
        self.i15_checkbox_1 = " "
        self.i15_checkbox_2 = " "
        self.i18_checkbox_1 = " "
        self.i18_checkbox_2 = " "
        self.i18_checkbox_3 = " "
        self.i20_checkbox_1 = " "
        self.i20_checkbox_2 = " "
        
        tval_var_to_format = None
        pval_var_to_format = None
        
        # set t_option and t_value defaults when there is a mixing calculation
        n_mixing_fluid_reactants=0
        for reactant in reactants:
            if isinstance(reactant, Mixing_Fluid):
                n_mixing_fluid_reactants += 1
                if t_option == None:
                    t_option = 3
                    self.t_option=t_option
                if t_value_1 == None:
                    t_value_1 = None # will be formatted with temp of fluid 1 later
                    self.t_value_1=t_value_1
                if t_value_2 == None:
                    t_value_2 = float(reactant.T) # temp of fluid 2
                    self.t_value_2=t_value_2
                if t_value_3 == None:
                    t_value_3 = reactant.mass_ratio # mass ratio factor
                    self.t_value_3=t_value_3
                
        if n_mixing_fluid_reactants == 0:
            # set t_option and t_value defaults when there is no mixing calculation
            if t_option == None:
                t_option = 0
                self.t_option = 0
#             if t_value_1 == None:
#                 t_value_1 = None
#                 self.t_value_1 = 0
#             if t_value_2 == None:
#                 t_value_2 = 0
#                 self.t_value_2 = 0
#             if t_value_3 == None:
#                 t_value_3 = 0
#                 self.t_value_3 = 0
        elif n_mixing_fluid_reactants == 1:
            pass
        else:
            self.err_handler.raise_exception((""
                    "There are {} mixing fluids ".format(n_mixing_fluid_reactants)+""
                    "in the list of reactants. There may only be one."))
        if t_option == 0:
            self.t_checkbox_1="x"
            if isinstance(t_value_1, numbers.Number):
                self.tval1=t_value_1
            else:
                tval_var_to_format = 1

        elif t_option == 1:
            self.t_checkbox_2="x"
            if isinstance(t_value_1, numbers.Number):
                self.tval2=t_value_1
                self.tval3=t_value_2
            else:
                tval_var_to_format = 2

        elif t_option == 2:
            self.t_checkbox_3="x"
            if isinstance(t_value_1, numbers.Number):
                self.tval4=t_value_1
                self.tval5=t_value_2
            else:
                tval_var_to_format = 4
                
        elif t_option == 3:
            self.t_checkbox_4="x"
            if isinstance(t_value_1, numbers.Number):
                self.tval6=t_value_1
                self.tval7=t_value_2
                self.tval8=t_value_3
            else:
                self.tval7=t_value_2
                self.tval8=t_value_3
                tval_var_to_format = 6
        else:
            raise Exception("t_option must be 0, 1, 2, or 3.")
            
        if p_option == 0:
            self.p_checkbox_1="x"

        elif p_option == 1:
            self.p_checkbox_2="x"

        elif p_option == 2:
            self.p_checkbox_3="x"
            if isinstance(p_value_1, numbers.Number):
                self.pval1=p_value_1
            else:
                pval_var_to_format = 1
                
        elif p_option == 3:
            self.p_checkbox_4="x"
            if isinstance(p_value_1, numbers.Number):
                self.pval2=p_value_1
                self.pval3=p_value_2
            else:
                pval_var_to_format = 2
                
        elif p_option == 4:
            self.p_checkbox_5="x"
            if isinstance(p_value_1, numbers.Number):
                self.pval4=p_value_1
                self.pval5=p_value_2
            else:
                pval_var_to_format = 4

        else:
            raise Exception("p_option must be 0, 1, 2, 3, or 4.")
        
        
        if physical_system_model == "closed":
            self.i1_checkbox_1 = "x"
        elif physical_system_model == "titration":
            self.i1_checkbox_2 = "x"
        elif physical_system_model == "fluid-centered flow-through open":
            self.i1_checkbox_3 = "x"
        else:
            msg = ("physical_system_model must either be 'closed', 'titration',"
                " or 'fluid-centered flow-through open'.")
            self.err_handler.raise_exception(msg)

        if kinetic_mode == "arbitrary":
            self.i2_checkbox_1 = "x"
        elif kinetic_mode == "true":
            self.i2_checkbox_2 = "x"
        else:
            msg = "kinetic_mode must either be 'arbitrary' or 'true'."
            self.err_handler.raise_exception(msg)
        
        if phase_boundary_search == 0:
            self.i3_checkbox_1 = "x"
        elif phase_boundary_search == 1:
            self.i3_checkbox_2 = "x"
        elif phase_boundary_search == 2:
            self.i3_checkbox_3 = "x"
        else:
            msg = "phase_boundary_search must be 0, 1, or 2."
            self.err_handler.raise_exception(msg)
            
        if permit_solid_solutions == False:
            self.i4_checkbox_1 = "x"
        elif permit_solid_solutions == True:
            self.i4_checkbox_2 = "x"
        else:
            msg = "permit_solid_solutions must be True or False."
            self.err_handler.raise_exception(msg)
            
        if clear_es_solids_read == False:
            self.i5_checkbox_1 = "x"
        elif clear_es_solids_read == True:
            self.i5_checkbox_2 = "x"
        else:
            msg = "clear_es_solids_read must be True or False."
            self.err_handler.raise_exception(msg)

        if clear_es_solids_initial == False:
            self.i6_checkbox_1 = "x"
        elif clear_es_solids_initial == True:
            self.i6_checkbox_2 = "x"
        else:
            msg = "clear_es_solids_initial must be True or False."
            self.err_handler.raise_exception(msg)

        if clear_es_solids_end == False:
            self.i7_checkbox_1 = "x"
        elif clear_es_solids_end == True:
            self.i7_checkbox_2 = "x"
        else:
            msg = "clear_es_solids_end must be True or False."
            self.err_handler.raise_exception(msg)
            
        # there is no iopt8
            
        if clear_prs_solids_read == False:
            self.i9_checkbox_1 = "x"
        elif clear_prs_solids_read == True:
            self.i9_checkbox_2 = "x"
        else:
            msg = "clear_prs_solids_read must be True or False."
            self.err_handler.raise_exception(msg)
            
        if clear_prs_solids_end == False:
            self.i10_checkbox_1 = "x"
        elif clear_prs_solids_end == True:
            self.i10_checkbox_2 = "x"
        else:
            msg = "clear_prs_solids_end must be True or False."
            self.err_handler.raise_exception(msg)
            
        if auto_basis_switching_pre_NR == False:
            self.i11_checkbox_1 = "x"
        elif auto_basis_switching_pre_NR == True:
            self.i11_checkbox_2 = "x"
        else:
            msg = "auto_basis_switching_pre_NR must be True or False."
            self.err_handler.raise_exception(msg)
            
        if auto_basis_switching_post_NR == False:
            self.i12_checkbox_1 = "x"
        elif auto_basis_switching_post_NR == True:
            self.i12_checkbox_2 = "x"
        else:
            msg = "auto_basis_switching_post_NR must be True or False."
            self.err_handler.raise_exception(msg)
            
        if calc_mode_selection == 0:
            self.i13_checkbox_1 = "x"
        elif calc_mode_selection == 1:
            self.i13_checkbox_2 = "x"
        elif calc_mode_selection == 2:
            self.i13_checkbox_3 = "x"
        else:
            msg = "calc_mode_selection must be 0, 1, or 2."
            self.err_handler.raise_exception(msg)

        if ODE_corrector_mode == 0:
            self.i14_checkbox_1 = "x"
        elif ODE_corrector_mode == 1:
            self.i14_checkbox_2 = "x"
        elif ODE_corrector_mode == 2:
            self.i14_checkbox_3 = "x"
        elif ODE_corrector_mode == 3:
            self.i14_checkbox_4 = "x"
        else:
            msg = "ODE_corrector_mode must be 0, 1, or 2."
            self.err_handler.raise_exception(msg)

        # The suppress redox option in EQ6 doesn't work with WORM data.
        # Use redox-isolated elements from the redox suppression
        # option in the AqEquil class instead.
        self.i15_checkbox_1 = "x"

        if write_tab == -1:
            self.i18_checkbox_1 = "x"
        elif write_tab == 0:
            self.i18_checkbox_2 = "x"
        elif write_tab == 1:
            self.i18_checkbox_3 = "x"
        else:
            self.err_handler.raise_exception("write_tab_option must be -1, 0, "
                    "or 1.")
            
        
        if fluid_mixing_setup == False:
            self.i20_checkbox_1 = "x"
        elif fluid_mixing_setup == True:
            self.i20_checkbox_2 = "x"
        else:
            msg = "fluid_mixing_setup must be True or False."
            self.err_handler.raise_exception(msg)

        now = datetime.now()
        self.date_created = now.strftime('%Y-%m-%d %I:%M %p')
        
        self.__format_reaction(tval_var_to_format, pval_var_to_format)
        
        if filename != None:
            with open(filename, 'w') as f:
                f.write(pr.formatted_reaction)
        
        
    def __format_reaction(self, tval_var_to_format=None, pval_var_to_format=None):
        reaction_options_formatted = dict(
            date_created=f"{self.date_created:<24}",
            start_xi=f"{'{:.5E}'.format(self.start_xi):>12}",
            max_xi=f"{'{:.5E}'.format(self.max_xi):>12}",
            start_time=f"{'{:.5E}'.format(self.start_time):>12}",
            max_time=f"{'{:.5E}'.format(self.max_time):>12}",
            min_pH=f"{'{:.5E}'.format(self.min_pH):>12}",
            max_pH=f"{'{:.5E}'.format(self.max_pH):>12}",
            min_Eh=f"{'{:.5E}'.format(self.min_Eh):>12}",
            max_Eh=f"{'{:.5E}'.format(self.max_Eh):>12}",
            min_fO2=f"{'{:.5E}'.format(self.min_fO2):>12}",
            max_fO2=f"{'{:.5E}'.format(self.max_fO2):>12}",
            min_aw=f"{'{:.5E}'.format(self.min_aw):>12}",
            max_aw=f"{'{:.5E}'.format(self.max_aw):>12}",
            max_n_steps=f"{self.max_n_steps:>12}",
            xi_print_int=f"{'{:.5E}'.format(self.xi_print_int):>12}",
            log_xi_print_int=f"{'{:.5E}'.format(self.log_xi_print_int):>12}",
            time_print_int=f"{'{:.5E}'.format(self.time_print_int):>12}",
            log_time_print_int=f"{'{:.5E}'.format(self.log_time_print_int):>12}",
            pH_print_interval=f"{'{:.5E}'.format(self.pH_print_interval):>12}",
            Eh_print_interval=f"{'{:.5E}'.format(self.Eh_print_interval):>12}",
            logfO2_print_interval=f"{'{:.5E}'.format(self.logfO2_print_interval):>12}",
            aw_print_interval=f"{'{:.5E}'.format(self.aw_print_interval):>12}",
            n_steps_print_interval=f"{self.n_steps_print_interval:>12}",
            max_finite_difference_order=f"{self.max_finite_difference_order:>12}",
            beta_convergence_tolerance=f"{self.beta_convergence_tolerance:>12}",
            del_convergence_tolerance=f"{self.del_convergence_tolerance:>12}",
            max_n_NR_iter=f"{self.max_n_NR_iter:>12}",
            search_find_convergeance_tolerance=f"{self.search_find_convergeance_tolerance:>12}",
            saturation_tolerance=f"{self.saturation_tolerance:>12}",
            max_n_phase_assemblage_tries=f"{self.max_n_phase_assemblage_tries:>12}",
            zero_order_step_size=f"{self.zero_order_step_size:>12}",
            max_interval_in_xi_between_PRS_transfers=f"{self.max_interval_in_xi_between_PRS_transfers:>12}",
            t_checkbox_1=self.t_checkbox_1,
            t_checkbox_2=self.t_checkbox_2,
            t_checkbox_3=self.t_checkbox_3,
            t_checkbox_4=self.t_checkbox_4,
            p_checkbox_1=self.p_checkbox_1,
            p_checkbox_2=self.p_checkbox_2,
            p_checkbox_3=self.p_checkbox_3,
            p_checkbox_4=self.p_checkbox_4,
            p_checkbox_5=self.p_checkbox_5,
            i1_checkbox_1=self.i1_checkbox_1,
            i1_checkbox_2=self.i1_checkbox_2,
            i1_checkbox_3=self.i1_checkbox_3,
            i2_checkbox_1=self.i2_checkbox_1,
            i2_checkbox_2=self.i2_checkbox_2,
            i3_checkbox_1=self.i3_checkbox_1,
            i3_checkbox_2=self.i3_checkbox_2,
            i3_checkbox_3=self.i3_checkbox_3,
            i4_checkbox_1=self.i4_checkbox_1,
            i4_checkbox_2=self.i4_checkbox_2,
            i5_checkbox_1=self.i5_checkbox_1,
            i5_checkbox_2=self.i5_checkbox_2,
            i6_checkbox_1=self.i6_checkbox_1,
            i6_checkbox_2=self.i6_checkbox_2,
            i7_checkbox_1=self.i7_checkbox_1,
            i7_checkbox_2=self.i7_checkbox_2,
            i9_checkbox_1=self.i9_checkbox_1,
            i9_checkbox_2=self.i9_checkbox_2,
            i10_checkbox_1=self.i10_checkbox_1,
            i10_checkbox_2=self.i10_checkbox_2,
            i11_checkbox_1=self.i11_checkbox_1,
            i11_checkbox_2=self.i11_checkbox_2,
            i12_checkbox_1=self.i12_checkbox_1,
            i12_checkbox_2=self.i12_checkbox_2,
            i13_checkbox_1=self.i13_checkbox_1,
            i13_checkbox_2=self.i13_checkbox_2,
            i13_checkbox_3=self.i13_checkbox_3,
            i14_checkbox_1=self.i14_checkbox_1,
            i14_checkbox_2=self.i14_checkbox_2,
            i14_checkbox_3=self.i14_checkbox_3,
            i14_checkbox_4=self.i14_checkbox_4,
            i15_checkbox_1=self.i15_checkbox_1,
            i15_checkbox_2=self.i15_checkbox_2,
            i18_checkbox_1=self.i18_checkbox_1,
            i18_checkbox_2=self.i18_checkbox_2,
            i18_checkbox_3=self.i18_checkbox_3,
            i20_checkbox_1=self.i20_checkbox_1,
            i20_checkbox_2=self.i20_checkbox_2,
        )
        
        # leave {tval1} (or {tval2}, {tval3}...) in pre_6i file so it can be updated after joining 3p
        if tval_var_to_format != None:
            reaction_options_formatted["tval"+str(tval_var_to_format)] = "{tval}"
            for i in range(1, 9):
                if i == tval_var_to_format:
                    continue
                else:
                    reaction_options_formatted["tval"+str(i)]=f"{'{:.5E}'.format(getattr(self, 'tval'+str(i))):>12}"
        else:
            reaction_options_formatted.update(dict(
                tval1=f"{'{:.5E}'.format(self.tval1):>12}",
                tval2=f"{'{:.5E}'.format(self.tval2):>12}",
                tval3=f"{'{:.5E}'.format(self.tval3):>12}",
                tval4=f"{'{:.5E}'.format(self.tval4):>12}",
                tval5=f"{'{:.5E}'.format(self.tval5):>12}",
                tval6=f"{'{:.5E}'.format(self.tval6):>12}",
                tval7=f"{'{:.5E}'.format(self.tval7):>12}",
                tval8=f"{'{:.5E}'.format(self.tval8):>12}",
            ))
            
        # leave {pval1} (or {pval2}, {pval3}...) in pre_6i file so it can be updated after joining 3p
        if pval_var_to_format != None:
            reaction_options_formatted["pval"+str(pval_var_to_format)] = "{pval}"
            for i in range(1, 6):
                if i == pval_var_to_format:
                    continue
                else:
                    reaction_options_formatted["pval"+str(i)]=f"{'{:.5E}'.format(getattr(self, 'pval'+str(i))):>12}"
        else:
            reaction_options_formatted.update(dict(
                pval1=f"{'{:.5E}'.format(self.pval1):>12}",
                pval2=f"{'{:.5E}'.format(self.pval2):>12}",
                pval3=f"{'{:.5E}'.format(self.pval3):>12}",
                pval4=f"{'{:.5E}'.format(self.pval4):>12}",
                pval5=f"{'{:.5E}'.format(self.pval5):>12}",
            ))

        
        reactant_blocks = "".join([r.formatted_block for r in self.reactants])
        gas_lines = "".join([r.formatted_line for r in self.gases])
        reaction_options_formatted.update(dict(reactant_blocks=reactant_blocks))
        reaction_options_formatted.update(dict(gas_lines=gas_lines))
        
        if self.mineral_suppression_option == "None":
            mineral_suppress_option_value = "None    "
        elif self.mineral_suppression_option == "All":
            mineral_suppress_option_value = "All     "
        else:
            msg = ("Error in Prepare_Reaction(): `mineral_suppression_option` "
                   "must be 'None' or 'All'.")
            self.err_handler.raise_exception(msg)
            
        min_supp_temp = copy.copy(mineral_suppress_template)
        min_supp_dict = dict(mineral_suppress_option=mineral_suppress_option_value)
        min_supp_temp = min_supp_temp.format(**min_supp_dict)
        
        reaction_options_formatted.update(
                dict(mineral_suppress_lines=min_supp_temp)
                )
        
        reaction_template = copy.copy(template)
        self.formatted_reaction = reaction_template.format(**reaction_options_formatted)