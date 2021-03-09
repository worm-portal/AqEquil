import os
import re
import sys
import shutil
import copy
import collections

import warnings
from subprocess import Popen
import pkg_resources
import pandas as pd
from plotnine import *

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()


def convert_to_RVector(value, force_Rvec=True):
    
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
        return ro.StrVector(value)


class Speciation(object):
    
    """
    Stores the output of a speciation calculation.
    
    Attributes
    ----------
    input : pd.Dataframe
        Pandas dataframe containing user-supplied sample chemistry data.
    
    aq_contribution : pd.Dataframe
        Pandas dataframe containing basis species contributions to mass balance
        of aqueous species.
    
    batch_3o : rpy2 ListVector
        An rpy2 ListVector (R object) containing speciation results, in case
        analysis in R is preferred.
    
    processed_input : pd.Dataframe
        Pandas dataframe containing user-supplied sample chemistry data that has
        been processed for `speciate`.
    
    report : pd.Dataframe
        Pandas dataframe reporting major results of speciation calculation in
        across all samples.
    
    report_divs : rpy2 ListVector
        An rpy2 ListVector of column names within the different sections of the
        speciation report.
    
    sample_data : dict
        Dictionary with sample names as keys and speciation results as values.
    
    """
    
    def __init__(self, args):
        for k in args:
            setattr(self, k, args[k])

    def __getitem__(self, item):
         return getattr(self, item)

    def col_lookup(self, column_list):
        
        """
        Look up desired columns in the speciation report.
        
        Parameters
        ----------
        column_list : list of str
            List of column names to look up.
            
        Returns
        ----------
        Pandas dataframe
            The speciation report with only the desired columns.
        """
        
        return self.report.iloc[:, self.report.columns.get_level_values(0).isin(set(column_list))]

    def viz_mass_contribution(self, basis):
        
        """
        Plot basis species contributions to mass balance of aqueous
        species.
        
        Parameters
        ----------
        basis : str
            Name of the basis species
            
        Returns
        ----------
        g : plotnine ggplot object
            A stacked bar plot.
        """
        
        df_spec = copy.deepcopy(
            self.aq_contribution.loc[self.aq_contribution['basis'] == basis])

        df_spec['percent'] = df_spec['percent'].astype(float)

        g = ggplot(df_spec, aes(fill="species", y="percent", x="sample")) + \
            geom_bar(stat="identity") + \
            ylab("%") + \
            ggtitle("Species Accounting for Mass Balance of " + basis) + \
            theme(axis_line=element_line(colour="black", size=0.25, linetype="solid"),
                  axis_text_x=element_text(angle=45, vjust=1, hjust=1),
                  axis_title_x=element_blank(),
                  panel_grid_major=element_blank(), panel_grid_minor=element_blank(),
                  panel_background=element_blank(),
                  legend_key=element_rect(fill=None, color=None),
                  legend_title=element_blank(),
                  plot_title=element_text(size=9, hjust=0.5)) + \
            guides(color=guide_legend(override_aes=None)) + \
            scale_y_continuous(limits=[0, 100], breaks=range(0, 125, 25))

        print(g)


class AqEquil():

    """
    Class containing functions to speciate aqueous water chemistry data using
    existing or custom thermodynamic datasets.
    
    Parameters
    ----------
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored. 
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.
    
    Attributes
    ----------
    eq36da : str
        Path to directory where data1 files are stored.
        
    eq36co : str
        Path to directory where EQ3 executables are stored.
        
    df_input_processed : pd.Dataframe
        Pandas dataframe containing user-supplied sample chemistry data that has
        been processed by `speciate`.
        
    out_dict : pd.Dataframe
        Pandas dataframe reporting results of last speciation calculation
        performed by `speciate`.
    
    messages : bool, default True
        Print messages during calculations?
        
    """

    def __init__(self,
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO')):

        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        self.out_dict = None
        self.messages = True

        os.environ['EQ36DA'] = self.eq36da  # set eq3 db directory
        os.environ['EQ36CO'] = self.eq36co  # set eq3 .exe directory
    
    def _check_sample_input_file(self, input_filename, exclude, db, custom_db):
        
        """
        Check for problems in sample input file.
        """
        
        if '.csv' in input_filename[-4:]:
            if os.path.exists(input_filename) and os.path.isfile(input_filename):
                df_in = pd.read_csv(input_filename, header=None) # no headers for now so colname dupes can be checked
            else:
                err = "Cannot locate input file {}.".format(input_filename)
                raise Exception(err)
        else:
            err = "Input file {}".format(input_filename) + \
                " must be in comma separated values (.csv) format."
            raise Exception(err)
        
        # are there any samples?
        if df_in.shape[0] <= 2:
            err_no_samples = "The file {} ".format(input_filename) + \
                "must contain at least three rows: the " + \
                "first for column names, the second for column subheaders, " + \
                "followed by one or more rows for sample data."
            raise Exception(err_no_samples)
        
        err_list = [] # for appending errors found in the sample input file
        
        # are there duplicate headers?
        col_list = list(df_in.iloc[0, 1:])
        dupe_cols = list(set([x for x in col_list if col_list.count(x) > 1]))
        if len(dupe_cols) > 0:
            err_dupe_cols = "Duplicate column names are not allowed. " + \
                "Duplicate column names were found for:\n{}".format(str(dupe_cols))
            err_list.append(err_dupe_cols)
        
        df_in.columns = df_in.iloc[0] # set column names
        df_in = df_in.drop(df_in.index[0], axis=0) # drop column name row
        df_in_headercheck = copy.deepcopy(df_in.iloc[:,1:]) # drop first column. Deepcopy slice because drop() doesn't work well with unnamed columns.
                                          
        try:
            df_in_headercheck = df_in_headercheck.drop(exclude, axis=1) # drop excluded columns
        except:
            err_bad_exclude = "err_bad_exclude"
            err_list.append(err_bad_exclude)
        
        # are there duplicate rows?
        row_list = list(df_in.iloc[1:, 0])
        dupe_rows = list(set([x for x in row_list if row_list.count(x) > 1]))
        if len(dupe_rows) > 0:
            err_dupe_rows = "Duplicate sample names are not allowed. " + \
                "Duplicate sample names were found for:\n{}".format(str(dupe_rows))
            err_list.append(err_dupe_rows)
        
        # are column names valid entries in the database?
        if custom_db:
            data0_path = "data0." + db
        else:
            data0_path = self.eq36da + "/data0." + db
        if os.path.exists(data0_path) and os.path.isfile(data0_path):
            with open(data0_path) as data0:
                data0_lines = data0.readlines()
                start_index = [i+1 for i, s in enumerate(data0_lines) if '*  species name' in s]
                end_index = [i-1 for i, s in enumerate(data0_lines) if 'elements' in s]
                db_species = [i.split()[0] for i in data0_lines[start_index[0]:end_index[0]]]
                for species in list(set(df_in_headercheck.columns)):
                    if species not in db_species and species != 'Temperature':
                        err_species_not_in_db = "The species '{}' ".format(species) + \
                            "was not found in {}".format(data0_path) + \
                            ". If the column contains data that should not be " + \
                            "included in the speciation calculation, add the " + \
                            "column name to the 'exclude' argument. Try " + \
                            "help(AqEquil.AqEquil.speciate) " + \
                            "for more information about 'exclude'."
                        err_list.append(err_species_not_in_db)
        else:
            err_no_data0 = "Could not locate {}. ".format(data0_path) + \
                "Unable to determine if column headers included in " + \
                "{} ".format(input_filename) + "match entries for species " + \
                "in the requested thermodynamic database '{}'.".format(db)
            err_list.append(err_no_data0)
        
        
        # are subheader units valid?
        subheaders = df_in_headercheck.iloc[0,]
        valid_subheaders = ["degC", "ppm", "ppb", "Suppressed", "Molality",
                            "Molarity", "mg/L", "mg/kg.sol", "Alk., eq/kg.H2O",
                            "Alk., eq/L", "Alk., eq/kg.sol", "Alk., mg/L CaCO3",
                            "Alk., mg/L HCO3-", "Log activity", "Log act combo",
                            "Log mean act", "pX", "pH", "pHCl", "pmH", "pmX",
                            "Hetero. equil.", "Homo. equil.", "Make non-basis"]
        for i, subheader in enumerate(subheaders):
            if subheader not in valid_subheaders:
                err_valid_sub = "The subheader '{}' ".format(subheader) + \
                    "for the column '{}'".format(df_in_headercheck.columns[i]) + \
                    " is not recognized. Valid subheaders are {}".format(str(valid_subheaders)) + \
                    ". If the column {} ".format(df_in_headercheck.columns[i]) + \
                    "contains data that is not meant for the " + \
                    "speciation calculation, add the column name " + \
                    "to the 'exclude' argument. Try help(AqEquil.AqEquil.speciate) " + \
                    "for more information about 'exclude'."
                err_list.append(err_valid_sub)
            
        # is a 'Temperature' column present?
        if "Temperature" not in df_in_headercheck.columns and "Temperature" not in exclude:
            err_temp = "The column 'Temperature' was not found in the input file. "+\
                "Please include a column with 'Temperature' in the first row, "+\
                "'degC' in the second row, and a temperature value for each "+\
                "sample in degrees Celsius."
            err_list.append(err_temp)
        
        # raise exception that outlines all errors found
        if len(err_list) > 0:
            errs = "\n\n*".join(err_list)
            errs = "The input file {}".format(input_filename)+" encountered" + \
                " errors:\n\n*" + errs
            raise Exception(errs)
        
        return
        
        
    def __clear_eqpt_extra_output(self):
        
        """
        Deletes all EQPT output except data1.
        """
        
        if os.path.exists("eqpt_log.txt") and os.path.isfile("eqpt_log.txt"):
            os.remove("eqpt_log.txt")
        if os.path.exists("data1f.txt") and os.path.isfile("data1f.txt"):
            os.remove("data1f.txt")
        if os.path.exists("slist.txt") and os.path.isfile("slist.txt"):
            os.remove("slist.txt")

    def runeqpt(self, db, extra_eqpt_output=False):
        
        """
        Convert a data0 into a data1 file with EQPT.
        
        Parameters
        ----------
        db : str
            Three letter code of database.
        
        extra_eqpt_output : bool, default False
            Keep additional output files from EQPT? These files include
            eqpt_log.txt, data1f.txt, and slist.txt.
        """

        if os.path.exists("data0."+db) and os.path.isfile("data0."+db):
            pass
        else:
            raise Exception("Error: could not locate custom database",
                            "data0.{} in {}.".format(db, os.getcwd()))

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            os.remove("data1."+db)

        self.__clear_eqpt_extra_output()

        os.environ['EQ36DA'] = os.getcwd()

        args = ['/bin/csh', self.eq36co+'/runeqpt', db]

        try:
            # run EQPT
            self.__run_script_and_wait(args)
        except:
            os.environ['EQ36DA'] = self.eq36da
            raise Exception(
                "Error: EQPT failed to run on {}.".format("data0."+db))

        if os.path.exists("data1") and os.path.isfile("data1"):
            os.rename("data1", "data1."+db)
        if os.path.exists("output") and os.path.isfile("output"):
            os.rename("output", "eqpt_log.txt")
        if os.path.exists("data1f") and os.path.isfile("data1f"):
            os.rename("data1f", "data1f.txt")
        if os.path.exists("slist") and os.path.isfile("slist"):
            os.rename("slist", "slist.txt")

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            if self.messages:
                print("Successfully created a data1."+db+" from data0."+db)
        else:
            raise Exception("EQPT could not create data1."+db+" from",
                            "data0."+db+". Check eqpt_log.txt for details.")

        if not extra_eqpt_output:
            self.__clear_eqpt_extra_output()

        os.environ['EQ36DA'] = self.eq36da  # reset default EQ3 db path

    def runeq3(self, filename_3i, db,
               path_3i=os.getcwd(),
               path_3o=os.getcwd(),
               path_3p=os.getcwd()):
        
        """
        Call EQ3 on a .3i input file.
        
        Parameters
        ----------
        filename_3i : str
            Name of 3i input file.
        
        db : str
            Three letter code of database.
        """

        # get current working dir
        cwd = os.getcwd()

        print('Calling EQ3 on ' + filename_3i + ' using ' + db)
        os.chdir(path_3i)  # step into 3i folder
        args = ['/bin/csh', self.eq36co+'/runeq3', db, filename_3i]

        # run EQ3
        self.__run_script_and_wait(args)

        # restore working dir
        os.chdir(cwd)

        filename_3o = filename_3i[:-1] + 'o'
        filename_3p = filename_3i[:-1] + 'p'

        try:
            # rename output
            os.rename(path_3i + '/output', path_3i + "/" + filename_3o)
        except:
            print('Error: EQ3 failed to produce output for ' + filename_3i)

        try:
            # move output
            shutil.move(path_3i + "/" + filename_3o,
                        path_3o + "/" + filename_3o)
        except:
            print('Error: Could not move', filename_3o, "to", path_3o)

        try:
            # rename pickup
            os.rename(path_3i + '/pickup', path_3i + "/" + filename_3p)
        except:
            print('Error: EQ3 failed to produce a pickup file for ' + filename_3i)

        try:
            # move pickup
            shutil.move(path_3i + "/" + filename_3p,
                        path_3p + "/" + filename_3p)
        except:
            print('Error: Could not move', filename_3p, "to", path_3p)

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

    def __read_inputs(self, file_type, location):
        
        """
        Finds all files of a filetype in all downstream folders.
        """
        
        file_name = []  # file names
        file_list = []  # file names with paths
        for root, dirs, files in os.walk(location):
            for file in files:
                if file.endswith(file_type):
                    if "-checkpoint" not in file:
                        file_name.append(file)
                        file_list.append(os.path.join(root, file))
        return file_name, file_list

    def __run_script_and_wait(self, args):
        
        """
        Runs shell commands.
        """
        
        with open(os.devnull, 'w') as fp:  # devnull supresses written output
            Popen(args, stdout=fp).wait()

    def _delete_rxn_folders(self):
        
        """
        Deletes folders storing raw EQ3 input and output.
        """
        
        if os.path.exists('rxn_3i') and os.path.isdir('rxn_3i'):
            shutil.rmtree('rxn_3i')
        if os.path.exists('rxn_3o') and os.path.isdir('rxn_3o'):
            shutil.rmtree('rxn_3o')
        if os.path.exists('rxn_3p') and os.path.isdir('rxn_3p'):
            shutil.rmtree('rxn_3p')

    def speciate(self,
                 input_filename,
                 db="jus",
                 redox_flag=0,
                 redox_aux="Fe+3",
                 default_logfO2=-6,
                 exclude=[],
                 suppress=[],
                 charge_balance_on="none",
                 suppress_missing=True,
                 verbose=1,
                 report_filename=None,
                 get_aq_dist=True,
                 aq_dist_type="log_activity",
                 get_aq_contrib=True,
                 aq_contrib_other=True,
                 get_mineral_sat=True,
                 mineral_sat_type="affinity",
                 get_redox=True,
                 redox_type="Eh",
                 get_affinity_energy=False,
                 rxn_filename=None,
                 not_limiting=["H+", "OH-", "H2O"],
                 get_charge_balance=True,
                 custom_db=False,
                 extra_eqpt_output=False,
                 batch_3o_filename=None,
                 delete_generated_folders=False):
        
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
        
        db : three letter str, default "jus"
            Three letter file extension for the desired thermodynamic database.
            If `custom_db` is False, this database must be named data1.xyz
            (where xyz is your desired three letter extension) and located
            in the EQ3/6 'EQ36DA' path. Otherwise, the database must be named
            data0.xyz and located in your current working directory. Note that
            data1 files are already compiled by EQPT, while data0 files will be
            automatically compiled for you if `custom_db` is True.
        
        redox_flag : int, default 0
            Values corresponding to redox options in the EQ3/6 v8.0 software.
            For more information see the 'Redox Option' section of the EQ3/6
            version 8.0 software user's manual. Set sample redox state with the
            following options:
            
            * -3 for O2(g)
            * -2 for pe (in pe units)
            * -1 for Eh (volts)
            *  0 for logfO2 (log bars), or dissolved O2 (see below)
            *  1 for defining a redox couple (see `redox_aux`)
             
            Note that if you are importing water chemistry data from a
            spreadsheet, a column must be supplied with data that corresponds to
            the redox option you chose. The column name is important:
            
            * -3 must have a column named: O2(g)
            * -2 must have a column named: pe
            * -1 must have a column named: Eh
            *  0 must have a column named: logfO2
            *  1 must have a column corresponding to the auxilliary basis species
              selected to form a redox couple with its linked species (see
              `redox_aux`). For example, the redox couple Fe+2/Fe+3 would need
              a column named: Fe+3
            
            If an appropriate header or redox data cannot be found to define
            redox state, `default_logfO2` is used to set sample logfO2.
            
            There is a special case where dissolved oxygen can be used to impose
            sample redox state if `redox_flag` is set to 0 and a column named
            logfO2 does not appear in the sample data sheet. If there is a
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
        
        get_aq_contrib : bool, default True
            Calculate basis species contributions to mass balance of aqueous
            species?
        
        aq_contrib_other : bool, default True
            Include an "other" species for the sake of summing percents of basis
            species contributions to 100%? Ignored if `get_aq_contrib` is False.
        
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
            
        get_affinity_energy : bool, default False
            Calculate affinities and energy supplies of reactions listed in a
            separate user-supplied file?
        
        rxn_filename : str, optional
            Name of file containing reactions used to calculate affinities and
            energy supplies. Ignored if `get_affinity_energy` is False.
        
        not_limiting : list, default ["H+", "OH-", "H2O"]
            List containing names of species that are not considered limiting
            when calculating energy supplies. Ignored if `get_affinity_energy`
            is False.
        
        get_charge_balance : bool, default True
            Calculate charge balance and ionic strength?
            
        custom_db : bool, default False
            Is the database defined by `db` a custom user-supplied database? If
            this is set to True, searches for a data0.xyz file in the current
            working directory, where 'xyz' corresponds to the three letter code
            assigned to `db`. This data0 file is automatically converted into a
            machine-readable file called data1 by software called EQPT. This
            data1 file is then used in speciation calculations.
        
        extra_eqpt_output : bool, default False
            Keep additional output files created by EQPT (see `custom_db`)?
            Ignored if `custom_db` is False.
        
        batch_3o_filename : str, optional
            Name of rds (R object) file exported after the speciation
            calculation? No file will be generated if this argument is not
            defined.
            
        delete_generated_folders : bool, default False
            Delete the 'rxn_3i', 'rxn_3o', and 'rxn_3p' folders containing raw
            EQ3NR input, output, and pickup files once the speciation
            calculation is complete?
        
        Returns
        -------
        speciation : object of class Speciation
            Contains the results of the speciation calculation.
        
        """
        
        # check input sample file for errors
        self._check_sample_input_file(input_filename, exclude, db, custom_db)
        
        # handle batch_3o naming
        if batch_3o_filename != None:
            if ".rds" in batch_3o_filename[-4:]:
                batch_3o_filename = batch_3o_filename
            else:
                batch_3o_filename = "batch_3o_{}.rds".format(db)
        else:
            batch_3o_filename = ro.r("NULL")

        if custom_db:
            # EQ3/6 cannot handle spaces in the 'EQ36DA' path name.
            if " " in os.getcwd():
                raise Exception(
                    "Error: the path to the custom database " + \
                    "cannot contain spaces. The current path " + \
                    "is: [ " + os.getcwd() + " ]. Remove or " + \
                    "replace spaces in folder names for this " + \
                    "feature. Example: [ " + \
                    os.getcwd().replace(" ", "-") + " ].")

            self.runeqpt(db, extra_eqpt_output)
            os.environ['EQ36DA'] = os.getcwd()

        if get_affinity_energy:
            if rxn_filename == None:
                warnings.warn(
                    "A reaction file was not specified. Affinities and " + \
                    "energy supplies will not be calculated.")
                get_affinity_energy = False
                rxn_filename = ""
            elif os.path.exists(rxn_filename) and os.path.isfile(rxn_filename):
                pass
            else:
                warnings.warn(
                    "Reaction file {} was not found. Affinities and " + \
                    "energy supplies will not be " + \
                    "calculated.".format(rxn_filename))
                get_affinity_energy = False
                rxn_filename = ""
        else:
            rxn_filename = ""

        # preprocess for eq3 using R scripts
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_prescript = pkg_resources.resource_string(
                __name__, 'preprocess_for_EQ3.r').decode("utf-8")
            ro.r(r_prescript)  # this will need to change when packaging
            df_input_processed = ro.r.preprocess(input_filename=input_filename,
                                                 exclude=convert_to_RVector(
                                                     exclude),
                                                 redox_flag=redox_flag,
                                                 default_logfO2=default_logfO2,
                                                 charge_balance_on=charge_balance_on,
                                                 suppress_missing=suppress_missing,
                                                 suppress=convert_to_RVector(
                                                     suppress),
                                                 verbose=verbose)

        for warning in w:
            print(warning.message)

        self.df_input_processed = pandas2ri.ri2py_dataframe(df_input_processed)

        # run EQ3 on each input file
        cwd = os.getcwd()

        self.__mk_check_del_directory('rxn_3o')
        self.__mk_check_del_directory('rxn_3p')
        three_i_files, three_i_file_paths = self.__read_inputs('3i', 'rxn_3i')

        input_dir = cwd + "/rxn_3i/"
        output_dir = cwd + "/rxn_3o/"
        pickup_dir = cwd + "/rxn_3p/"

        for file in three_i_files:
            self.runeq3(filename_3i=file, db=db, path_3i=input_dir,
                        path_3o=output_dir, path_3p=pickup_dir)

        if custom_db:
            os.environ['EQ36DA'] = self.eq36da

        # mine output
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_3o_mine = pkg_resources.resource_string(
                __name__, '3o_mine.r').decode("utf-8")
            ro.r(r_3o_mine)  # this will need to change when packaging
            batch_3o = ro.r.main_3o_mine(
                input_filename=input_filename,
                rxn_filename=rxn_filename,
                get_aq_dist=get_aq_dist,
                aq_dist_type=aq_dist_type,
                get_aq_contrib=get_aq_contrib,
                aq_contrib_other=aq_contrib_other,
                get_mineral_sat=get_mineral_sat,
                mineral_sat_type=mineral_sat_type,
                get_redox=get_redox,
                redox_type=redox_type,
                get_charge_balance=get_charge_balance,
                get_affinity_energy=get_affinity_energy,
                not_limiting=convert_to_RVector(not_limiting),
                batch_3o_filename=batch_3o_filename,
                df_input_processed=pandas2ri.py2ri(self.df_input_processed),
                # Needed for keeping symbols in column names after porting
                #   df_input_processed in the line above. Some kind of check.names
                #   option for pandas2ri.py2ri would be nice. Workaround:
                df_input_processed_names=convert_to_RVector(
                    list(self.df_input_processed.columns)),
            )
        for warning in w:
            print(warning.message)

        aq_contribution = pandas2ri.ri2py_dataframe(batch_3o[1])
        df_report = pandas2ri.ri2py_dataframe(batch_3o[2])
        df_input = pandas2ri.ri2py_dataframe(batch_3o[3])
        df_pinput = pandas2ri.ri2py_dataframe(batch_3o[4])
        report_divs = batch_3o.rx2('report_divs')

        input_cols = list(report_divs.rx2('input'))
        df_input = df_report.loc[:, input_cols]

        # handle headers and subheaders of input section
        headers = [col.split("_")[0] for col in list(df_input.columns)]
        headers = [header+"_(input)" for header in headers]
        subheaders = [subheader[1] if len(subheader) > 1 else "" for subheader in [
            col.split("_") for col in list(df_input.columns)]]
        multicolumns = pd.MultiIndex.from_arrays(
            [headers, subheaders], names=['Sample', ''])
        df_input.columns = multicolumns

        df_join = df_input

        if get_aq_dist:
            aq_distribution_cols = list(report_divs.rx2('aq_distribution'))
            df_aq_distribution = df_report.loc[:, aq_distribution_cols]

            # handle headers of aq_distribution section
            headers = df_aq_distribution.columns
            subheaders = [aq_dist_type]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_aq_distribution.columns = multicolumns
            df_join = df_join.join(df_aq_distribution)

        if get_mineral_sat:
            mineral_sat_cols = list(report_divs.rx2('mineral_sat'))
            df_mineral_sat = df_report.loc[:, mineral_sat_cols]

            # handle headers of df_mineral_sat section
            if mineral_sat_type == "affinity":
                mineral_sat_unit = "affinity_kcal"
            elif mineral_sat_type == "logQoverK":
                mineral_sat_unit = "logQ/K"
            else:
                raise Exception(
                    "mineral_sat_type must be either 'affinity' or 'logQoverK'")

            headers = df_mineral_sat.columns
            subheaders = [mineral_sat_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_mineral_sat.columns = multicolumns
            df_join = df_join.join(df_mineral_sat)

        if get_redox:
            redox_cols = list(report_divs.rx2('redox'))
            df_redox = df_report.loc[:, redox_cols]

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
                raise Exception(
                    "redox_type must be either 'Eh', 'pe', 'logfO2', or 'Ah'")

            headers = df_redox.columns
            subheaders = [redox_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_redox.columns = multicolumns
            df_join = df_join.join(df_redox)

        if get_charge_balance:
            charge_balance_cols = list(report_divs.rx2('charge_balance'))
            df_charge_balance = df_report.loc[:, charge_balance_cols]

            # handle headers of df_charge_balance section
            headers = df_charge_balance.columns
            subheaders = ["%"]*2 + ['eq/kg.H2O', 'molality'] + \
                ['eq/kg.H2O']*4 + ['molality']
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_charge_balance.columns = multicolumns
            df_join = df_join.join(df_charge_balance)

        if get_affinity_energy:
            affinity_cols = list(report_divs.rx2('affinity'))
            energy_cols = list(report_divs.rx2('energy'))
            df_affinity = df_report.loc[:, affinity_cols]
            df_energy = df_report.loc[:, energy_cols]

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
                    'aq_contribution': aq_contribution, 'report': df_join,
                    'input': df_input, 'processed_input': df_pinput, 'report_divs': report_divs}

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
                dict_sample_data.update(
                    {"aq_distribution": pandas2ri.ri2py_dataframe(sample.rx2('aq_distribution'))})

            if get_aq_contrib:
                sample_aq_contrib = sample.rx2('aq_contribution')
                dict_sample_data.update(
                    {"aq_contribution": pandas2ri.ri2py_dataframe(sample.rx2('aq_contribution'))})
                sample_aq_contrib_dict = {}
                for ii, species_df in enumerate(sample_aq_contrib):
                    species_name = sample_aq_contrib.names[ii]
                    sample_aq_contrib_dict.update(
                        {species_name: pandas2ri.ri2py_dataframe(species_df)})

            if get_mineral_sat:
                dict_sample_data.update(
                    {"mineral_sat": pandas2ri.ri2py_dataframe(sample.rx2('mineral_sat'))})

            if get_redox:
                dict_sample_data.update(
                    {"redox": pandas2ri.ri2py_dataframe(sample.rx2('redox'))})

            if get_charge_balance:
                cbal = sample.rx2('charge_balance')
                charge_balance_dict = {
                    'IS (molal)': float(cbal.rx2('IS (molal)')[0]),
                    'stoichiometric IS (molal)': float(cbal.rx2('stoichiometric IS (molal)')[0]),
                    'Sigma(mz) cations': float(cbal.rx2('Sigma(mz) cations')[0]),
                    'Sigma(mz) anions': float(cbal.rx2('Mean charge')[0]),
                    'Total charge': float(cbal.rx2('Total charge')[0]),
                    'Charge imbalance': float(cbal.rx2('Charge imbalance')[0]),
                    '%CI of total': float(cbal.rx2('%CI of total')[0]),
                    '%CI of mean': float(cbal.rx2('%CI of mean')[0]),
                    }
                dict_sample_data.update(
                    {"charge_balance": pandas2ri.ri2py_dataframe(sample.rx2('charge_balance'))})

            if get_affinity_energy:
                dict_sample_data.update({"affinity_energy_raw": pandas2ri.ri2py_dataframe(
                    sample.rx2('affinity_energy_raw'))})
                dict_sample_data.update(
                    {"affinity_energy": pandas2ri.ri2py_dataframe(sample.rx2('affinity_energy'))})

            out_dict["sample_data"].update(
                {sample_data.names[i]: dict_sample_data})

        out_dict.update({"batch_3o": batch_3o})
        speciation = Speciation(out_dict)

        if report_filename != None:
            if ".csv" in report_filename[-4:]:
                out_dict["report"].to_csv(report_filename)
            else:
                out_dict["report"].to_csv(report_filename+".csv")

        if delete_generated_folders:
            self._delete_rxn_folders()

        return speciation

    def create_data0(self,
                     filename,
                     supp_file,
                     supp_file_ss=None,
                     data0_formula_ox_name=None,
                     suppress_redox=[],
                     db="wrm",
                     exceed_Ttr=True,
                     grid_temps=[0.0100, 50.0000, 100.0000, 150.0000,
                                 200.0000, 250.0000, 300.0000, 350.0000],
                     grid_press="Psat",
                     infer_formula_ox=False,
                     # basis_pref_mods={}, # TODO: dict of replacements. e.g. iron oxidation example.
                     template_name=None,
                     ):
        """
        Create a data0 file from a custom thermodynamic dataset.
        
        Parameters
        ----------
        filename : str
            Name of csv file containing thermodynamic data in the OBIGT format.
            
        supp_file : str
            Path of file containing data0-specific parameters.
            
        supp_file_ss : str, optional
            Name of file containing solid solution parameters.
        
        data0_formula_ox_name : str, optional
            Name of supplementary file containing data0 parameters and inferred
            formula oxidation states. Ignored if `infer_formula_ox` is False.
            See `infer_formula_ox` for more detail.
        
        suppress_redox : list of str, default []
            Suppress equilibrium between oxidation states of listed elements
            (Cl, H, and O cannot be included).
        
        db : str, default "wrm"
            Desired three letter code of data0 output.

        exceed_Ttr : bool, default True
            Calculate Gibbs energies of mineral phases and other species
            beyond their transition temperatures?

        grid_temps : list of eight float, default [0.0100, 50.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000]
            Eight temperature values that make up the T-P grid.
        
        grid_press : list of float, default "Psat"
            Eight pressure values that make up the T-P grid. "Psat" for
            calculations along the liquid-vapor saturation curve.

        infer_formula_ox : bool, default False
            Create a supplementary file containing data0 parameters and
            inferred formula oxidation states? This option is useful for
            creating as many entries in the formula_ox column when creating a
            new supplementary file. Note that compounds like DySO4+ result in
            blank entries in formula_ox because the redox states of two
            elements, Dy and S, would have to be estimated together; S has many
            oxidation states and Dy's oxidation states are not hard-coded.
        
        template_name : str, optional
            Name of the sample input template file generated. If no name is
            supplied, defaults to 'sample_template_xyz.csv', where 'xyz' is
            the three letter code given to `db`.
        """
        
        template = pkg_resources.resource_string(
            __name__, 'data0.min').decode("utf-8")
        grid_temps = convert_to_RVector(grid_temps)
        grid_press = convert_to_RVector(grid_press)
        suppress_redox = convert_to_RVector(suppress_redox)

        if supp_file_ss == None:
            supp_file_ss = ro.r("NULL")
        if data0_formula_ox_name == None:
            data0_formula_ox_name = ro.r("NULL")
        if template_name == None:
            template_name = "sample_template_{}.csv".format(db)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_create_data0 = pkg_resources.resource_string(
                __name__, 'create_data0.r').decode("utf-8")
            ro.r(r_create_data0)  # this will need to change when packaging
            ro.r.main_create_data0(filename=filename,
                                   supp_file=supp_file,
                                   supp_file_ss=supp_file_ss,
                                   grid_temps=grid_temps,
                                   grid_press=grid_press,
                                   db=db,
                                   template=template,
                                   exceed_Ttr=exceed_Ttr,
                                   data0_formula_ox_name=data0_formula_ox_name,
                                   suppress_redox=suppress_redox,
                                   infer_formula_ox=infer_formula_ox,
                                   template_name=template_name,
                                   #basis_pref_mods=basis_pref_mods, # TODO: dict of replacements. e.g. iron oxidation example.
                                   )

        for warning in w:
            print(warning.message)
