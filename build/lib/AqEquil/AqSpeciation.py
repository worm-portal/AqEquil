import os, re, sys, shutil

import warnings
from subprocess import Popen
import pkg_resources
import pandas as pd

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
    def __init__(self, args):
        for k in args:
            setattr(self, k, args[k])

    def __getitem__(self, item):
         return getattr(self, item)
        
    def col_lookup(self, column_list):
        return self.report.iloc[:, self.report.columns.get_level_values(0).isin(set(column_list))]


class AqEquil():
    
    """
    AqEquil
    
    """
    # remember to add a blank line at the end of the AqEquil docstring
    
    def __init__(self,
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO')):
        
        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        self.out_dict = None
        self.messages = True
        
        os.environ['EQ36DA'] = self.eq36da # set eq3 db directory
        os.environ['EQ36CO'] = self.eq36co # set eq3 .exe directory

    def __clear_eqpt_extra_output(self):
        """Deletes non-data1 EQPT output."""
        if os.path.exists("eqpt_log.txt") and os.path.isfile("eqpt_log.txt"):
            os.remove("eqpt_log.txt")
        if os.path.exists("data1f.txt") and os.path.isfile("data1f.txt"):
            os.remove("data1f.txt")
        if os.path.exists("slist.txt") and os.path.isfile("slist.txt"):
            os.remove("slist.txt")
        
        
    def runeqpt(self, db, extra_eqpt_output=False):
        """Uses EQPT to convert a data0 into a data1 file readable by EQ3"""
        
        if os.path.exists("data0."+db) and os.path.isfile("data0."+db):
            pass
        else:
            raise Exception("Error: could not locate custom database",
                            "data0.{} in {}.".format(db, os.getcwd()))
        
        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            os.remove("data1."+db)
            
        self.__clear_eqpt_extra_output()
        
        os.environ['EQ36DA']=os.getcwd()
        
        args = ['/bin/csh', self.eq36co+'/runeqpt', db]
        
        try:
            # run EQPT
            self.__run_script_and_wait(args)
        except:
            os.environ['EQ36DA']=self.eq36da
            raise Exception("Error: EQPT failed to run on {}.".format("data0."+db))
        
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
        
        os.environ['EQ36DA']=self.eq36da # reset default EQ3 db path
    
    def runeq3(self, filename_3i, db,
                 path_3i=os.getcwd(),
                 path_3o=os.getcwd(),
                 path_3p=os.getcwd()):
        """
        Call EQ3 on a .3i input file.
        """

        # get current working dir
        cwd = os.getcwd()
        
        print('Calling EQ3 on ' + filename_3i + ' using ' + db)
        os.chdir(path_3i) # step into 3i folder
        args = ['/bin/csh', self.eq36co+'/runeq3', db, filename_3i]

        # run EQ3
        self.__run_script_and_wait(args)

        # restore working dir
        os.chdir(cwd)
        
        filename_3o =  filename_3i[:-1] + 'o'
        filename_3p =  filename_3i[:-1] + 'p'

        #try:
        # rename output
        os.rename(path_3i + '/output', path_3i + "/" + filename_3o)
        #except:
        #    print('Error: EQ3 failed to produce output for ' + filename_3i)
            
        try:
            # move output
            shutil.move(path_3i + "/" + filename_3o, path_3o + "/" + filename_3o)
        except:
            print('Error: Could not move', filename_3o, "to", path_3o)

            
        try:
            # rename pickup
            os.rename(path_3i + '/pickup', path_3i + "/" + filename_3p)
        except:
            print('Error: EQ3 failed to produce a pickup file for ' + filename_3i)
            
        try:
            # move pickup
            shutil.move(path_3i + "/" + filename_3p, path_3p + "/" + filename_3p)
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
        '''
        Finds all files of a filetype in all downstream folders.
        '''
        file_name = [] # file names
        file_list = [] # file names with paths
        for root, dirs, files in os.walk(location):
            for file in files:
                if file.endswith(file_type):
                    if "-checkpoint" not in file:
                        file_name.append(file)
                        file_list.append(os.path.join(root, file))
        return file_name, file_list
    
    def __run_script_and_wait(self, args):
        with open(os.devnull, 'w') as fp: # devnull supresses written output
            Popen(args, stdout=fp).wait()

    def _delete_rxn_folders(self):
        if os.path.exists('rxn_3i') and os.path.isdir('rxn_3i'):
            shutil.rmtree('rxn_3i')
        if os.path.exists('rxn_3o') and os.path.isdir('rxn_3o'):
            shutil.rmtree('rxn_3o')
        if os.path.exists('rxn_3p') and os.path.isdir('rxn_3p'):
            shutil.rmtree('rxn_3p')
    
    def speciate(self,
                 input_filename,
                 db = "jus",
                 redox_flag = 0,
                 redox_aux = "Fe+3",
                 default_logfO2 = -6,
                 exclude = [],
                 suppress = [],
                 charge_balance_on = "none", 
                 suppress_missing = True,
                 verbose = 1,
                 report_filename = None,
                 get_aq_dist  = True,
                 aq_dist_type = "log_activity",
                 get_aq_contrib  = True,
                 aq_contrib_other = True,
                 get_mineral_sat = True,
                 mineral_sat_type = "affinity",
                 get_redox = True,
                 redox_type = "Eh",
                 get_affinity_energy = False,
                 rxn_filename = None,
                 not_limiting = ["H+", "OH-", "H2O"],
                 get_charge_balance = True,
                 custom_db = False,
                 extra_eqpt_output = False,
                 batch_3o_filename = None,
                 delete_generated_folders = False):

        """
        Calculate the equilibrium distribution of chemical species in solution.
        Additionally, calculate chemical affinities and energy supplies for
        user-specified reactions.
        """
        
        if batch_3o_filename != None:
            if ".rds" in batch_3o_filename[-4:]:
                batch_3o_filename = batch_3o_filename
            else:
                batch_3o_filename = "batch_3o_{}.rds".format(db)
        else:
            batch_3o_filename = "" # no rds file generated
        
        if custom_db:
            self.runeqpt(db, extra_eqpt_output)
            os.environ['EQ36DA']=os.getcwd()
            
        
        # preprocess for eq3 using R scripts
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_prescript = pkg_resources.resource_string(__name__, 'preprocess_for_EQ3.r').decode("utf-8")
            ro.r(r_prescript) # this will need to change when packaging
            df_input_processed = ro.r.preprocess(input_filename=input_filename,
                            exclude=convert_to_RVector(exclude),
                            redox_flag=redox_flag,
                            default_logfO2=default_logfO2,
                            charge_balance_on=charge_balance_on,
                            suppress_missing=suppress_missing,
                            suppress=convert_to_RVector(suppress),
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
            os.environ['EQ36DA']=self.eq36da
        
        # mine output
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_3o_mine = pkg_resources.resource_string(__name__, '3o_mine.r').decode("utf-8")
            ro.r(r_3o_mine) # this will need to change when packaging
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
                              df_input_processed_names=convert_to_RVector(list(self.df_input_processed.columns)),
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
        subheaders = [subheader[1] if len(subheader)>1 else "" for subheader in [col.split("_") for col in list(df_input.columns)]]
        multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
        df_input.columns = multicolumns
        
        df_join = df_input
        
        if get_aq_dist:
            aq_distribution_cols = list(report_divs.rx2('aq_distribution'))
            df_aq_distribution = df_report.loc[:, aq_distribution_cols]
            
            # handle headers of aq_distribution section
            headers = df_aq_distribution.columns
            subheaders = [aq_dist_type]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
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
                raise Exception("mineral_sat_type must be either 'affinity' or 'logQoverK'")

            headers = df_mineral_sat.columns
            subheaders = [mineral_sat_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
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
                raise Exception("redox_type must be either 'Eh', 'pe', 'logfO2', or 'Ah'")

            headers = df_redox.columns
            subheaders = [redox_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
            df_redox.columns = multicolumns
            df_join = df_join.join(df_redox)

        if get_charge_balance:
            charge_balance_cols = list(report_divs.rx2('charge_balance'))
            df_charge_balance = df_report.loc[:, charge_balance_cols]
            
            # handle headers of df_charge_balance section
            headers = df_charge_balance.columns
            '''
            headers:
            ['%CI of mean', '%CI of total', 'Charge imbalance', 'IS (molal)',
           'Mean charge', 'Sigma(mz) anions', 'Sigma(mz) cations', 'Total charge',
           'stoichiometric IS (molal)']
            '''
            subheaders = ["%"]*2 + ['eq/kg.H2O', 'molality'] + ['eq/kg.H2O']*4 + ['molality']
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
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
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
            df_affinity.columns = multicolumns

            # handle headers of df_energy section
            headers = df_energy.columns
            subheaders = ['cal/kg.H2O']*len(headers)
            multicolumns = pd.MultiIndex.from_arrays([headers, subheaders], names=['Sample', ''])
            df_energy.columns = multicolumns
            df_join = df_join.join(df_affinity)
            df_join = df_join.join(df_energy)
            
        out_dict = {'sample_data':{},
                    'aq_contribution':aq_contribution, 'report':df_join,
                    'input':df_input, 'processed_input':df_pinput, 'report_divs':report_divs}
            
        sample_data = batch_3o.rx2('sample_data')
        
        # assemble sample data
        for i, sample in enumerate(sample_data):
            dict_sample_data = {"filename":str(sample.rx2('filename')[0]),
                    "name":str(sample.rx2('name')[0]),
                    "temperature":float(sample.rx2('temperature')[0]),
                    "pressure":float(sample.rx2('pressure')[0]),
                    "logact_H2O":float(sample.rx2('logact_H2O')[0]),
                    "H2O_density":float(sample.rx2('H2O_density')[0]),
                    "H2O_molality":float(sample.rx2('H2O_molality')[0]),
                    "H2O_log_molality":float(sample.rx2('H2O_log_molality')[0]),
                   }
            
            if get_aq_dist:
                dict_sample_data.update({"aq_distribution":pandas2ri.ri2py_dataframe(sample.rx2('aq_distribution'))})
            
            if get_aq_contrib:
                sample_aq_contrib = sample.rx2('aq_contribution')
                dict_sample_data.update({"aq_contribution":pandas2ri.ri2py_dataframe(sample.rx2('aq_contribution'))})
                sample_aq_contrib_dict = {}
                for ii, species_df in enumerate(sample_aq_contrib):
                    species_name = sample_aq_contrib.names[ii]
                    sample_aq_contrib_dict.update({species_name:pandas2ri.ri2py_dataframe(species_df)})
            
            if get_mineral_sat:
                dict_sample_data.update({"mineral_sat":pandas2ri.ri2py_dataframe(sample.rx2('mineral_sat'))})
                
            if get_redox:
                dict_sample_data.update({"redox":pandas2ri.ri2py_dataframe(sample.rx2('redox'))})
            
            if get_charge_balance:
                cbal = sample.rx2('charge_balance')
                charge_balance_dict = {'IS (molal)':float(cbal.rx2('IS (molal)')[0]),
                                       'stoichiometric IS (molal)':float(cbal.rx2('stoichiometric IS (molal)')[0]),
                                       'Sigma(mz) cations':float(cbal.rx2('Sigma(mz) cations')[0]),
                                       'Sigma(mz) anions':float(cbal.rx2('Mean charge')[0]),
                                       'Total charge':float(cbal.rx2('Total charge')[0]),
                                       'Charge imbalance':float(cbal.rx2('Charge imbalance')[0]),
                                       '%CI of total':float(cbal.rx2('%CI of total')[0]),
                                       '%CI of mean':float(cbal.rx2('%CI of mean')[0]),
                                      }
                dict_sample_data.update({"charge_balance":pandas2ri.ri2py_dataframe(sample.rx2('charge_balance'))})

            if get_affinity_energy:
                dict_sample_data.update({"affinity_energy_raw":pandas2ri.ri2py_dataframe(sample.rx2('affinity_energy_raw'))})
                dict_sample_data.update({"affinity_energy":pandas2ri.ri2py_dataframe(sample.rx2('affinity_energy'))})
                
            out_dict["sample_data"].update({sample_data.names[i]:dict_sample_data})
        
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
        suppress_redox = [],
        lettercode = "wrm",
        exceed_Ttr = True, 
        grid_temps = [0.0100, 50.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000],
        grid_press = "Psat",
        infer_formula_ox = False,
        basis_pref_mods = {},
        template_name=None,
        ):

        """
        Create a data0 file from a custom OBIGT table.
        
        Parameters
        ----------
        filename : str, # "example_input/wrm_main_data.csv"
            Path of csv file containing thermodynamic data.
            
        supp_file : str, optional
            Path of file containing data0-specific parameters.
            
        supp_file_ss : str, #"example_input/wrm_solid_solutions.csv"
            Path of file containing solid solution parameters.
            
        template : str
            data0 template
        
        data0_formula_ox_name : str, optional
            Name of supplementary file containing data0 parameters and inferred
            formula oxidation states. Ignored if `infer_formula_ox` is False.
            See `infer_formula_ox` for more detail.
        
        suppress_redox : list of str, default []
            Suppress equilibrium between oxidation states of listed elements
            (Cl, H, and O are not allowed).
        
        lettercode : str, default "wrm"
            Three letter code of data0 output

        exceed_Ttr : bool, default True
            Calculate Gibbs energies of mineral phases and other species
            beyond their transition temperatures?

        grid_temps : list of eight float values, default [0.0100, 50.0000,
        100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000]
            Eight temperature values that make up the T-P grid.
        
        grid_press = list of float or "Psat"
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
        
        basis_pref_mods : list of str
            ???
        
        Returns
        ----------
        
        """
        template = pkg_resources.resource_string(__name__, 'data0.min').decode("utf-8")
        grid_temps = convert_to_RVector(grid_temps)
        suppress_redox = convert_to_RVector(suppress_redox)
        
        if supp_file_ss == None:
            supp_file_ss = ro.r("NULL")
        if data0_formula_ox_name == None:
            data0_formula_ox_name = ro.r("NULL")
        if template_name == None:
            template_name = "sample_template_{}.csv".format(lettercode)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_create_data0 = pkg_resources.resource_string(__name__, 'create_data0.r').decode("utf-8")
            ro.r(r_create_data0) # this will need to change when packaging
            ro.r.main_create_data0(filename=filename,
                                   supp_file=supp_file,
                                   supp_file_ss=supp_file_ss,
                                   grid_temps=grid_temps,
                                   grid_press=grid_press,
                                   lettercode=lettercode,
                                   template=template,
                                   exceed_Ttr=exceed_Ttr,
                                   data0_formula_ox_name=data0_formula_ox_name,
                                   suppress_redox=suppress_redox,
                                   infer_formula_ox=infer_formula_ox,
                                   template_name=template_name,
                                   #basis_pref_mods=basis_pref_mods, #dict of replacements. e.g. iron oxidation example. TODO: implement.
                                  )
            
        for warning in w:
            print(warning.message)
            