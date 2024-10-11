Module AqEquil.AqSpeciation
===========================

Functions
---------

    
`check_balance(formulas, stoich)`
:   Check that a chemical reaction is balanced. If not, get missing composition.
    
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

    
`chemlabel(name, charge_sign_at_end=False)`
:   Format a chemical formula to display subscripts and superscripts in HTML
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

    
`compare(*args)`
:   Combine two or more speciations into a single speciation object for
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

    
`format_equation(species, stoich, charge_sign_at_end=False)`
:   Format a chemical equation to display in HTML
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

    
`load(filename, messages=True, hide_traceback=True)`
:   Load a speciation file.
    
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

Classes
-------

`AqEquil(eq36da=None, eq36co=None, db='WORM', elements=None, solid_solutions=None, logK=None, logK_S=None, logK_extrapolate='none', download_csv_files=False, exclude_organics=False, exclude_category=None, suppress_redox=None, input_template='none', water_model='SUPCRT92', exceed_Ttr=True, verbose=1, load_thermo=True, hide_traceback=True)`
:   Class containing functions to speciate aqueous water chemistry data using
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

    ### Class variables

    `Thermodata`
    :   Metaclass to store and load thermodynamic databases.
        Inherits attributes from its outer class, AqEquil.

    ### Methods

    `create_data0(self, db, filename_ss=None, activity_model='b-dot', exceed_Ttr=True, grid_temps=[0.01, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0], grid_press='Psat', P1=True, plot_poly_fit=False, logK_extrapolate='none', fill_data0=True, dynamic_db=False, dynamic_db_sample_temps=[], dynamic_db_sample_press=[], verbose=1)`
    :   Create a data0 file from a custom thermodynamic dataset.
        
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

    `make_redox_reactions(self, *args, **kwargs)`
    :   Deprecated

    `plot_logK_fit(self, name, plot_out=False, res=200, internal=True, logK_extrapolate=None, T_vals=[])`
    :   Plot the fit of logK values used in the speciation.
        
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

    `runeq3(self, filename_3i, db, samplename=None, path_3i='', path_3o='', path_3p='', data1_path='', dynamic_db_name=None)`
    :   Call EQ3 on a .3i input file.
        
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

    `runeq6(self, filename_6i, db, samplename=None, path_6i='', data1_path=None, dynamic_db_name=None)`
    :   Call EQ6 on a .6i input file.
        
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

    `runeqpt(self, db, dynamic_db=False)`
    :   Convert a data0 into a data1 file with EQPT.
        
        Parameters
        ----------
        db : str
            Three letter code of database.

    `show_redox_reactions(self, *args, **kwargs)`
    :   Deprecated

    `speciate(self, input_filename, db=None, db_solid_solution=None, db_logK=None, logK_extrapolate=None, activity_model='b-dot', redox_flag='logfO2', redox_aux='Fe+3', default_logfO2=-6, exclude=[], suppress=[], alter_options=[], charge_balance_on='none', suppress_missing=True, blanks_are_0=False, strict_minimum_pressure=True, aq_scale=1, verbose=1, report_filename=None, get_aq_dist=True, aq_dist_type='log_activity', get_mass_contribution=True, mass_contribution_other=True, get_mineral_sat=True, mineral_sat_type='affinity', get_redox=True, redox_type='Eh', get_ion_activity_ratios=True, get_fugacity=True, get_basis_totals=True, get_solid_solutions=True, get_affinity_energy=False, negative_energy_supplies=False, mineral_reactant_energy_supplies=False, rxn_filename=None, not_limiting=['H+', 'OH-', 'H2O'], get_charge_balance=True, custom_db=False, batch_3o_filename=None, delete_generated_folders=False, db_args={})`
    :   Calculate the equilibrium distribution of chemical species in solution.
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

`Error_Handler(clean=True)`
:   Handles how errors are printed in Jupyter notebooks. By default, errors that
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

    ### Static methods

    `hide_traceback(exc_tuple=None, filename=None, tb_offset=None, exception_only=False, running_compiled_code=False)`
    :   Return a modified ipython showtraceback function that does not display
        traceback when encountering an error.

    ### Methods

    `raise_exception(self, msg)`
    :   Raise an exception that displays the error message without traceback. This
        happens only when the exception is predicted by the AqEquil package
        (e.g., for common user errors).

`Speciation(args, hide_traceback=True)`
:   Stores the output of a speciation calculation.
    
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

    ### Static methods

    `format_reaction(coeffs, names, formatted=True, charge_sign_at_end=True, show=True)`
    :

    ### Methods

    `apply_redox_reactions(self, y_type='E', y_units='cal', limiting=None, grams_minerals=0, negative_energy_supplies=False, custom_grouping_filepath=None, append_report=True)`
    :   Calculate chemical affinities, energy supplies, and more in samples
        using the redox reactions generated by the function
        `make_redox_reactions`.
        
        Parameters
        ----------
        grams_minerals : float or dict, default 0
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

    `barplot(self, y='pH', title=None, convert_log=True, plot_zero=True, show_missing=True, plot_width=4, plot_height=3, ppi=122, colormap='WORM', save_as=None, save_format=None, save_scale=1, interactive=True, plot_out=False)`
    :   Show a bar plot to vizualize one or more variables across all samples.
        
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

    `calculate_energy(self, species, stoich, divisor=1, per_electron=False, grams_minerals=0, rxn_name='custom reaction', negative_energy_supplies=False, y_type='A', y_units='kcal', limiting=None, charge_sign_at_end=False, as_written=False, simple_df_output=False, append_report=True, custom_grouping_filepath=None, print_logK_messages=False, raise_nonlimiting_exception=True)`
    :   Calculate Gibbs free energy, logK, logQ, chemical affinity, or energy
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
        
        grams_minerals : float or dict, default 0
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

    `join_6i_p(self, filepath_6i, chain_mt)`
    :

    `lookup(self, col=None)`
    :   Look up desired columns in the speciation report.
        
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

    `make_redox_reactions(self, idx_list='all', show=True, formatted=True, charge_sign_at_end=False)`
    :   Generate an organized collection of redox reactions for calculating
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

    `mt(self, sample)`
    :   Retrieve mass transfer results for a sample.
        
        Parameters
        ----------
        sample : str
            Name of the sample for which to retrieve mass transfer results.
            
        Returns
        -------
        An object of class `AqEquil.MassTransfer.Mass_Transfer`.

    `plot_logK_fit(self, name, plot_out=False, res=200, internal=True, logK_extrapolate=None, T_vals=[])`
    :   Plot the fit of logK values used in the speciation.
        
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

    `plot_mass_contribution(self, basis, title=None, sort_by=None, ascending=True, sort_y_by=None, width=0.9, colormap='WORM', sample_label='sample', colors=None, plot_width=4, plot_height=3, ppi=122, save_as=None, save_format=None, save_scale=1, interactive=True, plot_out=False)`
    :   Plot basis species contributions to mass balance of aqueous species
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

    `plot_mineral_saturation(self, sample_name, title=None, mineral_sat_type='affinity', plot_width=4, plot_height=3, ppi=122, colors=['blue', 'orange'], save_as=None, save_format=None, save_scale=1, interactive=True, plot_out=False)`
    :   Vizualize mineral saturation states in a sample as a bar plot.
        
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

    `plot_solid_solutions(self, sample, title=None, width=0.9, colormap='WORM', affinity_plot=True, affinity_plot_colors=['blue', 'orange'], plot_width=4, plot_height=4, ppi=122, save_as=None, save_format=None, save_scale=1, interactive=True, plot_out=False)`
    :   Plot fractions of minerals of hypothetical solid solutions in a sample.
        
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

    `save(self, filename, messages=True)`
    :   Save the speciation as a '.speciation' file to your current working
        directory. This file can be loaded with `AqEquil.load(filename)`.
        
        Parameters
        ----------
        filename : str
            The desired name of the file.
            
        messages : str
            Print a message confirming the save?

    `scatterplot(self, x='pH', y='Temperature', samples=None, title=None, log_x=False, log_y=False, plot_zero=True, rxns_as_labels=True, charge_sign_at_end=False, plot_width=4, plot_height=3, ppi=122, fill_alpha=0.7, point_size=10, ylab=None, lineplot=False, linemarkers=True, colormap='WORM', save_as=None, save_format=None, save_scale=1, interactive=True, plot_out=False)`
    :   Vizualize two or more sample variables with a scatterplot.
        
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

    `show_redox_reactions(self, formatted=True, charge_sign_at_end=False, show=True)`
    :   Show a table of redox reactions generated with the function
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