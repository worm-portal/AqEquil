"""
preprocess_for_EQ3.py - Python conversion of preprocess_for_EQ3.r

Preprocessing functions for EQ3 input files.
"""

import pandas as pd
import numpy as np
import re
import pychnosz


def vprint(string, verbose=2):
    """Print messages depending on desired level of 'verbose'-ness."""
    if verbose > 0:
        print(string)
    elif verbose >= 0:
        if "warning" in string.lower() or "error" in string.lower():
            print(string)


def header_unit(header, keepcase=False):
    """
    Get units from a header.
    (unit should always be at end of header after an underscore)
    """
    split_header = header.split("_")
    if not keepcase:
        return split_header[-1].lower()
    else:
        return split_header[-1]


def uc_molal(value=1, chemical="H+", unit="ppm"):
    """Convert concentration units to molality."""
    # convert unit name to lowercase
    unit = unit.lower()

    # convert ppb to ppm
    if unit == "ppb":
        value = value / 1000
        unit = "ppm"

    # convert ppm or mg/L to molality
    if unit in ["ppm", "mg/l"]:
        value = value * 0.001 / pychnosz.mass(pychnosz.info(pychnosz.info(chemical, messages=False), messages=False)["formula"][0])
        return value
    elif unit in ["molality", "molal"]:
        return value
    else:
        raise ValueError("Error in uc_molal(): unit must be either ppb, ppm, or mg/L")


def header_species(header):
    """
    Get chemical species from header.
    (should always be at start of header before an underscore)
    """
    return header.split("_")[0]


# EQ3 jflags
EQ3_jflags = [
    "Suppressed", "Molality", "Molarity", "mg/L",
    "mg/kg.sol", "Alk., eq/kg.H2O", "Alk., eq/L",
    "Alk., eq/kg.sol", "Alk., mg/L CaCO3", "Alk., mg/L HCO3-",
    "Log activity", "Log act combo", "Log mean act", "pX",
    "pH", "pHCl", "pmH", "pmX", "Hetero. equil.",
    "Homo. equil.", "Make non-basis", "logfO2"
]

# function to convert ppb, ppm, or mg/L to molality
acceptable_units = ["ppb", "ppm", "mg/l", "molality", "molal"]


def preprocess(input_filename,
               exclude,
               grid_temps,
               grid_press,
               strict_minimum_pressure,
               dynamic_db,
               poly_coeffs_1,
               poly_coeffs_2,
               water_model,
               verbose=2):
    """
    Preprocess input file for EQ3 calculations.

    Parameters
    ----------
    input_filename : str
        Path to input CSV file
    exclude : list
        List of column names to exclude
    grid_temps : array-like
        Temperature grid
    grid_press : array-like
        Pressure grid
    strict_minimum_pressure : bool
        Whether to enforce minimum pressure
    dynamic_db : bool
        Whether using dynamic database
    poly_coeffs_1 : array-like or str
        First polynomial coefficients
    poly_coeffs_2 : array-like or str
        Second polynomial coefficients
    water_model : str
        Water model to use
    verbose : int
        Verbosity level

    Returns
    -------
    dict
        Dictionary containing processed dataframe and parameters
    """

    # Read input
    df = pd.read_csv(input_filename, header=0, keep_default_na=False, dtype=str)

    # Create a copy for numeric processing
    df_numeric = df.drop(columns=[col for col in exclude if col in df.columns])

    # Filter out columns with "Hetero. equil." in first row
    hetero_cols = [col for col in df_numeric.columns if df_numeric.iloc[0][col] == "Hetero. equil."]
    df_numeric = df_numeric.drop(columns=hetero_cols)

    num_cols = list(df_numeric.columns)
    # Remove first column (sample names) from numeric columns
    if len(num_cols) > 0 and num_cols[0] == df.columns[0]:
        num_cols = num_cols[1:]

    # Specify headers and subheaders
    header1 = list(df.columns)
    header2 = list(df.iloc[0])

    # Create concatenated header_subheader column names
    new_colnames = []
    for i in range(len(header1)):
        if pd.isna(header2[i]) or header2[i] == "":
            new_colnames.append(header1[i])
        else:
            new_colnames.append(f"{header1[i]}_{header2[i]}")

    df.columns = new_colnames
    df = df.iloc[1:]  # remove subheader row
    df = df.rename(columns={df.columns[0]: "Sample"})  # rename first column as sample name header

    # Convert data to numeric
    for col in num_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Convert special columns to numeric even if excluded from aqueous block
    special_numeric = ["rho_", "Pressure_", "redox_", "Temperature_", "pe_", "Eh_", "logfO2_"]

    for col_prefix in special_numeric:
        matching_cols = [col for col in df.columns if col_prefix in col]
        for col in matching_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Add columns to exclusion list
    exclude = list(exclude) + ["", "Sample", "rho", "Pressure", "redox", "Temperature", "pe", "Eh", "logfO2"]

    # Find temperature (in degK)
    if "Temperature_degC" in df.columns:
        temp_degC = df["Temperature_degC"].values
        temp_degK = temp_degC + 273.15
    elif "Temperature_degK" in df.columns:
        temp_degK = df["Temperature_degK"].values
        temp_degC = temp_degK - 273.15
    else:
        raise ValueError("Error: need a column for 'Temperature' with units in 'degC' or 'degK'")

    # Handle TP-dependent preprocessing
    if not dynamic_db:
        grid_temps = np.array(grid_temps, dtype=float)
        grid_press = np.array(grid_press, dtype=float)
        minimum_pressure = float(np.min(grid_press))

        # Interpolate pressure
        def f1(T):
            sum_val = 0
            for i in range(len(poly_coeffs_1)):
                sum_val = sum_val + poly_coeffs_1[i] * T**(i)
            return sum_val

        def f2(T):
            sum_val = 0
            for i in range(len(poly_coeffs_2)):
                sum_val = sum_val + poly_coeffs_2[i] * T**(i)
            return sum_val

        pressure_bar = []
        for T in temp_degC:
            if grid_temps[0] <= T <= grid_temps[len(poly_coeffs_1)-1]:
                pressure_bar.append(f1(T))
            elif grid_temps[len(poly_coeffs_1)-1] <= T <= grid_temps[-1]:
                pressure_bar.append(f2(T))
            else:
                raise ValueError(f"Error: one or more temperatures in this sample set is outside of the temperature range of this thermodynamic dataset ({grid_temps[0]} to {grid_temps[-1]} C).")

        pressure_bar = np.array(pressure_bar)

        if strict_minimum_pressure:
            for i in range(len(pressure_bar)):
                if pressure_bar[i] < minimum_pressure:
                    pressure_bar[i] = minimum_pressure

    else:
        if "Pressure_bar" in df.columns:
            if df["Pressure_bar"].isna().any():
                vprint("Warning: a pressure value is missing for one or more " +
                      "samples in the sample input file. Defaulting to water " +
                      "saturation pressure for these missing pressures...",
                      verbose=verbose)

                pychnosz.water(water_model, messages=False)
                psat_pressures = []
                for T in temp_degC + 273.15:
                    psat_result = pychnosz.water("Psat", T=T, messages=False)
                    psat_pressures.append(psat_result)
                psat_pressures = np.array(psat_pressures)

                # Fill missing pressures with PSAT
                pressure_bar = df["Pressure_bar"].values.copy()
                pressure_bar[pd.isna(df["Pressure_bar"])] = psat_pressures[pd.isna(df["Pressure_bar"])]
            else:
                pressure_bar = df["Pressure_bar"].values
        else:
            # Assume PSAT if no pressure column in sample input file
            vprint("Warning: a column for Pressure was not found in the " +
                  "sample input file. Defaulting to water saturation pressure...",
                  verbose=verbose)
            pychnosz.water(water_model, messages=False)
            pressure_bar = []
            for T in temp_degC + 273.15:
                psat_result = pychnosz.water("Psat", T=T, messages=False)
                pressure_bar.append(psat_result)
            pressure_bar = np.array(pressure_bar)

        minimum_pressure = 0

    # Create unique df row names that will become .3i filenames
    # Python equivalent of R's make.names with unique=TRUE
    row_order = []
    seen = {}
    for sample in df["Sample"]:
        # Replace invalid characters with dots
        name = re.sub(r'[^a-zA-Z0-9._]', '.', str(sample))
        # Ensure it doesn't start with a digit or dot followed by digit
        if re.match(r'^(\.|[0-9])', name):
            name = 'X' + name
        # Make unique
        if name in seen:
            seen[name] += 1
            name = f"{name}.{seen[name]}"
        else:
            seen[name] = 0
        row_order.append(name)

    df.index = row_order

    input_processed_list = {
        "df": df,
        "temp_degC": temp_degC,
        "pressure_bar": pressure_bar,
        "minimum_pressure": minimum_pressure,
        "exclude": exclude
    }

    return input_processed_list


def write_3i_file(df,
                  temp_degC,
                  pressure_bar,
                  minimum_pressure,
                  strict_minimum_pressure,
                  pressure_override,
                  suppress_missing,
                  exclude,
                  allowed_aq_block_species,
                  charge_balance_on,
                  suppress,
                  alter_options,
                  aq_scale,
                  get_solid_solutions,
                  input_dir,
                  redox_flag,
                  redox_aux,
                  default_logfO2,
                  water_model,
                  warned_about_redox_column,
                  activity_model,
                  verbose):
    """
    Write a .3i input file for EQ3.

    This is a large function that generates formatted text files for EQ3NR.

    Returns
    -------
    bool
        Updated warned_about_redox_column flag
    """

    # Make a copy to avoid SettingWithCopyWarning
    df = df.copy()

    # Set water model
    pychnosz.water(water_model, messages=False)

    row = 0  # Process first row (in R this would be row 1)

    # Calculate rho and append
    rho_result = pychnosz.water("rho", T=temp_degC+273.15, P=pressure_bar, messages=False)
    df.loc[df.index[row], "rho"] = rho_result / 1000  # convert rho to g/cm3

    # Handle redox block
    df["redox_flag"] = redox_flag
    df["redox_value"] = pd.NA
    df["redox_unit"] = pd.NA

    this_redox_flag = redox_flag
    assigned = False

    # use gaseous O2 in aqueous species block
    if redox_flag == -3:
        redox_col_index = [col for col in df.columns if "O2(g)" in col]
        redox_col_index = [col for col in redox_col_index if col.startswith("O2(g)")]

        if len(redox_col_index) > 0:
            redox_col = redox_col_index[0]
            if not pd.isna(df.iloc[row][redox_col]):
                unit = redox_col.split("_")[1] if "_" in redox_col else ""
                this_redox_value = unit  # subheader should always be "Hetero. equil."
                this_redox_unit = "using O2(g) in aqueous species block"
            else:
                vprint(f"Warning: non-numeric 'O2(g)' value in sample {df.iloc[row]['Sample']}. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True
        else:
            if not warned_about_redox_column:
                vprint(f"Warning: no 'O2(g)' column found. Resorting to using " +
                      f"Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                warned_about_redox_column = True
            this_redox_flag = 0
            this_redox_value = f"{default_logfO2:.4E}"
            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
            assigned = True

    # specify pe in pe units
    if redox_flag == -2:
        redox_col_index = [col for col in df.columns if col.startswith("pe_")]
        if len(redox_col_index) > 0:
            redox_col = redox_col_index[0]
            if not pd.isna(df.iloc[row][redox_col]):
                this_redox_value = f"{float(df.iloc[row][redox_col]):.4E}"
                this_redox_unit = "pe, pe units"
            else:
                vprint(f"Warning: non-numeric pe value in sample {df.iloc[row]['Sample']}. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True
        else:
            if not warned_about_redox_column:
                vprint(f"Warning: no 'pe' column found. Resorting to using " +
                      f"Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                warned_about_redox_column = True
            this_redox_flag = 0
            this_redox_value = f"{default_logfO2:.4E}"
            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
            assigned = True

    # specify Eh in volts
    if redox_flag == -1:
        redox_col_index = [col for col in df.columns if "Eh_volts" in col]
        if len(redox_col_index) > 0:
            redox_col = redox_col_index[0]
            if not pd.isna(df.iloc[row][redox_col]):
                this_redox_value = f"{float(df.iloc[row][redox_col]):.4E}"
                this_redox_unit = "Eh, volts"
            else:
                vprint(f"Warning: non-numeric Eh value in sample {df.iloc[row]['Sample']}. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True
        else:
            if not warned_about_redox_column:
                vprint(f"Warning: no 'Eh' column found with 'volts' units. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                warned_about_redox_column = True
            this_redox_flag = 0
            this_redox_value = f"{default_logfO2:.4E}"
            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
            assigned = True

    if redox_flag == 0 and assigned == False:
        redox_col_index = [col for col in df.columns if "logfO2" in col]
        if len(redox_col_index) > 0:
            redox_col = redox_col_index[0]
            if not pd.isna(df.iloc[row][redox_col]):
                this_redox_value = f"{float(df.iloc[row][redox_col]):.4E}"
                this_redox_unit = "Log fO2 (log bars) [from logfO2 column]"
            else:
                vprint(f"Warning: non-numeric logfO2 value in sample {df.iloc[row]['Sample']}. " +
                      f"Resorting to using a logfO2 value of {default_logfO2}",
                      verbose=verbose)
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True
        else:
            if not warned_about_redox_column:
                vprint("Warning: no 'logfO2' column found. Attempting to find a " +
                      "column for aqueous O2 to estimate logfO2 at sample temperature and pressure...",
                      verbose=verbose)
                warned_about_redox_column = True

            redox_col_index = [col for col in df.columns if re.match(r'^O2,AQ|^O2_', col)]
            if len(redox_col_index) > 0:
                redox_col = redox_col_index[0]
                # Treat empty string as NaN
                redox_val = df.iloc[row][redox_col]
                if redox_val == '':
                    redox_val = np.nan
                if not pd.isna(redox_val):
                    header_name = redox_col
                    this_header_unit = header_unit(header_name)

                    if this_header_unit in acceptable_units + EQ3_jflags:
                        O2_molal = uc_molal(value=float(redox_val),
                                            chemical="O2", unit=this_header_unit)

                        if pd.isna(temp_degC) or pd.isna(pressure_bar):
                            vprint(f"Warning: non-numeric temperature or pressure value in sample " +
                                  f"{df.iloc[row]['Sample']}. Resorting to using Log fO2 (log bars) " +
                                  f"with a value of {default_logfO2}",
                                  verbose=verbose)
                            this_redox_flag = 0
                            this_redox_value = f"{default_logfO2:.4E}"
                            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                            assigned = True
                        else:
                            logK = pychnosz.subcrt(
                                ["O2", "O2"],
                                [-1, 1],
                                ["aq", "gas"],
                                property="logK",
                                T=temp_degC,
                                P=pressure_bar,
                                exceed_rhomin=True,
                                exceed_Ttr=True,
                                messages=False,
                                show=False,
                            ).out['logK'][0]
                            logfO2 = np.log10(O2_molal * 10**logK)

                            if not np.isfinite(logfO2):
                                logfO2 = default_logfO2
                                this_redox_value = f"{default_logfO2:.4E}"
                                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                            else:
                                this_redox_unit = "Log fO2 (log bars) [calculated from O2(aq) = O2(g) at temperature and pressure of sample]"

                            assigned = True
                            this_redox_value = f"{logfO2:.4E}"
                    else:
                        vprint(f"Warning: column found for aqueous O2, but units are not recognized. " +
                              f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2} for all samples.",
                              verbose=verbose)
                        this_redox_flag = 0
                        this_redox_value = f"{default_logfO2:.4E}"
                        this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                        assigned = True
                else:
                    vprint(f"Warning: non-numeric aqueous O2 value in sample {df.iloc[row]['Sample']}. " +
                          f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                          verbose=verbose)
                    this_redox_flag = 0
                    this_redox_value = f"{default_logfO2:.4E}"
                    this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                    assigned = True
            else:
                if not warned_about_redox_column:
                    vprint(f"Warning: a column for aqueous O2 was not found. Resorting to using " +
                          f"Log fO2 (log bars) with a value of {default_logfO2}",
                          verbose=verbose)
                    warned_about_redox_column = True
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True

    # specify aux species redox couple
    if redox_flag == 1:
        redox_col_index = [col for col in df.columns if f"{redox_aux}_" in col]
        redox_col_index = [col for col in redox_col_index if col.startswith(redox_aux)]

        if len(redox_col_index) == 1:
            redox_col = redox_col_index[0]
            if not pd.isna(df.iloc[row][redox_col]):
                this_redox_value = f"{float(df.iloc[row][redox_col]):.4E}"
                this_redox_unit = f"{redox_aux} aux. sp."
            else:
                vprint(f"Warning: non-numeric {redox_aux} value in sample {df.iloc[row]['Sample']}. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                this_redox_flag = 0
                this_redox_value = f"{default_logfO2:.4E}"
                this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
                assigned = True
        elif len(redox_col_index) > 1:
            if not warned_about_redox_column:
                vprint(f"Warning: multiple matches for a {redox_aux} column found: {' '.join(redox_col_index)}. " +
                      f"Resorting to using Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                warned_about_redox_column = True
            this_redox_flag = 0
            this_redox_value = f"{default_logfO2:.4E}"
            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
            assigned = True
        else:
            if not warned_about_redox_column:
                vprint(f"Warning: no {redox_aux} column found. Resorting to using " +
                      f"Log fO2 (log bars) with a value of {default_logfO2}",
                      verbose=verbose)
                warned_about_redox_column = True
            this_redox_flag = 0
            this_redox_value = f"{default_logfO2:.4E}"
            this_redox_unit = f"Log fO2 (log bars) [default = {default_logfO2}]"
            assigned = True

    # append redox values
    df.loc[df.index[row], "redox_flag"] = this_redox_flag
    df.loc[df.index[row], "redox_value"] = this_redox_value
    df.loc[df.index[row], "redox_unit"] = this_redox_unit

    # check that the redox state of the sample is within the stability region of water
    logK_H2 = None
    if this_redox_flag != 1 and water_model != "DEW":
        try:
            logK_H2 = pychnosz.subcrt(
                ["H2O", "O2", "H2"],
                [-1, 0.5, 1],
                ["liq", "gas", "gas"],
                T=temp_degC+273.15,
                P=pressure_bar,
                property="logK",
                convert=False,
                messages=False,
                show=False,
            ).out['logK'][0]
        except Exception:
            # H2(g) may be absent from a reduced custom database (e.g., when
            # gases are excluded); the stability check is skipped in that case
            logK_H2 = None

    if logK_H2 is not None:
        T = temp_degC
        P = pressure_bar

        # reduction
        logaH2O = 0  # a good starting guess is 0 before actually doing the speciation...
        logfH2 = logaH2O
        logK = logK_H2
        logfO2_red = 2 * (logK - logfH2 + logaH2O)

        # oxidation
        logfO2_ox = 0

        if charge_balance_on != "H+" and (this_redox_flag == -1 or this_redox_flag == -2):
            pH = float(df.iloc[row]["H+_pH"])
            Eh_ox = pychnosz.convert(logfO2_ox, 'E0', T=T+273.15, P=P, pH=pH, logaH2O=logaH2O, messages=False)
            pe_ox = pychnosz.convert(pychnosz.convert(logfO2_ox, 'E0', T=T+273.15, P=P, pH=pH, logaH2O=logaH2O, messages=False), "pe", T=T+273.15, messages=False)
            Eh_red = pychnosz.convert(logfO2_red, 'E0', T=T+273.15, P=P, pH=pH, logaH2O=logaH2O, messages=False)
            pe_red = pychnosz.convert(pychnosz.convert(logfO2_red, 'E0', T=T+273.15, P=P, pH=pH, logaH2O=logaH2O, messages=False), "pe", T=T+273.15, messages=False)

            if this_redox_flag == -1:  # Eh
                if float(this_redox_value) <= Eh_red:
                    vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too reduced) " +
                          f"which may cause the speciation calculation to fail. " +
                          f"The sample has an Eh of {round(float(this_redox_value), 2)} volts, which is below the {round(Eh_red, 2)} " +
                          f"volt threshold for stable water at this temperature, pressure, and pH.")
                elif float(this_redox_value) >= Eh_ox:
                    vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too oxidized) " +
                          f"which may cause the speciation calculation to fail. " +
                          f"The sample has an Eh of {round(float(this_redox_value), 2)} volts, which is above the {round(Eh_ox, 2)} " +
                          f"volt threshold for stable water at this temperature, pressure, and pH.")
            else:  # pe
                if float(this_redox_value) <= pe_red:
                    vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too reduced) " +
                          f"which may cause the speciation calculation to fail. " +
                          f"The sample has a pe of {round(float(this_redox_value), 2)}, which is below the {round(pe_red, 2)} " +
                          f"pe threshold for stable water at this temperature, pressure, and pH.")
                elif float(this_redox_value) >= pe_ox:
                    vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too oxidized) " +
                          f"which may cause the speciation calculation to fail. " +
                          f"The sample has a pe of {round(float(this_redox_value), 2)}, which is above the {round(pe_ox, 2)} " +
                          f"pe threshold for stable water at this temperature, pressure, and pH.")
        elif this_redox_flag == 0:  # logfO2
            if float(this_redox_value) <= logfO2_red:
                vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too reduced) " +
                      f"which may cause the speciation calculation to fail. " +
                      f"The sample has a logfO2 of {round(float(this_redox_value), 2)}, which is below the {round(logfO2_red, 2)} " +
                      f"logfO2 threshold for stable water at this temperature and pressure.")
            elif float(this_redox_value) >= logfO2_ox:
                vprint(f"Warning: the sample {df.iloc[row]['Sample']} may be outside of the stability region of water (too oxidized) " +
                      f"which may cause the speciation calculation to fail. " +
                      f"The sample has a logfO2 of {round(float(this_redox_value), 2)}, which is above the {round(logfO2_ox, 2)} " +
                      f"logfO2 threshold for stable water at this temperature and pressure.")
        elif this_redox_flag == -3:  # O2(g)
            # there should be no reason to check this because O2(g) must be set to heterogeneous equilibrium
            # with minerals. Assumedly, this keeps samples from becoming too reduced or oxidized.
            pass

    if charge_balance_on == "none":
        eq3_cb_block = "\n|Electrical balancing option (iebal3):                                         |\n" + \
                      "|  [x] ( 0) No balancing is done                                               |\n" + \
                      "|  [ ] ( 1) Balance on species |None                    | (uebal)              |"
    else:
        eq3_cb_block = "\n|Electrical balancing option (iebal3):                                         |\n" + \
                      "|  [ ] ( 0) No balancing is done                                               |\n" + \
                      f"|  [x] ( 1) Balance on species |{charge_balance_on:<24}| (uebal)              |"

    eq3_filename = f"{df.index[row]}.3i"

    eq3_header1 = "|------------------------------------------------------------------------------|\n" + \
                 "| Title                  | (utitl(n))                                          |\n" + \
                 "|------------------------------------------------------------------------------|\n" + \
                 "|                                                                              |\n" + \
                 "|"

    eq3_samplename = f"Sample: {df.iloc[row]['Sample']:<70}"

    eq3_header2 = "|\n" + \
                 "|                                                                              |\n" + \
                 "|------------------------------------------------------------------------------|\n" + \
                 "|Special Basis Switches (for model definition only)       | (nsbswt)           |\n" + \
                 "|------------------------------------------------------------------------------|\n" + \
                 "|Replace |None                                            | (usbsw(1,n))       |\n" + \
                 "|   with |None                                            | (usbsw(2,n))       |\n" + \
                 "|------------------------------------------------------------------------------|\n" + \
                 "|Temperature (C)         | "

    eq3_temperature = f"{temp_degC:.5E}"

    if pressure_override:
        eq3_header3 = "| (tempc)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Pressure option (jpres3):                                                     |\n" + \
                     "|  [ ] ( 0) Data file reference curve value                                    |\n" + \
                     "|  [ ] ( 1) 1.013-bar/steam-saturation curve value                             |\n" + \
                     f"|  [x] ( 2) Value (bars) |{pressure_bar:>12.5E}| (press)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Density (g/cm3)         | "
    elif strict_minimum_pressure == False or pressure_bar > minimum_pressure:
        eq3_header3 = "| (tempc)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Pressure option (jpres3):                                                     |\n" + \
                     "|  [x] ( 0) Data file reference curve value                                    |\n" + \
                     "|  [ ] ( 1) 1.013-bar/steam-saturation curve value                             |\n" + \
                     "|  [ ] ( 2) Value (bars) | 1.00000E+00| (press)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Density (g/cm3)         | "
    else:
        eq3_header3 = "| (tempc)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Pressure option (jpres3):                                                     |\n" + \
                     "|  [ ] ( 0) Data file reference curve value                                    |\n" + \
                     "|  [ ] ( 1) 1.013-bar/steam-saturation curve value                             |\n" + \
                     f"|  [x] ( 2) Value (bars) |{minimum_pressure:>12.5E}| (press)                                |\n" + \
                     "|------------------------------------------------------------------------------|\n" + \
                     "|Density (g/cm3)         | "

    eq3_density = f"{df.iloc[row]['rho']:.5E}"

    eq3_header4 = "| (rho)                                  |\n" + \
                 "|------------------------------------------------------------------------------|\n" + \
                 "|Total dissolved solutes option (itdsf3):                                      |\n" + \
                 "|  [x] ( 0) Value (mg/kg.sol) | 0.00000E+00| (tdspkg)                          |\n" + \
                 "|  [ ] ( 1) Value (mg/L)      | 0.00000E+00| (tdspl)                           |\n" + \
                 "|------------------------------------------------------------------------------|"

    # cb_block will be pasted here (one for all samples)

    eq3_header5 = "\n|------------------------------------------------------------------------------|\n" + \
                 "|Default redox constraint (irdxc3):                                            |"

    default_redox_minus3 = "\n|  [ ] (-3) Use O2(g) line in the aqueous basis species block                  |"
    default_redox_minus2 = "|  [ ] (-2) pe (pe units)      | 0.00000E+00| (pei)                            |"
    default_redox_minus1 = "|  [ ] (-1) Eh (volts)         | 0.00000E+00| (ehi)                            |"
    default_redox_0 = "|  [ ] ( 0) Log fO2 (log bars) | 0.00000E+00| (fo2lgi)                         |"
    default_redox_1 = "|  [ ] ( 1) Couple (aux. sp.)  |None                    | (uredox)             |"

    if df.iloc[row]["redox_flag"] == -3:
        default_redox_minus3 = "\n|  [x] (-3) Use O2(g) line in the aqueous basis species block                  |"
    elif df.iloc[row]["redox_flag"] == -2:
        default_redox_minus2 = f"|  [x] (-2) pe (pe units)      |{df.iloc[row]['redox_value']:>12}| (pei)                            |"
    elif df.iloc[row]["redox_flag"] == -1:
        default_redox_minus1 = f"|  [x] (-1) Eh (volts)         |{df.iloc[row]['redox_value']:>12}| (ehi)                            |"
    elif df.iloc[row]["redox_flag"] == 0:
        default_redox_0 = f"|  [x] ( 0) Log fO2 (log bars) |{df.iloc[row]['redox_value']:>12}| (fo2lgi)                         |"
    elif df.iloc[row]["redox_flag"] == 1:
        default_redox_1 = f"|  [x] ( 1) Couple (aux. sp.)  | {redox_aux:<23}| (uredox)             |"
    else:
        raise ValueError(f"Error when writing .3i file for sample {df.iloc[row]['Sample']}. " +
                        "Redox flag was not recognized! Choose -3, -2, -1, 0, or 1.")

    redox_block = "\n".join([default_redox_minus3, default_redox_minus2,
                            default_redox_minus1, default_redox_0, default_redox_1])

    eq3_header6 = "\n|------------------------------------------------------------------------------|\n" + \
                 "|Aqueous Basis Species/Constraint Species        |Conc., etc. |Units/Constraint|\n" + \
                 "| (uspeci(n)/ucospi(n))                          | (covali(n))|(ujf3(jflgi(n)))|\n" + \
                 "|------------------------------------------------------------------------------|"

    if allowed_aq_block_species[0] == "all":
        allowed_aq_block_species = [header_species(col) for col in df.columns]

    # handle aqueous block
    aqueous_lines = []
    for column in df.columns:
        col_val = df.iloc[row][column]

        # Treat empty strings as NaN
        if col_val == '':
            col_val = np.nan

        if suppress_missing and pd.isna(col_val):
            # Use string "0" if column has string dtype, otherwise use numeric 0
            zero_val = "0" if pd.api.types.is_string_dtype(df[column]) else 0
            df.loc[df.index[row], column] = zero_val
            col_val = 0

        if not pd.isna(col_val):
            species_name = header_species(column)
            # Check if the column should be excluded - handle both the species name
            # and columns that start with an excluded name (to handle headers with underscores)
            should_exclude = False
            for exc in exclude:
                if exc and (species_name == exc or column.startswith(exc + "_")):
                    should_exclude = True
                    break

            if not should_exclude and (species_name in allowed_aq_block_species):
                species_value = df.iloc[row][column]

                # EQ3 won't balance on a species if its concentration is 0 so
                # change it to a very small non-zero value
                if charge_balance_on == species_name and float(species_value) == 0:
                    species_value = 1e-18

                species_unit = header_unit(column, keepcase=True)

                if species_unit not in EQ3_jflags:
                    if species_unit.lower() == "ppb":
                        species_value = float(species_value) / 1000
                        species_unit = "mg/L"
                    elif species_unit.lower() == "ppm":
                        species_unit = "mg/L"
                    else:
                        print(f"Error creating .3i file: {species_unit} " +
                             "is not a recognized aqueous block jflag. Try checking " +
                             "capitalization and spelling to match one of the following: " +
                             " ".join(EQ3_jflags))

                if species_unit == "Hetero. equil.":
                    species_value_split = str(species_value).split(" ")
                    if len(species_value_split) == 2:
                        # for gases
                        species_value = species_value_split[0]
                        hetero_equil_species = species_value_split[1]
                    else:
                        # for minerals
                        hetero_equil_species = species_value
                        species_value = 0

                species_value_str = f"{float(species_value):>12.5E}"
                this_aq_line = f"\n|{species_name:<48}|{species_value_str}|{species_unit:<16}|"

                # handle additional line for 'Hetero. equil.' jflag
                if species_unit == "Hetero. equil.":
                    this_aq_line += f"\n|->|{hetero_equil_species:<45}| {'(ucospi(n))':<28}|"

                aqueous_lines.append(this_aq_line)

    # suppressing species
    for species in suppress:
        this_aq_line = f"\n|{species:<48}| {0:.5E}|{'Suppressed':<16}|"
        aqueous_lines.append(this_aq_line)

    aqueous_block = "".join(aqueous_lines)

    eq3_ender1 = "\n|------------------------------------------------------------------------------|\n" + \
                "* Valid jflag strings (ujf3(jflgi(n))) are:                                    *\n" + \
                "*    Suppressed          Molality            Molarity                          *\n" + \
                "*    mg/L                mg/kg.sol           Alk., eq/kg.H2O                   *\n" + \
                "*    Alk., eq/L          Alk., eq/kg.sol     Alk., mg/L CaCO3                  *\n" + \
                "*    Alk., mg/L HCO3-    Log activity        Log act combo                     *\n" + \
                "*    Log mean act        pX                  pH                                *\n" + \
                "*    pHCl                pmH                 pmX                               *\n" + \
                "*    Hetero. equil.      Homo. equil.        Make non-basis                    *\n" + \
                "*------------------------------------------------------------------------------*\n" + \
                "|Create Ion Exchangers  | (net)                                                |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Advisory: no exchanger creation blocks follow on this file.                   |\n" + \
                "|Option: on further processing (writing a PICKUP file or running XCON3 on the  |\n" + \
                "|present file), force the inclusion of at least one such block (qgexsh):       |\n" + \
                "|  [ ] (.true.)                                                                |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Ion Exchanger Compositions      | (neti)                                      |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Exchanger phase |None                    | (ugexpi(n))                        |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|->|Moles/kg.H2O    |  0.0000    | (cgexpi(n))                                 |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|->|Exchange site   |None    | (ugexji(j,n))                                   |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|--->|Exchange species        |Eq. frac.   | (this is a table header)          |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|--->|None                    | 0.00000E+00| (ugexsi(i,j,n), egexsi(i,j,n))    |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Solid Solution Compositions     | (nxti)                                      |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Solid Solution          |None                    | (usoli(n))                 |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|->|Component               |Mole frac.  | (this is a table header)            |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|->|None                    | 0.00000E+00| (umemi(i,n), xbari(i,n))            |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Alter/Suppress Options  | (nxmod)                                             |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Species                                         |Option          |Alter value |\n" + \
                "| (uxmod(n))                                     |(ukxm(kxmod(n)))| (xlkmod(n))|\n" + \
                "|------------------------------------------------------------------------------|"

    alter_block_lines = []
    if len(alter_options) > 0:
        for species, option in alter_options.items():
            if len(option) == 2:
                alter_line = f"\n|{species:<48}| {option[0]:<15}|{float(option[1]):>12.5E}|"
                alter_block_lines.append(alter_line)
        alter_block = "".join(alter_block_lines)
    else:
        alter_block = "\n|None                                            |None            | 0.00000E+00|"

    if get_solid_solutions:
        ss_checkbox_ignore = " "
        ss_checkbox_permit = "x"
    else:
        ss_checkbox_ignore = "x"
        ss_checkbox_permit = " "

    eq3_ender2 = "\n|------------------------------------------------------------------------------|\n" + \
                "* Valid alter/suppress strings (ukxm(kxmod(n))) are:                           *\n" + \
                "*    Suppress            Replace             AugmentLogK                       *\n" + \
                "*    AugmentG                                                                  *\n" + \
                "*------------------------------------------------------------------------------*\n" + \
                "|Iopt Model Option Switches (\"( 0)\" marks default choices)                     |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopt(4) - Solid Solutions:                                                    |\n" + \
                f"|  [{ss_checkbox_ignore}] ( 0) Ignore                                                             |\n" + \
                f"|  [{ss_checkbox_permit}] ( 1) Permit                                                             |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopt(11) - Auto Basis Switching in pre-N-R Optimization:                      |\n" + \
                "|  [x] ( 0) Turn off                                                           |\n" + \
                "|  [ ] ( 1) Turn on                                                            |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopt(17) - PICKUP File Options:                                               |\n" + \
                "|  [ ] (-1) Don't write a PICKUP file                                          |\n" + \
                "|  [x] ( 0) Write a PICKUP file                                                |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopt(19) - Advanced EQ3NR PICKUP File Options:                                |\n" + \
                "|  [ ] ( 0) Write a normal EQ3NR PICKUP file                                   |\n" + \
                "|  [ ] ( 1) Write an EQ6 INPUT file with Quartz dissolving, relative rate law  |\n" + \
                "|  [ ] ( 2) Write an EQ6 INPUT file with Albite dissolving, TST rate law       |\n" + \
                "|  [x] ( 3) Write an EQ6 INPUT file with Fluid 1 set up for fluid mixing       |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Iopg Activity Coefficient Option Switches (\"( 0)\" marks default choices)      |\n" + \
                "|------------------------------------------------------------------------------|"

    davies_box = " "
    bdot_box = " "
    pitzer_box = " "
    if activity_model == "davies":
        davies_box = "x"
    elif activity_model == "b-dot":
        bdot_box = "x"
    elif activity_model == "pitzer":
        pitzer_box = "x"
    else:
        print("Error: activity model is not recognized. Try 'b-dot', 'davies', or 'pitzer'")

    eq3_ender3 = "\n|iopg(1) - Aqueous Species Activity Coefficient Model:                         |\n" + \
                f"|  [{davies_box}] (-1) The Davies equation                                                |\n" + \
                f"|  [{bdot_box}] ( 0) The B-dot equation                                                 |\n" + \
                f"|  [{pitzer_box}] ( 1) Pitzer's equations                                                 |\n" + \
                "|  [ ] ( 2) HC + DH equations                                                  |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopg(2) - Choice of pH Scale (Rescales Activity Coefficients):                |\n" + \
                "|  [ ] (-1) \"Internal\" pH scale (no rescaling)                                 |\n" + \
                "|  [x] ( 0) NBS pH scale (uses the Bates-Guggenheim equation)                  |\n" + \
                "|  [ ] ( 1) Mesmer pH scale (numerically, pH = -log m(H+))                     |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Iopr Print Option Switches (\"( 0)\" marks default choices)                     |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(1) - Print All Species Read from the Data File:                          |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print                                                              |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(2) - Print All Reactions:                                                |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print the reactions                                                |\n" + \
                "|  [ ] ( 2) Print the reactions and log K values                               |\n" + \
                "|  [ ] ( 3) Print the reactions, log K values, and associated data             |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(3) - Print the Aqueous Species Hard Core Diameters:                      |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print                                                              |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(4) - Print a Table of Aqueous Species Concentrations, Activities, etc.:  |\n" + \
                "|  [ ] (-3) Omit species with molalities < 1.e-8                               |\n" + \
                "|  [ ] (-2) Omit species with molalities < 1.e-12                              |\n" + \
                "|  [ ] (-1) Omit species with molalities < 1.e-20                              |\n" + \
                "|  [ ] ( 0) Omit species with molalities < 1.e-100                             |\n" + \
                "|  [x] ( 1) Include all species                                                |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(5) - Print a Table of Aqueous Species/H+ Activity Ratios:                |\n" + \
                "|  [ ] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print cation/H+ activity ratios only                               |\n" + \
                "|  [x] ( 2) Print cation/H+ and anion/H+ activity ratios                       |\n" + \
                "|  [ ] ( 3) Print ion/H+ activity ratios and neutral species activities        |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(6) - Print a Table of Aqueous Mass Balance Percentages:                  |\n" + \
                "|  [ ] (-1) Don't print                                                        |\n" + \
                "|  [x] ( 0) Print those species comprising at least 99% of each mass balance   |\n" + \
                "|  [ ] ( 1) Print all contributing species                                     |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(7) - Print Tables of Saturation Indices and Affinities:                  |\n" + \
                "|  [ ] (-1) Don't print                                                        |\n" + \
                "|  [x] ( 0) Print, omitting those phases undersaturated by more than 10 kcal   |\n" + \
                "|  [ ] ( 1) Print for all phases                                               |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(8) - Print a Table of Fugacities:                                        |\n" + \
                "|  [ ] (-1) Don't print                                                        |\n" + \
                "|  [x] ( 0) Print                                                              |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(9) - Print a Table of Mean Molal Activity Coefficients:                  |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print                                                              |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(10) - Print a Tabulation of the Pitzer Interaction Coefficients:         |\n" + \
                "|  [ ] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print a summary tabulation                                         |\n" + \
                "|  [x] ( 2) Print a more detailed tabulation                                   |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iopr(17) - PICKUP file format (\"W\" or \"D\"):                                   |\n" + \
                "|  [x] ( 0) Use the format of the INPUT file                                   |\n" + \
                "|  [ ] ( 1) Use \"W\" format                                                     |\n" + \
                "|  [ ] ( 2) Use \"D\" format                                                     |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Iodb Debugging Print Option Switches (\"( 0)\" marks default choices)           |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iodb(1) - Print General Diagnostic Messages:                                  |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print Level 1 diagnostic messages                                  |\n" + \
                "|  [ ] ( 2) Print Level 1 and Level 2 diagnostic messages                      |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iodb(3) - Print Pre-Newton-Raphson Optimization Information:                  |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print summary information                                          |\n" + \
                "|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |\n" + \
                "|  [ ] ( 3) Print more detailed information (including matrix equations)       |\n" + \
                "|  [ ] ( 4) Print most detailed information (including activity coefficients)  |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iodb(4) - Print Newton-Raphson Iteration Information:                         |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print summary information                                          |\n" + \
                "|  [ ] ( 2) Print detailed information (including the beta and del vectors)    |\n" + \
                "|  [ ] ( 3) Print more detailed information (including the Jacobian)           |\n" + \
                "|  [ ] ( 4) Print most detailed information (including activity coefficients)  |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|iodb(6) - Print Details of Hypothetical Affinity Calculations:                |\n" + \
                "|  [x] ( 0) Don't print                                                        |\n" + \
                "|  [ ] ( 1) Print summary information                                          |\n" + \
                "|  [ ] ( 2) Print detailed information                                         |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Numerical Parameters                                                          |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "| Beta convergence tolerance      | 0.00000E+00| (tolbt)                       |\n" + \
                "| Del convergence tolerance       | 0.00000E+00| (toldl)                       |\n" + \
                "| Max. Number of N-R Iterations   |   0        | (itermx)                      |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Ordinary Basis Switches (for numerical purposes only)    | (nobswt)           |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Replace |None                                            | (uobsw(1,n))       |\n" + \
                "|   with |None                                            | (uobsw(2,n))       |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|Sat. flag tolerance     | 0.00000E+00| (tolspf)                               |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                f"|Aq. Phase Scale Factor  |{float(aq_scale):>12.5E}| (scamas)                               |\n" + \
                "|------------------------------------------------------------------------------|\n" + \
                "|End of problem                                                                |\n" + \
                "|------------------------------------------------------------------------------|"

    this_file = "".join([eq3_header1, eq3_samplename, eq3_header2,
                        eq3_temperature, eq3_header3, eq3_density,
                        eq3_header4, eq3_cb_block, eq3_header5,
                        redox_block, eq3_header6, aqueous_block,
                        eq3_ender1, alter_block, eq3_ender2, eq3_ender3])

    with open(f"{input_dir}/{eq3_filename}", 'w') as f:
        f.write(this_file)

    return warned_about_redox_column
