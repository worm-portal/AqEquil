import re
import numpy as np
import pychnosz


def calc_bdot(T):
    """
    Calculate bdot parameter at a given temperature.

    GB notes:
    The equation used by dbcreate to approximate the curve in Fig 3
    of Helgeson 1969 results in numbers that are close to, but not
    exactly the same as those in data0.jus:

    Bdot parameter grid:
     0.0376   0.0443   0.0505   0.0529   0.0479   0.0322   0.0000   0.0000  # from dbcreate
     0.0374   0.0430   0.0460   0.0470   0.0470   0.0340   0.0000   0.0000  # from data0.jus

    Close but not exact! data0.jus is closer to what is depicted in Fig 3 of Helgeson 1969.
    Not sure what other equation to use, though. Will keep the dbcreate equation for now.
    TODO: look into alternative equations.

    Parameters
    ----------
    T : float or array-like
        Temperature in degrees Celsius

    Returns
    -------
    float or array-like
        Bdot parameter value(s)
    """
    b1 =  0.0374
    b2 =  1.3569e-4
    b3 =  2.6411e-7
    b4 = -4.6103e-9

    result = b1 + b2*T + b3*(T-25.0)**2 + b4*(T-25.0)**3

    # Return 0 if T >= 300, otherwise return result
    if isinstance(T, (list, tuple, np.ndarray)):
        return np.where(np.array(T) >= 300, 0, result)
    else:
        return 0 if T >= 300 else result


def fill_data0_head(data0_template, db, grid_temps, grid_press,
                    water_model, activity_model):
    """
    Fill the header section of a data0 file with calculated thermodynamic parameters.

    Parameters
    ----------
    data0_template : str or list
        Template for data0 file (can be string or list of lines)
    db : str
        Database name
    grid_temps : list or array-like
        Temperature grid in degrees Celsius
    grid_press : list or array-like
        Pressure grid in bars
    water_model : str
        Water model name (e.g., "SUPCRT92")
    activity_model : str
        Activity coefficient model ("b-dot", "davies", or "pitzer")

    Returns
    -------
    str
        Data0 template with filled header section
    """
    # Convert to single string if list
    if isinstance(data0_template, list):
        data0_template = "\n".join(data0_template)

    # Make sure the TP grid is in the correct form, esp. for single TP values
    grid_temps = np.atleast_1d(grid_temps).flatten()
    grid_press = np.atleast_1d(grid_press).flatten()
    grid_temps_original = grid_temps.copy()

    if len(grid_temps) == 1:
        # Only the first T value is valid, but this is needed for EQ3
        grid_temps = grid_temps[0] + 10 * np.arange(8)
        grid_press = np.repeat(grid_press[0], len(grid_temps))

    # Calculate debye huckel a and b parameters for the grid
    if len(grid_temps_original) >= 8:
        A_DH_grid = pychnosz.water("A_DH", T=273.15+grid_temps, P=grid_press, messages=False)
        B_DH_grid = pychnosz.water("B_DH", T=273.15+grid_temps, P=grid_press, messages=False) * 1e-8
    elif len(grid_temps_original) == 1:
        A_DH_grid = pychnosz.water("A_DH", T=273.15+grid_temps[0], P=grid_press[0], messages=False)
        A_DH_grid = np.concatenate([[A_DH_grid], np.zeros(7)])
        B_DH_grid = pychnosz.water("B_DH", T=273.15+grid_temps[0], P=grid_press[0], messages=False) * 1e-8
        B_DH_grid = np.concatenate([[B_DH_grid], np.zeros(7)])

    # Convert to numpy arrays if they're scalars
    A_DH_grid = np.atleast_1d(A_DH_grid)
    B_DH_grid = np.atleast_1d(B_DH_grid)

    if activity_model == "pitzer":
        # The Debye-Huckel Aphi parameter used by Pitzer's equations is
        # related to A_gamma (in log10 units) by Aphi = A_gamma * ln(10) / 3
        Aphi_DH_grid = (A_DH_grid * np.log(10)) / 3
        Aphi_DH_grid_f = [f"{val:.4f}" for val in Aphi_DH_grid]

    # Ensure that the number of characters expressing pressure allows for spaces between values in the header
    # e.g., avoid a line like: 50000.000050000.000050000.000050000.0000
    # This is important for DEW calculations
    if len(f"{grid_press[0]:.4f}") == 10:
        pdig = 3
    elif len(f"{grid_press[0]:.4f}") < 10:
        pdig = 4
    else:
        pdig = 4  # Default case

    # Format grid values
    grid_temps_f = [f"{val:.4f}" for val in grid_temps]
    grid_press_f = [f"{val:.{pdig}f}" for val in grid_press]
    A_DH_grid_f = [f"{val:.4f}" for val in A_DH_grid]
    B_DH_grid_f = [f"{val:.4f}" for val in B_DH_grid]

    # Calculate bdot parameter
    if len(grid_temps_original) == 1:
        bdot_val = calc_bdot(grid_temps[0])
        bdot_grid_f = [f"{bdot_val:.4f}"] + ["0.0000"] * 7
    else:
        bdot_grid = calc_bdot(grid_temps)
        bdot_grid_f = [f"{val:.4f}" for val in bdot_grid]

    # cco2 (coefficients for the drummond (1981) polynomial)
    # GB note: might not change with T or P?
    # Examination of various data0 files seems to indicate that DBcreate does not change these values.

    # Calculate the "log k for eh reaction" grid.
    # From eq. 9 in EQPT user manual (version 7.0) by Daveler and Wolery:
    if len(grid_temps_original) >= 8:
        subcrt_result = pychnosz.subcrt(
            ["H2O", "O2", "e-", "H+"],
            [-2, 1, 4, 4],
            ["liq", "gas", "aq", "aq"],
            T=grid_temps,
            P=np.round(grid_press, 9),
            exceed_rhomin=True,
            exceed_Ttr=True,
            messages=False,
            show=False,
        )
        logK_Eh_vals = subcrt_result.out["logK"]
        
    elif len(grid_temps_original) == 1:
        subcrt_result = pychnosz.subcrt(
            ["H2O", "O2", "e-", "H+"],
            [-2, 1, 4, 4],
            ["liq", "gas", "aq", "aq"],
            T=grid_temps[0],
            P=np.round(grid_press[0], 9),
            exceed_rhomin=True,
            exceed_Ttr=True,
            messages=False,
            show=False,
        )
        logK_Eh_vals = subcrt_result.out["logK"]
        # Ensure logK_Eh_vals is an array before concatenating
        logK_Eh_vals = np.atleast_1d(logK_Eh_vals)
        logK_Eh_vals = np.concatenate([logK_Eh_vals, np.zeros(7)])

    # Convert to array if scalar
    logK_Eh_vals = np.atleast_1d(logK_Eh_vals)
    logk_grid_f = [f"{val:.4f}" for val in logK_Eh_vals]
    
    # Build grids with proper spacing (10 characters per value)
    tempgrid = ["     "]
    presgrid = ["     "]
    A_DHgrid = ["     "]
    B_DHgrid = ["     "]
    bdotgrid = ["     "]
    logkgrid = ["     "]
    Aphi_DHgrid = ["     "]

    for i in range(len(grid_temps)):
        tempgrid.append(" " * (10 - len(grid_temps_f[i])) + grid_temps_f[i])
        presgrid.append(" " * (10 - len(grid_press_f[i])) + grid_press_f[i])
        A_DHgrid.append(" " * (10 - len(A_DH_grid_f[i])) + A_DH_grid_f[i])
        B_DHgrid.append(" " * (10 - len(B_DH_grid_f[i])) + B_DH_grid_f[i])
        bdotgrid.append(" " * (10 - len(bdot_grid_f[i])) + bdot_grid_f[i])
        logkgrid.append(" " * (10 - len(logk_grid_f[i])) + logk_grid_f[i])
        if activity_model == "pitzer":
            Aphi_DHgrid.append(" " * (10 - len(Aphi_DH_grid_f[i])) + Aphi_DH_grid_f[i])

        # Add newline every 6 values
        if (i + 1) % 6 == 0:
            tempgrid.append("\n     ")
            presgrid.append("\n     ")
            A_DHgrid.append("\n     ")
            B_DHgrid.append("\n     ")
            bdotgrid.append("\n     ")
            logkgrid.append("\n     ")
            if activity_model == "pitzer":
                Aphi_DHgrid.append("\n     ")

    tempgrid = "".join(tempgrid)
    presgrid = "".join(presgrid)
    A_DHgrid = "".join(A_DHgrid)
    B_DHgrid = "".join(B_DHgrid)
    bdotgrid = "".join(bdotgrid)
    logkgrid = "".join(logkgrid)
    if activity_model == "pitzer":
        Aphi_DHgrid = "".join(Aphi_DHgrid)
    
    # Insert minimum and maximum temperature values into data0 template
    temp_min_max_insertlines = r"\nTemperature limits \(degC\)\n.*?\ntemperatures\n"
    t_min = np.min(grid_temps)
    t_max = np.max(grid_temps)
    t_min_f = f"{t_min:.4f}"
    t_max_f = f"{t_max:.4f}"
    t_min_max = " " * (10 - len(t_min_f)) + t_min_f
    t_min_max = t_min_max + " " * (10 - len(t_max_f)) + t_max_f
    t_min_max = "     " + t_min_max
    data0_template = re.sub(temp_min_max_insertlines,
                           f"\nTemperature limits (degC)\n{t_min_max}\ntemperatures\n",
                           data0_template,
                           flags=re.DOTALL)

    # Insert temperature grid values into data0 template
    tempgrid_insertlines = r"\ntemperatures\n.*?\npressures\n"
    data0_template = re.sub(tempgrid_insertlines,
                           f"\ntemperatures\n{tempgrid}\npressures\n",
                           data0_template,
                           flags=re.DOTALL)

    # Insert pressure grid values into data0 template
    presgrid_insertlines = r"\npressures\n.*?\ndebye huckel a \(adh\)\n"
    data0_template = re.sub(presgrid_insertlines,
                           f"\npressures\n{presgrid}\ndebye huckel a (adh)\n",
                           data0_template,
                           flags=re.DOTALL)

    if activity_model == "b-dot" or activity_model == "davies":
        # Insert Debeye Huckel A and B parameter values into data0 template
        A_DHgrid_insertlines = r"\ndebye huckel a \(adh\)\n.*?\ndebye huckel b \(bdh\)\n"
        data0_template = re.sub(A_DHgrid_insertlines,
                               f"\ndebye huckel a (adh)\n{A_DHgrid}\ndebye huckel b (bdh)\n",
                               data0_template,
                               flags=re.DOTALL)

        B_DHgrid_insertlines = r"\ndebye huckel b \(bdh\)\n.*?\nbdot\n"
        data0_template = re.sub(B_DHgrid_insertlines,
                               f"\ndebye huckel b (bdh)\n{B_DHgrid}\nbdot\n",
                               data0_template,
                               flags=re.DOTALL)

        # Insert bdot grid values into data0 template
        bdotgrid_insertlines = r"\nbdot\n.*?\ncco2"
        data0_template = re.sub(bdotgrid_insertlines,
                               f"\nbdot\n{bdotgrid}\ncco2",
                               data0_template,
                               flags=re.DOTALL)

    elif activity_model == 'pitzer':
        # Insert Debeye Huckel A phi values into data0 template
        Aphi_DHgrid_insertlines = r"\ndebye huckel a \(adh\)\n.*?\nlog k for eh reaction\n"
        data0_template = re.sub(Aphi_DHgrid_insertlines,
                               f"\ndebye huckel aphi\n{Aphi_DHgrid}\nlog k for eh reaction\n",
                               data0_template,
                               flags=re.DOTALL)

    # Insert logk (eh) grid values into data0 template. The block that
    # follows is "bdot parameters" for B-dot/Davies data0 files and
    # "ca combinations" (the first Pitzer superblock) for Pitzer data0 files.
    logkgrid_insertlines = r"\nlog k for eh reaction\n.*?\n\+-+\n(bdot parameters|ca combinations)"
    logkgrid_end_insert = "\n+--------------------------------------------------------------------\n"
    data0_template = re.sub(logkgrid_insertlines,
                           lambda m: f"\nlog k for eh reaction\n{logkgrid}{logkgrid_end_insert}{m.group(1)}",
                           data0_template,
                           count=1,
                           flags=re.DOTALL)

    # Modify the data0 header lines
    desc = f"data0.{db}\nWater model: {water_model}\nTP points: {len(grid_temps_original)}"
    min_desc = "data0.min\nminimal working data0 file"
    data0_template = data0_template.replace(min_desc, desc)

    return data0_template  # returns a single string
