"""
Python implementation of create_data0.r
Converts thermodynamic data into EQ3/6 data0 format
"""

import re
import pandas as pd
import numpy as np
from pychnosz import makeup, mass, info, water

from .pitzer import build_pitzer_superblocks, pitzer_config_lines, DATA0_DELIMITER


# Helper functions

def vmessage(m, vlevel, verbose):
    """Print messages if 'verbose' setting >= vlevel of message."""
    if verbose >= vlevel:
        print(m)


def fillspace(string, nspaces, spaces_after=True):
    """
    Create a string of a certain length by adding additional spaces.
    e.g., "H2O" becomes "H2O    " if nspaces=7.
    Spaces can be added before the string by specifying spaces_after=False
    """
    try:
        if spaces_after:
            out = string + " " * (nspaces - len(string))
        else:
            out = " " * (nspaces - len(string)) + string
    except:
        out = ""
    return out


def s_d(x, k):
    """
    Specify how many decimals are printed.
    e.g. 12.433 becomes "12.4330" if k=4
    """
    return f"{float(x):.{k}f}".strip()


def create_data0(thermo_df,
                 element_df=None,
                 solid_solution_df=None,
                 db="",
                 water_model="SUPCRT92",
                 template="",
                 dissrxns=None,
                 basis_pref=None,
                 exceed_Ttr=False,
                 fixed_species=None,
                 redox_elem_states=None,
                 activity_model="b-dot",
                 pitzer_df=None,
                 verbose=1):
    """
    Main function to create a data0 file from thermodynamic data.

    Parameters
    ----------
    thermo_df : DataFrame
        Thermodynamic database
    element_df : DataFrame, optional
        Element database (not currently used, kept for compatibility)
    solid_solution_df : DataFrame, optional
        Solid solution database
    db : str
        Three-letter database code
    water_model : str
        Water model to use (SUPCRT92, IAPWS95, or DEW)
    template : str
        data0.min template content
    dissrxns : dict
        Dissociation reactions dictionary with 'basis_list' key
    basis_pref : list, optional
        Basis preference list (not currently used)
    exceed_Ttr : bool
        Whether to exceed transition temperatures
    fixed_species : list
        List of fixed species (H2O, H+, O2(g), water, Cl-, e-)
    activity_model : str
        Aqueous species activity coefficient model the data0 file will be
        used with: "b-dot" or "davies" produce a data0 file with a B-dot
        (hard core diameter) parameter block; "pitzer" produces a data0 file
        with Pitzer interaction parameter blocks.
    pitzer_df : DataFrame, optional
        Pitzer parameter database in the long CSV format described in
        `aqequil.pitzer`. Required when `activity_model` is "pitzer".
    verbose : int
        Verbosity level (0, 1, or 2)

    Returns
    -------
    str
        Formatted data0 file content
    """

    if fixed_species is None:
        fixed_species = ["H2O", "H+", "O2(g)", "water", "e-"]

    if dissrxns is None:
        dissrxns = {"basis_list": {}}

    # Set water model
    water(water_model, messages=False)

    # Initialize data structures
    azero_vec = {}
    neutral_ion_type_vec = {}
    dissociation_list = {}
    tag_vec = {}

    # Process each row in thermodynamic dataframe
    for i in range(len(thermo_df)):
        row = thermo_df.iloc[i]
        species_name = row["name"]

        # Look up azero bdot param
        azero_vec[species_name] = row["azero"]

        # Look up neutral ion type
        neutral_ion_type_vec[species_name] = row["neutral_ion_type"]

        # Look up dissociation reaction
        if pd.notna(row["dissrxn"]) and row["dissrxn"] != "" and species_name not in dissrxns:
            dissrxn = row["dissrxn"]
        elif species_name in dissrxns:
            # Use auto-balanced dissociation reaction
            dissrxn = dissrxns[species_name]
        elif row["tag"] == "basis":
            # Basis species don't need dissociation reactions
            dissrxn = None
        else:
            vmessage(f"Error: dissociation reaction could not be generated for:", 1, verbose)
            vmessage(species_name, 1, verbose)
            dissrxn = None

        # Parse dissociation reaction if not a basis species
        if row["tag"] != "basis" and dissrxn is not None:
            dissrxn_parts = dissrxn.split()
            dissrxn_names = dissrxn_parts[1::2]  # Get names (odd indices)
            dissrxn_coefs = [float(c) for c in dissrxn_parts[0::2]]  # Get coeffs (even indices)
            dissociation_list[species_name] = dict(zip(dissrxn_names, dissrxn_coefs))

        # Look up tag
        tag_vec[species_name] = row["tag"]

        # When redox is suppressed, auxiliary basis species representing pseudoelements
        # should be promoted to strict basis species
        if redox_elem_states and row["tag"] == "aux":
            # Check if this auxiliary species represents a redox-suppressed oxidation state
            # redox_elem_states format: {'Fe': {'Fe+2': 'Fejiip', 'Fe+3': 'Fejiiip'}}
            for element, ox_states_dict in redox_elem_states.items():
                # Check if species_name matches any of the oxidation state keys (e.g., 'Fe+2', 'Fe+3')
                if species_name in ox_states_dict.keys():
                    # This aux species represents a pseudoelement, promote to basis
                    tag_vec[species_name] = "basis"
                    vmessage(f"Promoting '{species_name}' from auxiliary to strict basis for redox suppression", 1, verbose)
                    # Remove dissociation reaction since this is now a strict basis species
                    if species_name in dissociation_list:
                        del dissociation_list[species_name]
                        vmessage(f"Removed dissociation reaction for promoted species '{species_name}'", 1, verbose)
                    break

    # Initialize vector of name differences between CHNOSZ and data0
    CHNOSZ_data0_name_diff = {}

    add_obigt_df = thermo_df.copy()

    # Remove raw line indicators '\r' in data0.min template
    data0_template = template.replace("\r", "")

    # Initialize list to store names of skipped species
    skipped_species = []

    # Loop through species in thermo_df
    for idx in range(len(thermo_df)):
        entry = thermo_df.iloc[idx]
        name = entry["name"]

        # Skip certain default basis species
        if name in ["O2(g)", "H2O", "H+"]:
            vmessage(f"{name} is included as a basis species by default. Moving to the next species...", 2, verbose)
            continue

        # Skip species that have been promoted to strict basis (for redox suppression)
        if name in tag_vec and tag_vec[name] == "basis":
            vmessage(f"'{name}' was promoted to strict basis species. Skipping auxiliary/non-basis processing...", 2, verbose)
            continue

        # Skip polymorphs
        if entry["state"] in [f"cr{i}" for i in range(2, 21)]:
            continue

        date = str(entry["date"])
        date = date[:9]  # Truncate date if greater than 9 characters

        elem = makeup(thermo_df.iloc[idx]["formula_modded"])  # Get elemental composition

        aux_basis = False

        # Check if this is a basis species
        # basis_list can be either a dict (from Python) or an R object
        if hasattr(dissrxns, 'rx2'):  # R object
            try:
                basis_list = dissrxns.rx2("basis_list")
                if hasattr(basis_list, 'names'):  # R StrVector with names
                    basis_list_names = [str(basis_list.rx2(k)[0]) for k in basis_list.names]
                else:
                    basis_list_names = []
            except:
                basis_list_names = []
        elif "basis_list" in dissrxns:  # Python dict
            basis_list = dissrxns["basis_list"]
            if hasattr(basis_list, 'items'):  # It's a dict-like object
                basis_list_names = [v[0] if isinstance(v, (list, tuple)) else v for v in basis_list.values()]
            elif hasattr(basis_list, '__iter__'):  # It's an iterable
                basis_list_names = [v[0] if isinstance(v, (list, tuple)) else v for v in basis_list]
            else:
                basis_list_names = []
        else:
            basis_list_names = []

        if name in basis_list_names:
            vmessage(f"'{name}' (basis species) processed successfully.", 2, verbose)
            continue
        elif entry["tag"] == "basis":
            vmessage(f"'{name}' (basis species) processed successfully.", 2, verbose)
            continue
        elif entry["tag"] == "aux":
            # Auxiliary basis species
            aux_basis = True

        # Format charge for this species' data0 entry
        if "Z" in elem:
            charge = elem["Z"]
            elem = {k: v for k, v in elem.items() if k != "Z"}
            formatted_charge = f"{charge:.1f}"
        else:
            formatted_charge = "0.0"

        # Format the element block of this species' data0 entry
        # Layout: 3 elements per line, each line starts with 7 spaces
        elem_list = []
        elem_names = list(elem.keys())

        for i, elem_name in enumerate(elem_names):
            elem_val = s_d(elem[elem_name], 4)

            # Conditional formatting based on position (3 elements per line)
            if i % 3 == 0:  # First entry of a line (indices 0, 3, 6, 9, ...)
                max_length = 8
                end_char = ""
                if len(elem_name) != 2:
                    elem_name = elem_name + " "
            elif (i + 1) % 3 == 0 and i != len(elem_names) - 1:  # Last entry of a line (not final element)
                max_length = 15
                end_char = "\n"
            else:  # Middle entry of a line
                max_length = 15
                end_char = ""
                if len(elem_name) != 2:
                    elem_name = elem_name + " "

            # Paste together value and element name
            pasted_entry = f"{elem_val} {elem_name}"

            # Get decimal position and format spaces accordingly
            decimal_position = pasted_entry.index(".")
            pasted_entry = " " * (max_length - decimal_position) + pasted_entry + end_char

            # Add entry to element list
            elem_list.append(pasted_entry)

        n_elements = str(len(elem_list))
        element_list = "".join(elem_list)

        # Format the dissociation reaction block of this species' data0 entry
        if name in dissociation_list:
            species_name_list = list(dissociation_list[name].keys())
            species_val_list = list(dissociation_list[name].values())
        else:
            species_name_list = []
            species_val_list = []

        n_species = len(species_name_list)
        nchar = 25  # spaces taken by entries in dissociation reaction (25 char)

        # Convert species names to their data0 counterparts
        for j, species in enumerate(species_name_list):
            if species in CHNOSZ_data0_name_diff:
                species_name_list[j] = CHNOSZ_data0_name_diff[species]

        logK_list = f"logK_grid_{name}"

        # Loop through species in dissociation reaction and format for data0
        spec_list = []

        for i in range(n_species):
            # Get species value and name
            species_val = f"{species_val_list[i]:.4f}"
            species_name = species_name_list[i]

            # Conditional formatting based on position
            # Always put 2 species per line (matching R behavior)
            # i=0, 2, 4, ... are first entries of lines
            # i=1, 3, 5, ... are second entries of lines (add newline after)
            if i % 2 == 0:  # First entry of a line (i=0, 2, 4, ...)
                max_length = 7
                end_char = ""
                species_name = fillspace(species_name, nchar)
            elif i % 2 == 1 and i != n_species - 1:  # Second entry of a line, not last (i=1, 3, 5, ...)
                max_length = 8
                end_char = "\n"
            else:  # Last entry (odd index, last species)
                max_length = 8
                end_char = ""
                species_name = fillspace(species_name, nchar)

            # Paste together coeff and element name
            pasted_entry = f"{species_val}  {species_name}"

            # Get decimal position and format spaces accordingly
            decimal_position = pasted_entry.index(".")
            pasted_entry = " " * (max_length - decimal_position) + pasted_entry + end_char

            # Add entry to element list
            spec_list.append(pasted_entry)

        # Instantiate template and begin formatting aq, cr, gas, liq entries
        species_list = "".join(spec_list)
        n_species = str(n_species)

        template_fmt = ("+--------------------------------------------------------------------\n%s\n"
                       "    date last revised = %s\n%s\n     charge  =   %s\n%s\n"
                       "     %s element(s):\n%s\n     %s species in aqueous dissociation reaction:      \n%s\n"
                       "**** logK grid [T, P @ Miscellaneous parameters]=%s\n%s")

        if entry["state"] == "aq":
            formatted_name = name
            volume = "*"
            if aux_basis:
                keys = " keys   = aux              active"
                insertline_regex = r"\+-+\naqueous species"
                insertline = "+--------------------------------------------------------------------\naqueous species"
            else:
                keys = " keys   = aqueous          active"
                insertline_regex = r"\+-+\nsolids"
                insertline = "+--------------------------------------------------------------------\nsolids"
        elif entry["state"] == "cr":
            formatted_name = fillspace(name, 24) + entry["formula"]
            tag = tag_vec[entry["name"]]
            formatted_tag = fillspace(tag, 17)
            keys = f" keys   = solid            {formatted_tag}active"
            V_entry = entry["V"] if "V" in entry.index else float("nan")
            if pd.isna(V_entry):
                V_entry = 0.0  # unknown molar volume; not used by EQ3NR
            formatted_V0PrTr = fillspace(str(V_entry), 9, spaces_after=False)
            volume = f"     V0PrTr = {formatted_V0PrTr} cm**3/mol"
            insertline_regex = r"\+-+\nliquids"
            insertline = "+--------------------------------------------------------------------\nliquids"
        elif entry["state"] == "gas":
            formatted_name = name
            tag = tag_vec[entry["name"]]
            formatted_tag = fillspace(tag, 17)
            keys = f" keys   = gas              {formatted_tag}active"
            volume = "     V0PrTr = 24465.000 cm**3/mol  (source = ideal gas law)"
            insertline_regex = r"\+-+\nsolid solutions"
            insertline = "+--------------------------------------------------------------------\nsolid solutions"
        elif entry["state"] == "liq":
            formatted_name = fillspace(name, 24) + entry["formula"]
            tag = tag_vec[entry["name"]]
            formatted_tag = fillspace(tag, 17)
            keys = f" keys   = liquid           {formatted_tag}active"
            formatted_V0PrTr = fillspace(str(entry["V"]), 9, spaces_after=False)
            volume = f"     V0PrTr = {formatted_V0PrTr} cm**3/mol"
            insertline_regex = r"\+-+\ngases"
            insertline = "+--------------------------------------------------------------------\ngases"
        else:
            raise ValueError(f"Error: in {entry['name']} ... {entry['state']} "
                           "is not a recognized state. Must be aq, gas, cr, liq, or ss.")

        # Append to aq, solid, gas, or liq portion of data0.min template
        output = template_fmt % (formatted_name, date, keys, formatted_charge, volume,
                                n_elements, element_list, n_species, species_list,
                                name, logK_list)

        if name == "O2":
            # Temporary placeholder values that are overwritten with calculated values later
            # Need DOTALL flag so . matches newlines
            O2_entry = r"\+-+\nO2\n.*?=O2\n         2.6560    3.0310    3.1080    3.0350\n         2.8740    2.6490    2.3540    1.8830"
            data0_template = re.sub(O2_entry, output, data0_template, flags=re.DOTALL)
        elif name == "H2O(g)":
            # Temporary placeholder values that are overwritten with calculated values later
            # Need DOTALL flag so . matches newlines
            steam_entry = r"\+-+\nH2O\(g\)\n.*?=H2O\(g\)\n         2.2990    0.9950   -0.0060   -0.6630\n        -1.1560   -1.5340   -1.8290   -2.0630"
            data0_template = re.sub(steam_entry, output, data0_template, flags=re.DOTALL)
        else:
            data0_template = re.sub(insertline_regex, output + "\n" + insertline, data0_template)

        vmessage(f"'{name}' processed successfully.", 2, verbose)

    # Handle basis species
    basis_entry_template = ("+--------------------------------------------------------------------\n%s\n"
                           "    date last revised =  %s\n keys   = basis            active\n"
                           "     charge  =   %s\n     %s element(s):\n%s")

    # Extract basis list
    if hasattr(dissrxns, 'rx2'):  # R object
        try:
            basis_list = dissrxns.rx2("basis_list")
            if hasattr(basis_list, 'names'):  # R StrVector with names
                basis_items = [str(basis_list.rx2(k)[0]) for k in basis_list.names]
            else:
                basis_items = []
        except:
            basis_items = []
    elif "basis_list" in dissrxns:  # Python dict
        basis_list = dissrxns["basis_list"]
        if hasattr(basis_list, 'items'):  # Dict-like
            basis_items = list(basis_list.values())
        elif hasattr(basis_list, '__iter__'):  # Iterable
            basis_items = list(basis_list)
        else:
            basis_items = []
    else:
        basis_items = []

    # Add promoted auxiliary basis species (for redox suppression)
    if redox_elem_states:
        for species_name, tag in tag_vec.items():
            if tag == "basis" and species_name not in basis_items:
                # Check if this was promoted from aux
                species_row = thermo_df[thermo_df['name'] == species_name]
                if len(species_row) > 0 and species_row.iloc[0]['tag'] == 'aux':
                    basis_items.append(species_name)
                    vmessage(f"Adding promoted species '{species_name}' to basis section", 1, verbose)

    for basis_item in basis_items:
        # Handle both single values and tuples/lists
        if isinstance(basis_item, (list, tuple)):
            basis = basis_item[0]
        else:
            basis = basis_item

        # Skip hard-coded EQ3 basis species
        if basis in fixed_species:
            continue

        # Get the date on the basis species entry
        try:
            date_basis = info(info(basis, messages=False), messages=False)["date"]
            # Handle pandas Series or scalar
            if hasattr(date_basis, 'iloc'):
                date_basis = str(date_basis.iloc[0])
            else:
                date_basis = str(date_basis)
        except:
            date_basis = "01-jan-1970"

        # Get the basis species formula
        basis_formula_matches = add_obigt_df[add_obigt_df["name"] == basis]["formula_modded"]
        if len(basis_formula_matches) > 0:
            basis_formula = basis_formula_matches.iloc[0]
        else:
            continue

        # Get the elemental makeup of the basis species
        elem_list_basis = makeup(basis_formula)

        # Extract charge from elemental composition
        if "Z" in elem_list_basis:
            charge_basis = elem_list_basis["Z"]
            elem_list_basis = {k: v for k, v in elem_list_basis.items() if k != "Z"}
            formatted_charge_basis = f"{charge_basis:.1f}"
        else:
            formatted_charge_basis = "0.0"

        # Begin formatting for data0 entry
        n_elem_basis = len(elem_list_basis)
        elem_list_basis_names = list(elem_list_basis.keys())
        elem_list_basis_coefs = [f"{elem_list_basis[k]:.4f}" for k in elem_list_basis_names]
        elem_list_basis_formatted = ""

        for i in range(len(elem_list_basis_names)):
            if i == 0:
                sp = "      "
            else:
                sp = "              "
            elem_list_basis_formatted += f"{sp}{elem_list_basis_coefs[i]} {elem_list_basis_names[i]}"

        basis_entry = basis_entry_template % (basis, date_basis, formatted_charge_basis,
                                              n_elem_basis, elem_list_basis_formatted)

        # Add basis entry to data0.min template
        # Use count=1 to only replace the first occurrence (like R's sub() function)
        basis_insertline_regex = r"\+--------------------------------------------------------------------\nO2\(g\)"
        basis_insertline = "+--------------------------------------------------------------------\nO2(g)"

        data0_template = re.sub(basis_insertline_regex, basis_entry + "\n" + basis_insertline,
                               data0_template, count=1)

    # Handle elements
    # CRITICAL: Elements must be ordered to match basis species order for EQPT
    # The n-th basis species must contain the n-th element

    # Build ordered element list from basis species
    ordered_elements = []

    # Hardcoded basis species in template: H2O, H+, plhold, O2(g)
    # H2O contains O and H
    ordered_elements.append('O')
    ordered_elements.append('H')
    # plhold contains pl
    ordered_elements.append('pl')

    # Now add elements from basis_items (the basis species we're adding)
    for basis_item in basis_items:
        # Handle both single values and tuples/lists
        if isinstance(basis_item, (list, tuple)):
            basis = basis_item[0]
        else:
            basis = basis_item

        # Skip if already in fixed_species (will be handled by template)
        if basis in fixed_species:
            continue

        # Get the formula for this basis species
        basis_formula_matches = add_obigt_df[add_obigt_df["name"] == basis]["formula_modded"]
        if len(basis_formula_matches) == 0:
            vmessage(f"Warning: Could not find formula for basis species {basis}", 1, verbose)
            continue

        basis_formula = basis_formula_matches.iloc[0]

        try:
            # Get elemental composition
            elem_makeup = makeup(basis_formula)

            # Find the non-H-O element (this is the element this basis species represents)
            elem_for_basis = None
            for elem in elem_makeup.keys():
                if elem not in ['H', 'O', 'Z']:
                    elem_for_basis = elem
                    break

            if elem_for_basis and elem_for_basis not in ordered_elements:
                ordered_elements.append(elem_for_basis)
        except Exception as e:
            vmessage(f"Warning: Could not parse formula for {basis}: {e}", 1, verbose)

    # O2(g) is last in template - it uses O which is already in the list
    # So we don't need to add anything for O2(g)

    # Now format the element block
    elem_temp_lines = []
    for elem in ordered_elements:
        try:
            weight = f"{mass(elem):.5f}".strip()
        except:
            # Special handling for placeholder element and pseudoelements
            if elem == 'pl':
                weight = "1.00000"
            elif elem.endswith('p') or elem.endswith('n') or elem.endswith('z'):
                # Pseudoelement - extract the base element name
                # e.g., Fejiip -> Fe, Fejiiip -> Fe
                base_elem_match = re.match(r'([A-Z][a-z]?)j', elem)
                if base_elem_match:
                    base_elem = base_elem_match.group(1)
                    try:
                        weight = f"{mass(base_elem):.5f}".strip()
                    except:
                        weight = "1.00000"
                else:
                    weight = "1.00000"
            else:
                weight = "1.00000"

        # Format: element name, spaces, atomic weight (total 18 chars)
        elem_temp_lines.append(elem + " " * (18 - (len(elem) + len(weight))) + weight)

    # Join the elements back together
    elem_block = "\n".join(elem_temp_lines)

    # Insert the full element block into the data0 template
    data0_template = re.sub(r"elements\n\+-+\n.*\n\+-+\nbasis species",
                           "elements\n+--------------------------------------------------------------------\n" + elem_block +
                           "\n+--------------------------------------------------------------------\nbasis species",
                           data0_template, flags=re.DOTALL)

    vmessage("Handling solid solutions...", 2, verbose)

    # Handle solid solutions
    if solid_solution_df is not None and len(solid_solution_df) > 0:
        for i in range(len(solid_solution_df)):
            entry = solid_solution_df.iloc[i]
            name = entry["name"]

            vmessage(f"Processing {name}", 2, verbose)

            date = entry["date"]
            species_name_list = entry["components"].split()
            nextflag = False

            # Check if all components are in the database
            for species in species_name_list:
                if species not in add_obigt_df["name"].values:
                    vmessage(f"Error encountered when processing the solid solution '{name}': "
                           f"'{species}' was not found in the data file as a pure mineral. Skipping it...",
                           1, verbose)
                    nextflag = True
                    break
                elif species in skipped_species:
                    vmessage(f"Error encountered when processing the solid solution '{name}': "
                           f"the dissociation reaction for '{species}' contained one or more NAs "
                           "in its logK grid. Skipping it...", 1, verbose)
                    nextflag = True
                    break

            if nextflag:
                continue

            n_species = len(species_name_list)
            species_val_list = [1] * n_species
            nchar = 23  # spaces taken by entries in ss endmember list

            # Loop through species and format for data0
            spec_list = []
            for j in range(n_species):
                # Get species value and name
                species_val = f"{species_val_list[j]:.4f}"
                species_name = species_name_list[j]

                # Conditional formatting based on position (0-indexed, so even=line start, odd=line end)
                if j % 2 == 0:  # Even index: first entry of a line (left column)
                    max_length = 7
                    end_char = ""
                    species_name = fillspace(species_name, nchar)
                elif j % 2 == 1 and j != n_species - 1:  # Odd index: second entry of a line (right column), not last
                    max_length = 8
                    end_char = "\n"
                else:  # Odd index but last entry
                    max_length = 8
                    end_char = ""
                    species_name = fillspace(species_name, nchar)

                # Paste together value and element name
                pasted_entry = f"{species_val}  {species_name}"

                # Get decimal position and format spaces accordingly
                # Note: index() returns 0-based position, but R's gregexpr returns 1-based
                # To match R behavior, add 1 to convert to 1-based indexing
                decimal_position = pasted_entry.index(".") + 1
                pasted_entry = " " * (max_length - decimal_position) + pasted_entry + end_char

                # Add entry to element list
                spec_list.append(pasted_entry)

            species_list = "".join(spec_list)
            n_species = str(n_species)

            ss_template = ("+--------------------------------------------------------------------\n%s\n"
                          "    date last revised = %s\n%s\n  %s components\n%s\n type   = %s\n%s\n%s")

            model_param_list = [entry[f"mp{i}"] for i in range(1, 7)]
            site_param_list = [entry[f"sp{i}"] for i in range(1, 7)]

            model_param_list = " ".join([s_d(p, 3) for p in model_param_list])
            site_param_list = " ".join([s_d(p, 3) for p in site_param_list])

            if float(entry["n_model_params"]) > 0:
                model_param_template = f"  {entry['n_model_params']} model parameters\n{model_param_list}"
            else:
                model_param_template = "  0 model parameters"

            if float(entry["n_site_params"]) > 0:
                site_param_template = f"  {entry['n_site_params']} site parameters\n{site_param_list}"
            else:
                site_param_template = "  0 site parameters"

            # Format solid solution entry
            formatted_name = fillspace(name, 24) + entry["formula"]
            tag = entry["tag"]
            formatted_tag = fillspace(tag, 17)
            keys = f" keys   = ss               {formatted_tag}active"

            insertline_regex = r"\+-+\nreferences"
            insertline = "+--------------------------------------------------------------------\nreferences"

            # Create output for solid solution entries
            output = ss_template % (formatted_name, date, keys, n_species, species_list,
                                   entry["type"], model_param_template, site_param_template)

            data0_template = re.sub(insertline_regex, output + "\n" + insertline, data0_template)
    else:
        vmessage("No solid solutions supplied. Moving on...", 2, verbose)

    def add_bdot_line(spec_name, spec_azero, neutral_ion_type):
        """Add a bdot parameter line to the data0 template"""
        nonlocal data0_template

        # Handle NaN values in neutral_ion_type
        if pd.isna(neutral_ion_type):
            neutral_ion_type_str = ""
        else:
            neutral_ion_type_str = str(int(neutral_ion_type))

        # Format: species_name + variable spaces + azero (total width 38) + 4 spaces + neutral_ion_type (right-aligned in 2 chars)
        # This accommodates azero values from 1-99 (e.g., "9.0000" is 6 chars, "11.0000" is 7 chars)
        # Calculate spacing: len(name) + spaces + len(azero) = 38
        spaces_needed = 38 - len(spec_name) - len(spec_azero)
        name_and_azero = spec_name + " " * spaces_needed + spec_azero
        neutral_ion_type_part = neutral_ion_type_str.rjust(2)

        bdot_entry = name_and_azero + "    " + neutral_ion_type_part

        # Add bdot entry to data0.min template
        bdot_insertline_regex = r"\+--------------------------------------------------------------------\nelements"
        bdot_insertline = "+--------------------------------------------------------------------\nelements"
        data0_template = re.sub(bdot_insertline_regex, bdot_entry + "\n" + bdot_insertline,
                               data0_template)

    if activity_model == "pitzer":
        # Replace the bdot parameter block with Pitzer interaction
        # parameter blocks. Pitzer data0 files do not have a bdot block.
        if pitzer_df is None:
            raise ValueError("A Pitzer parameter database is required to "
                             "create a data0 file for the Pitzer activity model.")

        # Charges of all aqueous species that appear in the data0 file
        species_charges = {"H2O": 0.0, "H+": 1.0}
        for i in range(len(add_obigt_df)):
            row = add_obigt_df.iloc[i]
            if row["state"] != "aq":
                continue
            try:
                z = makeup(row["formula_modded"]).get("Z", 0.0)
            except Exception:
                z = row["z.T"] if not pd.isna(row["z.T"]) else 0.0
            species_charges[row["name"]] = float(z)

        pitzer_out = build_pitzer_superblocks(pitzer_df, species_charges, verbose=verbose)

        vmessage(f"Added {pitzer_out['n_blocks']} Pitzer interaction parameter "
                 "blocks to the data0 file.", 2, verbose)

        pitzer_insertline_regex = r"\+-+\nbdot parameters\n\+-+\n.*?\n\+-+\nelements"
        data0_template = re.sub(pitzer_insertline_regex,
                                lambda m: DATA0_DELIMITER + "\n" + pitzer_out["text"] + "\nelements",
                                data0_template, count=1, flags=re.DOTALL)

        # Tell EQPT how the Pitzer data blocks are organized. These lines
        # are read from the title of the data0 file.
        config_anchor = "*   NO. OF POINTS IN RANGE 2  = 4\n"
        config_text = config_anchor + "\n".join(pitzer_config_lines(pitzer_out["n_terms"])) + "\n"
        if config_anchor in data0_template:
            data0_template = data0_template.replace(config_anchor, config_text, 1)
        else:
            raise ValueError("Could not find the data grid parameter lines in the "
                             "data0 template needed to insert Pitzer configuration lines.")

        return data0_template

    # Format basis and non-basis species for bdot parameter section
    if water_model in ["SUPCRT92", "IAPWS95"]:
        # azero and neutral ion type parameter defaults for fixed species
        Cl_azero = 3
        O2_azero = 3
        O2g_azero = 3
        OH_azero = 3
        H2O_azero = 3
        H_azero = 9
        Cl_neutral_ion_type = 0
        O2_neutral_ion_type = -1
        O2g_neutral_ion_type = 0
        H2O_neutral_ion_type = 0
        H_neutral_ion_type = 0
    elif water_model == "DEW":
        Cl_azero = 3.7
        O2_azero = -0.5
        O2g_azero = 0
        OH_azero = 3.7
        H2O_azero = 0
        H_azero = 3.7
        Cl_neutral_ion_type = 0
        O2_neutral_ion_type = 0
        O2g_neutral_ion_type = 0
        H2O_neutral_ion_type = 0
        H_neutral_ion_type = 0

    # Add fixed species to bdot section
    for sp in ["Cl-", "O2", "O2(g)", "H2O", "H+"]:
        if sp == "Cl-":
            spec_name = "Cl-"
            spec_azero = f"{Cl_azero:.4f}"
            neutral_ion_type = Cl_neutral_ion_type
        elif sp == "O2":
            spec_name = "O2"
            spec_azero = f"{O2_azero:.4f}"
            neutral_ion_type = O2_neutral_ion_type
        elif sp == "O2(g)":
            spec_name = "O2(g)"
            spec_azero = f"{O2g_azero:.4f}"
            neutral_ion_type = O2g_neutral_ion_type
        elif sp == "H2O":
            spec_name = "H2O"
            spec_azero = f"{H2O_azero:.4f}"
            neutral_ion_type = H2O_neutral_ion_type
        elif sp == "H+":
            spec_name = "H+"
            spec_azero = f"{H_azero:.4f}"
            neutral_ion_type = H_neutral_ion_type

        add_bdot_line(spec_name, spec_azero, neutral_ion_type)

    # Add all aqueous species to bdot section
    for i in range(len(add_obigt_df)):
        row = add_obigt_df.iloc[i]
        if row["name"] in ["Cl-", "O2", "H2O", "H+"]:
            continue
        elif row["state"] == "aq":
            spec_name = row["name"]
            spec_azero_val = azero_vec[row["name"]]

            if pd.isna(spec_azero_val) or np.isnan(spec_azero_val):
                if row["z.T"] != 0:
                    if water_model in ["SUPCRT92", "IAPWS95"]:
                        spec_azero = "4.0000"
                    elif water_model == "DEW":
                        spec_azero = "3.7000"
                else:
                    if water_model in ["SUPCRT92", "IAPWS95"]:
                        spec_azero = "0.0000"
                    elif water_model == "DEW":
                        spec_azero = "-0.5000"
            else:
                spec_azero = f"{spec_azero_val:.4f}"

            neutral_ion_type = neutral_ion_type_vec[row["name"]]
            add_bdot_line(spec_name, spec_azero, neutral_ion_type)

    return data0_template
