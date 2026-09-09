"""
Main module for processing dissociation reactions in thermodynamic databases.

This module handles:
1. Checking that dissociation reactions are balanced
2. Checking that dissociation reactions are written in terms of strict basis
   and auxiliary basis species
3. Generating balanced dissociation reactions for species that lack them or have unbalanced or invalid ones
4. Organizing basis species to prevent EQPT errors

This replaces the R script redox_and_dissrxns.r for the dissociation reaction
checking and generation tasks (redox suppression will be handled separately).
"""

import pandas as pd
import numpy as np
import pychnosz
from .dissrxn_balancer import check_dissrxn_balanced, format_dissrxn_string
from .dissrxn_generator import (
    find_species_needing_dissrxns,
    find_invalid_dissrxn_species,
    generate_dissrxn,
    get_simplest_basis,
    match_basis_comp
)
from .redox_suppression import (
    create_redox_element_states_dict,
    add_pseudoelements_to_chnosz,
    modify_formulas_for_redox_suppression
)


def order_thermo_df(thermo_df, fixed_species=None, verbose=1):
    """
    Order thermo_df to prevent EQPT errors.

    Strict basis species come first, then aux basis species that depend only
    on strict basis, then aux that depend on other aux, then everything else.

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database dataframe
    fixed_species : list, optional
        List of fixed species like H2O, H+, O2(g)
    verbose : int
        Verbosity level

    Returns
    -------
    pd.DataFrame
        Reordered thermodynamic database
    """
    if fixed_species is None:
        fixed_species = ["H2O", "H+", "O2(g)", "water", "e-"]

    # Separate by tag
    basis_entries = thermo_df[thermo_df['tag'] == 'basis'].copy()
    aux_entries = thermo_df[thermo_df['tag'] == 'aux'].copy()
    other_entries = thermo_df[(thermo_df['tag'] != 'basis') &
                               (thermo_df['tag'] != 'aux')].copy()

    # Helper function to check if all species in dissrxn are in a given list
    def all_in_basis(dissrxn, basis_names):
        if not dissrxn or dissrxn.strip() == '' or dissrxn == 'nan':
            return False
        try:
            parts = dissrxn.strip().split()
            species_in_rxn = [parts[i] for i in range(1, len(parts), 2)]
            # Remove the dissociating species itself (first one)
            species_in_rxn = species_in_rxn[1:]
            # Remove fixed species
            species_in_rxn = [s for s in species_in_rxn if s not in fixed_species]
            return all(s in basis_names for s in species_in_rxn)
        except:
            return False

    # Get basis names
    basis_names = list(basis_entries['name'])

    # Categorize aux entries
    aux_with_only_strict = aux_entries[
        aux_entries['dissrxn'].apply(lambda x: all_in_basis(x, basis_names))
    ].copy()

    aux_with_other = aux_entries[
        ~aux_entries['dissrxn'].apply(lambda x: all_in_basis(x, basis_names))
    ].copy()

    # Further categorize aux_with_other
    combined_basis = basis_names + list(aux_with_only_strict['name'])
    aux_with_strict_aux = aux_with_other[
        aux_with_other['dissrxn'].apply(lambda x: all_in_basis(x, combined_basis))
    ].copy()

    aux_with_nonstrict_aux = aux_with_other[
        ~aux_with_other['dissrxn'].apply(lambda x: all_in_basis(x, combined_basis))
    ].copy()

    # Sort aux_with_nonstrict_aux to prevent circular dependencies
    if len(aux_with_nonstrict_aux) > 0:
        # Simple topological sort attempt
        max_iterations = 10
        for iteration in range(max_iterations):
            changed = False
            ordered_list = list(aux_with_nonstrict_aux['name'])

            for i, row in aux_with_nonstrict_aux.iterrows():
                dissrxn = row['dissrxn']
                if not dissrxn or dissrxn.strip() == '':
                    continue

                # Get species in this dissrxn
                try:
                    parts = dissrxn.strip().split()
                    species_in_rxn = [parts[j] for j in range(1, len(parts), 2)][1:]
                    species_in_rxn = [s for s in species_in_rxn
                                       if s not in fixed_species and s != row['name']]

                    # Check if any species in dissrxn come after this species
                    current_pos = ordered_list.index(row['name'])
                    for sp in species_in_rxn:
                        if sp in ordered_list:
                            sp_pos = ordered_list.index(sp)
                            if sp_pos > current_pos:
                                # Move current species after sp
                                ordered_list.remove(row['name'])
                                ordered_list.insert(sp_pos, row['name'])
                                changed = True
                                break
                except:
                    continue

            if not changed:
                break

            # Reorder dataframe according to ordered_list
            aux_with_nonstrict_aux = aux_with_nonstrict_aux.set_index('name').loc[ordered_list].reset_index()

    # Concatenate all parts in correct order
    result = pd.concat([
        basis_entries,
        aux_with_only_strict,
        aux_with_strict_aux,
        aux_with_nonstrict_aux,
        other_entries
    ], ignore_index=True)

    return result


def process_dissrxns(thermo_df, water_model='SUPCRT92', exceed_Ttr=True,
                      suppress_redox=None, exclude_category=None,
                      element_df=None, fixed_species=None, verbose=1,
                      already_rejected_species=None):
    """
    Process dissociation reactions for a thermodynamic database.

    This function:
    1. Handles redox suppression by creating pseudoelements
    2. Checks which species have missing or unbalanced dissociation reactions
    3. Checks which species have dissociation reactions that are written in
       terms of species that are not strict or auxiliary basis species
    4. Generates balanced dissociation reactions for those species
    5. Orders the database to prevent EQPT errors
    6. Performs cascading rejection of species dependent on rejected species

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database dataframe
    water_model : str
        Water model to use (currently not used in Python version)
    exceed_Ttr : bool
        Whether to exceed transition temperatures (currently not used)
    suppress_redox : list or None
        List of elements to suppress redox for (e.g., ["Fe", "S"])
    exclude_category : dict or None
        Categories to exclude (currently not used)
    element_df : pd.DataFrame or None
        Element database (for future use)
    fixed_species : list or None
        List of fixed species (H2O, H+, O2(g), etc.)
    verbose : int
        Verbosity level (0=silent, 1=normal, 2=detailed)

    Returns
    -------
    dict
        Dictionary with keys:
        - 'thermo_df': Updated thermodynamic database
        - 'dissrxns': Dictionary of generated dissociation reactions
        - 'basis_pref': Dictionary of preferred basis species
        - 'redox_elem_states': Dictionary of pseudoelement mappings (if redox suppressed)
    """
    if fixed_species is None:
        fixed_species = ["H2O", "H+", "O2(g)", "water", "e-"]

    if suppress_redox is None:
        suppress_redox = []

    if already_rejected_species is None:
        already_rejected_species = []

    # Make a copy to avoid modifying original
    thermo_df = thermo_df.copy()

    # Initialize formula_modded and formula_ox_modded columns
    if 'formula_modded' not in thermo_df.columns:
        thermo_df['formula_modded'] = thermo_df['formula']

    if 'formula_ox_modded' not in thermo_df.columns:
        thermo_df['formula_ox_modded'] = thermo_df.get('formula_ox', 'nan')

    # Handle redox suppression if requested
    redox_elem_states = {}
    redox_rejected_species = []
    if len(suppress_redox) > 0:
        if verbose >= 1:
            print(f"Suppressing redox for elements: {', '.join(suppress_redox)}")

        # Create pseudoelement mappings and get rejected species
        redox_elem_states, redox_rejected_species = create_redox_element_states_dict(thermo_df, suppress_redox)

        if verbose >= 2:
            print("Pseudoelement mappings:")
            for elem, states in redox_elem_states.items():
                print(f"  {elem}:")
                for ox_state, pseudoelem in states.items():
                    print(f"    {ox_state} -> {pseudoelem}")

        if redox_rejected_species and verbose >= 1:
            print(f"Rejecting {len(redox_rejected_species)} species without aqueous basis for their oxidation states")
            if verbose >= 2:
                print(f"  Rejected: {', '.join(redox_rejected_species[:10])}")
                if len(redox_rejected_species) > 10:
                    print(f"  ... and {len(redox_rejected_species) - 10} more")

        # Remove rejected species from thermo_df
        thermo_df = thermo_df[~thermo_df['name'].isin(redox_rejected_species)].copy()

        # Add pseudoelements to CHNOSZ element database
        add_pseudoelements_to_chnosz(redox_elem_states, element_df)

        # Modify formulas to use pseudoelements
        thermo_df = modify_formulas_for_redox_suppression(thermo_df, redox_elem_states, verbose)

    # Add regenerate_dissrxn column if it doesn't exist
    # This column is used by _HKF_cgl.py
    if 'regenerate_dissrxn' not in thermo_df.columns:
        thermo_df['regenerate_dissrxn'] = False

    # Clean up whitespace in dissrxn column
    if 'dissrxn' in thermo_df.columns:
        thermo_df['dissrxn'] = thermo_df['dissrxn'].apply(
            lambda x: x.strip() if isinstance(x, str) else x
        )

    # Promote auxiliary basis species to strict basis when they represent pseudoelements
    # This must happen BEFORE building basis_pref so the pseudoelements are included
    # This replicates lines 1151-1165 of the R script redox_and_dissrxns.r
    if len(suppress_redox) > 0 and redox_elem_states:
        # First, get all elements represented in formulas (pseudoelements)
        all_elements = set()
        for idx, row in thermo_df.iterrows():
            try:
                formula_to_check = row.get('formula_modded', row['formula'])
                if pd.notna(formula_to_check) and formula_to_check != '' and formula_to_check != 'nan':
                    mkp = pychnosz.makeup(str(formula_to_check))
                    for elem in mkp.keys():
                        if elem not in ['H', 'O', 'Z']:
                            all_elements.add(elem)
            except:
                continue

        # Check which pseudoelements lack basis species
        basis_df = thermo_df[thermo_df['tag'] == 'basis'].copy()
        basis_elements = set()
        for idx, row in basis_df.iterrows():
            try:
                formula_to_check = row.get('formula_modded', row['formula'])
                if pd.notna(formula_to_check) and formula_to_check != '' and formula_to_check != 'nan':
                    mkp = pychnosz.makeup(str(formula_to_check))
                    for elem in mkp.keys():
                        if elem not in ['H', 'O', 'Z']:
                            basis_elements.add(elem)
            except:
                continue

        # Find elements that need basis representation
        elem_need_basis = all_elements - basis_elements

        # Create reverse mapping from pseudoelement names to oxidation state strings
        # This is needed because formula_modded uses pseudoelement names (e.g., "Sjz")
        # but redox_elem_states uses oxidation state strings as keys (e.g., "S")
        pseudoelem_to_oxstate = {}
        for element, ox_states_dict in redox_elem_states.items():
            for ox_state_str, pseudoelem_name in ox_states_dict.items():
                pseudoelem_to_oxstate[pseudoelem_name] = (element, ox_state_str)

        # Try to promote auxiliary basis species or other suitable species to represent these elements
        # This replicates the R code behavior where get_dissrxn searches all aqueous/gas species
        # (lines 400, 421-425 of redox_and_dissrxns.r)
        # Track elements that couldn't be promoted for later rejection
        elem_without_basis = []
        if elem_need_basis:
            for elem in elem_need_basis:
                # First, look for aux species
                aux_df = thermo_df[thermo_df['tag'] == 'aux'].copy()
                promoted = False

                for idx, row in aux_df.iterrows():
                    try:
                        formula_to_check = row.get('formula_modded', row['formula'])
                        if pd.notna(formula_to_check) and formula_to_check != '' and formula_to_check != 'nan':
                            mkp = pychnosz.makeup(str(formula_to_check))
                            aux_elements = [e for e in mkp.keys() if e not in ['H', 'O', 'Z']]

                            # Check if this aux species can represent the element (pseudoelement)
                            if len(aux_elements) == 1 and aux_elements[0] == elem:
                                # Promote this aux species to basis
                                thermo_df.at[idx, 'tag'] = 'basis'
                                if verbose >= 1:
                                    # Get oxidation state string for clearer output
                                    if elem in pseudoelem_to_oxstate:
                                        orig_elem, ox_state_str = pseudoelem_to_oxstate[elem]
                                        print(f"Promoted '{row['name']}' from aux to basis to represent {ox_state_str}")
                                    else:
                                        print(f"Promoted '{row['name']}' from aux to basis to represent element {elem}")
                                promoted = True
                                break
                    except:
                        continue

                # If no aux species found, look for any aqueous or gaseous species
                # EQ3/6 basis species must be aqueous or gaseous (e.g., O2(g), S2(g))
                if not promoted:
                    # Filter for aqueous and gaseous species (excluding solids and liquids)
                    candidate_df = thermo_df[thermo_df['state'].isin(['aq', 'gas'])].copy()

                    for idx, row in candidate_df.iterrows():
                        try:
                            formula_to_check = row.get('formula_modded', row['formula'])
                            if pd.notna(formula_to_check) and formula_to_check != '' and formula_to_check != 'nan':
                                mkp = pychnosz.makeup(str(formula_to_check))
                                elements = [e for e in mkp.keys() if e not in ['H', 'O', 'Z']]

                                # Check if this species can represent the element (pseudoelement)
                                # Prefer species with only H, O, Z, and the target element
                                if len(elements) == 1 and elements[0] == elem:
                                    # Promote this species to basis
                                    thermo_df.at[idx, 'tag'] = 'basis'
                                    if verbose >= 1:
                                        old_tag = row.get('tag', 'non-basis')
                                        # Handle NaN values
                                        if pd.isna(old_tag) or old_tag == 'nan':
                                            old_tag = 'non-basis'
                                        # Get oxidation state string for clearer output
                                        if elem in pseudoelem_to_oxstate:
                                            orig_elem, ox_state_str = pseudoelem_to_oxstate[elem]
                                            print(f"Promoted '{row['name']}' from {old_tag} to basis to represent {ox_state_str}")
                                        else:
                                            print(f"Promoted '{row['name']}' from {old_tag} to basis to represent element {elem}")
                                    promoted = True
                                    break
                        except:
                            continue

                    if not promoted:
                        # Track this element for species rejection
                        elem_without_basis.append(elem)
                        if verbose >= 1:
                            # Convert pseudoelement name to scientific notation for clearer output
                            if elem in pseudoelem_to_oxstate:
                                orig_elem, ox_state_str = pseudoelem_to_oxstate[elem]
                                print(f"Warning: Could not find a suitable basis species for {ox_state_str}")
                            else:
                                print(f"Warning: Could not find a suitable basis species for element {elem}")

        # Reject species containing pseudoelements without basis representation
        # Also reject species with blank formula_ox when they contain suppressed redox elements
        if elem_without_basis or suppress_redox:
            additional_rejected = []
            blank_formula_ox_species = []
            missing_basis_species = {}  # Map oxidation state to list of species

            for idx, row in thermo_df.iterrows():
                species_name = row['name']
                formula = row['formula']
                formula_ox = row.get('formula_ox', 'nan')
                formula_modded = row.get('formula_modded', formula)

                # Check if species contains any of the suppressed redox elements in its formula
                contains_suppressed_elem = False
                if suppress_redox:
                    try:
                        mkp = pychnosz.makeup(str(formula))
                        for elem in mkp.keys():
                            # Check if this is a base element (not pseudoelement) that's being suppressed
                            for supp_elem in suppress_redox:
                                if elem == supp_elem or elem.startswith(supp_elem + 'j'):
                                    contains_suppressed_elem = True
                                    break
                            if contains_suppressed_elem:
                                break
                    except:
                        pass

                # If species contains suppressed element and has blank formula_ox, reject it
                if contains_suppressed_elem and (pd.isna(formula_ox) or formula_ox == 'nan' or formula_ox == ''):
                    if species_name not in redox_rejected_species and species_name not in additional_rejected:
                        additional_rejected.append(species_name)
                        blank_formula_ox_species.append(species_name)
                    continue

                # Check if species contains pseudoelements without basis
                if elem_without_basis:
                    try:
                        mkp = pychnosz.makeup(str(formula_modded))
                        for elem in mkp.keys():
                            if elem in elem_without_basis:
                                if species_name not in redox_rejected_species and species_name not in additional_rejected:
                                    additional_rejected.append(species_name)
                                    # Get oxidation state string for grouping
                                    if elem in pseudoelem_to_oxstate:
                                        orig_elem, ox_state_str = pseudoelem_to_oxstate[elem]
                                        if ox_state_str not in missing_basis_species:
                                            missing_basis_species[ox_state_str] = []
                                        missing_basis_species[ox_state_str].append(species_name)
                                    else:
                                        if elem not in missing_basis_species:
                                            missing_basis_species[elem] = []
                                        missing_basis_species[elem].append(species_name)
                                break
                    except:
                        pass

            # Print consolidated rejection messages
            if verbose >= 1:
                if blank_formula_ox_species:
                    print(f"\nThe following {len(blank_formula_ox_species)} species will be excluded because they contain a redox-suppressed element in their formulae and do not have element oxidation states defined in the \"formula_ox\" column of the thermodynamic database CSV:")
                    # Print species list in a formatted way
                    species_str = ", ".join([f"'{sp}'" for sp in blank_formula_ox_species[:10]])
                    if len(blank_formula_ox_species) > 10:
                        species_str += f", ... and {len(blank_formula_ox_species) - 10} more"
                    print(f"  {species_str}")

                if missing_basis_species:
                    for ox_state_str, species_list in missing_basis_species.items():
                        print(f"\nRejecting {len(species_list)} species containing {ox_state_str} which lacks a basis species:")
                        species_str = ", ".join([f"'{sp}'" for sp in species_list[:10]])
                        if len(species_list) > 10:
                            species_str += f", ... and {len(species_list) - 10} more"
                        print(f"  {species_str}")

            # Remove rejected species from thermo_df
            if additional_rejected:
                redox_rejected_species.extend(additional_rejected)
                thermo_df = thermo_df[~thermo_df['name'].isin(additional_rejected)].copy()
                if verbose >= 1:
                    print(f"\nTotal: Rejected {len(additional_rejected)} additional species")

    # Cascading rejection: reject species that depend on rejected species
    # This handles cases where basis/auxiliary species are rejected and species
    # depending on them in their dissociation reactions must also be rejected
    all_rejected = list(set(redox_rejected_species + already_rejected_species))
    if all_rejected:
        cascaded_rejected = []
        max_iterations = 10  # Prevent infinite loops

        for iteration in range(max_iterations):
            iteration_rejected = []
            rejected_set = set(all_rejected + cascaded_rejected)

            for idx, row in thermo_df.iterrows():
                species_name = row['name']
                dissrxn = row.get('dissrxn', 'nan')

                # Skip if already rejected or no dissociation reaction
                if species_name in rejected_set:
                    continue
                if pd.isna(dissrxn) or dissrxn == 'nan' or dissrxn == '':
                    continue

                # Parse dissociation reaction to find dependencies
                try:
                    parts = dissrxn.strip().split()
                    # Species in dissrxn are at odd indices (1, 3, 5, ...)
                    species_in_rxn = [parts[i] for i in range(1, len(parts), 2)]
                    # Remove the dissociating species itself (first one)
                    species_in_rxn = species_in_rxn[1:]

                    # Check if any required species have been rejected
                    for req_species in species_in_rxn:
                        if req_species in rejected_set:
                            if species_name not in iteration_rejected and species_name not in cascaded_rejected:
                                iteration_rejected.append(species_name)
                                if verbose >= 2:
                                    print(f"  Rejecting '{species_name}': depends on rejected species '{req_species}'")
                            break
                except:
                    continue

            # If no new rejections this iteration, we're done
            if not iteration_rejected:
                break

            cascaded_rejected.extend(iteration_rejected)

        # Apply cascading rejections
        if cascaded_rejected:
            redox_rejected_species.extend(cascaded_rejected)
            thermo_df = thermo_df[~thermo_df['name'].isin(cascaded_rejected)].copy()
            if verbose >= 1:
                print(f"\nCascading rejection: {len(cascaded_rejected)} additional species rejected due to dependencies on rejected species")
                if verbose >= 2:
                    # Show first few examples
                    examples = ", ".join([f"'{sp}'" for sp in cascaded_rejected[:10]])
                    if len(cascaded_rejected) > 10:
                        examples += f", ... and {len(cascaded_rejected) - 10} more"
                    print(f"  Examples: {examples}")

    # Build basis_pref dictionary from basis species in thermo_df
    # Use formula_modded for redox-suppressed elements
    basis_df = thermo_df[thermo_df['tag'] == 'basis'].copy()
    basis_pref = {}

    if verbose >= 2 and len(suppress_redox) > 0:
        print(f"\nBuilding basis_pref from {len(basis_df)} basis species...")

    for idx, row in basis_df.iterrows():
        try:
            # Use formula_modded if redox suppression is active
            formula_to_use = row['formula_modded'] if len(suppress_redox) > 0 else row['formula']

            if verbose >= 3 and len(suppress_redox) > 0 and 'C' in row['name']:
                print(f"  Processing basis species {row['name']}: formula_to_use={formula_to_use}")

            mkp = pychnosz.makeup(formula_to_use)
            # Find the non-H, non-O, non-Z element (including pseudoelements)
            basis_elem = [e for e in mkp.keys() if e not in ['H', 'O', 'Z']]

            if len(basis_elem) == 1:
                basis_pref[basis_elem[0]] = row['name']
                if verbose >= 2 and len(suppress_redox) > 0:
                    print(f"    {basis_elem[0]} -> {row['name']}")
            elif len(basis_elem) > 1:
                # Basis species has multiple non-H-O elements
                if verbose >= 1:
                    print(f"Warning: Basis species '{row['name']}' has multiple "
                          f"non-H-O elements: {basis_elem}. This should be an "
                          f"auxiliary basis species instead.")
        except Exception as e:
            if verbose >= 2:
                print(f"Error processing basis species {row['name']}: {e}")

    # EQ3 has H2O and O2(g) hard-coded
    basis_pref['H'] = 'H2O'
    basis_pref['O'] = 'O2(g)'

    # Build aux_pref dictionary from auxiliary basis species
    # Use formula_modded for redox-suppressed elements
    aux_df = thermo_df[thermo_df['tag'] == 'aux'].copy()
    aux_pref = {}

    for idx, row in aux_df.iterrows():
        try:
            # Use formula_modded if redox suppression is active
            formula_to_use = row['formula_modded'] if len(suppress_redox) > 0 else row['formula']
            mkp = pychnosz.makeup(formula_to_use)
            aux_elem = [e for e in mkp.keys() if e not in ['H', 'O', 'Z']]

            # Only include aux species with one non-H-O element and one atom of it
            if len(aux_elem) == 1 and mkp[aux_elem[0]] == 1:
                elem = aux_elem[0]
                if elem not in aux_pref:
                    aux_pref[elem] = []
                aux_pref[elem].append(row['name'])
        except Exception as e:
            if verbose >= 2:
                print(f"Error processing aux species {row['name']}: {e}")

    # Generate dissociation reactions
    generated_dissrxns = {}

    # NOTE: Do NOT reset or reload CHNOSZ database here!
    # The custom thermodynamic database has already been loaded during AqEquil() initialization
    # Resetting would clear that custom data and break dissociation reaction generation

    # However, if redox suppression is active, we need to update the OBIGT database
    # to use the modified formulas with pseudoelements
    if len(suppress_redox) > 0:
        # Update CHNOSZ OBIGT database with modified formulas
        # This is necessary for proper dissociation reaction generation
        # Get current OBIGT database
        obigt = pychnosz.thermo().OBIGT

        # Update formulas for species in our database
        for idx, row in thermo_df.iterrows():
            sp_name = row['name']
            formula_modded = row['formula_modded']

            # Skip if formula_modded is not a valid string
            if pd.isna(formula_modded) or formula_modded == '' or formula_modded == 'nan':
                continue

            # Ensure it's a string
            formula_modded = str(formula_modded)

            # Find this species in OBIGT
            obigt_idx = obigt[obigt['name'] == sp_name].index
            if len(obigt_idx) > 0:
                # Update the formula (ensure we're setting a string value)
                obigt.loc[obigt_idx, 'formula'] = formula_modded

        # Ensure all object columns are converted to string dtype
        # This prevents ".str accessor" errors when chnosz processes the data
        for col in obigt.columns:
            if obigt[col].dtype == 'object':
                obigt[col] = obigt[col].astype(str)

        # Reload the OBIGT database
        pychnosz.thermo(OBIGT=obigt)

        if verbose >= 2:
            print("Updated CHNOSZ OBIGT database with modified formulas for redox suppression")

    # Find species needing dissociation reactions
    # IMPORTANT: This must be called AFTER CHNOSZ OBIGT update for redox suppression
    # so that balance checks use the pseudoelement formulas
    needs_dissrxn = find_species_needing_dissrxns(thermo_df, verbose=verbose)

    # Remove duplicates while preserving order (polymorphs have same name)
    # Use dict.fromkeys() to preserve order while removing duplicates
    needs_dissrxn_unique = list(dict.fromkeys(needs_dissrxn))

    if len(needs_dissrxn_unique) > 0:
        if verbose >= 1:
            print(f"\nBalanced dissociation reactions are missing for {len(needs_dissrxn_unique)} species:")
            # Format the list of species nicely
            species_str = ", ".join(needs_dissrxn_unique[:20])
            if len(needs_dissrxn_unique) > 20:
                species_str += f", ... and {len(needs_dissrxn_unique) - 20} more"
            print(f"  {species_str}")

    # A balanced dissociation reaction is still invalid if it is written in
    # terms of species that are not strict basis or auxiliary basis species,
    # because EQPT cannot compile it. Mark these reactions so they get
    # regenerated, too.
    invalid_dissrxns = find_invalid_dissrxn_species(thermo_df, fixed_species=fixed_species)

    if len(invalid_dissrxns) > 0:
        if verbose >= 1:
            db_names = set(thermo_df['name'])
            print(f"\nWarning: dissociation reactions for the following "
                  f"{len(invalid_dissrxns)} species are written in terms of "
                  "species that are not strict basis or auxiliary basis "
                  "species:")
            for sp_name in list(invalid_dissrxns.keys())[:10]:
                required = ["'{}'".format(sp) if sp in db_names
                            else "'{}' (not in the database)".format(sp)
                            for sp in invalid_dissrxns[sp_name]]
                print(f"  '{sp_name}' requires {', '.join(required)}")
            if len(invalid_dissrxns) > 10:
                print(f"  ... and {len(invalid_dissrxns) - 10} more")
            print("These dissociation reactions are invalid and will be "
                  "regenerated.")

    # Species with missing, unbalanced, or invalid dissociation reactions all
    # get a freshly generated dissociation reaction
    regenerate_dissrxn_species = needs_dissrxn_unique + [
        sp for sp in invalid_dissrxns if sp not in needs_dissrxn_unique]

    if len(regenerate_dissrxn_species) > 0:
        thermo_df.loc[thermo_df['name'].isin(regenerate_dissrxn_species),
                      'regenerate_dissrxn'] = True
        if verbose >= 1:
            print("\nGenerating dissociation reactions for these species using "
                  "strict and auxiliary basis species containing a maximum of "
                  "one atom of one element besides O and H...\n")

    # Generate dissociation reactions for species that need them
    failed_dissrxn_species = []
    for sp_name in regenerate_dissrxn_species:
        dissrxn = generate_dissrxn(
            sp_name,
            thermo_df,
            basis_pref,
            aux_pref=aux_pref,
            verbose=verbose,
            redox_elem_states=redox_elem_states
        )

        if dissrxn:
            generated_dissrxns[sp_name] = dissrxn
            # Update thermo_df - this updates ALL polymorphs with the same name
            thermo_df.loc[thermo_df['name'] == sp_name, 'dissrxn'] = dissrxn

            if verbose >= 1:
                print(f"{sp_name} : {dissrxn}")
        else:
            # Failed to generate dissociation reaction - mark for rejection
            failed_dissrxn_species.append(sp_name)
            if verbose >= 1:
                print(f"{sp_name} : Failed to generate dissociation reaction")

    # Reject species that failed dissociation reaction generation
    if failed_dissrxn_species:
        if verbose >= 1:
            print(f"\nRejecting {len(failed_dissrxn_species)} species that failed dissociation reaction generation:")
            species_str = ", ".join([f"'{sp}'" for sp in failed_dissrxn_species[:10]])
            if len(failed_dissrxn_species) > 10:
                species_str += f", ... and {len(failed_dissrxn_species) - 10} more"
            print(f"  {species_str}")

        # Remove from thermo_df
        thermo_df = thermo_df[~thermo_df['name'].isin(failed_dissrxn_species)].copy()

        # Add to rejected species list
        redox_rejected_species.extend(failed_dissrxn_species)

    # Order thermo_df to prevent EQPT errors
    thermo_df = order_thermo_df(thermo_df, fixed_species, verbose)

    # Prepare output
    out_list = {
        'thermo_df': thermo_df,
        'dissrxns': {
            **generated_dissrxns,
            'basis_list': basis_pref,
            'rejected_species': redox_rejected_species
        },
        'basis_pref': basis_pref,
        'redox_elem_states': redox_elem_states  # Include pseudoelement mappings
    }

    return out_list
