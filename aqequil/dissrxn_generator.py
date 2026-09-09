"""
Module for generating dissociation reactions for species that lack them.

This module uses CHNOSZ to automatically generate balanced dissociation reactions
using strict and auxiliary basis species.
"""

import pychnosz
import pandas as pd
import numpy as np
from fractions import Fraction
import math
from .formula_parser import (
    parse_formula_ox_string,
    get_element_oxidation_states_from_formula_ox,
    calculate_average_oxidation_state
)


def gcd(a, b):
    """Greatest common divisor."""
    while b:
        a, b = b, a % b
    return abs(a)


def lcm(a, b):
    """Least common multiple."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)


def get_lcm_list(numbers):
    """Get the least common multiple of a list of numbers."""
    if len(numbers) == 0:
        return 1
    result = numbers[0]
    for num in numbers[1:]:
        result = lcm(result, num)
    return result


def coefficients_to_integers(coeffs, tolerance=1e-6):
    """
    Convert fractional coefficients to integers by finding LCM of denominators.

    Parameters
    ----------
    coeffs : list
        List of coefficient values (can be floats)
    tolerance : float
        Tolerance for considering a coefficient as zero

    Returns
    -------
    list
        List of integer coefficients
    """
    # Convert to fractions and find LCM of denominators
    fractions = []
    for c in coeffs:
        if abs(c) < tolerance:
            fractions.append(Fraction(0))
        else:
            frac = Fraction(c).limit_denominator(10000)
            fractions.append(frac)

    # Get LCM of all denominators
    denominators = [f.denominator for f in fractions]
    lcm_val = get_lcm_list(denominators) if denominators else 1

    # Multiply all coefficients by LCM
    integer_coeffs = [float(f.numerator * lcm_val // f.denominator) for f in fractions]

    return integer_coeffs


def reduce_coefficients(coeffs, tolerance=1e-6):
    """
    Reduce integer coefficients to their simplest form by dividing by their GCD.

    Parameters
    ----------
    coeffs : list
        List of integer coefficient values
    tolerance : float
        Tolerance for considering a coefficient as zero

    Returns
    -------
    list
        List of reduced integer coefficients
    """
    # Filter out zero coefficients for GCD calculation
    non_zero_coeffs = [abs(int(round(c))) for c in coeffs if abs(c) > tolerance]

    if not non_zero_coeffs:
        return coeffs

    # Calculate GCD of all non-zero coefficients
    result_gcd = non_zero_coeffs[0]
    for coeff in non_zero_coeffs[1:]:
        result_gcd = gcd(result_gcd, coeff)

    # If GCD is 1, coefficients are already reduced
    if result_gcd <= 1:
        return coeffs

    # Divide all coefficients by GCD
    reduced_coeffs = [c / result_gcd for c in coeffs]

    return reduced_coeffs


def match_basis_comp(sp_elements, elem):
    """
    Check if a species contains only elem, H, O, and charge.

    Parameters
    ----------
    sp_elements : list
        List of elements in the species
    elem : str
        Element to check for

    Returns
    -------
    bool
        True if species matches basis composition requirements
    """
    if elem in sp_elements:
        valid_elements = set([elem, 'H', 'O', 'Z'])
        return all(e in valid_elements for e in sp_elements)
    return False


def get_simplest_basis(basis_candidates, elem):
    """
    Get the simplest basis species for an element (fewest total atoms).

    Parameters
    ----------
    basis_candidates : list
        List of candidate basis species names
    elem : str
        Element to represent

    Returns
    -------
    str
        Name of simplest basis species, or None if none found
    """
    min_complexity = float('inf')
    simplest = None

    for candidate in basis_candidates:
        try:
            mkp = chnosz.makeup(candidate)
            # Calculate complexity as sum of absolute values
            complexity = sum(abs(v) for v in mkp.values())

            # Prefer species with 1 atom of the element
            if elem in mkp and mkp[elem] == 1:
                complexity -= 0.5  # Bonus for having 1 atom

            if complexity < min_complexity:
                min_complexity = complexity
                simplest = candidate

        except:
            continue

    return simplest


# These functions are now imported from formula_parser module:
# - parse_formula_ox_string (replaces get_oxidation_state)
# - get_element_oxidation_states_from_formula_ox (replaces get_element_oxidation_states)
# - calculate_average_oxidation_state (replaces get_average_oxidation_state)

# Keep aliases for backward compatibility
get_element_oxidation_states = get_element_oxidation_states_from_formula_ox
get_average_oxidation_state = calculate_average_oxidation_state


def generate_dissrxn(sp_name, thermo_df, basis_pref, aux_pref=None,
                     HOZ_balancers=None, verbose=1, redox_elem_states=None):
    """
    Generate a dissociation reaction for a species.

    This function replicates the logic from the R function spec_diss() in
    redox_and_dissrxns.r (lines 536-775), which uses oxidation state information
    to choose the most appropriate basis species.

    Parameters
    ----------
    sp_name : str
        Name of species to generate dissociation reaction for
    thermo_df : pd.DataFrame
        Thermodynamic database dataframe with formula, formula_ox, tag columns
    basis_pref : dict
        Dictionary mapping elements (including pseudoelements) to preferred basis species names
    aux_pref : dict, optional
        Dictionary mapping elements to lists of preferred auxiliary basis species
    HOZ_balancers : list, optional
        List of species to use for balancing H, O, and charge
    verbose : int
        Verbosity level
    redox_elem_states : dict, optional
        Dictionary mapping elements to their pseudoelement names for each oxidation state
        Format: {element: {ox_state_string: pseudoelement_name}}

    Returns
    -------
    str
        Dissociation reaction string, or None if generation fails
    """
    if HOZ_balancers is None:
        HOZ_balancers = ["H+", "O2(g)", "H2O"]

    if aux_pref is None:
        aux_pref = {}

    if redox_elem_states is None:
        redox_elem_states = {}

    try:
        # Get formula for this species
        sp_row = thermo_df[thermo_df['name'] == sp_name]
        if sp_row.empty:
            return None

        # Use formula_modded if available (for redox-suppressed elements), else use formula
        sp_formula = sp_row.iloc[0].get('formula_modded', sp_row.iloc[0]['formula'])
        # If formula_modded is NaN or empty, fall back to formula
        if pd.isna(sp_formula) or sp_formula == '' or sp_formula == 'nan':
            sp_formula = sp_row.iloc[0]['formula']

        # Get formula_ox_modded for oxidation state information
        sp_formula_ox_modded = sp_row.iloc[0].get('formula_ox_modded', 'nan')
        if pd.isna(sp_formula_ox_modded) or sp_formula_ox_modded == 'nan' or sp_formula_ox_modded == '':
            sp_formula_ox_modded = sp_row.iloc[0].get('formula_ox', 'nan')

        # Get elemental composition
        comp = pychnosz.makeup(sp_formula)

        # Get elements that are not H, O, Z (including pseudoelements)
        basis_elem = [e for e in comp.keys() if e not in ['H', 'O', 'Z']]

        # Special case: if this is a gas species and there's an aqueous species with the same formula_modded
        # that is a basis species, create a simple phase-change reaction
        # e.g., CH4(g) -> CH4(aq), CO(g) -> CO(aq), S2(g) -> S2(aq)
        if sp_row.iloc[0].get('state') == 'gas':
            # Look for an aqueous basis species with matching formula_modded
            matching_aq_basis = thermo_df[
                (thermo_df['tag'] == 'basis') &
                (thermo_df['state'] == 'aq') &
                (thermo_df['formula_modded'] == sp_formula)
            ]
            if not matching_aq_basis.empty:
                aq_species_name = matching_aq_basis.iloc[0]['name']
                # Simple phase change reaction
                dissrxn_str = f"-1.0000 {sp_name} 1.0000 {aq_species_name}"
                if verbose >= 2:
                    print(f"  Phase change reaction: {dissrxn_str}")
                return dissrxn_str

        # Determine which basis species to use for each element
        # This implements the oxidation-state-aware selection logic from R's spec_diss()
        chosen_basis_species = []
        chosen_basis_coeffs = []

        # If we don't have oxidation state information, use simple method
        if pd.isna(sp_formula_ox_modded) or sp_formula_ox_modded == 'nan' or sp_formula_ox_modded == '':
            # Simple method: just use basis_pref
            for elem in basis_elem:
                if elem in basis_pref:
                    chosen_basis_species.append(basis_pref[elem])
                    chosen_basis_coeffs.append(comp[elem])
        else:
            # Oxidation-state-aware method (R code lines 583-706)
            # For each element, match oxidation states to basis/aux species

            for elem in basis_elem:
                # Get count of this element in the species
                elem_count = comp[elem]

                # Check if this is a pseudoelement (redox-suppressed element)
                # Pseudoelements have names like "Fejiiip", "Sjin", "Cjivn"
                is_pseudoelement = 'j' in elem and any(c in elem for c in ['i', 'v', 'x', 'p', 'n', 'z'])

                if is_pseudoelement:
                    # This is a pseudoelement from redox suppression
                    # Use the preferred basis species directly
                    if elem in basis_pref:
                        chosen_basis_species.append(basis_pref[elem])
                        chosen_basis_coeffs.append(elem_count)
                    continue

                # Element is NOT redox-suppressed
                # Use oxidation-state matching to find best basis/aux species
                # This replicates R code lines 605-706

                # Check if we have aux_pref for this element
                if elem not in aux_pref or not aux_pref[elem]:
                    # No aux preference - just use simplest basis
                    if elem in basis_pref:
                        chosen_basis_species.append(basis_pref[elem])
                        chosen_basis_coeffs.append(elem_count)
                    continue

                # Get oxidation states of this element in the species
                elem_ox_states = get_element_oxidation_states(sp_formula_ox_modded, elem)

                if not elem_ox_states:
                    # Couldn't parse oxidation states, fall back to simple method
                    if elem in basis_pref:
                        chosen_basis_species.append(basis_pref[elem])
                        chosen_basis_coeffs.append(elem_count)
                    continue

                # Get available basis and aux species for this element
                # Build a list of (species_name, average_ox_state, formula_ox)
                basis_aux_candidates = []

                # Get the strict basis species for this element
                if elem in basis_pref:
                    basis_sp_name = basis_pref[elem]
                    basis_row = thermo_df[thermo_df['name'] == basis_sp_name]
                    if not basis_row.empty:
                        basis_formula_ox = basis_row.iloc[0].get('formula_ox_modded',
                                                                   basis_row.iloc[0].get('formula_ox', ''))
                        avg_ox = get_average_oxidation_state(basis_formula_ox, elem)
                        basis_aux_candidates.append((basis_sp_name, avg_ox, basis_formula_ox))

                # Get aux species for this element
                for aux_sp_name in aux_pref[elem]:
                    # Don't allow a species to dissociate into itself
                    if aux_sp_name == sp_name:
                        continue
                    aux_row = thermo_df[thermo_df['name'] == aux_sp_name]
                    if not aux_row.empty:
                        aux_formula_ox = aux_row.iloc[0].get('formula_ox_modded',
                                                               aux_row.iloc[0].get('formula_ox', ''))
                        avg_ox = get_average_oxidation_state(aux_formula_ox, elem)
                        basis_aux_candidates.append((aux_sp_name, avg_ox, aux_formula_ox))

                if not basis_aux_candidates:
                    # No candidates found, fall back
                    if elem in basis_pref:
                        chosen_basis_species.append(basis_pref[elem])
                        chosen_basis_coeffs.append(elem_count)
                    continue

                # For each oxidation state instance in the species, find the closest matching basis/aux
                # This implements R code lines 651-705
                if len(elem_ox_states) == 1:
                    # Simple case: element has one oxidation state
                    ox_state, count = elem_ox_states[0]

                    # Find closest matching basis/aux species
                    closest_sp = min(basis_aux_candidates, key=lambda x: abs(x[1] - ox_state))
                    chosen_sp_name = closest_sp[0]

                    # Get how many atoms of this element are in the chosen basis species
                    chosen_row = thermo_df[thermo_df['name'] == chosen_sp_name]
                    chosen_formula = chosen_row.iloc[0].get('formula_modded', chosen_row.iloc[0]['formula'])
                    if pd.isna(chosen_formula) or chosen_formula == '' or chosen_formula == 'nan':
                        chosen_formula = chosen_row.iloc[0]['formula']
                    chosen_makeup = pychnosz.makeup(chosen_formula)
                    n_elem_in_chosen = chosen_makeup.get(elem, 1)

                    # Coefficient = (count of elem in species) / (count of elem in basis)
                    coeff = count / n_elem_in_chosen

                    chosen_basis_species.append(chosen_sp_name)
                    chosen_basis_coeffs.append(coeff)
                else:
                    # Complex case: element has multiple oxidation states
                    # e.g., linnaeite has Co+2 and Co+3
                    # Match each oxidation state to the closest basis/aux species
                    chosen_dict = {}  # species_name -> coefficient

                    for ox_state, count in elem_ox_states:
                        # Find closest matching basis/aux species for this oxidation state
                        closest_sp = min(basis_aux_candidates, key=lambda x: abs(x[1] - ox_state))
                        chosen_sp_name = closest_sp[0]

                        # Get how many atoms of this element are in the chosen basis species
                        chosen_row = thermo_df[thermo_df['name'] == chosen_sp_name]
                        chosen_formula = chosen_row.iloc[0].get('formula_modded', chosen_row.iloc[0]['formula'])
                        if pd.isna(chosen_formula) or chosen_formula == '' or chosen_formula == 'nan':
                            chosen_formula = chosen_row.iloc[0]['formula']
                        chosen_makeup = pychnosz.makeup(chosen_formula)
                        n_elem_in_chosen = chosen_makeup.get(elem, 1)

                        # Coefficient for this instance
                        coeff = count / n_elem_in_chosen

                        # Accumulate coefficients for the same species
                        if chosen_sp_name in chosen_dict:
                            chosen_dict[chosen_sp_name] += coeff
                        else:
                            chosen_dict[chosen_sp_name] = coeff

                    # Add all chosen species and coefficients
                    for sp_name_inner, coeff in chosen_dict.items():
                        chosen_basis_species.append(sp_name_inner)
                        chosen_basis_coeffs.append(coeff)

        # New approach: use formula_ox_modded to directly identify which basis/aux species to use
        # This replicates the workflow that successfully works in testing

        # Build basis_to_elem mapping (basis species -> element/pseudoelement)
        basis_to_elem = {}
        for idx, row in thermo_df[thermo_df['tag'] == 'basis'].iterrows():
            try:
                formula_modded = row.get('formula_modded', row['formula'])
                if pd.isna(formula_modded) or formula_modded == '':
                    continue
                mkp = pychnosz.makeup(str(formula_modded))
                elems = [e for e in mkp.keys() if e not in ['H', 'O', 'Z']]
                if len(elems) == 1:
                    basis_to_elem[row['name']] = elems[0]
            except:
                continue

        # Reverse mapping: element/pseudoelement -> basis species
        elem_to_basis = {v: k for k, v in basis_to_elem.items()}

        # Parse formula_ox_modded to get element+oxidation combinations
        # Use shared function that handles pseudoelements
        element_ox_modded = parse_formula_ox_string(sp_formula_ox_modded, include_pseudoelements=True)

        # For each element+oxidation in formula_ox_modded, find the appropriate basis/aux species
        basis_to_add = []
        for el_with_charge in element_ox_modded.keys():
            # Extract just the element/pseudoelement part (remove charge and digits)
            el = el_with_charge.rstrip('0123456789+-')

            # Extract oxidation state from el_with_charge (e.g., "Cu+1" -> 1, "S-2" -> -2)
            ox_state_str = el_with_charge[len(el):]  # e.g., "+1", "-2"
            try:
                ox_state_val = int(ox_state_str) if ox_state_str else 0
            except:
                ox_state_val = 0

            # Try to find basis/aux species with matching oxidation state
            # Match by oxidation state in formula_ox, not just by species name
            if ox_state_val != 0:
                # Look for basis/aux species that contain this element in this oxidation state
                # Check aux_pref first, then basis_pref
                found = False

                # Check if we have aux species for this element
                if el in aux_pref and aux_pref[el]:
                    # Find aux species with closest matching oxidation state
                    best_match = None
                    min_ox_diff = float('inf')

                    for aux_sp_name in aux_pref[el]:
                        aux_row = thermo_df[thermo_df['name'] == aux_sp_name]
                        if not aux_row.empty:
                            aux_formula_ox = aux_row.iloc[0].get('formula_ox_modded',
                                                                   aux_row.iloc[0].get('formula_ox', ''))
                            # Get average oxidation state of element in this aux species
                            avg_ox = get_average_oxidation_state(aux_formula_ox, el)
                            ox_diff = abs(avg_ox - ox_state_val)

                            if ox_diff < min_ox_diff:
                                min_ox_diff = ox_diff
                                best_match = aux_sp_name

                    if best_match and min_ox_diff < 1.0:  # Close match
                        basis_to_add.append(best_match)
                        if verbose >= 2:
                            print(f"  Found {best_match} for {el_with_charge} (oxidation state match)")
                        found = True

                # If no aux match, try basis species
                if not found and el in basis_pref:
                    basis_sp_name = basis_pref[el]
                    basis_row = thermo_df[thermo_df['name'] == basis_sp_name]
                    if not basis_row.empty:
                        basis_to_add.append(basis_sp_name)
                        if verbose >= 2:
                            print(f"  Using basis {basis_sp_name} for {el_with_charge}")
                        found = True

                if found:
                    continue

            # Check if el is a pseudoelement and map to basis
            if el in elem_to_basis:
                basis_sp = elem_to_basis[el]
                basis_to_add.append(basis_sp)
                if verbose >= 2:
                    print(f"  Mapped {el} -> {basis_sp}")
                continue

        # Remove duplicates while preserving order
        seen = set()
        basis_to_add = [x for x in basis_to_add if not (x in seen or seen.add(x))]

        # Oxidation states are not always available. When formula_ox is blank
        # there is nothing to match basis species against, so fall back on the
        # basis species chosen from elemental composition alone.
        if len(basis_to_add) == 0 and len(chosen_basis_species) > 0:
            seen = set()
            basis_to_add = [x for x in chosen_basis_species
                            if not (x in seen or seen.add(x))]
            if verbose >= 2:
                print(f"  No oxidation states are available for {sp_name}, so "
                      "basis species were chosen from its elemental "
                      f"composition: {basis_to_add}")

        if verbose >= 2:
            print(f"  Basis species to add: {basis_to_add}")

        # Try with just H+ first (works for linnaeite)
        # If that fails, try with partial HOZ [H+, H2O] (needed for ilvaite)
        # If that fails, try with full HOZ_balancers (needed for laurite)
        hoz_for_basis = [hb.replace('(g)', '') if hb == 'O2(g)' else hb for hb in HOZ_balancers]

        try:
            # Try with just H+
            basis_list = basis_to_add + ['H+']
            seen = set()
            basis_list = [x for x in basis_list if not (x in seen or seen.add(x))]

            if verbose >= 2:
                print(f"  Trying basis with H+: {basis_list}")

            b = pychnosz.basis(basis_list, messages=False, global_state=False)

            # Use balance_reaction instead of subcrt for efficiency
            # This only calculates the stoichiometry without thermodynamic properties
            balance_result = pychnosz.balance_reaction([sp_name], [-1], basis=b, messages=False)

            if balance_result is None:
                raise Exception("balance_reaction returned None")

            balanced_species, coeffs = balance_result

        except Exception as e:
            # Try with partial HOZ (H+ and H2O only)
            if verbose >= 2:
                print(f"  H+ only failed ({e}), trying partial HOZ [H+, H2O]")

            try:
                basis_list = basis_to_add + ['H+', 'H2O']
                seen = set()
                basis_list = [x for x in basis_list if not (x in seen or seen.add(x))]

                if verbose >= 2:
                    print(f"  Trying basis with H+ and H2O: {basis_list}")

                b = pychnosz.basis(basis_list, messages=False, global_state=False)

                # Use balance_reaction instead of subcrt for efficiency
                balance_result = pychnosz.balance_reaction([sp_name], [-1], basis=b, messages=False)

                if balance_result is None:
                    raise Exception("balance_reaction returned None")

                balanced_species, coeffs = balance_result

            except Exception as e2:
                # Try with full HOZ_balancers
                if verbose >= 2:
                    print(f"  Partial HOZ failed ({e2}), trying full HOZ")

                basis_list = basis_to_add + hoz_for_basis
                seen = set()
                basis_list = [x for x in basis_list if not (x in seen or seen.add(x))]

                if verbose >= 2:
                    print(f"  Trying basis with full HOZ: {basis_list}")

                try:
                    b = pychnosz.basis(basis_list, messages=False, global_state=False)

                    # Use balance_reaction instead of subcrt for efficiency
                    balance_result = pychnosz.balance_reaction([sp_name], [-1], basis=b, messages=False)

                    if balance_result is None:
                        return None

                    balanced_species, coeffs = balance_result

                except Exception as e3:
                    if verbose >= 2:
                        print(f"  Error: {e3}")
                    return None

        # Convert species to names if they are indices
        # CRITICAL: indices from balance_reaction() are 1-based CHNOSZ OBIGT indices
        # These must be looked up using .loc[] not .iloc[]
        species_names = []
        thermo_sys = pychnosz.thermo()
        for sp in balanced_species:
            if isinstance(sp, (int, np.int64, np.int32)):
                # Species index from CHNOSZ - use .loc[] to look up (1-based label)
                if sp in thermo_sys.obigt.index:
                    species_names.append(thermo_sys.obigt.loc[sp]['name'])
                else:
                    # Fallback: species might not be in CHNOSZ but in thermo_df
                    sp_name_lookup = thermo_df[thermo_df['name'].notna()].query(f'name == "{sp}"')
                    if len(sp_name_lookup) > 0:
                        species_names.append(sp_name_lookup.iloc[0]['name'])
                    else:
                        species_names.append(str(sp))
            else:
                species_names.append(sp)

        # Filter out polymorphs of the dissociating species
        filtered_coeffs = []
        filtered_names = []
        found_primary = False
        for coeff, name in zip(coeffs, species_names):
            if name == sp_name:
                if not found_primary and coeff < 0:
                    filtered_coeffs.append(coeff)
                    filtered_names.append(name)
                    found_primary = True
            else:
                filtered_coeffs.append(coeff)
                filtered_names.append(name)

        coeffs = filtered_coeffs
        species_names = filtered_names

        # Replace 'water' with 'H2O' for EQ3 compatibility
        species_names = ['H2O' if n == 'water' else n for n in species_names]

        # Restore (g) suffix for gas species that were in HOZ_balancers
        # e.g., convert 'O2' back to 'O2(g)' if 'O2(g)' was in HOZ_balancers
        for i, name in enumerate(species_names):
            for hb in HOZ_balancers:
                if '(g)' in hb:
                    # Extract base name without (g) suffix
                    base_name = hb.replace('(g)', '')
                    if name == base_name:
                        species_names[i] = hb
                        break

        # Convert coefficients to integers if needed (R code lines 728-761)
        if any(abs(c % 1) > 1e-6 for c in coeffs):
            coeffs = coefficients_to_integers(coeffs)

        # Reduce coefficients to their simplest form by dividing by their GCD
        coeffs = reduce_coefficients(coeffs)

        # Format as string
        parts = []
        for coeff, name in zip(coeffs, species_names):
            parts.append(f"{coeff:.4f}")
            parts.append(name)

        dissrxn_str = " ".join(parts)

        if verbose >= 2:
            print(f"Generated dissociation reaction for {sp_name}:")
            print(f"  {dissrxn_str}")

        return dissrxn_str

    except Exception as e:
        if verbose >= 1:
            print(f"Error generating dissociation reaction for {sp_name}: {e}")
        import traceback
        if verbose >= 2:
            traceback.print_exc()
        return None


def find_species_needing_dissrxns(thermo_df, verbose=1):
    """
    Find all species that need dissociation reactions (missing or unbalanced).

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database dataframe
    verbose : int
        Verbosity level

    Returns
    -------
    list
        List of species names that need dissociation reactions
    """
    from .dissrxn_balancer import check_dissrxn_balanced

    needs_dissrxn = []

    for idx, row in thermo_df.iterrows():
        # Skip basis species
        if row.get('tag') == 'basis':
            continue

        sp_name = row['name']
        dissrxn = row.get('dissrxn', '')

        # Check if dissociation reaction is missing or unbalanced
        # Handle both string and NaN/None cases
        if pd.isna(dissrxn) or dissrxn == '' or (isinstance(dissrxn, str) and dissrxn.strip() == ''):
            needs_dissrxn.append(sp_name)
        elif not check_dissrxn_balanced(str(dissrxn), sp_name, verbose=verbose >= 2):
            needs_dissrxn.append(sp_name)

    return needs_dissrxn


def find_invalid_dissrxn_species(thermo_df, fixed_species=None):
    """
    Find species whose dissociation reactions are written in terms of species
    that are not strict basis or auxiliary basis species.

    A dissociation reaction may only be written in terms of strict basis
    species ('basis' in the 'tag' column) and auxiliary basis species ('aux' in
    the 'tag' column). A reaction that requires any other species is balanced
    but cannot be compiled, so EQPT fails with an error like:

        * Error - (EQPT/pcraq) The reaction which destroys non-basis species
              NaCO3- is written in terms of an unrecognized basis species
              "CO3-2"

    Note that a balance check will not catch these reactions because the
    offending species can still have a valid chemical formula.

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database dataframe with 'name', 'tag', and 'dissrxn'
        columns
    fixed_species : list, optional
        List of fixed species (H2O, H+, O2(g), etc.). These are always valid
        because EQ3/6 has them hard-coded.

    Returns
    -------
    dict
        Dictionary mapping species names to the species in their dissociation
        reactions that are neither strict nor auxiliary basis species. Species
        that are missing from the database entirely are included, too.
    """
    if fixed_species is None:
        fixed_species = ["H2O", "H+", "O2(g)", "water", "e-"]

    if 'dissrxn' not in thermo_df.columns or 'tag' not in thermo_df.columns:
        return {}

    valid_names = set(thermo_df.loc[thermo_df['tag'].isin(['basis', 'aux']), 'name'])
    valid_names.update(fixed_species)

    invalid_dissrxns = {}

    for idx, row in thermo_df.iterrows():
        sp_name = row['name']
        dissrxn = row.get('dissrxn', '')

        if not isinstance(dissrxn, str) or dissrxn.strip() == '' or dissrxn == 'nan':
            continue

        # species names are at odd positions in a dissociation reaction, and
        # the first one is the dissociating species itself
        parts = dissrxn.strip().split()
        species_in_rxn = [parts[i] for i in range(1, len(parts), 2)][1:]

        invalid = [s for s in species_in_rxn
                   if s not in valid_names and s != sp_name]

        if len(invalid) > 0:
            # polymorphs share a name, so combine what each row requires
            # while removing duplicates and preserving order
            invalid = invalid_dissrxns.get(sp_name, []) + invalid
            invalid_dissrxns[sp_name] = list(dict.fromkeys(invalid))

    return invalid_dissrxns
