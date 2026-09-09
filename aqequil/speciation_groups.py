"""
Speciation groups and limiting reactant selection shared by energy supply
calculations.

Before an energy supply can be reported, the code has to decide how much of
each reactant is actually available. That decision involves speciation
groups: a microbe oxidizing sulfide can draw on every dissolved form of
sulfide (H2S, HS-, and the sulfide complexes), not just the one form that
happens to appear in the reaction the user wrote.

Both `Speciation` (AqSpeciation.py) and `Mass_Transfer` (MassTransfer.py)
need this, so the logic lives here as a mixin. Keeping it in one place is
what prevents the two classes from disagreeing about the energy supply of
the same reaction written two different ways.
"""

import collections
import copy

import pandas as pd
import pychnosz

from wormutils import chemlabel, format_equation, import_package_file, isnotebook


class SpeciationGroups(object):

    """
    Mixin providing speciation group lookups and limiting reactant selection.

    Classes mixing this in are expected to provide:
      - `self.thermo.csv_db` : WORM-style thermodynamic database dataframe
      - `self.err_handler`   : a wormutils Error_Handler
      - `self.verbose`       : verbosity level
    """

    # Class-level defaults so classes that never set these in their __init__
    # (e.g. Mass_Transfer) still work. Instance attributes shadow these.
    custom_grouping_filepath = None
    reactant_dict_scalar = None
    speciation_group_dict_unpacked = None

    # the grouping file the current dictionaries were built from
    _speciation_group_source = None

    def _set_speciation_groups(self, custom_grouping_filepath=None, reset=False):

        """
        Make sure the speciation group dictionaries exist and are up to date.

        They are rebuilt if they have never been built, or if the grouping file
        they were built from is not the one being asked for now. Without the
        latter check, a `custom_grouping_filepath` passed to a second energy
        calculation would be silently ignored.

        Parameters
        ----------
        custom_grouping_filepath : str, optional
            Filepath of a TXT file of customized speciation groups.

        reset : bool, default False
            If True, `custom_grouping_filepath=None` means "go back to the
            built-in speciation groups" instead of "keep whatever grouping file
            was set before".
        """

        if isinstance(custom_grouping_filepath, str) or reset:
            self.custom_grouping_filepath = custom_grouping_filepath

        if (not isinstance(self.reactant_dict_scalar, dict)
                or self._speciation_group_source != self.custom_grouping_filepath):
            self._make_speciation_group_dict()
            self._speciation_group_source = self.custom_grouping_filepath

    def _make_speciation_group_dict(self):

        if isinstance(self.custom_grouping_filepath, str):
            with open(self.custom_grouping_filepath) as file:
                lines = [line.rstrip() for line in file]
        else:
            content = import_package_file('aqequil.databases', 'speciation_groups_WORM.txt')
            content = content.split("\n")
            lines = [line.rstrip() for line in content]

        lines = [l for l in lines if l != ""] # remove blank lines

        # Creates a dictionary with this format:
        #
        # {...
        #  "sulfide 1": ["H2S", "HS-"],
        #  "sulfide 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "sulfide 3": ["Pb(HS)3-"],
        #  "iron(II) 1": [...],
        #  ...
        # }
        #
        # Used by calculate_energy() to calculate concentrations of limiting reactants.
        try:
            reactant_dict_scalar = {l.split(":")[0]:l.split(":")[1].strip() for l in lines}
        except:
            bad_line_indices = []
            bad_lines = []
            for i,l in enumerate(lines):
                if ":" not in l:
                    bad_line_indices.append(i)
                    bad_lines.append(l)

            bad_line_indices = [i+1 for i in bad_line_indices]
            self.err_handler.raise_exception("Line(s) in the speciation group file do not contain a colon ':'. Fix the formatting on line(s):\n"+str(bad_line_indices)+"\nLines to blame:\n"+str(bad_lines))


        self.reactant_dict_scalar = {k:v.split(" ") for k,v in zip(reactant_dict_scalar.keys(), reactant_dict_scalar.values())}

        # Creates a dictionary with this format:
        # {...
        #  "H2S": ["H2S", "HS-"],
        #  "HS-": ["H2S", "HS-"],
        #  "Pb(HS)2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Ag(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Au(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  ...
        # }
        #
        # Used to switch a user-specified limiting reactant to one within the same scalar group.
        speciation_group_dict = {}
        for i,line in enumerate(lines):
            line = line.strip().split(":")
            assert len(line) == 2 # will fail if : in species names
            group_name = line[0].strip()
            group_species = line[1].strip().split(" ")
            group_species = list(collections.OrderedDict.fromkeys(group_species)) # remove duplicates
            speciation_group_dict[group_name] = group_species

        speciation_group_dict_unpacked = {}
        for key in list(speciation_group_dict.keys()):
            for sp in speciation_group_dict[key]:
                speciation_group_dict_unpacked[sp] = speciation_group_dict[key]

        self.speciation_group_dict_unpacked = speciation_group_dict_unpacked

        # Creates a dictionary with this format:
        # {...
        #  "sulfate": ['SO4-2','HSO4-','PdSO4','RhSO4', ..., 'Ru(SO4)3-4'],
        #  "formate": ['formic-acid', 'formate', 'Am(For)+2', ..., 'Zn(For)2'],
        #  ...
        # }
        # Used when the user wants to exclude all organics except for acetate (e.g.)
        # and all species in the acetate speciation group
        group_categories = []
        for group_name in self.reactant_dict_scalar.keys():
            group_name_no_number = group_name.split(" ")
            group_name_no_number = " ".join(group_name_no_number[:-1])
            group_categories.append(group_name_no_number)
        group_categories = list(set(group_categories))

        sp_dict_groups = {}
        for g1 in group_categories:
            for i,g2 in enumerate(self.reactant_dict_scalar.keys()):
                g2_cat_name = g2.split(" ")
                g2_cat_name = " ".join(g2_cat_name[:-1])
                if g1 == g2_cat_name:
                    if g1 not in sp_dict_groups.keys():
                        sp_dict_groups[g1] = copy.deepcopy(self.reactant_dict_scalar[list(self.reactant_dict_scalar.keys())[i]])
                    else:
                        sp_dict_groups[g1] += copy.deepcopy(self.reactant_dict_scalar[list(self.reactant_dict_scalar.keys())[i]])
        self.sp_dict_groups = sp_dict_groups


    def _match_grouped_species(self, s):
        """
        Match whether a species is in a speciation group, and get a list of relevant groups and their scalars.

        e.g.,

        self.reactant_dict_scalar = {...
                                     "sulfides 1": ["H2S", "HS-"],
                                     "sulfides 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
                                     "sulfides 3": ["Pb(HS)3-"],
                                     "ferrous iron 1": [...],
                                     ...
                                     }
        If `s` = "Ag(HS)2-"
        then `scalars` = [1, 2, 3]
        and `groups` = ["sulfides 1", "sulfides 2", "sulfides 3"]

        """
        scalars = []
        groups = []

        for i,grp_list in enumerate(list(self.reactant_dict_scalar.values())):
            if s in grp_list:
                s_key = list(self.reactant_dict_scalar.keys())[i] #e.g., s_key can be "sulfates 1"
                s_key_split = s_key.split(" ")
                s_key_grp = " ".join(s_key_split[:-1])
                possible_scalars = [k.split(s_key_grp)[-1].strip() for k in list(self.reactant_dict_scalar.keys()) if s_key_grp in k]
                for k in possible_scalars:
                    try:
                        float(k) # test whether the scalar is a number. If so, append.
                        scalars.append(k)
                    except:
                        continue
                groups = [s_key_grp+" "+str(s) for s in scalars]
                break

        scalars = [float(s) for s in scalars]

        # prevent duplicate groups and scalars from appearing
        scalars_final = []
        groups_final = []
        for i,group in enumerate(groups):
            if group not in groups_final:
                scalars_final.append(scalars[i])
                groups_final.append(groups[i])

        groups_final = [self.reactant_dict_scalar[g] for g in groups_final]
        return scalars_final, groups_final


    def _switch_limiting(self, limiting, stoich, species, lenient=False):
        if limiting != None:
            reactant_idx = [1 if i<0 else 0 for i in stoich]
            reactants = [species[i] for i,idx in enumerate(reactant_idx) if idx == 1]

            if limiting not in reactants and limiting in list(self.speciation_group_dict_unpacked.keys()):
                # if the user specifies a limiting reactant like "HCO3-" but the
                # reaction has "CO2" as a reactant, check the dict that
                # contains speciated groups and switch the limiting reactant to
                # the relevant limiting reactant to appear in reports,
                # e.g., "HCO3-" -> "CO2"
                lim_species_group_list = self.speciation_group_dict_unpacked[limiting]
                for s in lim_species_group_list:
                    if s in reactants:
                        if self.verbose > 0:
                            print("The specified limiting reactant", str(limiting),
                                  "has been switched to", str(s), "because the latter",
                                  "appears as a reactant in the reaction:")
                            if isnotebook() and hasattr(self, "format_reaction"):
                                _ = self.format_reaction(coeffs=stoich,
                                   names=species,
                                   formatted=True,
                                   charge_sign_at_end=True,
                                   show=True)
                            else:
                                equation = format_equation(species, stoich,
                                                      charge_sign_at_end=True)
                                try:
                                    print(equation)
                                except UnicodeEncodeError:
                                    # e.g., a Windows console that cannot encode
                                    # the rightwards arrow
                                    print(equation.replace("→", "->"))
                        limiting = s
                        break

            if limiting not in reactants and lenient:
                # if the user specifies "CO2" as the limiting reactant during a
                # batch calculation of many different reactions, some reactions
                # won't actually have "CO2" as a reactant. In this case, set
                # limiting to None so that CO2 will be limiting when applicable.
                limiting = None

        return limiting


    def _species_state(self, s):
        """
        Physical state ('aq', 'cr', 'liq', 'gas') of a species according to the
        thermodynamic database CSV.
        """
        csv_db = self.thermo.csv_db
        return list(csv_db[csv_db["name"] == s]["state"])[0]


    @staticmethod
    def _mineral_amounts_known(grams_minerals):
        """
        Whether `grams_minerals` says enough about how much mineral or liquid
        reactant is present for those reactants to be treated as potentially
        limiting.
        """
        return isinstance(grams_minerals, (pd.DataFrame, dict, float, int))


    def _mineral_molality(self, s, n_rows, grams_minerals):

        """
        Moles of a mineral or liquid reactant available at each sample or step
        of reaction progress.

        Unlike an aqueous species, a mineral has no molality in the EQ3/6
        output; how much of it a microbe has to work with is something the user
        has to state. `grams_minerals` is that statement, and it is converted to
        moles here using the formula mass from the thermodynamic database.

        Parameters
        ----------
        s : str
            Name of the mineral or liquid species.

        n_rows : int
            Number of samples (`Speciation`) or steps of reaction progress
            (`Mass_Transfer`) being evaluated.

        grams_minerals : float, int, dict, pandas.DataFrame, or None
            Grams of mineral present. A float or int applies to every mineral
            reactant; a dict gives a value per mineral; a dataframe with a
            column per mineral gives a value for each row. Anything else (e.g.
            None) means the amount is unknown, in which case NaN is returned and
            the mineral cannot be limiting. Note that this is not the same as a
            mass of 0 grams, which makes the mineral limiting straight away.

        Returns
        ----------
        list of float
            Moles of `s` available, one value per row.
        """

        if not self._mineral_amounts_known(grams_minerals):
            return [float("NaN")]*n_rows

        if pychnosz.thermo().element is None:
            # touching the thermodynamic system initializes it, and that is what
            # loads the element masses that pychnosz.mass() needs. Nothing can
            # be lost by doing this here because an uninitialized system has no
            # thermodynamic data in it yet.
            _ = pychnosz.thermo().OBIGT

        csv_db = self.thermo.csv_db
        sp_formula = list(csv_db[csv_db["name"] == s]["formula"])[0]
        sp_mass = pychnosz.mass(sp_formula)

        if isinstance(grams_minerals, pd.DataFrame):
            if s not in grams_minerals.columns:
                return [0]*n_rows
            sp_grams_list = list(grams_minerals[s])
            if len(sp_grams_list) != n_rows:
                self.err_handler.raise_exception("The '"+str(s)+"' column of "
                        "grams_minerals has "+str(len(sp_grams_list))+" values "
                        "but "+str(n_rows)+" are needed, one for each row of "
                        "results.")
            return [sp_grams/sp_mass for sp_grams in sp_grams_list]

        if isinstance(grams_minerals, dict):
            sp_grams = grams_minerals.get(s, 0)
        else:
            sp_grams = grams_minerals

        return [sp_grams/sp_mass]*n_rows


    def _grouped_molality(self, s, molal_df, as_written=False):

        """
        Total molality available to a reactant, summed over its speciation
        group and weighted by how many atoms of the grouped element each
        species contributes.

        For example, the sulfide group is
        1*(H2S, HS-) + 2*(Pb(HS)2, Ag(HS)2-, Au(HS)2-) + 3*(Pb(HS)3-),
        so a reaction written with either "H2S" or "HS-" draws on the same
        total inventory of dissolved sulfide.

        Parameters
        ----------
        s : str
            Name of the species as written in the reaction.

        molal_df : pandas.DataFrame
            Table of molalities with species names as columns. One row per
            sample (`Speciation`) or per step of reaction progress
            (`Mass_Transfer`).

        as_written : bool, default False
            If True, do not sum over the speciation group; report only the
            molality of `s` itself, scaled by its group's scalar. Use this
            when the reaction is meant to be limited by one particular form
            of a species rather than by the total.

        Returns
        ----------
        list of float
            Molality available to `s`, one value per row of `molal_df`.
        """

        self._set_speciation_groups()

        # check whether the species matches any of the groups in
        # reactant_dict_scalar and retrieve scalars and groups
        scalars, groups = self._match_grouped_species(s)

        if as_written:
            # if energy supplies are to be calculated as written with no
            # grouping, then use the current species and its scalar
            for i, sc in enumerate(scalars):
                if s in groups[i]:
                    scalars = [sc]
                    groups = [[s]]
                    break

        if len(scalars) == 0:
            # `s` does not belong to any speciation group (e.g. O2, H2, CH4)
            return list(molal_df[s])

        total_summed_scaled = [0]*molal_df.shape[0]

        for i,scalar in enumerate(scalars):
            col_subset = [col for col in groups[i] if col in molal_df.columns]
            scaled_df = molal_df[col_subset].apply(lambda x: x*scalar)
            summed_scaled = list(scaled_df.sum(axis=1, numeric_only=True))
            total_summed_scaled = [ii+iii for ii,iii in zip(total_summed_scaled, summed_scaled)]

        return total_summed_scaled


    def _select_limiting_reactant(self, species, stoich, s_molal_dict, step,
                                        limiting=None, allow_minerals=False,
                                        charge_sign_at_end=False):

        """
        Identify the limiting reactant of a reaction at one sample or one step
        of reaction progress.

        Parameters
        ----------
        species, stoich : lists
            The reaction, as passed to the energy calculation.

        s_molal_dict : dict
            Maps each species name to a list of molalities (already summed
            over speciation groups by `_grouped_molality`).

        step : int
            Index into the lists in `s_molal_dict`; the sample or step of
            reaction progress being evaluated.

        limiting : str, optional
            Force this species to be the limiting reactant instead of
            choosing one by concentration and stoichiometry.

        allow_minerals : bool, default False
            Permit 'cr' and 'liq' reactants to be limiting. Only meaningful
            when mineral abundances are known (see the `grams_minerals`
            parameter of `Speciation.calculate_energy`).

        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-) in the
            reported limiting reactant name.

        Returns
        ----------
        tuple
            (lr_name, lr_reported, lr_concentration, lr_stoich). All four are
            None/NaN when the reaction has no valid limiting reactant, or when
            a limiting reactant concentration is unknown (NaN).
        """

        none_result = (None, float("nan"), float("nan"), float("nan"))

        if isinstance(limiting, str):
            lrc_dict = {}
            lr_list = [limiting]
        else:
            lrc_dict = {}
            for i_s,s in enumerate(species):
                # identify valid limiting reactants and record concentrations
                # 1. negative coefficient (reactant)
                # 2. can't be OH-, H+, H2O
                # 3. can't be cr or liq unless mineral abundances are known
                if stoich[i_s] >= 0:
                    continue
                if s in ["H2O", "H+", "OH-"]:
                    continue
                if not allow_minerals and self._species_state(s) in ["cr", "liq"]:
                    continue
                lrc_dict[s] = s_molal_dict[s][step]/abs(stoich[i_s])

            if len(lrc_dict) == 0:
                return none_result

            lr_val = min(lrc_dict.values())

            # handle situations where there might be multiple limiting reactants
            lr_list = [str(k) for k,v in lrc_dict.items() if v == lr_val]

        # a reactant of unknown concentration makes the energy supply undefined
        if any([pd.isna(v) for v in lrc_dict.values()]):
            return none_result

        lr_list_formatted = [chemlabel(lr_name, charge_sign_at_end=charge_sign_at_end) for lr_name in lr_list]
        if len(lr_list_formatted) > 1:
            lr_reported = ", ".join(lr_list_formatted)
        else:
            lr_reported = lr_list[0]

        lr_name = lr_list[0] # doesn't matter which lr is used to calculate
        lr_concentration = s_molal_dict[lr_name][step]
        lr_stoich = -stoich[species.index(lr_name)]

        return lr_name, lr_reported, lr_concentration, lr_stoich
