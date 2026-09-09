"""
Support for Pitzer ion-interaction parameters in aqequil.

The Pitzer parameter database is a CSV file in "long" format, with one row
per interaction parameter:

    species1,species2,species3,param,a1,a2,a3,a4,alpha,ref1,ref2,date,note

Columns
-------
species1, species2, species3 : str
    Names of aqueous species exactly as they appear in the thermodynamic
    database used to build the data0 file (e.g., "Na+", "SO4-2", "CO2").
    `species3` is left blank for pair parameters. For a lambda or mu
    parameter involving a single neutral species (the EQ3/6 "nn" block),
    `species2` may be left blank or set equal to `species1`.
param : str
    One of "beta0", "beta1", "beta2", "cphi" (cation-anion pairs),
    "theta" (cation-cation or anion-anion pairs), "lambda" (pairs involving
    a neutral species), "psi" (cation-cation-anion or anion-anion-cation
    triplets), "zeta" (neutral-cation-anion triplets), or "mu"
    (neutral-neutral-neutral triplets).
a1, a2, ... aN : float
    Coefficients of the EQ3/6 "Livermore" temperature function
    (T in Kelvin, Tr = 298.15 K):

        x(T) = a1 + a2*(1/T - 1/Tr) + a3*ln(T/Tr) + a4*(T - Tr) + a5*(T^2 - Tr^2)

    Up to five coefficient columns are supported. A database with only an
    `a1` column holds 25 °C values that are treated as constant with
    temperature. Blank coefficients are treated as zero.
alpha : float, optional
    Ionic strength dependence parameter. Only used on "beta1" rows (where
    it is alpha(1)) and "beta2" rows (where it is alpha(2)). When blank,
    EQ3/6 assigns the usual defaults (alpha(1) = 2.0 for pairs that include
    a monovalent ion; alpha(1) = 1.4 and alpha(2) = 12.0 otherwise).
ref1, ref2, date, note : str, optional
    Provenance information. Not used in calculations.
"""

import re
import copy

import numpy as np
import pandas as pd


PITZER_PAIR_PARAMS = ["beta0", "beta1", "beta2", "cphi", "theta", "lambda"]
PITZER_TRIPLET_PARAMS = ["psi", "zeta", "mu"]
PITZER_PARAMS = PITZER_PAIR_PARAMS + PITZER_TRIPLET_PARAMS
PITZER_CA_PARAMS = ["beta0", "beta1", "beta2", "cphi"]

PITZER_REQUIRED_COLUMNS = ["species1", "species2", "param", "a1"]
PITZER_OPTIONAL_COLUMNS = ["species3", "a2", "a3", "a4", "a5", "alpha",
                           "ref1", "ref2", "date", "note"]

# The maximum number of temperature function terms supported by the
# EQ3/6 "Livermore" Pitzer temperature function.
PITZER_MAX_TERMS = 5

# Superblock headers, in the order required by EQPT. The first 24
# characters of each header are matched by EQPT, so they must not change.
PITZER_SUPERBLOCK_HEADERS = {
    "ca": "ca combinations: beta(n)(ca) and Cphi(ca) [optional: alpha(n)(ca)]",
    "theta": "cc' and aa' combinations: theta(cc') and theta(aa')",
    "lambda_ion": "nc and na combinations: lambda(nc) and lambda(na)",
    "nn": "nn combinations: lambda(nn) and mu(nnn)",
    "nn_prime": "nn' combinations: lambda(nn')",
    "psi": "cc'a and aa'c combinations: psi(cc'a) and psi(aa'c)",
    "zeta": "nca combinations: zeta(nca)",
    "mu": "nnn' combinations: mu(nnn')",
}
PITZER_SUPERBLOCK_ORDER = ["ca", "theta", "lambda_ion", "nn", "nn_prime",
                           "psi", "zeta", "mu"]

DATA0_DELIMITER = "+--------------------------------------------------------------------"


class PitzerDatabaseError(ValueError):
    """Raised when a Pitzer parameter database is malformed."""
    pass


def _is_blank(x):
    if x is None:
        return True
    if isinstance(x, float) and np.isnan(x):
        return True
    return str(x).strip() == ""


def _fmt_coeff(x):
    """
    Format a temperature function coefficient for a data0 file.

    EQPT reads these values with a Fortran E25.18 edit descriptor, which
    treats a number without a decimal point as having 18 implied decimal
    places. Always writing an exponential form with a decimal point avoids
    that pitfall.
    """
    x = float(x)
    if x == 0:
        return "0.0"
    s = repr(x)  # shortest string that round-trips to the same float
    if "e" in s:
        mant, exp = s.split("e")
        if "." not in mant:
            mant += ".0"
        s = f"{mant}E{int(exp):+03d}"
    elif "." not in s:
        s += ".0"
    if len(s) > 25:
        s = f"{x:.16E}"
    return s


def get_pitzer_coeff_columns(df):
    """
    Return the list of temperature-function coefficient columns (a1, a2, ...)
    present in a Pitzer parameter DataFrame, in order.
    """
    cols = []
    for i in range(1, PITZER_MAX_TERMS + 1):
        col = f"a{i}"
        if col in df.columns:
            cols.append(col)
        else:
            break
    return cols


def pitzer_param_keys(df):
    """
    Return a list with one key per row of a validated Pitzer parameter
    database. A key identifies which parameter a row holds: the species the
    parameter applies to (order does not matter) paired with the parameter
    name. Two rows with the same key are two definitions of the same parameter.
    """
    keys = []
    for _, row in df.iterrows():
        names = [row["species1"], row["species2"]]
        if not _is_blank(row["species3"]):
            names.append(row["species3"])
        keys.append((tuple(sorted(names)), row["param"]))
    return keys


def format_pitzer_param_key(key):
    """
    Format a key returned by `pitzer_param_keys` for display, e.g.,
    "beta0 for Ca+2, Cl-".
    """
    names, param = key
    return f"{param} for {', '.join(names)}"


def validate_pitzer_db(df, filename="Pitzer parameter database"):
    """
    Validate and normalize a Pitzer parameter database.

    Parameters
    ----------
    df : pd.DataFrame
        Pitzer parameter database in the long CSV format described in the
        module docstring.
    filename : str
        Name used in error messages.

    Returns
    -------
    pd.DataFrame
        A normalized copy of the database (stripped species names,
        lower-case parameter names, numeric coefficient columns, and an
        `species3` column that is present even if it was missing).
    """

    if not isinstance(df, pd.DataFrame):
        raise PitzerDatabaseError(f"{filename} must be a pandas DataFrame.")

    df = copy.deepcopy(df)
    df.columns = [str(c).strip() for c in df.columns]

    missing = [c for c in PITZER_REQUIRED_COLUMNS if c not in df.columns]
    if len(missing) > 0:
        raise PitzerDatabaseError(
            f"{filename} is missing the required column(s) {missing}. "
            "A Pitzer parameter database must have the columns "
            f"{PITZER_REQUIRED_COLUMNS} and may have the optional columns "
            f"{PITZER_OPTIONAL_COLUMNS}.")

    coeff_cols = get_pitzer_coeff_columns(df)
    unexpected_a = [c for c in df.columns
                    if re.fullmatch(r"a\d+", c) and c not in coeff_cols]
    if len(unexpected_a) > 0:
        raise PitzerDatabaseError(
            f"{filename} has coefficient column(s) {unexpected_a} that are "
            f"not supported. Coefficient columns must be named a1 through "
            f"a{PITZER_MAX_TERMS} without gaps.")

    if "species3" not in df.columns:
        df["species3"] = np.nan
    if "alpha" not in df.columns:
        df["alpha"] = np.nan

    for col in ["species1", "species2", "species3", "param"]:
        df[col] = df[col].apply(lambda x: np.nan if _is_blank(x) else str(x).strip())
    df["param"] = df["param"].str.lower()

    for col in coeff_cols + ["alpha"]:
        try:
            df[col] = pd.to_numeric(df[col].apply(lambda x: np.nan if _is_blank(x) else x))
        except Exception as e:
            raise PitzerDatabaseError(
                f"Column '{col}' of {filename} contains non-numeric values: {e}")

    for col in coeff_cols:
        df[col] = df[col].fillna(0.0)

    bad_params = sorted(set([p for p in df["param"].dropna() if p not in PITZER_PARAMS]))
    if len(bad_params) > 0:
        raise PitzerDatabaseError(
            f"{filename} contains unrecognized parameter name(s) {bad_params} "
            f"in the 'param' column. Valid names are {PITZER_PARAMS}.")

    if df["param"].isnull().any() or df["species1"].isnull().any():
        raise PitzerDatabaseError(
            f"Every row of {filename} must have a 'species1' and a 'param' value.")

    # Pair parameters need two species (the nn lambda/mu case allows a blank
    # species2). Triplet parameters need three species, except mu for a
    # single neutral species.
    for i, row in df.iterrows():
        p = row["param"]
        s1, s2, s3 = row["species1"], row["species2"], row["species3"]
        if p in ["lambda", "mu"]:
            if _is_blank(s2):
                df.at[i, "species2"] = s1
            if p == "mu" and _is_blank(s3):
                df.at[i, "species3"] = s1
        elif p in PITZER_PAIR_PARAMS:
            if _is_blank(s2):
                raise PitzerDatabaseError(
                    f"Row {i} of {filename} ({p} for {s1}) needs a 'species2' value.")
            if not _is_blank(s3):
                raise PitzerDatabaseError(
                    f"Row {i} of {filename} ({p} for {s1}, {s2}) is a pair "
                    f"parameter but has a 'species3' value ('{s3}').")
        else:  # psi, zeta
            if _is_blank(s2) or _is_blank(s3):
                raise PitzerDatabaseError(
                    f"Row {i} of {filename} ({p} for {s1}) is a triplet "
                    "parameter and needs 'species2' and 'species3' values.")

    # check for duplicate parameter entries
    keys = pitzer_param_keys(df)
    dupes = sorted(set([k for k in keys if keys.count(k) > 1]))
    if len(dupes) > 0:
        dupe_str = "; ".join([format_pitzer_param_key(k) for k in dupes])
        raise PitzerDatabaseError(
            f"{filename} contains duplicate entries: {dupe_str}")

    return df


def _classify_pair(z1, z2, s1, s2):
    """Return the superblock key for a pair of species with charges z1, z2."""
    if z1 * z2 < 0:
        return "ca"
    if z1 != 0 and z2 != 0:
        return "theta"
    if z1 == 0 and z2 == 0:
        return "nn" if s1 == s2 else "nn_prime"
    return "lambda_ion"


def build_pitzer_superblocks(pitzer_df, species_charges, verbose=1):
    """
    Build the text of the Pitzer parameter superblocks of a data0 file
    (everything between the miscellaneous parameters and the elements
    block).

    Parameters
    ----------
    pitzer_df : pd.DataFrame
        Validated Pitzer parameter database (see `validate_pitzer_db`).
    species_charges : dict
        Dictionary of {species name: charge} for the aqueous species that
        will be present in the data0 file. Parameters involving species that
        are not in this dictionary are left out of the data0 file.
    verbose : int
        Verbosity level.

    Returns
    -------
    dict
        {"text": str, "n_terms": int, "skipped_species": list,
         "n_blocks": int}
        `text` starts with the "ca combinations" header line and ends with
        the delimiter line that precedes "elements".
    """

    pitzer_df = validate_pitzer_db(pitzer_df)
    coeff_cols = get_pitzer_coeff_columns(pitzer_df)
    n_terms = max(len(coeff_cols), 1)

    def a_lines(row):
        lines = []
        for j, col in enumerate(coeff_cols):
            lines.append(f"    a{j+1} = {_fmt_coeff(row[col])}")
        return lines

    def zeros():
        return [f"    a{j+1} = 0.0" for j in range(n_terms)]

    def names_line(*names):
        return "  ".join([f"{n:<24}" for n in names]).rstrip()

    superblocks = {k: [] for k in PITZER_SUPERBLOCK_ORDER}
    skipped_species = set()
    skipped_rows = []
    ca_blocks = {}    # (cation, anion) -> {param: row}
    nn_blocks = {}    # neutral -> {"lambda": row, "mu": row}
    lambda_pairs = set()   # (neutral, ion) pairs with lambda data
    nn_prime_pairs = []    # (neutral, neutral') pairs with lambda data
    psi_triplets = []      # (like ion, like ion', counter ion)
    zeta_triplets = []     # (neutral, cation, anion)

    for i, row in pitzer_df.iterrows():
        p = row["param"]
        names = [row["species1"], row["species2"]]
        if not _is_blank(row["species3"]):
            names.append(row["species3"])

        absent = [n for n in names if n not in species_charges]
        if len(absent) > 0:
            skipped_species.update(absent)
            skipped_rows.append(i)
            continue

        z = [float(species_charges[n]) for n in names]

        if p in PITZER_CA_PARAMS:
            if len(names) != 2 or z[0] * z[1] >= 0:
                raise PitzerDatabaseError(
                    f"The parameter {p} for {', '.join(names)} requires one "
                    "cation and one anion.")
            cation, anion = (names[0], names[1]) if z[0] > 0 else (names[1], names[0])
            ca_blocks.setdefault((cation, anion), {})[p] = row

        elif p == "theta":
            if z[0] == 0 or z[1] == 0 or z[0] * z[1] < 0 or names[0] == names[1]:
                raise PitzerDatabaseError(
                    f"The parameter theta for {', '.join(names)} requires two "
                    "different cations or two different anions.")
            superblocks["theta"].append(
                [names_line(*names), "  theta:"] + a_lines(row))

        elif p == "lambda":
            kind = _classify_pair(z[0], z[1], names[0], names[1])
            if kind == "ca" or kind == "theta":
                raise PitzerDatabaseError(
                    f"The parameter lambda for {', '.join(names)} requires at "
                    "least one neutral species.")
            if kind == "nn":
                nn_blocks.setdefault(names[0], {})["lambda"] = row
            elif kind == "lambda_ion":
                # neutral species first
                neutral, ion = (names[0], names[1]) if z[0] == 0 else (names[1], names[0])
                lambda_pairs.add((neutral, ion))
                superblocks["lambda_ion"].append(
                    [names_line(neutral, ion), "  lambda:"] + a_lines(row))
            else:
                nn_prime_pairs.append((names[0], names[1]))
                superblocks["nn_prime"].append(
                    [names_line(*names), "  lambda:"] + a_lines(row))

        elif p == "psi":
            if any([zi == 0 for zi in z]):
                raise PitzerDatabaseError(
                    f"The parameter psi for {', '.join(names)} requires three ions.")
            pos = [n for n, zi in zip(names, z) if zi > 0]
            neg = [n for n, zi in zip(names, z) if zi < 0]
            if len(pos) == 2 and len(neg) == 1:
                ordered = pos + neg
            elif len(neg) == 2 and len(pos) == 1:
                ordered = neg + pos
            else:
                raise PitzerDatabaseError(
                    f"The parameter psi for {', '.join(names)} requires two "
                    "ions of like charge and one ion of opposite charge.")
            if ordered[0] == ordered[1]:
                raise PitzerDatabaseError(
                    f"The parameter psi for {', '.join(names)} requires two "
                    "different ions of like charge.")
            psi_triplets.append(tuple(ordered))
            superblocks["psi"].append(
                [names_line(*ordered), "  psi:"] + a_lines(row))

        elif p == "zeta":
            neutral = [n for n, zi in zip(names, z) if zi == 0]
            pos = [n for n, zi in zip(names, z) if zi > 0]
            neg = [n for n, zi in zip(names, z) if zi < 0]
            if not (len(neutral) == 1 and len(pos) == 1 and len(neg) == 1):
                raise PitzerDatabaseError(
                    f"The parameter zeta for {', '.join(names)} requires one "
                    "neutral species, one cation, and one anion.")
            zeta_triplets.append((neutral[0], pos[0], neg[0]))
            superblocks["zeta"].append(
                [names_line(neutral[0], pos[0], neg[0]), "  zeta:"] + a_lines(row))

        elif p == "mu":
            if any([zi != 0 for zi in z]):
                raise PitzerDatabaseError(
                    f"The parameter mu for {', '.join(names)} requires neutral species.")
            unique = list(dict.fromkeys(names))
            if len(unique) == 1:
                nn_blocks.setdefault(names[0], {})["mu"] = row
            elif len(unique) == 2:
                # the species that appears twice goes first
                doubled = [n for n in unique if names.count(n) == 2][0]
                other = [n for n in unique if n != doubled][0]
                superblocks["mu"].append(
                    [names_line(doubled, doubled, other), "  mu:"] + a_lines(row))
            else:
                raise PitzerDatabaseError(
                    f"The parameter mu for {', '.join(names)} is not supported "
                    "by EQ3/6: at most two different neutral species may "
                    "appear in a mu triplet.")

    # EQPT requires supporting pair data for every triplet parameter: a psi
    # (cc'a or aa'c) triplet needs cation-anion blocks for both like-charged
    # ions with the counter ion, a zeta (nca) triplet needs lambda blocks for
    # the neutral species with both ions, and a lambda(nn') pair needs nn
    # blocks for both neutral species. Pairs that have no data in the
    # database are given blocks with zero-valued parameters, which is the
    # convention used in EQ3/6 Pitzer data files (e.g., lambda(CO2, H+) = 0).
    added_zero_blocks = []
    for like1, like2, opp in psi_triplets:
        for like in (like1, like2):
            pair = (like, opp) if species_charges[like] > 0 else (opp, like)
            if pair not in ca_blocks:
                ca_blocks[pair] = {}
                added_zero_blocks.append(f"beta/Cphi {pair[0]}, {pair[1]}")
    for neutral, cation, anion in zeta_triplets:
        for ion in (cation, anion):
            if (neutral, ion) not in lambda_pairs:
                lambda_pairs.add((neutral, ion))
                superblocks["lambda_ion"].append(
                    [names_line(neutral, ion), "  lambda:"] + zeros())
                added_zero_blocks.append(f"lambda {neutral}, {ion}")
    for n1, n2 in nn_prime_pairs:
        for n in (n1, n2):
            if n not in nn_blocks:
                nn_blocks[n] = {}
                added_zero_blocks.append(f"lambda/mu {n}")
    if verbose > 1 and len(added_zero_blocks) > 0:
        print("Added zero-valued Pitzer blocks required by EQPT for:", "; ".join(added_zero_blocks))

    # assemble cation-anion blocks
    for (cation, anion), params in ca_blocks.items():
        lines = [names_line(cation, anion)]
        alpha1 = params["beta1"]["alpha"] if "beta1" in params else np.nan
        alpha2 = params["beta2"]["alpha"] if "beta2" in params else np.nan
        if not (_is_blank(alpha1) and _is_blank(alpha2)):
            # EQPT requires both alpha lines if either is given
            if _is_blank(alpha1):
                alpha1 = 2.0 if (abs(species_charges[cation]) == 1 or abs(species_charges[anion]) == 1) else 1.4
            if _is_blank(alpha2):
                alpha2 = 12.0
            lines.append(f"  alpha(1) = {_fmt_coeff(alpha1)}")
            lines.append(f"  alpha(2) = {_fmt_coeff(alpha2)}")
        for p, label in zip(PITZER_CA_PARAMS, ["beta(0)", "beta(1)", "beta(2)", "Cphi"]):
            lines.append(f"  {label}:")
            if p in params:
                lines += a_lines(params[p])
            else:
                lines += zeros()
        superblocks["ca"].append(lines)

    # assemble single-neutral blocks
    for neutral, params in nn_blocks.items():
        lines = [names_line(neutral)]
        for p in ["lambda", "mu"]:
            lines.append(f"  {p}:")
            if p in params:
                lines += a_lines(params[p])
            else:
                lines += zeros()
        superblocks["nn"].append(lines)

    # write out the superblocks
    out = []
    n_blocks = 0
    for key in PITZER_SUPERBLOCK_ORDER:
        out.append(PITZER_SUPERBLOCK_HEADERS[key])
        out.append(DATA0_DELIMITER)
        for block in superblocks[key]:
            out += block
            out.append(DATA0_DELIMITER)
            n_blocks += 1

    if verbose > 0 and len(skipped_species) > 0:
        skipped_sorted = sorted(skipped_species)
        msg = ("Pitzer parameters for the following species were not "
               "included in the data0 file because the species are not "
               "in the thermodynamic database: ")
        if len(skipped_sorted) > 15:
            msg += ", ".join(skipped_sorted[:15]) + f", and {len(skipped_sorted)-15} more"
        else:
            msg += ", ".join(skipped_sorted)
        print(msg)

    return {"text": "\n".join(out),
            "n_terms": n_terms,
            "skipped_species": sorted(skipped_species),
            "n_blocks": n_blocks}


def pitzer_config_lines(n_terms=4):
    """
    Lines that EQPT reads from the title of a Pitzer data0 file to determine
    how the Pitzer data blocks are organized. These lines must not begin
    with an asterisk because EQPT strips comment lines before reading the
    title.
    """
    return [
        "BEGIN CONFIGURATION DATA BLOCK",
        "Pitzer data parameters:",
        "   PITZER DATA BLOCK ORG.= NEW",
        "   PITZER TEMP FUNCTION= LIVERMORE",
        f"   NO. OF PITZER TEMP FUNC TERMS= {int(n_terms)}",
        "END CONFIGURATION DATA",
    ]


def data0_is_pitzer(data0_text):
    """
    Determine whether a data0 file (given as its text) is formatted for
    Pitzer's equations rather than the B-dot/Davies (SEDH) equations.
    """
    lines = data0_text.split("\n")
    for line in lines:
        s = line.strip().lower()
        if s.startswith("ca combinations:") or s.startswith("single-salt parameters") \
                or s.startswith("debye huckel aphi"):
            return True
        if s.startswith("bdot parameters"):
            return False
    return False


def data1_is_pitzer(data1_bytes):
    """
    Determine whether a data1 file (given as bytes) was produced from a
    Pitzer data0 file. EQPT writes the keystring "Pitzer" or "SEDH" in the
    header record of the data1 file.
    """
    head = data1_bytes[:256]
    if b"Pitzer" in head:
        return True
    return False


def eq36_to_worm_name(name):
    """
    Convert an EQ3/6-style species name (e.g., "Ca++", "SO4--", "CO2(aq)",
    "NpO2(CO3)3(5-)") to the WORM database style ("Ca+2", "SO4-2", "CO2",
    "NpO2(CO3)3-5"). Names that do not end in a charge designation are
    returned unchanged apart from the removal of an "(aq)" suffix.
    """
    name = name.strip()
    if name.endswith("(aq)"):
        name = name[:-4]
    m = re.fullmatch(r"(.*?)\((\d+)([+-])\)", name)
    if m:
        return f"{m.group(1)}{m.group(3)}{m.group(2)}"
    m = re.fullmatch(r"(.*?)(\+{2,}|-{2,})", name)
    if m:
        sign = m.group(2)[0]
        return f"{m.group(1)}{sign}{len(m.group(2))}"
    return name


def pitzer_data0_to_csv(data0_filename, csv_filename=None, worm_names=False,
                        name_map=None, verbose=1):
    """
    Extract the Pitzer interaction parameters from an EQ3/6 Pitzer data0
    file (with the "NEW" Pitzer data block organization, e.g. data0.ypf)
    into the aqequil Pitzer parameter CSV format.

    Parameters
    ----------
    data0_filename : str
        Path of the Pitzer data0 file.
    csv_filename : str, optional
        If given, the extracted parameters are written to this CSV file.
    worm_names : bool, default False
        Convert EQ3/6-style species names to WORM-style names (e.g.,
        "Ca++" to "Ca+2" and "CO2(aq)" to "CO2")? See `eq36_to_worm_name`.
    name_map : dict, optional
        Dictionary of {name in data0: desired name} applied after
        `worm_names` conversion.
    verbose : int
        Verbosity level.

    Returns
    -------
    pd.DataFrame
        Pitzer parameters in the long CSV format.
    """

    if name_map is None:
        name_map = {}

    with open(data0_filename, "r") as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f.readlines()]

    n_terms = 4
    for line in lines:
        m = re.search(r"NO\. OF PITZER TEMP FUNC TERMS\s*=\s*(\d+)", line)
        if m:
            n_terms = int(m.group(1))
            break
        if line.strip().startswith("+---"):
            break
    for line in lines:
        if "PITZER DATA BLOCK ORG." in line and "NEW" not in line.upper():
            raise PitzerDatabaseError(
                f"{data0_filename} uses the old ('single-salt') Pitzer data "
                "block organization, which is not supported by this converter.")

    label_to_param = {"beta(0)": "beta0", "beta(1)": "beta1", "beta(2)": "beta2",
                      "cphi": "cphi", "theta": "theta", "lambda": "lambda",
                      "psi": "psi", "zeta": "zeta", "mu": "mu"}

    # locate the Pitzer superblocks
    header_prefixes = [h[:24] for h in PITZER_SUPERBLOCK_HEADERS.values()]
    start = None
    end = None
    for i, line in enumerate(lines):
        if start is None and line[:24] == PITZER_SUPERBLOCK_HEADERS["ca"][:24]:
            start = i
        elif start is not None and line.strip() == "elements":
            end = i
            break
    if start is None or end is None:
        raise PitzerDatabaseError(
            f"Could not find Pitzer data blocks in {data0_filename}.")

    def rename(n):
        if worm_names:
            n = eq36_to_worm_name(n)
        return name_map.get(n, n)

    rows = []
    block = []

    def flush(block):
        if len(block) == 0:
            return
        names_line = block[0]
        names = [names_line[0:24].strip(), names_line[26:50].strip(), names_line[52:76].strip()]
        names = [n for n in names if n != ""]
        # handle names lines that are not aligned on 24-character fields
        if len(names) <= 1 and len(names_line.split()) > 1:
            names = names_line.split()
        names = [rename(n) for n in names]
        params = {}
        alphas = {}
        current = None
        notes = []
        for l in block[1:]:
            s = l.strip()
            if s.startswith("*"):
                notes.append(s.lstrip("*").strip())
                continue
            m = re.fullmatch(r"alpha\((\d)\)\s*=\s*(\S+)", s)
            if m:
                alphas[int(m.group(1))] = float(m.group(2))
                continue
            m = re.fullmatch(r"([A-Za-z()0-9]+)\s*:", s)
            if m:
                label = m.group(1)
                current = label_to_param.get(label.lower().replace("c(phi)", "cphi"), None)
                if current is None:
                    current = label_to_param.get(label, None)
                if current is None:
                    raise PitzerDatabaseError(
                        f"Unrecognized parameter label '{label}' in the block "
                        f"for {', '.join(names)} of {data0_filename}.")
                params[current] = {}
                continue
            m = re.fullmatch(r"a(\d+)\s*=\s*(\S+)", s)
            if m and current is not None:
                params[current][f"a{m.group(1)}"] = float(m.group(2))
                continue
        for p, coeffs in params.items():
            row = {"species1": names[0] if len(names) > 0 else "",
                   "species2": names[1] if len(names) > 1 else "",
                   "species3": names[2] if len(names) > 2 else "",
                   "param": p}
            if p == "mu" and len(names) == 1:
                row["species2"] = names[0]
                row["species3"] = names[0]
            if p == "lambda" and len(names) == 1:
                row["species2"] = names[0]
            for j in range(1, n_terms + 1):
                row[f"a{j}"] = coeffs.get(f"a{j}", 0.0)
            if p == "beta1" and 1 in alphas:
                row["alpha"] = alphas[1]
            elif p == "beta2" and 2 in alphas:
                row["alpha"] = alphas[2]
            else:
                row["alpha"] = np.nan
            row["ref1"] = ""
            row["ref2"] = ""
            row["date"] = ""
            row["note"] = " ".join(notes)
            rows.append(row)

    for line in lines[start:end]:
        if line[:24] in header_prefixes:
            flush(block)
            block = []
            continue
        if line.startswith("+---"):
            flush(block)
            block = []
            continue
        if line.strip() == "":
            continue
        if line.startswith("*") and len(block) == 0:
            # a comment between blocks, not part of a block
            continue
        block.append(line)
    flush(block)

    cols = ["species1", "species2", "species3", "param"] + \
           [f"a{j}" for j in range(1, n_terms + 1)] + \
           ["alpha", "ref1", "ref2", "date", "note"]
    df = pd.DataFrame(rows, columns=cols)

    if verbose > 0:
        print(f"Extracted {df.shape[0]} Pitzer parameters from {data0_filename}.")

    if csv_filename is not None:
        df.to_csv(csv_filename, index=False)
        if verbose > 0:
            print(f"Wrote {csv_filename}.")

    return df
