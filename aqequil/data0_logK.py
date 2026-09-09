"""
Convert the species blocks of an EQ3/6 data0 file into an aqequil logK CSV.

The logK CSV format (see `wrm_data_logK.csv` in the WORM database) tabulates
the log K of a species' dissociation reaction at up to eight temperatures and
pressures. This module reads the auxiliary basis species, aqueous species,
solids, liquids and gases of a data0 file (e.g., `data0.ypf`), converts the
species names to WORM conventions, rewrites the dissociation reactions in
terms of the strict basis species of a WORM-style CSV database, and writes
the result as a logK CSV that can be loaded with
`aqequil.AqEquil(logK=...)`.

The main use is to pair the log K values of a Pitzer data0 file with the
Pitzer parameters extracted from the same file by
`aqequil.pitzer_data0_to_csv`, so that the aqueous complexes, mineral
solubilities and interaction parameters all come from one internally
consistent model.

Reactions are converted by linear algebra: every species is first expressed
as a formation vector over the data0 file's own basis, and those vectors are
then re-expressed over the target basis. Log K grids are combined the same
way, so a species whose data0 reaction is written in terms of, e.g., Cr+++
gets a reaction written in terms of CrO4-2, O2(g), H+ and H2O when the
target database uses CrO4-2 as its chromium basis species. Grid points that
the data0 file marks as "no data" (500) are dropped, so a species with data
only at 25 °C gets a one-point grid.
"""

import os
import re
import datetime

import numpy as np
import pandas as pd

from .pitzer import eq36_to_worm_name


# Elements whose oxidation state is taken as fixed when building the
# formula_ox column.
FIXED_OXIDATION_STATES = {
    "H": 1, "O": -2, "F": -1,
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Fr": 1,
    "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Ra": 2,
    "B": 3, "Al": 3, "Si": 4, "Zr": 4, "Hf": 4, "Th": 4,
    "Zn": 2, "Cd": 2, "Ni": 2, "Ga": 3, "In": 3, "Sc": 3, "Y": 3, "La": 3,
    "Gd": 3, "Nd": 3, "Cm": 3,
}

# Default oxidation states for elements that can take several; if the
# defaults do not balance the charge and exactly one such element is present,
# its oxidation state is inferred from the charge.
DEFAULT_VARIABLE_OXIDATION_STATES = {
    "C": 4, "N": 5, "P": 5, "S": 6, "Cl": -1, "Br": -1, "I": -1,
    "Fe": 2, "Mn": 2, "Cr": 3, "Co": 2, "Cu": 2, "Mo": 6, "W": 6, "Tc": 7,
    "V": 5, "Ti": 4, "U": 6, "Np": 5, "Pu": 4, "Am": 3, "Sn": 2, "Pb": 2,
    "Hg": 2, "Ag": 1, "Au": 1, "As": 5, "Sb": 3, "Se": 6, "Ce": 3, "Eu": 3,
    "Sm": 3, "Ru": 3, "Rh": 3, "Pd": 2, "Pt": 2, "Tl": 1, "Bi": 3,
}

SECTION_NAMES = ["basis species", "auxiliary basis species", "aqueous species",
                 "solids", "liquids", "gases", "solid solutions", "references"]

SECTION_STATE = {"basis species": "aq", "auxiliary basis species": "aq",
                 "aqueous species": "aq", "solids": "cr", "liquids": "liq",
                 "gases": "gas"}

NO_DATA = 500.0

# Cation prefixes used by the WORM database for exchangeable-cation clays
# and zeolites, e.g. data0 "Montmorillonite-Ca" is WORM "Ca-montmorillonite".
CATION_SUFFIXES = ["Na", "K", "Ca", "Mg", "H", "Cs", "Sr", "NH4", "Fe", "Li", "Ba"]


class Data0ParseError(ValueError):
    pass


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _read_numbers(lines, i):
    """Read consecutive lines of numbers starting at index i."""
    vals = []
    while i < len(lines):
        s = lines[i].strip()
        if s == "" or s.startswith("*") or s.startswith("+"):
            break
        toks = s.split()
        if not all(re.fullmatch(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?", t) or
                   re.fullmatch(r"(?i)no_?data", t) for t in toks):
            break
        for t in toks:
            try:
                vals.append(float(t))
            except ValueError:
                # "No_Data" placeholders in some data0 files
                vals.append(np.nan)
        i += 1
    return vals, i


def _split_name_line(line):
    """Name (and optional formula) from the first line of a species block.
    The name occupies the first 24 characters and may contain spaces
    (e.g. "Soda Niter              NaNO3")."""
    if len(line) > 24 and line[:24].strip() != "" and line[24:].strip() != "":
        return line[:24].strip(), line[24:].strip()
    return line.strip(), None


def _split_reaction_tokens(text):
    """Tokens of a reaction block: coefficients and names are separated by at
    least two spaces, names may contain single spaces."""
    return [t for t in re.split(r"\s{2,}", text.strip()) if t != ""]


def _parse_block(block, section):
    """Parse one species block (list of lines) into a dict."""
    name, formula = _split_name_line(block[0])
    entry = {"name": name, "formula": formula, "section": section,
             "state": SECTION_STATE.get(section, "aq"), "charge": 0.0,
             "elements": {}, "reaction": {}, "logK": None, "V": None,
             "azero": None, "sp_type": None, "source": []}

    i = 1
    while i < len(block):
        line = block[i]
        s = line.strip()
        m = re.match(r"charge\s*=\s*(\S+)", s)
        if m:
            entry["charge"] = float(m.group(1))
            i += 1
            continue
        m = re.match(r"sp\.type\s*=\s*(\S+)", s)
        if m:
            entry["sp_type"] = m.group(1)
            i += 1
            continue
        m = re.match(r"V0PrTr\s*=\s*(\S+)", s)
        if m:
            try:
                v = float(m.group(1))
                entry["V"] = v if (v != 0 and abs(v - NO_DATA) > 1e-9) else None
            except ValueError:
                pass
            i += 1
            continue
        m = re.match(r"\*\s*DHazero\s*=\s*(\S+)", s)
        if m:
            entry["azero"] = float(m.group(1))
            i += 1
            continue
        m = re.match(r"(\d+)\s+element\(s\):", s)
        if m:
            n = int(m.group(1))
            i += 1
            toks = []
            while len(toks) < 2 * n and i < len(block):
                toks += block[i].split()
                i += 1
            for k in range(0, 2 * n, 2):
                entry["elements"][toks[k + 1]] = float(toks[k])
            continue
        m = re.match(r"(\d+)\s+species in reaction:", s)
        if m:
            n = int(m.group(1))
            i += 1
            toks = []
            while len(toks) < 2 * n and i < len(block):
                toks += _split_reaction_tokens(block[i])
                i += 1
            for k in range(0, 2 * n, 2):
                entry["reaction"][toks[k + 1]] = float(toks[k])
            continue
        if s.startswith("**** logK grid") or s.lower().startswith("**** log k grid"):
            vals, i = _read_numbers(block, i + 1)
            # keep only the first grid; a later "logK grid" header may
            # introduce commented-out alternate data
            if entry["logK"] is None or (len(entry["logK"]) == 0 and len(vals) > 0):
                entry["logK"] = vals
            continue
        m = re.match(r"\*\s*Source:\s*(.*)", s)
        if m:
            entry["source"].append(m.group(1).strip())
            i += 1
            continue
        i += 1
    return entry


def parse_data0_species(data0_filename):
    """
    Parse the temperature grid and the species blocks of an EQ3/6 data0 file.

    Returns
    -------
    dict with keys "temperatures", "pressures", "basis" (list of entries),
    "species" (list of entries for auxiliary basis, aqueous, solid, liquid
    and gas species, in file order).
    """
    with open(data0_filename, "r") as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f.readlines()]

    temps, press = [], []
    for i, line in enumerate(lines):
        if line.strip() == "temperatures":
            temps, _ = _read_numbers(lines, i + 1)
        elif line.strip() == "pressures":
            press, _ = _read_numbers(lines, i + 1)
        if len(temps) > 0 and len(press) > 0:
            break
    if len(temps) == 0:
        raise Data0ParseError(f"Could not find the temperature grid in {data0_filename}")

    # locate the start of the basis species section
    start = None
    for i, line in enumerate(lines):
        if line.strip() == "basis species":
            start = i
            break
    if start is None:
        raise Data0ParseError(f"Could not find 'basis species' in {data0_filename}")

    basis, species = [], []
    section = None
    block = []

    def flush(block):
        if len(block) == 0:
            return
        head = block[0].strip()
        if head in SECTION_NAMES:
            return head
        if head == "stop." or section is None or section in ["solid solutions", "references"]:
            return None
        e = _parse_block(block, section)
        if section == "basis species":
            basis.append(e)
        else:
            species.append(e)
        return None

    for line in lines[start:]:
        if line.startswith("+---"):
            new_section = flush(block)
            if new_section is not None:
                section = new_section
            block = []
            continue
        if line.strip() == "stop.":
            flush(block)
            block = []
            break
        if line.strip() == "":
            continue
        if len(block) == 0 and line.strip() in SECTION_NAMES:
            section = line.strip()
            continue
        block.append(line)
    flush(block)

    return {"temperatures": temps, "pressures": press, "basis": basis, "species": species}


# ---------------------------------------------------------------------------
# Naming helpers
# ---------------------------------------------------------------------------

def data0_to_worm_name(name, state, worm_names=None):
    """
    Convert a data0 species name to WORM conventions. Aqueous species use
    `eq36_to_worm_name` ("Ca++" -> "Ca+2", "CO2(aq)" -> "CO2"). Minerals are
    lower-cased unless the name is a chemical formula, ":" in hydrate names
    becomes "*", and "Mineral-Ca" style names become "Ca-mineral". If
    `worm_names` (an iterable of names in the target database) contains a
    case-insensitive match, that spelling is used.
    """
    name = name.strip()
    if state in ["aq", "gas"]:
        return eq36_to_worm_name(name)

    lookup = {}
    if worm_names is not None:
        lookup = {n.lower(): n for n in worm_names}

    n = re.sub(r"\s+", "_", name).replace(":", "*")
    if n.lower() in lookup:
        return lookup[n.lower()]
    m = re.fullmatch(r"([A-Za-z]+)-(" + "|".join(CATION_SUFFIXES) + r")", n)
    if m:
        swapped = f"{m.group(2)}-{m.group(1).lower()}"
        if swapped.lower() in lookup:
            return lookup[swapped.lower()]
        return swapped
    if re.search(r"\d", n) or "*" in n or "(" in n or "." in n:
        # chemical formulas (hydrates, "Fe(OH)3", "SiO2(am)") keep their case
        return n
    return n.lower()


def _formula_from_elements(elements):
    """A parseable chemical formula built from an element composition."""
    parts = []
    order = [el for el in elements if el not in ["O", "H"]] + \
            [el for el in ["O", "H"] if el in elements]
    for el in order:
        n = elements[el]
        if abs(n) < 1e-9:
            continue
        n_str = "" if abs(n - 1) < 1e-9 else (f"{n:g}")
        parts.append(f"{el}{n_str}")
    return "".join(parts)


def _solid_formula(entry, elements):
    """Formula of a solid: the formula given in the data0 block, or the
    species name when it is itself a formula (e.g. "Fe(OH)3", "KB5O8:4H2O"),
    provided it reproduces the element composition; otherwise a formula
    built from the element composition."""
    try:
        from pychnosz import makeup
    except Exception:
        makeup = None
    candidates = [entry["formula"], entry["name"]]
    for cand in candidates:
        if not cand:
            continue
        f = cand.replace(":", "*")
        if not re.fullmatch(r"[A-Za-z0-9().*]+", f):
            continue
        if makeup is None:
            return f
        try:
            mk = {k: float(v) for k, v in dict(makeup(f)).items() if k != "Z"}
        except Exception:
            continue
        ref = {k: float(v) for k, v in elements.items() if abs(v) > 1e-12}
        if set(mk) == set(ref) and all(abs(mk[k] - ref[k]) < 1e-6 for k in ref):
            return f
    return _formula_from_elements(elements)


def formula_ox_from_elements(elements, charge):
    """
    Build a WORM-style formula_ox string ("Fe+3 4O-2 ..."). Elements with a
    fixed oxidation state use it; if exactly one element with a variable
    oxidation state is present, its state is inferred from the charge.
    Returns "" when the oxidation states cannot be assigned.
    """
    els = {el: n for el, n in elements.items() if abs(n) > 1e-12 and el != "Z"}
    ox = {}
    variable = []
    for el in els:
        if el in FIXED_OXIDATION_STATES:
            ox[el] = float(FIXED_OXIDATION_STATES[el])
        elif el in DEFAULT_VARIABLE_OXIDATION_STATES:
            ox[el] = float(DEFAULT_VARIABLE_OXIDATION_STATES[el])
            variable.append(el)
        else:
            return ""
    # a single element in a neutral species (H2, O2, native metals) has
    # oxidation state zero, written without a charge as in the WORM database
    if len(els) == 1 and abs(charge) < 1e-9:
        el, n = next(iter(els.items()))
        n_str = "" if abs(n - 1) < 1e-9 else f"{n:g}"
        return f"{n_str}{el}"

    def _part(el, n, o):
        n_str = "" if abs(n - 1) < 1e-9 else f"{n:g}"
        sign = "+" if o > 0 else "-"
        mag = "" if abs(o) == 1 else f"{abs(o):g}"
        return f"{n_str}{el}{sign}{mag}"

    split = {}  # element -> [(count, state), (count, state)] for mixed valence
    total = sum(els[el] * ox[el] for el in els)
    if abs(total - charge) > 1e-6:
        # let one variable element absorb the difference; metals are tried
        # before non-metals (U(SO4)2 is U+4 with S+6, not U+6 with S+5)
        nonmetals = ["C", "N", "P", "S", "Se", "As", "Cl", "Br", "I"]
        order = [el for el in variable if el not in nonmetals] +                 [el for el in variable if el in nonmetals]
        for el in order:
            trial = ox[el] + (charge - total) / els[el]
            if abs(trial - round(trial)) < 1e-6:
                ox[el] = float(round(trial))
                break
            # mixed valence (magnetite: Fe+2 + 2 Fe+3): split an integer
            # number of atoms between two adjacent integer states
            n = els[el]
            if abs(n - round(n)) < 1e-6 and el not in nonmetals:
                need = ox[el] * n + (charge - total)
                lo = np.floor(need / n)
                n_hi = need - lo * n
                n_lo = n - n_hi
                if abs(n_hi - round(n_hi)) < 1e-6 and n_hi > 0 and n_lo > 0:
                    split[el] = [(n_lo, lo), (n_hi, lo + 1)]
                    break
        else:
            return ""

    parts = []
    for el, n in els.items():
        if el in split:
            for cnt, o in split[el]:
                parts.append(_part(el, cnt, o))
        else:
            parts.append(_part(el, n, ox[el]))
    return " ".join(parts)


def _default_worm_db_path():
    return os.path.join(os.path.dirname(__file__), "databases", "wrm_data_latest.csv")


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def data0_to_logK_csv(data0_filename,
                      csv_filename=None,
                      worm_db=None,
                      basis_species=None,
                      name_map=None,
                      dataset_name=None,
                      ref=None,
                      include_states=("aq", "cr", "gas", "liq"),
                      verbose=1):
    """
    Convert the species blocks of an EQ3/6 data0 file into an aqequil logK
    CSV with WORM species names and dissociation reactions written in terms
    of the strict basis species of a WORM-style CSV database.

    Parameters
    ----------
    data0_filename : str
        Path of the data0 file (e.g., "data0.ypf").
    csv_filename : str, optional
        If given, the logK database is written to this CSV file.
    worm_db : str or pd.DataFrame, optional
        The WORM-style CSV thermodynamic database whose strict basis species
        and species names are the target. Defaults to the bundled
        wrm_data_latest.csv.
    basis_species : list of str, optional
        Names of the target strict basis species. Defaults to H2O, H+, O2(g)
        and the species tagged "basis" in `worm_db`. Each must exist in the
        data0 file (as a basis, auxiliary basis or aqueous species) to be
        used; species of the data0 file containing an element with no
        target basis species are skipped and reported.
    name_map : dict, optional
        {data0 name: desired name} overrides applied after the automatic
        conversion to WORM naming conventions.
    dataset_name : str, optional
        Short name used for the category_1 column (e.g., "ypf" gives
        "ypf_aq", "ypf_cr", "ypf_gas"). Defaults to the data0 suffix.
    ref : str, optional
        Text prepended to the ref2 column. The data0 "Source:" comment of
        each species follows it.
    include_states : tuple of str
        Which states to convert. Default: aqueous, solids, gases, liquids.
    verbose : int
        Verbosity level.

    Returns
    -------
    dict
        {"logK": DataFrame, "report": dict}. The report lists skipped
        species with reasons, renamed species, and the basis mapping.
    """
    if name_map is None:
        name_map = {}
    if dataset_name is None:
        base = os.path.basename(data0_filename)
        dataset_name = re.sub(r"^data0\.", "", base) or base
    if ref is None:
        ref = ""

    # ---- target database
    if worm_db is None:
        worm_db = _default_worm_db_path()
    if isinstance(worm_db, pd.DataFrame):
        wdf = worm_db.astype(str)
        worm_db_name = "the supplied WORM database"
    else:
        wdf = pd.read_csv(worm_db, dtype=str)
        worm_db_name = os.path.basename(worm_db)
    worm_names_all = list(wdf["name"])
    worm_state = dict(zip(wdf["name"], wdf["state"]))
    worm_aq_names = set(wdf.loc[wdf["state"] == "aq", "name"])
    if basis_species is None:
        basis_species = ["H2O", "H+", "O2(g)"] + \
            [n for n in wdf.loc[wdf["tag"] == "basis", "name"] if n not in ["H2O", "H+", "O2(g)"]]

    ds = parse_data0_species(data0_filename)
    temps = ds["temperatures"]
    if len(temps) > 8:
        raise Data0ParseError("The aqequil logK CSV format supports at most eight "
                              f"temperatures, but {data0_filename} has {len(temps)}.")
    nT = len(temps)
    report = {"skipped": {}, "renamed": {}, "basis_old": [], "basis_new": [],
              "temperatures": temps, "collisions": {}}

    old_basis = [b["name"] for b in ds["basis"]]
    # a gas that is also a basis species (O2(g)) is the same species; drop the
    # gas block so that names stay unique
    ds["species"] = [e for e in ds["species"] if e["name"] not in old_basis]
    entries = {e["name"]: e for e in ds["species"]}
    n_old = len(old_basis)
    old_index = {n: i for i, n in enumerate(old_basis)}

    charge = {b["name"]: b["charge"] for b in ds["basis"]}
    comp = {b["name"]: dict(b["elements"]) for b in ds["basis"]}
    for e in ds["species"]:
        charge[e["name"]] = e["charge"] if e["state"] == "aq" else 0.0
        comp[e["name"]] = dict(e["elements"])

    # ---- WORM names
    worm = {}
    for n in old_basis:
        worm[n] = name_map.get(n, eq36_to_worm_name(n))
    for e in ds["species"]:
        w = data0_to_worm_name(e["name"], e["state"], worm_names_all)
        worm[e["name"]] = name_map.get(e["name"], w)

    # resolve name collisions: aqequil matches logK entries to database
    # species by name alone, so a solid whose converted name equals an
    # aqueous species name (in this file or in the target database) gets a
    # "(cr)" suffix.
    aq_names_here = set(worm[n] for n in old_basis) | \
        set(worm[e["name"]] for e in ds["species"] if e["state"] == "aq")
    for e in ds["species"]:
        if e["state"] != "aq" and e["name"] not in name_map:
            w = worm[e["name"]]
            if w in aq_names_here or (w in worm_aq_names and worm_state.get(w) == "aq"):
                new = w + "(cr)" if e["state"] == "cr" else w + f"({e['state']})"
                report["collisions"][w] = new
                worm[e["name"]] = new
    seen = {}
    for n, w in worm.items():
        seen.setdefault(w, []).append(n)
    dupes = {w: ns for w, ns in seen.items() if len(ns) > 1}
    if len(dupes) > 0:
        raise Data0ParseError(f"Name conversion produced duplicate names: {dupes}. "
                              "Use name_map to resolve them.")
    for n, w in worm.items():
        if n != w:
            report["renamed"][n] = w

    # ---- formation vectors and formation log K grids over the old basis
    form_vec = {}
    form_logK = {}
    for b in old_basis:
        v = np.zeros(n_old)
        v[old_index[b]] = 1.0
        form_vec[b] = v
        form_logK[b] = np.zeros(nT)

    def grid_of(e):
        g = e["logK"]
        if g is None or len(g) == 0:
            return None
        arr = np.full(nT, np.nan)
        for k in range(min(nT, len(g))):
            arr[k] = np.nan if (np.isnan(g[k]) or abs(g[k] - NO_DATA) < 1e-9) else g[k]
        return arr

    def normalized_reaction(e):
        """{species: coef} of the dissociation of one mole of the species
        (the species itself excluded) and the matching log K grid."""
        c0 = e["reaction"].get(e["name"], -1.0)
        if c0 == 0:
            c0 = -1.0
        scale = 1.0 / (-c0)
        coefs = {sp: coef * scale for sp, coef in e["reaction"].items() if sp != e["name"]}
        g = grid_of(e)
        return coefs, (None if g is None else g * scale)

    def resolve(name, stack=()):
        if name in form_vec:
            return True
        if name not in entries:
            return False
        if name in stack:
            raise Data0ParseError(f"Circular reaction definitions involving {name}")
        e = entries[name]
        if len(e["reaction"]) == 0 or e["logK"] is None:
            report["skipped"][worm.get(name, name)] = "no dissociation reaction or log K grid"
            return False
        coefs, grid = normalized_reaction(e)
        v = np.zeros(n_old)
        lk = -grid  # formation log K = -(dissociation log K)
        for sp, coef in coefs.items():
            if not resolve(sp, stack + (name,)):
                report["skipped"][worm.get(name, name)] = f"reaction involves unknown species {sp}"
                return False
            v += coef * form_vec[sp]
            lk = lk + coef * form_logK[sp]
        form_vec[name] = v
        form_logK[name] = lk
        return True

    for name in list(entries.keys()):
        resolve(name)

    # ---- new basis: target basis species present in the data0 file
    inv_worm = {w: n for n, w in worm.items()}
    new_basis = []
    for bs in basis_species:
        n = inv_worm.get(bs)
        if n is not None and n in form_vec and (n in old_basis or entries[n]["state"] == "aq"):
            new_basis.append(n)
    new_basis_set = set(new_basis)
    M = np.array([form_vec[n] for n in new_basis])  # rows: new basis in old basis
    if np.linalg.matrix_rank(M) < len(new_basis):
        raise Data0ParseError("The target basis species are not linearly independent "
                              f"in terms of the data0 basis: {[worm[n] for n in new_basis]}")
    report["basis_old"] = [worm[n] for n in old_basis]
    report["basis_new"] = [worm[n] for n in new_basis]

    def to_new_basis(name):
        c, _, _, _ = np.linalg.lstsq(M.T, form_vec[name], rcond=None)
        if np.linalg.norm(M.T @ c - form_vec[name]) > 1e-8:
            return None
        c[np.abs(c) < 1e-9] = 0.0
        return c

    report["basis_unmapped"] = [worm[n] for n in old_basis
                                if n not in new_basis_set and to_new_basis(n) is None]

    def strict_basis_reaction(name):
        """Dissociation reaction and log K grid of a species written in terms
        of the new strict basis (None if not expressible)."""
        c = to_new_basis(name)
        if c is None:
            return None, None
        coefs = {nb: c[j] for j, nb in enumerate(new_basis) if c[j] != 0.0}
        logK_form_new = form_logK[name] - sum(c[j] * form_logK[nb] for j, nb in enumerate(new_basis))
        return coefs, -logK_form_new

    # ---- species that can serve as auxiliary basis species in the target:
    # the data0 auxiliary basis species, plus data0 basis species that are
    # not basis species in the target (e.g. Cr+3 when the target uses CrO4-2)
    aux_available = set()
    for n in old_basis:
        if n not in new_basis_set:
            coefs, grid = strict_basis_reaction(n)
            if coefs is not None and grid is not None and (~np.isnan(grid)).any():
                aux_available.add(n)
    for e in ds["species"]:
        n = e["name"]
        if (e["section"] == "auxiliary basis species" and e["state"] == "aq"
                and n in form_vec and n not in new_basis_set and to_new_basis(n) is not None):
            g = grid_of(e)
            if g is not None and (~np.isnan(g)).any():
                aux_available.add(n)
    available = new_basis_set | aux_available

    def fmt_coef(x):
        # four decimals unless more are needed to keep the reaction balanced
        if abs(x - round(x, 4)) < 1e-9:
            return f"{x:.4f}"
        return f"{x:.10f}".rstrip("0")

    # ---- rows
    today = datetime.date.today().strftime("%m/%d/%Y")
    rows = []

    def balanced_at_four_decimals(name, coefs, elements):
        """aqequil writes reaction coefficients with four decimals and EQPT
        rejects reactions with element imbalances around 1e-4, so check the
        balance of the rounded reaction against the element compositions."""
        tot = {el: -n for el, n in elements.items()}
        for sp, coef in coefs.items():
            for el, n in comp.get(sp, {}).items():
                tot[el] = tot.get(el, 0.0) + round(coef, 4) * n
        return all(abs(v) < 1e-6 for v in tot.values())

    def make_row(name, state, coefs, logK_diss, elements, z, formula, azero, source, V):
        ok = ~np.isnan(logK_diss)
        if not ok.any():
            report["skipped"][worm[name]] = "no log K data at any grid temperature"
            return
        if not balanced_at_four_decimals(name, coefs, elements):
            report["skipped"][worm[name]] = ("reaction cannot be balanced with four-decimal "
                                            "coefficients (EQPT would reject it)")
            return
        T_sp = [temps[k] for k in range(nT) if ok[k]]
        lk_sp = [float(logK_diss[k]) for k in range(nT) if ok[k]]
        dissrxn = f"-1.0000 {worm[name]}"
        for sp, coef in coefs.items():
            if coef != 0.0:
                dissrxn += f" {fmt_coef(coef)} {worm[sp]}"
        ref2 = "; ".join([t for t in [ref] + source if t])
        tag = "aux" if (state == "aq" and name in aux_available) else ""
        row = {"name": worm[name], "abbrv": "", "formula": formula, "state": state,
               "ref1": os.path.basename(data0_filename), "ref2": ref2, "date": today}
        for k in range(8):
            row[f"logK{k+1}"] = round(lk_sp[k], 5) if k < len(T_sp) else np.nan
        for k in range(8):
            row[f"T{k+1}"] = T_sp[k] if k < len(T_sp) else np.nan
        for k in range(8):
            row[f"P{k+1}"] = "psat" if k < len(T_sp) else np.nan
        row["azero"] = azero
        row["dissrxn"] = dissrxn
        row["tag"] = tag
        row["formula_ox"] = formula_ox_from_elements(elements, z)
        row["category_1"] = f"{dataset_name}_{state}"
        row["V"] = V if (state == "cr" and V is not None) else np.nan
        rows.append(row)

    # data0 basis species that become auxiliary species in the target
    for b in ds["basis"]:
        n = b["name"]
        if n in new_basis_set or n not in aux_available or "aq" not in include_states:
            continue
        coefs, grid = strict_basis_reaction(n)
        make_row(n, "aq", coefs, grid, dict(b["elements"]), b["charge"],
                 worm[n], 4.0 if b["charge"] != 0 else 3.0, [], None)

    for e in ds["species"]:
        name = e["name"]
        if name in new_basis_set or name not in form_vec:
            continue
        if e["state"] not in include_states:
            continue
        if worm[name] in ["O2(g)", "H2O", "H+"]:
            continue

        coefs, grid = normalized_reaction(e)
        if not all(sp in available for sp in coefs):
            # a participant is neither a target basis species nor an
            # available auxiliary species: rewrite in the strict basis
            coefs, grid = strict_basis_reaction(name)
            if coefs is None:
                bad = [el for el in comp.get(name, {})
                       if not any(el in comp.get(nb, {}) for nb in new_basis)]
                report["skipped"][worm[name]] = ("no target basis species for element(s) "
                                                f"{bad}" if bad else "cannot be expressed in the target basis")
                continue
            report.setdefault("rewritten", []).append(worm[name])

        elements = comp.get(name, {})
        z = charge.get(name, 0.0)
        state = e["state"]
        if state == "aq":
            formula = worm[name]
        elif state == "gas":
            formula = re.sub(r"\(g\)$", "", worm[name])
        else:
            formula = _solid_formula(e, elements)
        azero = e["azero"] if e["azero"] is not None else (4.0 if z != 0 else 3.0)
        make_row(name, state, coefs, grid, elements, z, formula, azero, e["source"], e["V"])

    logK_cols = ["name", "abbrv", "formula", "state", "ref1", "ref2", "date"] + \
                [f"logK{k}" for k in range(1, 9)] + [f"T{k}" for k in range(1, 9)] + \
                [f"P{k}" for k in range(1, 9)] + ["azero", "dissrxn", "tag", "formula_ox",
                                                  "category_1", "V"]
    logK_df = pd.DataFrame(rows, columns=logK_cols)

    if verbose > 0:
        n_full = int((~logK_df["T8"].isna()).sum()) if len(logK_df) else 0
        print(f"{os.path.basename(data0_filename)}: {len(logK_df)} logK entries "
              f"({(logK_df['state']=='aq').sum()} aqueous, {(logK_df['state']=='cr').sum()} "
              f"solids, {(logK_df['state']=='gas').sum()} gases, "
              f"{(logK_df['state']=='liq').sum()} liquids); {n_full} have data at all "
              f"{nT} grid temperatures ({', '.join(f'{T:g}' for T in temps)} C).")
        print("Basis species of", os.path.basename(data0_filename), "->", worm_db_name, "basis:")
        print("  ", ", ".join(report["basis_old"]))
        print("  ", ", ".join(report["basis_new"]))
        if len(report["basis_unmapped"]) > 0:
            print("Data0 basis species with no counterpart in the target basis:",
                  ", ".join(report["basis_unmapped"]))
        if len(report.get("rewritten", [])) > 0:
            print(f"{len(report['rewritten'])} reactions rewritten in terms of the target basis "
                  "species (their data0 reactions involve a basis species that is not a basis "
                  "species in the target):", ", ".join(report["rewritten"]))
        if len(report["collisions"]) > 0:
            print("Solids renamed to avoid a clash with an aqueous species name:",
                  ", ".join(f"{a} -> {b}" for a, b in report["collisions"].items()))
        if len(report["skipped"]) > 0:
            reasons = {}
            for k, v in report["skipped"].items():
                reasons.setdefault(v, []).append(k)
            for v, ks in reasons.items():
                print(f"Skipped {len(ks)} species ({v}):", ", ".join(ks))

    if csv_filename is not None:
        logK_df.to_csv(csv_filename, index=False)
        if verbose > 0:
            print("Wrote", csv_filename)

    return {"logK": logK_df, "report": report}
