"""
Convert Geochemist's Workbench (GWB) thermodynamic datasets (.tdat files)
that use the Harvie-Moller-Weare (Pitzer) activity model, such as
thermo_frezchem.tdat and thermo_coldchem.tdat, into aqequil's CSV formats:

- a logK database CSV in the format of wrm_data_logK.csv, with log K values
  tabulated at the dataset's principal temperatures and dissociation
  reactions rewritten in terms of the WORM basis species, and
- a Pitzer parameter CSV in the format described in `aqequil.pitzer`.

GWB datasets describe the temperature dependence of log K values and Pitzer
parameters with the six-term polynomial

    x(T) = a + b*(T - Tr) + c*(T^2 - Tr^2) + d*(1/T - 1/Tr) + e*(1/T^2 - 1/Tr^2)
           + f*ln(T/Tr)

(T in Kelvin, Tr = 298.15 K). EQ3/6's "Livermore" Pitzer temperature
function has the terms a1 + a2*(1/T - 1/Tr) + a3*ln(T/Tr) + a4*(T - Tr) +
a5*(T^2 - Tr^2), so every GWB term except e*(1/T^2 - 1/Tr^2) maps directly
(a1=a, a2=d, a3=f, a4=b, a5=c). Parameters with a non-zero e coefficient are
refit by least squares to the five Livermore terms over the temperature
range of the dataset; the maximum deviation of the refit is recorded in the
`note` column of the Pitzer CSV. Log K values are tabulated directly, so no
refitting is needed for them.
"""

import os
import re
import datetime

import numpy as np
import pandas as pd

from .pitzer import eq36_to_worm_name

TR_K = 298.15

# Oxidation states used to build the 'formula_ox' column for the simple
# electrolyte chemistry that low-temperature Pitzer datasets cover.
DEFAULT_OXIDATION_STATES = {
    "H": 1, "O": -2, "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1,
    "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "F": -1, "Cl": -1, "Br": -1, "I": -1,
    "S": 6, "C": 4, "N": 5, "B": 3, "Si": 4, "Al": 3, "P": 5,
}

GWB_SECTIONS = ["elements", "basis species", "redox couples", "aqueous species",
                "free electron", "minerals", "solid solutions", "gases", "oxides"]

GWB_VIRIAL_PARAMS = ["beta0", "beta1", "beta2", "cphi", "theta", "lambda", "psi",
                     "zeta", "mu"]


class GWBDatasetError(ValueError):
    pass


def gwb_poly(coeffs, T_degC):
    """
    Evaluate the GWB six-term temperature polynomial at one or more
    temperatures in degrees Celsius. `coeffs` is a dict with any of the keys
    a, b, c, d, e, f.
    """
    Tk = np.asarray(T_degC, dtype=float) + 273.15
    g = lambda k: float(coeffs.get(k, 0.0) or 0.0)
    return (g("a") + g("b") * (Tk - TR_K) + g("c") * (Tk**2 - TR_K**2)
            + g("d") * (1.0 / Tk - 1.0 / TR_K) + g("e") * (1.0 / Tk**2 - 1.0 / TR_K**2)
            + g("f") * np.log(Tk / TR_K))


def livermore_poly(a, T_degC):
    """
    Evaluate EQ3/6's Livermore Pitzer temperature function
    a1 + a2*(1/T - 1/Tr) + a3*ln(T/Tr) + a4*(T - Tr) + a5*(T^2 - Tr^2).
    `a` is a sequence of up to five coefficients.
    """
    Tk = np.asarray(T_degC, dtype=float) + 273.15
    a = list(a) + [0.0] * (5 - len(a))
    return (a[0] + a[1] * (1.0 / Tk - 1.0 / TR_K) + a[2] * np.log(Tk / TR_K)
            + a[3] * (Tk - TR_K) + a[4] * (Tk**2 - TR_K**2))


def gwb_to_livermore(coeffs, T_min_degC, T_max_degC, n_fit=201):
    """
    Convert GWB polynomial coefficients to the five-term Livermore function.

    Returns
    -------
    (list of 5 floats, float)
        The Livermore coefficients a1..a5 and the maximum absolute deviation
        between the two functions over the fitting temperature range (zero
        when the conversion is exact, i.e., when the GWB e coefficient is
        zero).
    """
    g = lambda k: float(coeffs.get(k, 0.0) or 0.0)
    a = [g("a"), g("d"), g("f"), g("b"), g("c")]
    if g("e") == 0.0:
        return a, 0.0

    if T_max_degC <= T_min_degC:
        # no range to fit over: keep the exactly-mappable terms and drop e
        return a, abs(g("e")) * abs(1.0 / (T_min_degC + 273.15)**2 - 1.0 / TR_K**2)

    # Keep a1 = a so that the value at 25 C (where every temperature term
    # vanishes) is reproduced exactly, and fit the remaining four terms to
    # the temperature dependence.
    T = np.linspace(T_min_degC, T_max_degC, n_fit)
    Tk = T + 273.15
    y = gwb_poly(coeffs, T) - g("a")
    X = np.column_stack([1.0 / Tk - 1.0 / TR_K, np.log(Tk / TR_K),
                         Tk - TR_K, Tk**2 - TR_K**2])
    fit, *_ = np.linalg.lstsq(X, y, rcond=None)
    dev = float(np.max(np.abs(X @ fit - y)))
    return [g("a")] + [float(v) for v in fit], dev


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_coeff_line(line):
    """Parse 'a= 1.0 b= 2.0 ...' pairs into a dict."""
    out = {}
    for k, v in re.findall(r"\b([a-f])\s*=\s*([-+0-9.eE]+)", line):
        out[k] = float(v)
    return out


def _parse_pairs(tokens, n):
    """Parse n (coefficient, name) pairs from a flat token list."""
    pairs = {}
    if len(tokens) < 2 * n:
        raise GWBDatasetError(f"Expected {n} coefficient/name pairs but found tokens {tokens}")
    for i in range(n):
        coef = float(tokens[2 * i])
        name = tokens[2 * i + 1]
        pairs[name] = pairs.get(name, 0.0) + coef
    return pairs


def parse_gwb_tdat(filename):
    """
    Parse a GWB thermodynamic dataset (.tdat) that uses the H-M-W (Pitzer)
    activity model.

    Returns
    -------
    dict with keys
        'header': {'temperatures': [...], 'pressures': [...], 'adh': {a..f},
                   'bdh': {a..f}, 'title': str, 'activity_model': str}
        'elements': {symbol: molwt}
        'basis': list of {'name', 'charge', 'elements': {sym: coef}}
        'aqueous': list of {'name', 'charge', 'reaction': {name: coef},
                            'logK': {a..f}, 'Tmin', 'Tmax'}
        'minerals': list of {'name', 'type', 'formula', 'V', 'reaction',
                             'logK', 'Tmin', 'Tmax'}
        'gases': list of {'name', 'reaction', 'logK', 'Tmin', 'Tmax'}
        'virial': list of {'species': [...], 'params': {param: {a..f}},
                           'alpha1', 'alpha2', 'set': str, 'Tmin', 'Tmax'}
    """

    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        raw_lines = [l.rstrip("\r\n") for l in f.readlines()]

    out = {"header": {"temperatures": [], "pressures": [], "adh": {}, "bdh": {},
                      "title": "", "activity_model": ""},
           "elements": {}, "basis": [], "aqueous": [], "minerals": [], "gases": [],
           "virial": []}

    # ---- header: everything before the elements section
    section_re = re.compile(r"^\s*(\d+)\s+(" + "|".join(GWB_SECTIONS) + r")\s*$")
    i_elements = None
    for i, l in enumerate(raw_lines):
        m = section_re.match(l)
        if m and m.group(2) == "elements":
            i_elements = i
            break
    if i_elements is None:
        raise GWBDatasetError(f"Could not find the elements section in {filename}")

    header_lines = raw_lines[:i_elements]
    if len(header_lines) > 0:
        out["header"]["title"] = header_lines[0].strip()
    numbers = []
    coeff_groups = []
    for l in header_lines:
        s = l.strip()
        if s.startswith("*") or s == "":
            continue
        m = re.match(r"^activity model:\s*(.*)$", s)
        if m:
            out["header"]["activity_model"] = m.group(1).strip()
            continue
        if re.fullmatch(r"[-+0-9.eE\s]+", s):
            numbers += [float(x) for x in s.split()]
            continue
        if re.search(r"\b[a-f]\s*=", s):
            c = _parse_coeff_line(s)
            if "a" in c or len(coeff_groups) == 0:
                coeff_groups.append(c)
            else:
                coeff_groups[-1].update(c)
    n_t = len(numbers) // 2
    out["header"]["temperatures"] = numbers[:n_t]
    out["header"]["pressures"] = numbers[n_t:2 * n_t]
    if len(coeff_groups) > 0:
        out["header"]["adh"] = coeff_groups[0]
    if len(coeff_groups) > 1:
        out["header"]["bdh"] = coeff_groups[1]

    if "h-m-w" not in out["header"]["activity_model"].lower() and \
            "pitzer" not in out["header"]["activity_model"].lower():
        raise GWBDatasetError(f"{filename} does not use the H-M-W (Pitzer) activity model "
                              f"(activity model: '{out['header']['activity_model']}').")

    # ---- sections
    lines = raw_lines[i_elements:]
    i = 0
    current = None
    end_of_sections = len(lines)
    while i < len(lines):
        l = lines[i]
        m = section_re.match(l)
        if m:
            current = m.group(2)
            i += 1
            block_lines = []
            # collect until -end-
            while i < len(lines) and not lines[i].strip().startswith("-end-"):
                block_lines.append(lines[i])
                i += 1
            _parse_section(current, block_lines, out)
            if current == "oxides":
                end_of_sections = i + 1
                break
        i += 1

    # ---- virial coefficients follow the oxides section
    _parse_virial(lines[end_of_sections:], out)

    return out


def _parse_section(section, block_lines, out):
    if section in ["free electron", "solid solutions", "oxides"]:
        return

    if section == "elements":
        for l in block_lines:
            m = re.match(r"^\s*\S+\s+\((\w+)\s*\)\s+mole wt\.=\s*([-+0-9.eE]+)", l)
            if m:
                out["elements"][m.group(1)] = float(m.group(2))
        return

    # split into species blocks: a block starts with a non-indented line
    blocks = []
    for l in block_lines:
        if l.strip() == "" or l.lstrip().startswith("*"):
            continue
        if not l[0].isspace():
            blocks.append([l])
        elif len(blocks) > 0:
            blocks[-1].append(l)

    for b in blocks:
        first = b[0].strip()
        m = re.match(r"^(\S+)(?:\s+type=\s*(.*))?$", first)
        name = m.group(1)
        entry = {"name": name, "type": (m.group(2) or "").strip(), "charge": 0.0,
                 "elements": {}, "reaction": {}, "logK": {}, "Tmin": None, "Tmax": None,
                 "formula": None, "V": None}
        j = 1
        while j < len(b):
            s = b[j].strip()
            j += 1
            mc = re.search(r"charge=\s*([-+0-9.]+)", s)
            if mc:
                entry["charge"] = float(mc.group(1))
            mf = re.search(r"formula=\s*(\S+)", s)
            if mf:
                entry["formula"] = mf.group(1)
            mv = re.search(r"mole vol\.=\s*([-+0-9.eE]+)", s)
            if mv:
                entry["V"] = float(mv.group(1))
            mt = re.search(r"TminK=\s*([-+0-9.eE]+)\s+TmaxK=\s*([-+0-9.eE]+)", s)
            if mt:
                entry["Tmin"] = float(mt.group(1)) - 273.15
                entry["Tmax"] = float(mt.group(2)) - 273.15
                continue
            me = re.match(r"^(\d+)\s+elements in species", s)
            mr = re.match(r"^(\d+)\s+species in reaction", s)
            if me or mr:
                n = int((me or mr).group(1))
                tokens = []
                while len(tokens) < 2 * n and j < len(b):
                    tokens += b[j].split()
                    j += 1
                pairs = _parse_pairs(tokens, n)
                if me:
                    entry["elements"] = pairs
                else:
                    entry["reaction"] = pairs
                continue
            if re.search(r"\b[a-f]\s*=", s) and not s.startswith("chi") and not s.startswith("Pcrit"):
                entry["logK"].update(_parse_coeff_line(s))

        if section == "basis species":
            out["basis"].append(entry)
        elif section in ["aqueous species", "redox couples"]:
            out["aqueous"].append(entry)
        elif section == "minerals":
            out["minerals"].append(entry)
        elif section == "gases":
            out["gases"].append(entry)


def _parse_virial(lines, out):
    set_names = ["beta", "theta", "lambda", "psi", "extra"]
    set_i = 0
    current = None
    for l in lines:
        s = l.strip()
        if s == "" or s.startswith("*"):
            continue
        if s.startswith("-end-"):
            if current is not None:
                out["virial"].append(current)
                current = None
            set_i += 1
            continue
        if not l[0].isspace():
            if current is not None:
                out["virial"].append(current)
            names = s.split()
            current = {"species": names, "params": {}, "alpha1": None, "alpha2": None,
                       "set": set_names[min(set_i, len(set_names) - 1)],
                       "Tmin": None, "Tmax": None}
            continue
        if current is None:
            continue
        ma = re.match(r"^alpha([12])\s*=\s*([-+0-9.eE]+)", s)
        if ma:
            current["alpha" + ma.group(1)] = float(ma.group(2))
            continue
        mt = re.search(r"TminK=\s*([-+0-9.eE]+)\s+TmaxK=\s*([-+0-9.eE]+)", s)
        if mt:
            current["Tmin"] = float(mt.group(1)) - 273.15
            current["Tmax"] = float(mt.group(2)) - 273.15
            continue
        mp = re.match(r"^(beta0|beta1|beta2|cphi|theta|lambda|psi|zeta|mu)\s+(.*)$", s)
        if mp:
            current["params"][mp.group(1)] = _parse_coeff_line(mp.group(2))
    if current is not None:
        out["virial"].append(current)


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def gwb_to_worm_name(name, mineral=False):
    """
    Convert a GWB species name to WORM conventions: "Ca++" -> "Ca+2",
    "CO2(aq)" -> "CO2", and for minerals "Ice(s)" -> "ice",
    "MgCl2:8H2O" -> "MgCl2*8H2O", "Halite" -> "halite".
    """
    name = name.strip()
    if mineral:
        if name.endswith("(s)"):
            name = name[:-3]
        name = name.replace(":", "*")
        # Names that are chemical formulas keep their capitalization; mineral
        # names are lower case in the WORM database.
        if re.search(r"\d", name) or "*" in name or "(" in name:
            return name
        return name.lower()
    return eq36_to_worm_name(name)


def _worm_formula(formula):
    return formula.replace(":", "*")


def _elements_of(species_name, comp_map):
    """Element composition {sym: coef} of a species from the GWB composition map."""
    return comp_map.get(species_name)


def formula_ox_string(elements, charge, oxidation_states=None):
    """
    Build a WORM-style 'formula_ox' string (e.g., "Ca+2 S+6 6O-2 4H+") from an
    element composition dictionary and a species charge, using fixed
    oxidation states. Returns "" if the oxidation states do not sum to the
    charge.
    """
    if oxidation_states is None:
        oxidation_states = DEFAULT_OXIDATION_STATES
    total = 0.0
    parts = []
    for el, n in elements.items():
        if el == "Z" or abs(n) < 1e-12:
            continue
        if el not in oxidation_states:
            return ""
        ox = oxidation_states[el]
        total += n * ox
        n_str = "" if abs(n - 1) < 1e-9 else (f"{n:g}")
        sign = "+" if ox > 0 else "-"
        mag = "" if abs(ox) == 1 else str(abs(ox))
        parts.append(f"{n_str}{el}{sign}{mag}")
    if abs(total - charge) > 1e-6:
        return ""
    return " ".join(parts)


def _fmt_coef(x):
    return f"{x:.4f}"


def _default_worm_basis():
    """Names of strict basis species in the bundled WORM database."""
    path = os.path.join(os.path.dirname(__file__), "databases", "wrm_data_latest.csv")
    basis = ["H2O", "H+", "O2(g)"]
    try:
        df = pd.read_csv(path, dtype=str)
        basis += list(df.loc[df["tag"] == "basis", "name"])
    except Exception:
        basis += ["Na+", "K+", "Ca+2", "Mg+2", "Cl-", "SO4-2", "HCO3-", "F-", "Br-",
                  "NO3-", "HPO4-2", "SiO2", "Li+", "Sr+2", "Ba+2", "Fe+2", "Al+3"]
    return basis


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def gwb_tdat_to_csv(tdat_filename,
                    logK_csv=None,
                    pitzer_csv=None,
                    basis_species=None,
                    name_map=None,
                    dataset_name=None,
                    ref=None,
                    oxidation_states=None,
                    verbose=1):
    """
    Convert a GWB Pitzer (H-M-W) thermodynamic dataset into an aqequil logK
    CSV and an aqequil Pitzer parameter CSV.

    Parameters
    ----------
    tdat_filename : str
        Path of the GWB .tdat file (e.g., "thermo_frezchem.tdat").
    logK_csv : str, optional
        If given, write the logK database to this CSV file.
    pitzer_csv : str, optional
        If given, write the Pitzer parameter database to this CSV file.
    basis_species : list of str, optional
        Names (WORM style) of the basis species that dissociation reactions
        should be written in terms of. Defaults to the strict basis species of
        the bundled WORM database (e.g., HCO3- rather than CO3-2). The dataset
        must contain exactly as many of these species as it has basis species
        of its own, and they must be independent.
    name_map : dict, optional
        {GWB name: desired name} overrides applied after the automatic
        conversion to WORM naming conventions.
    dataset_name : str, optional
        Short name used for the category_1 column (e.g., "frezchem" gives
        "frezchem_aq" and "frezchem_cr"). Defaults to the file name without
        "thermo_" and ".tdat".
    ref : str, optional
        Text for the ref2 column of the logK CSV. ref1 is the .tdat file name.
    oxidation_states : dict, optional
        {element: oxidation state} used to build the formula_ox column.
    verbose : int
        Verbosity level.

    Returns
    -------
    dict
        {"logK": DataFrame, "pitzer": DataFrame, "report": dict}
    """

    if name_map is None:
        name_map = {}
    if basis_species is None:
        basis_species = _default_worm_basis()
    if dataset_name is None:
        dataset_name = os.path.basename(tdat_filename)
        dataset_name = re.sub(r"^thermo[_-]", "", dataset_name)
        dataset_name = re.sub(r"\.tdat$", "", dataset_name, flags=re.I)
    if ref is None:
        ref = ""

    ds = parse_gwb_tdat(tdat_filename)
    temps = list(ds["header"]["temperatures"])
    if len(temps) == 0:
        raise GWBDatasetError(f"No principal temperatures were found in {tdat_filename}")
    if len(temps) > 8:
        raise GWBDatasetError("The aqequil logK CSV format supports at most eight "
                              f"temperatures, but {tdat_filename} has {len(temps)}.")
    T_min, T_max = min(temps), max(temps)
    report = {"skipped": {}, "refit_max_dev": {}, "basis_old": [], "basis_new": [],
              "temperatures": temps}

    # ---- names
    old_basis = [b["name"] for b in ds["basis"]]
    gwb_names = old_basis + [s["name"] for s in ds["aqueous"]] + \
                [m["name"] for m in ds["minerals"]] + [g["name"] for g in ds["gases"]]
    mineral_names = set(m["name"] for m in ds["minerals"])

    def wname(n):
        w = gwb_to_worm_name(n, mineral=(n in mineral_names))
        return name_map.get(n, name_map.get(w, w))

    worm = {n: wname(n) for n in gwb_names}
    dupes = [n for n in set(worm.values()) if list(worm.values()).count(n) > 1]
    if len(dupes) > 0:
        raise GWBDatasetError(f"Name conversion produced duplicate names: {dupes}. "
                              "Use name_map to resolve them.")

    # ---- element compositions (from GWB elements blocks, or from reactions)
    comp = {}
    charge = {}
    for b in ds["basis"]:
        comp[b["name"]] = dict(b["elements"])
        charge[b["name"]] = b["charge"]

    entries = {}  # gwb name -> entry dict (aqueous, mineral, gas)
    for s in ds["aqueous"]:
        s = dict(s); s["state"] = "aq"; entries[s["name"]] = s
        charge[s["name"]] = s["charge"]
    for m in ds["minerals"]:
        m = dict(m); m["state"] = "cr"; entries[m["name"]] = m
        charge[m["name"]] = 0.0
    for g in ds["gases"]:
        g = dict(g); g["state"] = "gas"; entries[g["name"]] = g
        charge[g["name"]] = 0.0

    # ---- formation vectors in the old (GWB) basis, resolved recursively
    n_old = len(old_basis)
    old_index = {n: i for i, n in enumerate(old_basis)}
    form_vec = {}
    form_logK = {}  # callable T -> logK of formation from the old basis

    for b in old_basis:
        v = np.zeros(n_old); v[old_index[b]] = 1.0
        form_vec[b] = v
        form_logK[b] = (lambda T: np.zeros_like(np.asarray(T, dtype=float)))

    def resolve(name, stack=()):
        if name in form_vec:
            return True
        if name not in entries:
            return False
        if name in stack:
            raise GWBDatasetError(f"Circular reaction definitions involving {name}")
        e = entries[name]
        v = np.zeros(n_old)
        terms = []
        for sp, coef in e["reaction"].items():
            if not resolve(sp, stack + (name,)):
                report["skipped"][worm.get(name, name)] = f"reaction involves unknown species {sp}"
                return False
            v += coef * form_vec[sp]
            terms.append((coef, sp))
        form_vec[name] = v
        logK_coeffs = dict(e["logK"])

        def fl(T, logK_coeffs=logK_coeffs, terms=terms):
            # species -> products has log K_diss; formation log K = -log K_diss
            # plus the formation log K of any non-basis products
            val = -gwb_poly(logK_coeffs, T)
            for coef, sp in terms:
                val = val + coef * form_logK[sp](T)
            return val
        form_logK[name] = fl

        # element composition from the reaction if not given
        if name not in comp:
            c = {}
            for coef, sp in terms:
                for el, n in comp[sp].items():
                    c[el] = c.get(el, 0.0) + coef * n
            comp[name] = {el: n for el, n in c.items() if abs(n) > 1e-9}
        return True

    for name in list(entries.keys()):
        resolve(name)

    # ---- new basis
    new_basis = [n for n in gwb_names if worm[n] in basis_species and n in form_vec
                 and entries.get(n, {}).get("state", "aq") == "aq"]
    if len(new_basis) != n_old:
        raise GWBDatasetError(
            f"The dataset has {n_old} basis species but {len(new_basis)} of its species "
            f"are basis species in the target database ({[worm[n] for n in new_basis]}). "
            "Provide a `basis_species` list with exactly one basis species per element "
            "(plus H2O and H+).")
    M = np.array([form_vec[n] for n in new_basis])  # rows: new basis in old basis
    if abs(np.linalg.det(M)) < 1e-10:
        raise GWBDatasetError("The chosen basis species are not linearly independent: "
                              f"{[worm[n] for n in new_basis]}")
    report["basis_old"] = [worm[n] for n in old_basis]
    report["basis_new"] = [worm[n] for n in new_basis]

    def to_new_basis(name):
        # solve c such that c @ M = form_vec[name]
        c = np.linalg.solve(M.T, form_vec[name])
        c[np.abs(c) < 1e-9] = 0.0
        return c

    # ---- logK rows
    today = datetime.date.today().strftime("%m/%d/%Y")
    rows = []
    for name in gwb_names:
        if name in new_basis:
            continue
        if name not in form_vec:
            continue
        e = entries.get(name)
        state = e["state"] if e is not None else "aq"
        c = to_new_basis(name)

        # temperature grid for this species
        Tmin_sp = e["Tmin"] if (e is not None and e["Tmin"] is not None) else -np.inf
        Tmax_sp = e["Tmax"] if (e is not None and e["Tmax"] is not None) else np.inf
        if np.isfinite(Tmin_sp) and np.isfinite(Tmax_sp) and abs(Tmax_sp - Tmin_sp) < 1e-6:
            # A validity range collapsed to a single temperature marks a
            # placeholder value (e.g., OH- in thermo_coldchem.tdat). Such an
            # entry would be rejected by aqequil for any other temperature, so
            # it is left out; the main thermodynamic database supplies it.
            report["skipped"][worm[name]] = (f"valid only at {Tmin_sp:g} C (placeholder value); "
                                             "not converted")
            continue
        T_sp = [T for T in temps if Tmin_sp - 1e-6 <= T <= Tmax_sp + 1e-6]
        if len(T_sp) == 0:
            # the validity range does not contain a principal temperature;
            # tabulate at the range itself
            T_sp = sorted(set([max(Tmin_sp, T_min), min(Tmax_sp, T_max)]))
        T_arr = np.array(T_sp, dtype=float)
        logK_form_new = form_logK[name](T_arr) - sum(
            c[j] * form_logK[nb](T_arr) for j, nb in enumerate(new_basis))
        logK_diss = -logK_form_new

        dissrxn = f"-1.0000 {worm[name]}"
        for j, nb in enumerate(new_basis):
            if c[j] != 0.0:
                dissrxn += f" {_fmt_coef(c[j])} {worm[nb]}"

        if state == "cr":
            formula = _worm_formula(e["formula"]) if e["formula"] else worm[name]
        elif state == "gas":
            # WORM gas formulas omit the "(g)" suffix, e.g. CO2(g) has formula CO2
            formula = re.sub(r"\(g\)$", "", worm[name])
        else:
            formula = worm[name]

        elements = comp.get(name, {})
        z = charge.get(name, 0.0)
        row = {"name": worm[name], "abbrv": "", "formula": formula, "state": state,
               "ref1": os.path.basename(tdat_filename), "ref2": ref, "date": today}
        for k in range(8):
            row[f"logK{k+1}"] = round(float(logK_diss[k]), 5) if k < len(T_sp) else np.nan
        for k in range(8):
            row[f"T{k+1}"] = T_sp[k] if k < len(T_sp) else np.nan
        for k in range(8):
            row[f"P{k+1}"] = "psat" if k < len(T_sp) else np.nan
        row["azero"] = 4.0 if z != 0 else 3.0
        row["dissrxn"] = dissrxn
        row["tag"] = ""
        row["formula_ox"] = formula_ox_string(elements, z, oxidation_states)
        row["category_1"] = f"{dataset_name}_{state}"
        row["V"] = e["V"] if (state in ["cr"] and e is not None and e["V"] is not None) else np.nan
        rows.append(row)

    logK_cols = ["name", "abbrv", "formula", "state", "ref1", "ref2", "date"] + \
                [f"logK{k}" for k in range(1, 9)] + [f"T{k}" for k in range(1, 9)] + \
                [f"P{k}" for k in range(1, 9)] + ["azero", "dissrxn", "tag", "formula_ox",
                                                  "category_1", "V"]
    logK_df = pd.DataFrame(rows, columns=logK_cols)

    # ---- Pitzer rows
    prows = []
    for blk in ds["virial"]:
        names = blk["species"]
        if any(n not in worm for n in names):
            unknown = [n for n in names if n not in worm]
            report["skipped"][" ".join(names)] = f"Pitzer block involves unknown species {unknown}"
            continue
        wnames = [worm[n] for n in names]
        zs = [charge.get(n, 0.0) for n in names]
        Tlo = blk["Tmin"] if blk["Tmin"] is not None else T_min
        Thi = blk["Tmax"] if blk["Tmax"] is not None else T_max
        note_range = f"{Tlo:g} to {Thi:g} C"

        params = dict(blk["params"])
        # GWB writes neutral-cation-anion (zeta) triplets in its psi set
        if "psi" in params and len(names) == 3 and any(z == 0 for z in zs):
            params["zeta"] = params.pop("psi")
        if "lambda" in params and len(names) == 1:
            names = names * 2; wnames = wnames * 2

        # alpha values for beta(1) and beta(2); a zero alpha is only
        # meaningful when the corresponding beta is identically zero, in
        # which case the usual EQ3/6 defaults are written instead
        alpha1 = blk["alpha1"]
        alpha2 = blk["alpha2"]
        beta1_zero = all(v == 0 for v in params.get("beta1", {"a": 0}).values())
        beta2_zero = all(v == 0 for v in params.get("beta2", {"a": 0}).values())
        if alpha1 is not None and alpha1 == 0 and beta1_zero:
            alpha1 = 2.0
        if alpha2 is not None and alpha2 == 0 and beta2_zero:
            alpha2 = 12.0

        for p, coeffs in params.items():
            if p not in GWB_VIRIAL_PARAMS:
                continue
            a, dev = gwb_to_livermore(coeffs, Tlo, Thi)
            note = f"GWB {blk['set']} set; valid {note_range}"
            if dev > 0:
                note += f"; 1/T^2 term refit to Livermore function, max dev {dev:.2e}"
                report["refit_max_dev"][f"{p} {' '.join(wnames)}"] = dev
            prow = {"species1": wnames[0], "species2": wnames[1] if len(wnames) > 1 else "",
                    "species3": wnames[2] if len(wnames) > 2 else "", "param": p,
                    "a1": a[0], "a2": a[1], "a3": a[2], "a4": a[3], "a5": a[4],
                    "alpha": np.nan, "ref1": os.path.basename(tdat_filename), "ref2": ref,
                    "date": today, "note": note}
            if p == "beta1" and alpha1 is not None:
                prow["alpha"] = alpha1
            if p == "beta2" and alpha2 is not None:
                prow["alpha"] = alpha2
            prows.append(prow)

    pitzer_cols = ["species1", "species2", "species3", "param", "a1", "a2", "a3", "a4", "a5",
                   "alpha", "ref1", "ref2", "date", "note"]
    pitzer_df = pd.DataFrame(prows, columns=pitzer_cols)

    # drop the a5 column if no parameter uses it
    if len(pitzer_df) > 0 and (pitzer_df["a5"] == 0).all():
        pitzer_df = pitzer_df.drop(columns=["a5"])

    if verbose > 0:
        print(f"{os.path.basename(tdat_filename)}: {len(logK_df)} logK entries "
              f"({(logK_df['state']=='aq').sum()} aqueous, {(logK_df['state']=='cr').sum()} "
              f"minerals, {(logK_df['state']=='gas').sum()} gases) tabulated at "
              f"{', '.join(f'{T:g}' for T in temps)} C; {len(pitzer_df)} Pitzer parameters.")
        print("Basis species:", ", ".join(report["basis_old"]), "->", ", ".join(report["basis_new"]))
        if len(report["refit_max_dev"]) > 0:
            worst = max(report["refit_max_dev"].values())
            print(f"{len(report['refit_max_dev'])} Pitzer parameters with a 1/T^2 term were "
                  f"refit to the EQ3/6 temperature function (max deviation {worst:.2e}).")
        for k, v in report["skipped"].items():
            print("Skipped", k, ":", v)

    if logK_csv is not None:
        logK_df.to_csv(logK_csv, index=False)
        if verbose > 0:
            print("Wrote", logK_csv)
    if pitzer_csv is not None:
        pitzer_df.to_csv(pitzer_csv, index=False)
        if verbose > 0:
            print("Wrote", pitzer_csv)

    return {"logK": logK_df, "pitzer": pitzer_df, "report": report, "dataset": ds}
