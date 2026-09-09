# AqEquil

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5534831.svg)](https://doi.org/10.5281/zenodo.5534831)

Boyer, G., Robare, J., Park, N., Ely, T., Shock, E.L.

## About

AqEquil is a Python 3 package that enables users to rapidly perform aqueous speciation calculations of water chemistry data for multiple samples by interfacing with [geochemical speciation software EQ3/6](https://github.com/LLNL/EQ3_6) (Wolery 2013, [Wolery 1979](https://inis.iaea.org/collection/NCLCollectionStore/_Public/10/474/10474294.pdf)). AqEquil uses [a modified version of EQ3/6](https://github.com/39alpha/eq3_6/tree/main) created by the 39Alpha research team for easy local installation.

Water sample data in CSV format is automatically converted to a format readable by EQ3 and then speciated. Distributions of aqueous species, mineral saturation indices, oxidation reduction potentials, and more are data-mined and returned as Pandas tables and interactive Plotly visualizations.

Speciated fluids can be further reacted with minerals or other fluids in mass transfer calculations to produce tables and interactive diagrams of reaction paths and composition changes as a function of reaction progress.

Development of AqEquil was made possible by National Science Foundation (NSF) grants EAR-1949030 and EAR-2149016.

## Requirements

AqEquil works on Linux, macOS, and Windows.

**Python Requirements:**
- Python >= 3.10
- pandas, numpy, matplotlib, plotly, and other dependencies (automatically installed)

**Note:** As of version 1.0.0, EQ3/6 executables are bundled with aqequil and no longer need to be installed separately. The package includes pre-compiled binaries for Linux, macOS, and Windows.

## Installation

Install AqEquil using pip:

```bash
pip install aqequil
```

The bundled EQ3/6 executables will be automatically installed with the package. No additional configuration or environment variables are needed.

## Usage

See this [demo notebook](https://nbviewer.jupyter.org/github/worm-portal/WORM-Library/blob/master/3-Aqueous-Speciation/1-Introduction-to-Aq-Speciation/2-Intro-to-Multi-Aq-Speciation.ipynb) for usage examples.

### Activity models

Speciation calculations use the B-dot equation by default (`activity_model="b-dot"`), or the Davies equation (`activity_model="davies"`). High ionic strength solutions can be modeled with Pitzer's equations (`activity_model="pitzer"`), which require a thermodynamic database that contains Pitzer ion-interaction parameters. EQ3/6 requires the activity model to match the database, so `speciate()` chooses the model automatically (`activity_model="auto"`) based on the database that was loaded.

A Pitzer calculation can be run with an EQ3/6 Pitzer data0 file, e.g. `data0.ypf`, placed in the working directory:

```python
import aqequil
ae = aqequil.AqEquil(db="ypf")
speciation = ae.speciate(input_filename="sylhal.csv", charge_balance_on="Cl-")
```

or by loading a CSV of Pitzer parameters alongside a WORM-style CSV thermodynamic database:

```python
ae = aqequil.AqEquil(db="wrm_data_latest.csv", pitzer="pitzer_params.csv")
speciation = ae.speciate(input_filename="sylhal.csv", charge_balance_on="Cl-")
```

The Pitzer parameter CSV has one row per parameter, with the columns `species1, species2, species3, param, a1, a2, a3, a4, alpha, ref1, ref2, date, note`, where `param` is one of `beta0`, `beta1`, `beta2`, `cphi`, `theta`, `lambda`, `psi`, `zeta`, or `mu`, and `a1` through `a4` are coefficients of the EQ3/6 temperature function `x(T) = a1 + a2*(1/T - 1/298.15) + a3*ln(T/298.15) + a4*(T - 298.15)`. See the docstring of `aqequil.pitzer` for details. The parameters of an existing Pitzer data0 file can be extracted into this format with `aqequil.pitzer_data0_to_csv("data0.ypf", "pitzer_params.csv", worm_names=True)`. An example file derived from `data0.ypf.R2` is bundled as `aqequil/test_data/pitzer_params_ypf.csv`.

EQPT checks every possible ion pair and triplet for parameters, so compiling a Pitzer data0 built from a large CSV database takes a few minutes per data0 file. Reduce the database with `exclude_organics` or `exclude_category` to speed this up.

#### Pitzer parameters and aqueous complexes belong together

A set of Pitzer parameters and the list of aqueous complexes it was fit with form one model. Most Pitzer datasets treat electrolytes such as NaCl, KCl and MgCl2 as fully dissociated and absorb the short-range Na-Cl interaction into the beta and C-phi parameters; a dataset that does include a complex (for example the CaCl+ and CaCl2(aq) of `data0.ypf`) fitted its Ca-Cl parameters with those complexes present. Combining Pitzer parameters with complexes from a different source counts the same interaction twice. With `data0.ypf` parameters and the WORM aqueous complexes, the halite and sylvite saturated brine of `sylhal.csv` comes out with 17.5 molal Na instead of the correct 4.6 molal, because WORM's NaCl(aq) complex sits on top of Na-Cl parameters that already describe that association.

To use the parameters of a Pitzer data0 file with a WORM-style CSV database, convert the species blocks of the data0 file (auxiliary basis species, aqueous complexes, minerals and gases with their log K grids) into a logK CSV as well, and exclude the WORM aqueous species and minerals so that only the data0 file's own species remain:

```python
import aqequil
aqequil.pitzer_data0_to_csv("data0.ypf", "ypf_pitzer.csv", worm_names=True)
aqequil.data0_to_logK_csv("data0.ypf", "ypf_logK.csv")

ae = aqequil.AqEquil(db="wrm_data_latest.csv",
                     logK="ypf_logK.csv",
                     pitzer="ypf_pitzer.csv",
                     exclude_category={"category_1": ["inorganic_aq", "organic_aq", "inorganic_cr", "organic_cr"]})
speciation = ae.speciate(input_filename="sylhal_worm.csv", charge_balance_on="Cl-")
```

`data0_to_logK_csv` converts species names to WORM conventions (`Ca++` to `Ca+2`, `Halite` to `halite`, `MgCl2:4H2O` to `MgCl2*4H2O`), rewrites reactions in terms of the WORM strict basis species where the data0 file uses a different basis (for example `Cr+++` becomes an auxiliary species written in terms of `CrO4-2`), drops grid entries marked as no data so that a species with a 25 °C value only gets a one-point grid, and skips species whose elements have no basis species in the WORM database (Np, Pu and Cm in `data0.ypf`). Solids whose name coincides with an aqueous species get a `(cr)` suffix (`CaCl2(cr)`). The bundled example `aqequil/test_data/logK_ypf.csv` was made this way. Speciation results obtained with these CSVs agree with `data0.ypf` used directly to within about 0.01 in log gamma and mineral saturation indices; the small differences come from the Debye-Hückel A-phi values, which aqequil computes from the water model rather than reading from the data0 file.

Note that a logK species is only used for a sample when the sample temperature lies within the species' grid, so a species with data at 25 °C only participates in 25 °C calculations. If one input file mixes temperatures, the check is made against all sample temperatures in the file at once; put samples at different temperatures in separate input files to keep species with limited temperature ranges.

#### Using a Geochemist's Workbench Pitzer dataset (FREZCHEM, ColdChem)

GWB thermodynamic datasets that use the H-M-W (Pitzer) activity model, such as `thermo_frezchem.tdat` and `thermo_coldchem.tdat`, can be converted into a logK CSV (log K values tabulated at the dataset's temperatures, reactions rewritten in terms of the WORM basis species) and a Pitzer parameter CSV with WORM species names:

```python
import aqequil
aqequil.gwb_tdat_to_csv("thermo_coldchem.tdat", "coldchem_logK.csv", "coldchem_pitzer_params.csv")

ae = aqequil.AqEquil(logK="coldchem_logK.csv",
                     pitzer="coldchem_pitzer_params.csv",
                     exclude_category={"category_1": ["inorganic_aq", "organic_aq", "inorganic_cr", "organic_cr"]})
speciation = ae.speciate(input_filename="coldchem_worm.csv", charge_balance_on="Cl-")
```

Species in a logK CSV replace species of the same name in the main thermodynamic database, so the dataset's own log K values are used for minerals such as halite or gypsum. The `exclude_category` argument above removes the rest of the WORM aqueous species and minerals (strict basis species are always kept), which restricts the model to the species the Pitzer parameters were fit for and makes the EQPT step fast. GWB temperature functions include a 1/T² term that EQ3/6 lacks; parameters that use it are refit to the EQ3/6 function over the dataset's temperature range and the refit error is recorded in the `note` column. Water properties come from pyCHNOSZ. With the default SUPCRT92 water model, samples must be at or above 0.01 °C. To use these datasets at sub-zero temperatures, select the IAPWS-95 water model, which extrapolates into the supercooled region down to about -30 °C at 1 bar (`AqEquil(..., water_model="IAPWS95")`); pyCHNOSZ then supplies density, the dielectric constant, Debye-Hückel parameters and Born functions from IAPWS-95 and Archer and Wang (1990), exactly as R CHNOSZ does. Below about -30 °C the equation of state has no liquid root, so that range is not accessible.

## Bundled Software

This package includes pre-compiled binaries from [EQ3/6 v8.0a](https://github.com/39alpha/eq3_6), a software package for geochemical modeling developed by Thomas Wolery at Lawrence Livermore National Laboratory and updated by [39 Alpha](https://github.com/39alpha/eq3_6).

**EQ3/6 License:** BSD 3-Clause License
**Copyright:** (c) 1987, 1990-1993, 1995, 1997, 2002, 2013 The Regents of the University of California, Lawrence Livermore National Laboratory.

See `THIRD_PARTY_LICENSES.txt` for the full EQ3/6 license text.

**References:**
- Wolery, T. J., and USDOE. EQ3/6 A Software Package for Geochemical Modeling. Computer software. December 13, 2010. https://www.osti.gov//servlets/purl/1231666. doi:https://doi.org/10.11578/dc.20210416.44.
- Wolery, T. J. and R. L. Jarek. Software User's Manual EQ36, Version 8.0. U.S. Tech. Rep. 2003. Department of Energy, Office of Civilian Radioactive Waste Management, Office of Repository Development. 10813-UM-8.0-00.
