# AqEquil

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7601102.svg)](https://doi.org/10.5281/zenodo.7601102)

Boyer, G., Robare, J., Ely, T., Shock, E.L.

## About

AqEquil is a Python 3 package that enables users to rapidly perform aqueous speciation calculations of water chemistry data for multiple samples by interfacing with [geochemical speciation software EQ3/6](https://github.com/LLNL/EQ3_6) (Wolery 2013, [Wolery 1979](https://inis.iaea.org/collection/NCLCollectionStore/_Public/10/474/10474294.pdf)).

Water sample data in CSV format is automatically converted to a format readable by EQ3 and then speciated. Distributions of aqueous species, mineral saturation indices, oxidation reduction potentials, and more are data-mined and returned as Pandas tables and interactive Plotly visualizations.

Speciated fluids can be further reacted with minerals in mass transfer calculations to produce tables and interactive diagrams of reaction paths and composition changes as a function of reaction progress.

## Requirements

AqEquil has only been tested with Ubuntu LTS 20.04.

This installation requires the Linux version of EQ3/6 v8.0a, which can downloaded [here](https://github.com/LLNL/EQ3_6). Installation instructions are provided there.

AqEquil must be installed into an environment with an R installation. See [these instructions](https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/) for installing R with Anaconda.

Additionally, the CHNOSZ package must be installed in R (see instructions below).

## Installation

### Installing EQ3/6 for Linux

I recommend using [this github version of EQ3/6 v.8.0a adapted by the 39Alpha team](https://github.com/39alpha/eq3_6/tree/main). Installation instructions are found there.

Create an environment variable called `EQ36DO` and set it to wherever you installed EQ3/6. (`/usr/local/bin` by default). Set another environment variable called `EQ36DA` to the directory containing your data1 thermodynamic database files (if you have one).

### Installing CHNOSZ version 1.4.3

Open an R session. Install CHNOSZ version 1.4.3 package with:

```install.packages('http://cran.r-project.org/src/contrib/Archive/CHNOSZ/CHNOSZ_1.4.3.tar.gz', repos=NULL, type='source')```

Once CHNOSZ is installed you can quit the R session.

Compatibility with CHNOSZ v.2.0.0 is forthcoming.

### Installing AqEquil

Install AqEquil using pip:

```pip install AqEquil```

### Usage

See this [demo notebook](https://nbviewer.jupyter.org/github/worm-portal/WORM-Library/blob/master/3-Aqueous-Speciation/1-Introduction-to-Aq-Speciation/2-Intro-to-Multi-Aq-Speciation.ipynb) for usage examples.
