# AqEquil

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6382932.svg)](https://doi.org/10.5281/zenodo.6382932)

Boyer, G., Robare, J., Ely, T., Shock, E.L.

## About

AqEquil is a Python 3 package that enables users to rapidly perform aqueous speciation calculations of water chemistry data for multiple samples by interfacing with [geochemical speciation software EQ3/6](https://github.com/LLNL/EQ3_6) (Wolery 2013, [Wolery 1979](https://inis.iaea.org/collection/NCLCollectionStore/_Public/10/474/10474294.pdf)). Water sample data in CSV format is automatically converted to a format readable by EQ3 and then speciated. Distributions of aqueous species, mineral saturation indices, oxidation reduction potentials, and more are data-mined and returned as Pandas tables and Matplotlib visualizations. Additionally, EQ6 can be called to speciate .6i input files.

## Requirements

AqEquil has only been tested with Ubuntu LTS 20.04.

This installation requires the Linux version of EQ3/6 v8.0a, which can downloaded [here](https://www-gs.llnl.gov/energy-homeland-security/geochemistry). Installation instructions are provided there.

AqEquil must be installed into an environment with an R installation. See [these instructions](https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/) for installing R with Anaconda.

Additionally, the CHNOSZ package must be installed in R (see instructions below).

## Installation

### Installing EQ3/6 for Linux

Installation instructions are packaged along with the Linux distribution of [EQ3/6 v8.0a](https://www-gs.llnl.gov/energy-homeland-security/geochemistry).

Set the environment variable EQ36DA to the EQ3/6 database directory containing data1 files (the 'db' folder by default). Likewise, set EQ36CO to the directory with the eq3 executable (the 'bin' folder by default).

### Installing CHNOSZ

Open an R session. Install the CHNOSZ package with:

```install.packages("CHNOSZ")```

Once CHNOSZ is installed you can quit the R session.

### Installing AqEquil

Install AqEquil using pip:

```pip install AqEquil```

### Usage

See this [demo notebook](https://nbviewer.jupyter.org/github/worm-portal/WORM-Library/blob/master/3-Aqueous-Speciation/1-Introduction-to-Aq-Speciation/2-Intro-to-Multi-Aq-Speciation.ipynb) for usage examples.
