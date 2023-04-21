# khipu: generalized tree structure to annotate untargeted metabolomics and stable isotope tracing data

[![Documentation Status](https://readthedocs.org/projects/khipu/badge/?version=latest)](https://khipu.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://img.shields.io/badge/DOI-doi%2F10.1021%2Facs.analchem.2c05810-blue)](https://pubs.acs.org/doi/10.1021/acs.analchem.2c05810)


Pre-annotation tool to annotate degenerate ions in relationships to the original compound and infer neutral mass. 

This applies to regular LC-MS data, but also enables easy analysis of isotope tracing and chemical labeling data.

![khipugram](doc/khipugram.png)

## Implementation overview
Khipu is developed as an open source Python 3 package, and available to install from the standard PyPi repository via the pip tool. It is freely available on GitHub (https://github.com/shuzhao-li/khipu) under a BSD 3-Clause License. The graph operations are supported by the networkx library, tree visualization aided by the treelib library. Khipu uses our package mass2chem for search functions. The data model of “empirical compound” is described in the metDataModel package. The package is designed in a modular way to encourage reuse.

The classes of Weavor and Khipu contain main algorithms, supported by numerous utility functions. All functions are documented in the source via docstrings. Examples of reuse are given in wrapper functions and in Jupyter notebooks. It can be run as a standalone command line tool. Users can use a feature table from any preprocessing tool as input and get annotated empirical compounds in JSON and tab delimited formats.

## Installation and Use
Install as a package (some systems may require pip3):

    pip install khipu-metabolomics

Run as a command line tool after installation:

    khipu -i testdata/ecoli_pos.tsv -o this_test

This will output pre-annotation to two files of JSON and tab delimited formats, this_test.json and this_test.tsv.

Run from source code:

    python3 -m khipu.main -i testdata/ecoli_pos.tsv -o this_test

Run test:

    python3 -m khipu.test
    (This downloads and uses test data from GitHub.)

Best used as a library for software development or in a Jupyter Notebook for data analysis. 

## Demo notebooks
We have provided multiple demo notebooks under

    notebooks/

They include algorithm demostrations, data analysis examples, use of custom isotope and adduct patterns.

## Algorithm overview 
1. Start with an initial list of isotope patterns and adduct patterns (see khipu grid below). Search feature list to get all pairs that match any of the pattern. The initial adduct patterns are trimmed to reduce ambiguity. 

2. Connect all pattern-matched feature pairs to an overall network, which is further partitioned into connected subnetworks.

3. Each subnetwork becomes a khipu instance. The subnetwork is inspected, redundant nodes removed, and converted to an optimal tree structure (see below). A khipu is essentially an 'empirical compound' that is used for downstream annotation and analysis.

4. This library supports tree and grid visualization in plain text. Once imported to a Jupyter Notebook, one can use enhanced visualization schemes. 

5. The library can also be used by others for extended tools. Our data processing tool, asari, uses khipu for preannotation. Additional documentation, more for developers, is provided under `doc/`.

## Assignment of ion species in a khipu to grid
1. Separate isotope edges and adduct edges.
2. The isotope edges form their own groups by shared nodes, each group belong to one adduct type. Each group of connected isotope edges is treated as one "branch".
3. Establish a "trunk" of adducts with a root and a path for adducts, by optimizing the number of nodes explained.
4. Assign each isotopic branch to the adduct trunk.
5. Re-align isotopes in all branches to establish optimal match to the khipu grid. 
6. Based on available ions and the theoretical "khipu grid", the neutral mass can be obtained via linear regression. 

Some ions may come into the initial network by mistakes or unresolved signals.
The are removed from the established khipu, and sent off to form a new khipu.

## The khipu grid
Initial grid may look like this:

                   M+H[+]   M+NH4[+]   M+Na[+]   M+HCl+H[+]   M+K[+]   M+ACN+H[+]
    M0           1.007276  18.033826  22.989276  36.983976  38.963158  42.033825
    13C/12C      2.010631  19.037181  23.992631  37.987331  39.966513  43.037180
    13C/12C*2    3.013986  20.040536  24.995986  38.990686  40.969868  44.040535
    13C/12C*3    4.017341  21.043891  25.999341  39.994041  41.973223  45.043890
    13C/12C*4    5.020696  22.047246  27.002696  40.997396  42.976578  46.047245
    13C/12C*5    6.024051  23.050601  28.006051  42.000751  43.979933  47.050600
    13C/12C*6    7.027406  24.053956  29.009406  43.004106  44.983288  48.053955

This can be extended by searching for additional ions. But the core construction should be done first.

## Applicable to isotope tracing
The search pattern for isotopes is often dependent on the biochemical experiment.
Users can overwrite the default by supplying their search patterns (see demo notebooks).
Search patterns are separate from search functions, lending flexibility to data analysis.

The next step is to apply Khipu to chemical derivatization experiments.
In chemical derivatization experiments, the origin compound and derivatized compound can be both measured in the LC-MS data.
We have separate khipu trees for each, then link them by the m/z shift from derivatization.
Because derivatization is a reaction that occurs before LC-MS, and
LC-MS measures whatever compounds that are present in samples.

## Test data
Three datasets are included under testdata/. All three tables were generated by asari v1.9.2.
The automated khipu.test downloads ecoli_pos.tsv from GitHub remotely.
- The ecoli_pos.tsv was generated by Li lab using the credentialed E. coli sample from Cambridge Isotopes.
- The yeast datasets were from the NetID paper by Rabinowitz lab. The yeast_neg table is features that are filted by SNR > 100 to serve as a cleaner demo.

Input tables are tab delimited text files.
The first columns are feature ID, m/z, rtime, followed by intensities.
Users can specify the start column and end column of intensity data.

## Detailed use of command and parameters

    >>> khipu -h

    usage: main.py [-h] [-v] [-m MODE] [--ppm PPM] [--rtol RTOL] [-i INPUT]
                [-s START] [-e END] [-o OUTPUT]

    khipu, annotating metabolomics features to empCpds

    optional arguments:
    -h, --help            show this help message and exit
    -v, --version         print version and exit
    -m MODE, --mode MODE  mode of ionization, pos or neg
    --ppm PPM             mass precision in ppm (part per million), same as
                            mz_tolerance_ppm
    --rtol RTOL           tolerance of retention time match, arbitrary unit
                            dependent on preprocessing tool
    -i INPUT, --input INPUT
                            input file as feature table
    -s START, --start START
                            start column for intensity in input table
    -e END, --end END     end column for intensity in input table
    -o OUTPUT, --output OUTPUT
                            prefix of output files


## What's "khipu"?
Khipu is a recording device using knots, often 2-level of strings,
historically used by people in Andean South America, includign Inca (https://en.wikipedia.org/wiki/Quipu).
The format is similar to how we represent isotopes and adducts in our data.
We chose "khipu" over the spelling of "quipu", to pay respect to the indigenous people.
