# Overview

This README describes how to replicate the evolutionary analysis reported in Riback, Katanski et al. 2017.

Outline of steps:

## For poly(A)-binding protein alignment.
0. Begin with the poly(A)-binding protein sequences obtained from [SMART][http://smart.embl-heidelberg.de], in `../data/smart-pabp-orthologs.fa`
1. Filter the sequences to eliminate fragments and sequences noncanonical amino acid
2. Align the sequences
3. Trim alignment to remove sequences from the same species and sequences with >0.95 sequence identity
4. Fetch the species and their taxonomic relationships
5. Filter out sequences from the alignment for which species can't be found. Sort by phylogeny.
6. Extract the aligned regions corresponding to the proline-rich domain and everything else.
7. Calculate residue frequencies in these regions.

<!--## For DisProt-->

<!--## For the yeast proteome.-->

The makefile `build-alignment.mak` contains commands for each of these steps.

# Required packages

Several utility libraries and scripts are in the dad/base package.
```
cd <your-dir> 
git clone git@github.com:dad/base.git
```
and add `<your-dir>/base/src` to your PYTHONPATH. `<your-dir>` should be the directory in which `pab1-phase-2017` is also contained, if you wish to use the Makefile.

Also required: 
1. [R](http://r-project.org) with the `taxize` package, for tree-building.
1. [MUSCLE][http://www.drive5.com/muscle/] alignment software 

