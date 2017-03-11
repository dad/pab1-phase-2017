# Overview

This README describes how to replicate the evolutionary analysis reported in Riback, Katanski et al. 2017.

The makefile `build-alignment.mak` automates most of the steps.

Outline of steps:

## For poly(A)-binding protein alignment.
0. Begin with the poly(A)-binding protein sequences obtained from [SMART][http://smart.embl-heidelberg.de], in `../data/smart-pabp-orthologs.fa`
1. Filter the sequences to eliminate fragments and sequences noncanonical amino acid
2. Align the sequences
3. Trim alignment to remove sequences from the same species and sequences with >0.95 sequence identity
4. Fetch the species and their taxonomic relationships
5. Filter out sequences from the alignment for which species can't be found. Sort by phylogeny.
6. Extract the aligned regions corresponding to the proline-rich P domain and everything else.
7. Calculate residue frequencies in these regions.

## For DisProt
Calculate residue frequencies in regions/proteins included in data/disprot-6.02.fasta. A command to do this is included in `build-alignment.mak`.

## For the yeast proteome
Use the frequencies already calculated in data/scer-proteome-aa-freqs.txt

## Making figures
All figures are made using the `figures-pabp-evol.R` [R] script using data either provided in the repository or calculated as outlined above and shown explicitly in `build-alignment.mak`.


# Required packages

Several utility libraries and scripts are in the dad/base package.
```
cd <your-dir> 
git clone git@github.com:dad/base.git
```
and add `<your-dir>/base/src` to your PYTHONPATH. `<your-dir>` should be the directory in which `pab1-phase-2017` is also contained, if you wish to use the makefile.

Also required: 

1. [Python](python.org) 3.5+
1. [R] with the `taxize` package, for tree-building, `ggplot2` and `cowplot` for plotting, `data.table` and `Hmisc` for miscellaneous analysis functions.
1. [MUSCLE](http://www.drive5.com/muscle/) sequence alignment software.

[R]: http://r-project.org