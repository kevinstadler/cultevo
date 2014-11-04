## Simple Mantel tests for R

This package provides functions for computing distance matrices and performing Mantel tests on them.

## Installation

    devtools::install_github("kevinstadler/mantel")

If this doesn't work you probably need to install the `devtools` package first:

    install.packages("devtools")

## Usage

A nice HTML documentation is available at https://rawgit.com/kevinstadler/mantel/master/doc/index.html

See [mantel.file()](https://rawgit.com/kevinstadler/mantel/master/doc/mantel.file.html) for straightforward computation of Mantel tests from an input file. While the command is customisable, it generally helps if the column order of your `.csv` files list the meaning specifications *before* the corresponding strings.

## Links

* [A guide to the Mantel test for linguists](http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists)
* [Mantel (1967) The Detection of Disease Clustering and a Generalized Regression Approach](http://cancerres.aacrjournals.org/content/27/2_Part_1/209.short)

