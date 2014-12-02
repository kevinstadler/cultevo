## Simple Mantel tests for R

This package provides functions for computing distance matrices and performing Mantel tests on them.

## Installation

    devtools::install_github("kevinstadler/mantel")

If this doesn't work you probably need to install the [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package first: `install.packages("devtools")`

## Usage

Call

    library(mantel)

and then fire away some Mantel tests! A nice HTML documentation is available at https://rawgit.com/kevinstadler/mantel/master/doc/index.html

See [mantel.file()](https://rawgit.com/kevinstadler/mantel/master/doc/mantel.file.html) for straightforward computation of Mantel tests from an input file. While the command is highly customisable, if you don't want to mess about with all its parameters it is simplest if you have *tab*-separated `.csv` files with the strings in the first column, followed by meaning specifications in the remaining columns.

If you already have the distance matrices computed and you just want to run the actual tests, have a look at [mantel.test()](https://rawgit.com/kevinstadler/mantel/master/doc/mantel.test.html). To visualise the outcome of one or more tests, simply call [plot()](https://rawgit.com/kevinstadler/mantel/master/doc/plot.mantel.html) on the object returned by the tests. If you want to inspect the raw data returned by the test more closely, have a look at the [mantel class](http://rawgit.com/kevinstadler/mantel/master/doc/mantel.html).

## Links

* [A guide to the Mantel test for linguists](http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists)
* [Mantel (1967) The Detection of Disease Clustering and a Generalized Regression Approach](http://cancerres.aacrjournals.org/content/27/2_Part_1/209.short)

## License

This project is licensed under the terms of the [MIT license](http://opensource.org/licenses/MIT), Copyright (c) 2014 Kevin Stadler.
