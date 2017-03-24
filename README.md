# Cultural evolution tests for R

This package provides implementations of functions and statistical tests often used in the study of cultural evolution. It provides tools for:

* computing edit distances for strings and atrices (providing convenient wrappers around R's efficient `dist` functions)
* performing Mantel tests on them and visualising the results of those tests
* an improved implementation of Page's L test for monotonicity (with exact p values up to k=22)

Full function documentation is available at <http://kevinstadler.github.io/cultevo/reference/>

## Installation

In order to install the latest version you first need the [devtools](https://CRAN.R-project.org/package=devtools) package

    install.packages("devtools")

then install the latest code from the github repository via

    devtools::install_github("kevinstadler/cultevo")


## Links

* [A guide to the Mantel test for linguists](http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists)
* [Mantel (1967) The Detection of Disease Clustering and a Generalized Regression Approach](http://cancerres.aacrjournals.org/content/27/2_Part_1/209.short)

## License

This project is licensed under the terms of the [MIT license](http://opensource.org/licenses/MIT), Copyright (c) 2014-2017 Kevin Stadler.
