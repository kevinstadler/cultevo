#!/bin/sh
Rscript -e 'roxygen2::roxygenise();pkgdown::build_site(run_dont_run=TRUE, preview=FALSE);devtools::check();devtools::build()'
