# LeaveOutKSS-R

## Overview

This repository provides example scripts, local data notes, and replication
materials for the `LeaveOutKSS` R package.

The package is available on CRAN:  
https://cran.r-project.org/package=LeaveOutKSS

`LeaveOutKSS` implements leave-out correction methods for variance components
in two-way fixed effects models, following:

> Kline, P., Saggio, R., and Sølvsten, M. (2020), *Leave-Out Estimation of Variance Components*, Econometrica.

The R implementation is conceptually aligned with the original MATLAB package:
https://github.com/rsaggio87/LeaveOutTwoWay


---

## Repository Contents

- `examples/`  
  Runnable R scripts demonstrating core workflows of `LeaveOutKSS`.  
  Several examples are adapted from, or inspired by, the original MATLAB
  repository to facilitate comparison across implementations.
- `data/`: local data notes and auxiliary data for selected examples
- `package/`: package workspace for `LeaveOutKSS` development

The examples are written for the CRAN release of `LeaveOutKSS`, so results can
be tied to the same package version used in analysis and citation.

Installation
------------
Install the package from CRAN:

```r
install.packages("LeaveOutKSS")
```

Run the examples from the repository root so relative paths resolve correctly.
For example:

```r
source("examples/01_basic_no_controls.R")
```

Examples
--------
- `01_basic_no_controls.R`
    Minimal leave-out KSS run using the small panel bundled with the package.
- `02_with_controls_year_FE.R`
    Adds year dummies as controls.
- `03_lincom_example.R`
    Projects estimated firm effects on a region indicator.
- `04_large_no_controls.R`
    Large local dataset runner using JLA and parallel computation.
- `05_with_absorbed_FE_function.R`
    Uses `leave_out_KSS_fe()` to absorb categorical controls internally.
- `06_rsquared_basic.R`
    Reports TWFE and saturated-model fit with `rsquared_comp()`.

Usage Notes
-----------
- Leave-out granularity: `"matches"` is the default for robustness to
  within-match serial correlation.
- For large datasets, use `type_algorithm = "JLA"` and tune
  `simulations_JLA`.
- Local outputs should be written to temporary folders or ignored paths.

Provenance
----------
Method: Kline, Saggio, and Solvsten (Econometrica, 2020).

Original MATLAB implementation:
https://github.com/rsaggio87/LeaveOutTwoWay

R package:
https://cran.r-project.org/web/packages/LeaveOutKSS/index.html

No original methodological contributions are introduced here; all credit for
the estimation framework belongs to the original authors.

License
-------
MIT for this repository's R examples and supporting files. See `LICENSE`.
