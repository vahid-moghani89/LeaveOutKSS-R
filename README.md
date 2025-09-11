LeaveOutKSS-R
========================

What this is
------------
R translation of the leave-out correction (KSS) for variance components in
two-way fixed effects (AKM). This repo mirrors my MATLAB workflow but in R.
Goal: make it easy to run on large admin-style panels.

What you get (main pieces)
--------------------------
- functions/: all core functions (I source them, not an R package build)
- examples/: runnable .R scripts that source functions then run end-to-end
- data/: see data/README.txt

How to run (quick)
------------------
1) Put your CSV(s) in ./data
2) Open an example from ./examples and run it top to bottom
   (each example sources ./functions first)

Examples
--------
01_basic_no_controls.R
    Minimal run (no controls), prints decomposition, writes a CSV of pe/fe.
02_with_controls_year_FE.R
    Adds year dummies as controls (partialling-out), then runs KSS.
03_lincom_example.R
    After estimating FE/PE, regresses firm FE on a Z (e.g., region).
04_large_no_controls.R
    Large dataset runner. Downloads or reads a big CSV, uses JLA and parallel.

Notes (for myself)
------------------
- Leave-out granularity: "matches" is my default (robust to within-match serial correlation).
- For big n: use type_algorithm="JLA" with simulations_JLA ~ 50–200.
- Parallel: set paral=TRUE in leave_out_KSS. foreach/doParallel must be available.
- Outputs: CSV file with y_original, y_control_adjusted (if controls), id, firm_id, pe, fe.

Provenance
----------
Method: Kline–Saggio–Sølvsten (Econometrica, 2020).
This R code is my own translation; link to the original MATLAB package in CITATION.txt.

License
-------
MIT for my R code. See LICENSE.txt.
