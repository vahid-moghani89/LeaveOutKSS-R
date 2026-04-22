Data folder guide
=================

The small package examples use the dataset bundled with `LeaveOutKSS`:

  system.file("extdata", "test.csv", package = "LeaveOutKSS")

Local files used by examples:

1) lincom.csv
   Expected columns:
     V1 = id
     V2 = firmid
     V3 = year (not used)
     V4 = region indicator (-1 or 1)
     V5 = y (outcome)
   Used by: examples/03_lincom_example.R

2) large_fake_data.csv (optional, not tracked)
   Expected columns:
     col1 = id
     col2 = firmid
     col3 = year
     col4 = y
   Used by: examples/04_large_no_controls.R
