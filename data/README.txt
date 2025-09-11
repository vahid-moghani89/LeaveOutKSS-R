Data folder guide
=================

Files used by examples:

1) test.csv  (small)
   Expected columns:
     V1 = id (worker)
     V2 = firmid
     V3 = year (for controls example)
     V4 = y (outcome)
   Used by: 01_basic_no_controls.R, 02_with_controls_year_FE.R

2) lincom.csv (small)
   Expected columns:
     V1 = id
     V2 = firmid
     V3 = year (not used)
     V4 = region indicator (-1 or 1)
     V5 = y (outcome)
   Used by: 03_lincom_example.R

3) large_fake_data.csv (large) -- Let examples/04_large_no_controls.R download it.
   Expected columns:
     col1 = id
     col2 = firmid
     col3 = year
     col4 = y
