#!/bin/csh
cd /work/hallc/c-lad/ehingerl/analysis/lad_cosmic/general/scripts

# Loop through run numbers 343 to 356
foreach run_num (344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364)
  # Call hcana with the specified script and run number
  root -l -q "cosmic_histos_hall.C(${run_num})"
end