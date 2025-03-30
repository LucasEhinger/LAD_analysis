#!/bin/csh

# Call the alias golad
golad

# Change directory to ../lad_replay
cd ../lad_replay

# Loop through run numbers 343 to 356
foreach run_num (343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364)
  # Call hcana with the specified script and run number
  hcana -l -q "SCRIPTS/LAD/PRODUCTION/replay_production_lad.C(${run_num},-1)"
end