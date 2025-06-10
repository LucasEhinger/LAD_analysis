#!/bin/bash

# Array of run numbers or file identifiers
run_numbers=(22609 22610 22611 22613 22614 22615)

# Loop over each run number and call the ROOT macro
for run in "${run_numbers[@]}"; do
  root -q "gem_tracking_by_clust.C(${run})"
done