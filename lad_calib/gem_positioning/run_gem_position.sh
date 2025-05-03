#!/bin/bash

# Call 'golad' command
golad

# Define an array of positions
positions=("pos1" "pos2" "pos3" "pos4") # Replace with actual position values

# Initialize the starting value for the second parameter
start_value=20000

# Loop over the positions array
for position in "${positions[@]}"; do
  # Call the hcana command with incremented parameter
  hcana -l -q "SCRIPTS/LAD_COIN/PRODUCTION/replay_production_lad_spec.C(22281,$start_value,0)"
  
  # Increment the start_value by 1
  ((start_value++))
  # Pause for 30 seconds
  sleep 30
done