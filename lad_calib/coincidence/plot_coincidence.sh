#!/bin/bash

# Check if run numbers are provided as arguments
if [ "$#" -eq 0 ]; then
  # Default list of run numbers if none are provided
  run_numbers=("1001" "1002" "1003") # Replace with your desired default run numbers
else
  # Use provided run numbers
  run_numbers=("$@")
fi

# Loop through each run number
for run_number in "${run_numbers[@]}"; do
  # Assign arguments to variables
  inputFileName="/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/CALIB/LAD_COIN_${run_number}_-1.root"
  outputPdfName="coin_histos_${run_number}.pdf"
  plotsPerPage="1" # Default to 4 if not provided

  # Run the ROOT command to execute the plot_coincidence function
  root -l -b -q "plot_coincidence.C(\"${inputFileName}\", \"${outputPdfName}\", ${plotsPerPage})"
done