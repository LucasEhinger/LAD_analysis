# THICKNEESS=5mm
# INPUTFILE=g4bg_out_${THICKNEESS}.root
# TREE=T
# OUTPUTFILE=survival_${THICKNEESS}.root

# root -q -b 'survival.cc("'${INPUTFILE}'","'${OUTPUTFILE}'","'${TREE}'")'

TREE=T
fileloc=../../trees/proton/full_detector/proton_only/
for file in ${fileloc}*; do
    INPUTFILE="$file"
    OUTPUTFILE="../../hists/proton/full_detector/proton_only/$(basename "$file" | sed 's/g4bg/hist/')"
    if [ ! -d "$(dirname "$OUTPUTFILE")" ]; then
        mkdir -p "$(dirname "$OUTPUTFILE")"
    fi
    root -q -b 'proton_momentum.cc("'${INPUTFILE}'","'${OUTPUTFILE}'","'${TREE}'")'
done
# INPUTFILE=../../trees/proton/full_detector/proton_only/g4bg_out_protons_TC_1mm_10_6.root
# TREE=T
# OUTPUTFILE=../../hists/proton/full_detector/proton_only/hists_g4bg_out_protons_1mm_TC_10_6.root

# root -q -b 'proton_momentum.cc("'${INPUTFILE}'","'${OUTPUTFILE}'","'${TREE}'")'