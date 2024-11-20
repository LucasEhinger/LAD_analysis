TREE=T
fileloc=../../trees/proton/full_detector/all_particles/
for file in ${fileloc}*; do
    INPUTFILE="$file"
    OUTPUTFOLDER=$(dirname "$file" | sed 's/trees/hists/')
    OUTPUTFILE="$OUTPUTFOLDER/$(basename "$file" | sed 's/g4bg/hist/')"
    if [ ! -d "$(dirname "$OUTPUTFILE")" ]; then
        mkdir -p "$(dirname "$OUTPUTFILE")"
    fi
    root -q -b 'full_sim_protons.cc("'${INPUTFILE}'","'${OUTPUTFILE}'","'${TREE}'")'
done