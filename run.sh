#!/bin/sh

# Run using:

# singularity exec /home/iwatson/Images/CCMadgraph.img make
# singularity exec /home/iwatson/Images/CCMadgraph.img ./run.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/code/MG5_aMC_v2_6_0/Delphes


#for f in /home/scratch/tsW/*t*; do
for f in /home/jang00747/tsW/*k.root*; do
    OUT="wj_out_`basename $f`"
    if [[ $f -nt result/`basename $f` ]]; then
        echo "Processing " $f
        ./wj_analysis $f result/$OUT
    else
        echo "Not processing" $f
    fi
done

#./wj_analysis /home/jang00747/tsW/tsWbW_1k.root /home/jang00747/tsW/wj_out_tsWbW_1k.root

