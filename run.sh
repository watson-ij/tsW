#!/bin/sh

# Run using:

# singularity exec /home/iwatson/Images/CCMadgraph.img make
# singularity exec /home/iwatson/Images/CCMadgraph.img ./run.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/code/MG5_aMC_v2_6_0/Delphes

#for f in /home/scratch/tsW/; do
#    if [[ $f -nt result/`basename $f` ]]; then
#        echo "Processing " $f
#        ./wj_test $f result/`basename $f`
#    else
#        echo "Not processing" $f
#    fi
#done

./wj_analysis /home/scratch/tsW/tsWbW_100k.root

