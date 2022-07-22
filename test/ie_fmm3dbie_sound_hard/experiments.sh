#!/bin/bash

# Constant parameters
ZK=1
RTOL=0.00005
M=50

# Experiment matrix
NPUS=(10 20 30 40)
NORDERS=(4)
OCCS=(500)

# Construct and submit a job from experiment matrix
for NORDER in ${NORDERS[@]}
do
    for NPU in ${NPUS[@]}
    do
        for OCC in ${OCCS[@]}
        do
        scommand="sbatch -J NORDER=$NORDER.NPU=${NPU}.OCC=${OCC}.run "
        scommand+="--export=NORDER=${NORDER},NPU=${NPU},OCC=${OCC},ZK=${ZK},RTOL=${RTOL},M=${M} "
        scommand+="--output=NORDER_${NORDER}_NPU_${NPU}_OCC_${OCC}.out "
        scommand+="ie_fmm3dbie_sound_hard.slurm"
        echo "submitting: $scommand"
        echo ""
        $scommand
        done
    done
done