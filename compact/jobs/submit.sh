#!/bin/bash

for i in 10 15 20 25 30 35 40 45 50 100
do
    echo "$i GeV job submission"
    #qsub -N pbs-executable-FSCEPonly-dial-${i}_GeV-e.sh   submit.pbs
    #qsub -N pbs-executable-FSCEPonly-dial-${i}_GeV-pi.sh  submit.pbs
    #qsub -N pbs-executable-DualTestBeam-dial-${i}_GeV-e.sh   submit.pbs
    #qsub -N pbs-executable-DualTestBeam-dial-${i}_GeV-pi.sh  submit.pbs
    qsub -N pbs-executable-Conly-dial-${i}_GeV-e.sh   submit.pbs
    #qsub -N pbs-executable-Conly-dial-${i}_GeV-pi.sh  submit.pbs
done
