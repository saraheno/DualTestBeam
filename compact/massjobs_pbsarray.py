from array import *
import argparse

# python massjobs.py -g FSCEPonly

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry)

# Check your working area
basedir="/cms/data/hatake/ana/CalVision/DD4hep_102b/DD4hep/"

# Various directories
compactdir=basedir+"/examples/DualTestBeam/compact/"
outputarea=basedir+"/examples/DualTestBeam/compact/output/"
hostarea=basedir+"/examples/DualTestBeam/compact/jobs/"

# Directions & position
#direction="0 0.176 1."
#position="0.,0*mm,-1*cm"
direction="0 0.05 0.99875"
position="0.,-7*mm,-1*cm"
#direction="0 0.0 1."
#position="0.,0*mm,-1*cm"

# LCG setup script
LCGsetup="/cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos8-gcc11-opt/setup.sh"

nenergy=10
energies=[10,15,20,25,30,35,40,45,50,100]
name="pbs-executable-"+args.geometry+"-dial-"

# create the .sh files for electrons
i=0
while (i<nenergy):
    print(i)
    shfile = open(hostarea+name+str(energies[i])+'_GeV-e.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd '+compactdir+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source '+LCGsetup+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  '+basedir+'/install/bin/thisdd4hep.sh'+'\n')
    shfile.write('module load gl_fix' + '\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    shfile.write('mkdir -p '+outputarea+'\n')
    shfile.write('ddsim --compactFile='+compactdir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+compactdir+'/SCEPCALsteering.py --outputFile='+outputarea+'/out_'+args.geometry+'-dial_'+str(energies[i])+'GeV_e-.${PBS_ARRAY_INDEX}.root --part.userParticleHandler='' -G --gun.position="'+position+'" --gun.direction "'+direction+'" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="e-" -N 50 >& '+outputarea+'sce_e_'+args.geometry+'-dial_'+str(energies[i])+'.${PBS_ARRAY_INDEX}.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")

# create the .sh files for pions
i=0
while (i<nenergy):
    print(i)
    shfile = open(hostarea+name+str(energies[i])+'_GeV-pi.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd '+compactdir+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source '+LCGsetup+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  '+basedir+'/install/bin/thisdd4hep.sh'+'\n')
    shfile.write('module load gl_fix' + '\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    shfile.write('mkdir -p '+outputarea+'\n')
    shfile.write('ddsim --compactFile='+compactdir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+compactdir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'-dial_'+str(energies[i])+'GeV_pi-.${PBS_ARRAY_INDEX}.root --part.userParticleHandler='' -G --gun.position="'+position+'" --gun.direction "'+direction+'" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="pi-" -N 50 >& '+outputarea+'sce_pi_'+args.geometry+'-dial_'+str(energies[i])+'.${PBS_ARRAY_INDEX}.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")
