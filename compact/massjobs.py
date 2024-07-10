import os
from array import *
import argparse
import numpy as np


# python massjobs.py -g FSCEPonly -p e-

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")
argParser.add_argument("-p", "--particle", help="particle gun")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry)


parent_dir = os.getcwd()
source_dir = os.getcwd()+'/../../..'
print(source_dir)
os.mkdir('output_1')
os.mkdir('jobs_1')
outputarea = os.getcwd() + '/output/'
hostarea = os.getcwd() + '/jobs/'

nenergy=3
energies=[10,15,20,25,30,35,40,45,50,100]
name="condor-executable-"+args.geometry+"-dial-"


# create the .sh files for electrons
shfile = open(hostarea+name+args.particle+'.sh',"w")
shfile.write('#!/bin/bash'+'\n')
shfile.write('cd '+ parent_dir + '\n')
shfile.write('START_TIME=`/bin/date`'+'\n')
shfile.write('echo "started at $START_TIME"'+'\n')
shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh'+'\n')
shfile.write('source '+source_dir+'/bin/thisdd4hep.sh')
shfile.write('echo "ran setup"'+'\n')
shfile.write('source ' + parent_dir + '\n')
shfile.write('echo "ran thisdd4hep"'+'\n')

for i in energies[:3]:
  shfile.write('ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'-dial_'+args.particle+'.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 0.0 1." --gun.energy "'+str(i)+'*GeV" --gun.particle="'+args.particle+'" -N 50 >& '+outputarea+'sce_e_'+args.geometry+'-dial_'+str(i)+'.log'+'\n')
shfile.write('exitcode=$?'+'\n')
shfile.write('echo ""'+'\n')
shfile.write('END_TIME=`/bin/date`'+'\n')
shfile.write('echo "finished at $END_TIME"'+'\n')
shfile.write('exit $exitcode'+'\n')
shfile.close()
print("file closed")


# create the .jdl files for electrons
jdlfile = open(hostarea+name+args.particle+'.jdl',"w")
jdlfile.write("universe = vanilla"+'\n')
jdlfile.write("Executable ="+hostarea+name+args.particle+".sh"+'\n')
jdlfile.write("should_transfer_files = NO"+'\n')
jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
jdlfile.write("Output = "+hostarea+name+args.particle+"_sce_$(cluster)_$(process).stdout"+'\n')
jdlfile.write("Error = "+hostarea+name+args.particle+"_sce_$(cluster)_$(process).stderr"+'\n')
jdlfile.write("Log = "+hostarea+name+args.particle+"_sce_$(cluster)_$(process).condor"+'\n')
jdlfile.write("Arguments = $(process)"+'\n')
jdlfile.write('Queue '+str(nenergy) + '\n')
#for i in energies[:3]: jdlfile.write(str(i)+ ' ')
#print("file closed")
jdlfile.close()


# create the submitter file
f = open("massjobs.sh",'w')
f.write('chmod +x '+hostarea+'*'+'\n')
f.write("condor_submit "+hostarea+name+args.particle+'.jdl'+'\n')
f.write("condor_q")
f.close()



