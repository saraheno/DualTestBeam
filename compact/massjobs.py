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
if not os.path.exists('output'):
   os.makedirs('output/' + args.geometry)
   os.makedirs('jobs/' + args.geometry)
if not os.path.exists('output/' + args.geometry):
   os.makedirs('output/' + args.geometry)
   os.makedirs('jobs/' + args.geometry)

outputarea = os.getcwd() + '/output/' + args.geometry + '/'
hostarea = os.getcwd() + '/jobs/' + args.geometry + '/'

energies=[10,15,20,25,30,35,40,45,50,100]
name="condor-executable-"+args.geometry+"-dial-"

# create the .sh files
shfile = open(hostarea+name+args.particle+'.sh',"w")
shfile.write('#!/bin/bash'+'\n')
shfile.write('cd '+ parent_dir + '\n')
shfile.write('START_TIME=`/bin/date`'+'\n')
shfile.write('echo "started at $START_TIME"'+'\n')
shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
# getting centos version
shfile.write('. /etc/os-release' + '\n')
shfile.write('echo "machine is centos${VERSION_ID%.*}"' + '\n')
shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos${VERSION_ID%.*}-gcc11-opt/setup.sh' + '\n')
shfile.write('source '+source_dir+'/install/bin/thisdd4hep.sh'+'\n')
shfile.write('echo "ran setup"'+'\n')
for i in energies[:3]:
  shfile.write('ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'-dial_'+args.particle+str(i)+'_.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 0.0 1." --gun.energy "'+str(i)+'*GeV" --gun.particle="'+args.particle+'" -N 100 >& '+outputarea+'Log_'+args.geometry+'-dial_'+str(i)+'_'+args.particle+'.log'+'\n')
shfile.write('exitcode=$?'+'\n')
shfile.write('echo ""'+'\n')
shfile.write('END_TIME=`/bin/date`'+'\n')
shfile.write('echo "finished at $END_TIME"'+'\n')
shfile.write('exit $exitcode'+'\n')
shfile.close()
print("sh file closed")

# create the .jdl files
jdlfile = open(hostarea+name+args.particle+'.jdl',"w")
jdlfile.write("universe = vanilla"+'\n')
jdlfile.write("Executable ="+hostarea+name+args.particle+".sh"+'\n')
jdlfile.write("should_transfer_files = NO"+'\n')
jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
jdlfile.write("Output = "+hostarea+args.particle+"$(cluster)_$(process).stdout"+'\n')
jdlfile.write("Error = "+hostarea+args.particle+"$(cluster)_$(process).stderr"+'\n')
jdlfile.write("Log = "+hostarea+name+args.particle+"$(cluster)_$(process).condor"+'\n')
jdlfile.write("Arguments = $(process)"+'\n')
jdlfile.write('Queue '+str(len(energies[:3])) + '\n')
print("jdl file closed")
jdlfile.close()


# create the submitter file
f = open("massjobs.sh",'w')
f.write('chmod +x '+hostarea+'*'+'\n')
f.write("condor_submit "+hostarea+name+args.particle+'.jdl'+'\n')
f.write("condor_q"+'\n')
f.close()
