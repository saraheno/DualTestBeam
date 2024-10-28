import os
from array import *
import argparse
import numpy as np


# python massjobs_ddsim.py -g FSCEPonly -p e- -d 1

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")
argParser.add_argument("-p", "--particle", help="particle gun")
argParser.add_argument("-d", "--direction", help="0 is straight 1 is fiber angle")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry)


direct="0. 0.0 1."

print(args.direction)
if args.direction=="1" :
    direct="0. 0.05 0.99875"
poss="0. 0.*mm 0*cm"
if args.direction=="1" :
    poss="0.,-7*mm,-1*mm"


parent_dir = os.getcwd()
source_dir = os.getcwd()+'/../../..'

if not os.path.exists('output'):
   os.makedirs('output/' + args.geometry)
   os.makedirs('jobs/' + args.geometry)
if not os.path.exists('output/' + args.geometry):
   os.makedirs('output/' + args.geometry)
   os.makedirs('jobs/' + args.geometry)
if not os.path.exists('jobs/' + args.geometry+'/exefiles/'):
   os.makedirs('jobs/' + args.geometry+'/exefiles/')
   os.makedirs('jobs/' + args.geometry+'/stdfiles/')

outputarea = os.getcwd() + '/output/' + args.geometry + '/'
exearea = os.getcwd() + '/jobs/' + args.geometry + '/exefiles/'
stdarea = os.getcwd() + '/jobs/' + args.geometry + '/stdfiles/'

energies=[10,15,20,25,30,35,40,45,50,100]
#energies=[15,20,25,30,35,40,45,50,100]
name="condor-executable-"+args.geometry+'_'

# create the .sh files
for i in energies:
    shfile = open(exearea+name+args.particle+str(i)+'gev.sh',"w")
    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd '+ parent_dir + '\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('echo "{args.geometry}, {args.particle} gun, {i} GeV"'+'\n')
    # getting centos version
    shfile.write('. /etc/os-release' + '\n')
    shfile.write('echo "machine is centos${VERSION_ID%.*}"' + '\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos${VERSION_ID%.*}-gcc11-opt/setup.sh' + '\n')
    shfile.write('source '+source_dir+'/install/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('process_id=$1'+'\n')
    shfile.write('echo "process $process_id"'+'\n')
    shfile.write('echo "ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 0.176 1." --gun.energy "'+str(i)+'*GeV" --gun.particle="'+args.particle+'" -N 10 >& '+outputarea+'Log_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.log"'+'\n')
    #shfile.write('ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 0.176 1." --gun.energy "'+str(i)+'*GeV" --gun.particle="'+args.particle+'" -N 10 >& '+outputarea+'Log_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.log'+'\n')
    shfile.write('ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.root --part.userParticleHandler='' -G --gun.position="'+poss+'" --gun.direction "'+direct+'" --gun.energy "'+str(i)+'*GeV" --gun.particle="'+args.particle+'" -N 10 >& '+outputarea+'Log_'+args.geometry+'_'+args.particle+str(i)+'gev_$process_id.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')
    shfile.close()
print("sh file closed")

# create the .jdl files
for i in energies:
    jdlfile = open(exearea+name+args.particle+str(i)+'gev.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+exearea+name+args.particle+str(i)+"gev.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    #jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    #jdlfile.write("Requirements = (machine == \"hepcms-namenode.privnet\") || (machine == \"r540-0-20.privnet\") || (machine == \"r720-0-1.privnet\")"+'\n')
    jdlfile.write("Requirements = (machine == \"r540-0-20.privnet\") || (machine == \"r720-0-1.privnet\")"+'\n')
    #jdlfile.write("Requirements = (machine == \"r720-0-1.privnet\") || (machine == \"hepcms-namenode.privnet\")"+'\n') #alternative req. for hepcms cluster
    jdlfile.write("Output = "+stdarea+args.particle+"$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+stdarea+args.particle+"$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+stdarea+args.particle+"$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = $(process)"+'\n')
    jdlfile.write("request_memory = 4GB"+'\n')
    jdlfile.write('Queue 50' + '\n') # run jobs in parallel
    jdlfile.close()
print("jdl file closed")

# create the submitter file
f = open("massjobs_ddsim.sh",'w')
f.write('chmod +x '+exearea+'*'+'\n')
for i in energies:
    f.write("condor_submit "+exearea+name+args.particle+str(i)+'gev.jdl'+'\n')
f.write("condor_q"+'\n')
f.close()
