import os
from array import *
import argparse

import numpy as np


# python massjobs_ddsim_angle.py -g FSCEPonly -p e-

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")
argParser.add_argument("-p", "--particle", help="particle gun")
#argParser.add_argument("-d", "--gdir",     help="gun direction")


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
if not os.path.exists('jobs/' + args.geometry+'/exefiles/'):
   os.makedirs('jobs/' + args.geometry+'/exefiles/')
   os.makedirs('jobs/' + args.geometry+'/stdfiles/')

if not os.path.exists('output/' + args.geometry+'/angles/'):
    os.makedirs('output/' + args.geometry+'/angles/')
    os.makedirs('jobs/' + args.geometry+'/angles/')

outputarea = os.getcwd() + '/output/' + args.geometry + '/angles/'
exearea = os.getcwd() + '/jobs/' + args.geometry + '/exefiles/'
stdarea = os.getcwd() + '/jobs/' + args.geometry + '/stdfiles/'

name="condor-executable-"+args.geometry+'_'
angles = [10,15,20,25,30,35,40,45,50]
# create the .sh files
for i in angles:
    gdir = np.tan(np.radians(i))
    shfile = open(exearea+name+args.particle+str(i)+'degrees.sh',"w")
    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd '+ parent_dir + '\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('echo "{args.geometry}, {args.particle} gun, {i} direction"'+'\n')
    # getting centos version
    shfile.write('. /etc/os-release' + '\n')
    shfile.write('echo "machine is centos${VERSION_ID%.*}"' + '\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos${VERSION_ID%.*}-gcc11-opt/setup.sh' + '\n')
    shfile.write('source '+source_dir+'/install/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('process_id=$1'+'\n')
    shfile.write('echo "process $process_id"'+'\n')
    shfile.write('echo "ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_'+args.particle+str(i)+'degrees_$process_id.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 '+str(gdir)+' 1." --gun.energy "50*GeV" --gun.particle="'+args.particle+'" -N 10 >& '+outputarea+'Log_'+args.geometry+'_'+args.particle+str(i)+'degrees_$process_id.log"'+'\n')
    shfile.write('ddsim --compactFile='+parent_dir+'/DR'+args.geometry+'.xml --runType=batch -G --steeringFile '+parent_dir+'/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_'+args.particle+str(i)+'degrees_$process_id.root --part.userParticleHandler='' -G --gun.position="0.,0*mm,-1*cm" --gun.direction "0 '+ str(gdir) + ' 1." --gun.energy "50GeV" --gun.particle="'+args.particle+'" -N 10 >& '+outputarea+'Log_'+args.geometry+'_'+args.particle+str(i)+'degrees_$process_id.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')
    shfile.close()
print("sh file closed")

# create the .jdl files
for i in angles:
    jdlfile = open(exearea+name+args.particle+str(i)+'degrees.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+exearea+name+args.particle+str(i)+"gev.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    #jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Requirements = (machine == \"r720-0-1.privnet\") || (machine == \"hepcms-namenode.privnet\")"+'\n') #alternative req. for hepcms cluster
    jdlfile.write("Output = "+stdarea+args.particle+"$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+stdarea+args.particle+"$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+stdarea+args.particle+"$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = $(process)"+'\n')
    jdlfile.write("request_memory = 4GB"+'\n')
    jdlfile.write('Queue 12' + '\n') # run jobs in parallel
    jdlfile.close()
print("jdl file closed")

# create the submitter file
f = open("massjobs_ddsim_angles.sh",'w')
f.write('chmod +x '+exearea+'*'+'\n')
for i in angles:
    f.write("condor_submit "+exearea+name+args.particle+str(i)+'degrees.jdl'+'\n')
f.write("condor_q"+'\n')
f.close()
