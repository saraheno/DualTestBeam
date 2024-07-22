import os
from array import *
import argparse

# python massjobs_s2.py -g1 DualTestBeam -g2 FSCEPonly -p e-

argParser = argparse.ArgumentParser()
argParser.add_argument("-g1", "--geometry1", help="main geometry ")
argParser.add_argument("-g2", "--geometry2", help="hcal calibration geometry ")
argParser.add_argument("-p", "--particle", help="particle gun")

args = argParser.parse_args()
print("args=%s" % args)
print("args.name=%s" % args.geometry1)
print("args.name=%s" % args.geometry2)

parent_dir = os.getcwd()
source_dir = os.getcwd()+'/../../..'
lis = os.listdir()
if len([i for i in lis if i == 'output_s2'])==0: #check if the output and jobs directory are made before
  os.mkdir('output_s2')
  os.mkdir('jobs_s2')
outputarea = os.getcwd() + '/output_s2/'
hostarea = os.getcwd() + '/jobs_s2/'

nenergy=3
energies=[10,15,20,25,30,35,40,45,50,100]
name="condor-executable-"+args.geometry1+"_"+args.geometry2+"-"

# create the .sh files 
shfile = open(hostarea+name+args.particle+'_GeV.sh',"w")

shfile.write('#!/bin/bash'+'\n')
shfile.write('cd '+parent_dir+'\n')
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
   shfile.write('root -b -l -q \'Resolution.C(500,"./output/out_'+args.geometry1+"-dial_"+str(energies[i])+'GeV_'+args.particle+'.root","./output/out_'+args.geometry1+"_"+str(energies[i])+'GeV_pi-.root","./output/out_'+args.geometry2+"_"+str(energies[i])+'GeV_e-.root",'+str(energies[i])+',0,1,1,1,3,"./output/hists_'+args.geometry1+"_"+str(energies[i])+'GeV_3.root","DRSNoSegment","DRSNoSegment")\' >& ./output/s2_'+str(energies[i])+'GeV.log \n' );
shfile.write('exitcode=$?'+'\n')
shfile.write('echo ""'+'\n')
shfile.write('END_TIME=`/bin/date`'+'\n')
shfile.write('echo "finished at $END_TIME"'+'\n')
shfile.write('exit $exitcode'+'\n')
shfile.close()

# create the .jdl files 
jdlfile = open(hostarea+name+str(energies[i])+'_GeV.jdl',"w")
jdlfile.write("universe = vanilla"+'\n')
jdlfile.write("Executable ="+hostarea+name+args.particle+"_GeV.sh"+'\n')
jdlfile.write("should_transfer_files = NO"+'\n')
jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
jdlfile.write("Output = "+hostarea+name+args.particle+"_$(cluster)_$(process).stdout"+'\n')
jdlfile.write("Error = "+hostarea+name+args.particle+"_$(cluster)_$(process).stderr"+'\n')
jdlfile.write("Log = "+hostarea+name+args.particle+"_$(cluster)_$(process).condor"+'\n')
jdlfile.write("Arguments = SCE"+'\n')
jdlfile.write("Queue 1"+'\n')
jdlfile.close()


# create the submitter file
f = open("massjobs_s2.sh",'w')
f.write('chmod 777 '+hostarea+'*'+'\n')
f.write("condor_submit "+hostarea+name+args.particle+'_GeV.jdl'+'\n')
f.write("condor_q")
f.close()
