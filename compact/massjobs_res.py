import os
from array import *
import argparse

'''
 wraper code to run Resolution.C in condor; code set to use hcal fiber type
 to run the code with DualTestBeam, ecal+hcal with different energies:
#   python massjobs_res.py -g DualTestBeam -ho=1 -hc=1 -ec=1 -ed=1 -gd=3 -e=1
'''

argParser = argparse.ArgumentParser()
argParser.add_argument("-g",  "--geometry", help="main geo")
argParser.add_argument("-ho", "--hcalonly", help="e-file hcal-cal", type=int, default=0)
argParser.add_argument("-hc",  "--hcal",    help="hcal-leaf",       type=int, default=0)
argParser.add_argument("-ec",  "--ecal",    help="ecal-leaf",       type=int, default=0)
argParser.add_argument("-ed", "--edge",     help="edge detector",   type=int, default=1)
argParser.add_argument("-gd", "--gendet",   help="gen type",        type=int, default=3)
argParser.add_argument("-a", "--angle",     help="run angles",      type=int, default=0)
argParser.add_argument("-e", "--energy",    help="run energies",    type=int, default=0)
argParser.add_argument("-ht", "--hcaltype", help="hcal-type",       type=int, default=0) # fiber: hc=0, sampl: hc=1

args = argParser.parse_args()
print("args=%s" % args)
print("args.name=%s" % args.geometry)

parent_dir = os.getcwd()
source_dir = os.getcwd()+'/../../..'
if not os.path.exists('jobs/' + args.geometry+'/res'):
   os.makedirs('output/' + args.geometry + '/res')
   os.makedirs('jobs/' + args.geometry + '/res')

inputfilearea = os.getcwd() + '/output/' + args.geometry + '/hadd'
outputarea    = os.getcwd() + '/output/' + args.geometry + '/res/'
hostarea      = os.getcwd() + '/jobs/' + args.geometry + '/res/'

energies=[10,15,20,25,30,35,40,45,50,100]
#angles=[10,15,20,25,30,35,40,45,50]
angles=[5,7,9,12]
nenergy=len(energies)
nangles=len(angles)

einputfile = inputfilearea + '/out_' + args.geometry + '_e-'
pinputfile = inputfilearea + '/out_' + args.geometry + '_pi-'

hinputfile = ""

outfile = outputarea+'res_'+args.geometry+'_'
outlogfile = outputarea+'log_'+args.geometry+'_'
if args.ecal: ecaleaf = "DRCNoSegment"
else: ecaleaf = ""
if args.hcal: hcaleaf = "DRSNoSegment"
else: hcaleaf = ""

name=args.geometry+"_"

# create the .sh files 
if args.energy: shfile = open(hostarea+name+'GeV.sh',"w")
if args.angle:  shfile = open(hostarea+name+'degrees.sh',"w")
shfile.write('#!/bin/bash'+'\n')
shfile.write('cd '+parent_dir+'\n')
shfile.write('START_TIME=`/bin/date`'+'\n')
shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
# getting centos version
shfile.write('. /etc/os-release' + '\n')
shfile.write('echo "machine is centos${VERSION_ID%.*}"' + '\n')
shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos${VERSION_ID%.*}-gcc11-opt/setup.sh' + '\n')
shfile.write('source '+source_dir+'/install/bin/thisdd4hep.sh'+'\n')
shfile.write('echo "ran setup"'+'\n')
if args.energy==1:
    for i in energies:
        if args.hcalonly==1: hinputfile = outputarea+'out_'+args.geometry+'_'+str(i)+'gev.root'
        shfile.write('echo root -b -l -q \'Resolution.C(500,"'+einputfile+str(i)+'gev.root","'+pinputfile+str(i)+'gev.root","'+hinputfile+'",'+str(i)+','+str(args.ecal)+","+str(args.hcal)+","+str(args.hcaltype)+","+str(args.edge)+","+str(args.gendet)+',"'+outfile+str(i)+'GeV.root","'+ecaleaf+'","'+hcaleaf+'")\' >&'+ outlogfile+str(i)+"GeV.log"+'\n')
        shfile.write('root -b -l -q \'Resolution.C(500,"'+einputfile+str(i)+'gev.root","'+pinputfile+str(i)+'gev.root","'+hinputfile+'",'+str(i)+','+str(args.ecal)+","+str(args.hcal)+","+str(args.hcaltype)+","+str(args.edge)+","+str(args.gendet)+',"'+outfile+str(i)+'GeV.root","'+ecaleaf+'","'+hcaleaf+'")\' >&'+ outlogfile+str(i)+"GeV.log"+'\n')
if args.angle==1:
    for a in angles:
        shfile.write('root -b -l -q \'Resolution.C(120,"'+einputfile+str(a)+'degrees.root","'+pinputfile+str(a)+'degrees.root","'+hinputfile+'",50,'+str(args.ecal)+","+str(args.hcal)+","+str(args.hcaltype)+","+str(args.edge)+","+str(args.gendet)+',"'+outfile+str(a)+'degrees.root","'+ecaleaf+'","'+hcaleaf+'")\' >&'+ outlogfile+str(a)+"degrees.log"+'\n')

shfile.write('exitcode=$?'+'\n')
shfile.write('echo ""'+'\n')
shfile.write('END_TIME=`/bin/date`'+'\n')
shfile.write('echo "finished at $END_TIME"'+'\n')
shfile.write('exit $exitcode'+'\n')
shfile.close()

# create the .jdl files 
if args.energy: jdlfile = open(hostarea+name+'GeV.jdl',"w")
if args.angle:  jdlfile = open(hostarea+name+'degrees.jdl',"w")
jdlfile.write("universe = vanilla"+'\n')
if args.energy: jdlfile.write("Executable ="+hostarea+name+"GeV.sh"+'\n')
if args.angle:  jdlfile.write("Executable ="+hostarea+name+"degrees.sh"+'\n')
jdlfile.write("should_transfer_files = NO"+'\n')
jdlfile.write("Requirements = (machine == \"r720-0-1.privnet\") || (machine == \"hepcms-namenode.privnet\")"+'\n') #alternative req. for hepcms cluster
jdlfile.write("Output = "+hostarea+name+"angles_$(cluster)_$(process).stdout"+'\n')
jdlfile.write("Error = "+hostarea+name+"angles_$(cluster)_$(process).stderr"+'\n')
jdlfile.write("Log = "+hostarea+name+"angles_$(cluster)_$(process).condor"+'\n')
jdlfile.write("Arguments = $(process)"+'\n')
if args.energy: jdlfile.write("Queue "+str(nenergy)+'\n')
if args.angle:  jdlfile.write("Queue "+str(nangles)+'\n')
jdlfile.close()


# create the submitter file
f = open("massjobs_res.sh",'w')
f.write('chmod 777 '+hostarea+'*'+'\n')
if args.energy: f.write("condor_submit "+hostarea+name+'GeV.jdl'+'\n')
if args.angle:  f.write("condor_submit "+hostarea+name+'degrees.jdl'+'\n')
f.write("condor_q")
f.close()
