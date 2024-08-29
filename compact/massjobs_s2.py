from array import *
import argparse


# python massjobs_s2.py -g1 DualTestBeam -g2 FSCEPonly

argParser = argparse.ArgumentParser()
argParser.add_argument("-g1", "--geometry1", help="main geometry ")
argParser.add_argument("-g2", "--geometry2", help="hcal calibration geometry ")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry1)
print("args.name=%s" % args.geometry2)



outputarea="/data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/output/"
#hostarea="/data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/jobs/"
hostarea = "/data/users/snabili/"



nenergy=10
energies=[10,15,20,25,30,35,40,45,50,100]
name="s2-condor-executable-"+args.geometry1+"_"+args.geometry2+"-"

# create the .sh files 
i=0
while (i<nenergy):
    print(i)
    shfile = open(hostarea+name+str(energies[i])+'_GeV.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/'+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    shfile.write('root -b -l -q \'Resolution.C(500,"./output/out_'+args.geometry1+"_"+str(energies[i])+'GeV_e-.root","./output/out_'+args.geometry1+"_"+str(energies[i])+'GeV_pi-.root","./output/out_'+args.geometry2+"_"+str(energies[i])+'GeV_e-.root",'+str(energies[i])+',0,1,1,1,3,"./output/hists_'+args.geometry1+"_"+str(energies[i])+'GeV_3.root","DRSNoSegment","DRSNoSegment")\' >& ./output/s2_'+str(energies[i])+'GeV.log \n' );
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")



# create the .jdl files 
i=0
while (i<nenergy):
    print(i)
    jdlfile = open(hostarea+name+str(energies[i])+'_GeV.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+hostarea+name+str(energies[i])+"_GeV.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Output = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")


# create the submitter file
f = open("massjobs_s2.sh",'w')
f.write('chmod 777 '+hostarea+'*'+'\n')
i=0
while (i<nenergy):
    f.write("condor_submit "+hostarea+name+str(energies[i])+'_GeV.jdl'+'\n')
    i=i+1
f.write("condor_q")
f.close()




