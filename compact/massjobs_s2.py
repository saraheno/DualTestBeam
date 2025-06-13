from array import *
import argparse


# python massjobs_s2.py -g1 DualTestBeam -g2 FSCEPonly -doecal 1 -dohcal 1 -hcaltype 0 -doedge 1 -gendet 3 -dodual 0 -NNN 2 -dotwocalcor 1

argParser = argparse.ArgumentParser()
argParser.add_argument("-g1", "--geometry1", help="main geometry ")
argParser.add_argument("-g2", "--geometry2", help="hcal calibration geometry ")
argParser.add_argument("-doecal","--doecal", help="0 no 1 yes")
argParser.add_argument("-dohcal","--dohcal", help="0 no 1 yes")
argParser.add_argument("-hcaltype","--hcaltype", help="0 fiber 1 sampling")
argParser.add_argument("-doedge","--doedge", help="0 no 1 yes")
argParser.add_argument("-dodual","--dodual", help="0 no 1 yes")
argParser.add_argument("-dotwocalcor","--dotwocalcor", help="0 no 1 yes")
argParser.add_argument("-NNN","--NNN", help="number to process")
argParser.add_argument("-gendet","--gendet", help="1 active media photons 2 photodetector  3 deposited energies")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry1)
print("args.name=%s" % args.geometry2)



outputarea="/data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/output/"
hostarea="/data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/jobs/"




#nenergy=9
#energies=[10,15,20,25,30,35,40,45,50]
nenergy=1
energies=[20]
name="s2-condor-executable-"+str(args.geometry1)+"_"+str(args.geometry2)+'_g'+str(args.gendet)+"-"

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
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos8-gcc11-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/DD4hep/install/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    if args.hcaltype=='0':
       shfile.write('root -b -l -q \'Resolution.C('+str(args.NNN)+',"./output/out_'+args.geometry1+"_"+str(energies[i])+'GeV_e-.root","./output/out_'+str(args.geometry1)+"_"+str(energies[i])+'GeV_pi-.root","./output/out_'+str(args.geometry2)+"_"+str(energies[i])+'GeV_e-.root","./output/out_'+str(args.geometry2)+"_"+str(energies[i])+'GeV_pi-.root",'+str(energies[i])+','+str(args.doecal)+','+str(args.dohcal)+','+str(args.hcaltype)+','+str(args.doedge)+',0,0,0.,'+str(args.gendet)+',"./output/hists_'+str(energies[i])+'GeV_'+str(args.geometry1)+'_g'+str(args.gendet)+'.root","DRCNoSegment","DRFNoSegment",0,0,'+str(args.dodual)+','+str(args.dotwocalcor)+')\' >& ./Figures/s2_'+str(energies[i])+'GeV'+'_'+str(args.geometry1)+'_'+str(args.geometry2)+'_g'+str(args.gendet)+'.log \n' );
    if args.hcaltype=='1':
       shfile.write('root -b -l -q \'Resolution.C('+str(args.NNN)+',"./output/out_'+args.geometry1+"_"+str(energies[i])+'GeV_e-.root","./output/out_'+str(args.geometry1)+"_"+str(energies[i])+'GeV_pi-.root","./output/out_'+str(args.geometry2)+"_"+str(energies[i])+'GeV_e-.root","./output/out_'+str(args.geometry2)+"_"+str(energies[i])+'GeV_pi-.root",'+str(energies[i])+','+str(args.doecal)+','+str(args.dohcal)+','+str(args.hcaltype)+','+str(args.doedge)+',0,0,0.,'+str(args.gendet)+',"./output/hists_'+str(energies[i])+'GeV_'+str(args.geometry1)+'_g'+str(args.gendet)+'.root","DRCNoSegment","DRSNoSegment",0,0,'+str(args.dodual)+','+str(args.dotwocalcor)+')\' >& ./Figures/s2_'+str(energies[i])+'GeV'+'_'+str(args.geometry1)+'_'+str(args.geometry2)+'_g'+str(args.gendet)+'.log \n' );
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
#    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"r720\""+'\n')
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



