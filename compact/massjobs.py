from array import *
import argparse


# python massjobs.py -g DualTestBeam -N 500 -d 0 -o 2

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")


argParser.add_argument("-N", "--number", help="number to make")
argParser.add_argument("-d", "--direction", help="0 is straight 1 is fiber angle")
argParser.add_argument("-o", "--origin", help="0 is 0,0,-1  1 is 0,-7,-1  2 is 0,0,-80")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry)



outputarea="/data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/output/"
hostarea="/data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/jobs/"





#nenergy=9
#energies=[10,15,20,25,30,35,40,45,50]
nenergy=1
energies=[20]
name="condor-executable-"+str(args.geometry)+"-"
direct="0. 0.0 1."
poss=" 0. 0. -80.*cm"

print(args.direction)
if args.direction=="0" :
    direct="0. 0. 1."
if args.direction=="1" :
    direct="0. 0.05 0.99875"


print(args.origin)

if args.origin=="0" :
    poss="0. 0.*mm -1*cm"
if args.origin=="1" :
    poss="0.,-7*mm,-1*mm"
if args.origin=="2" :
    poss="0. 0.*mm -80*cm"


print(direct)
print(poss)

# create the .sh files for electrons
i=0
while (i<nenergy):
    print(i)
    shfile = open(hostarea+name+str(energies[i])+'_GeV-e.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd /data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/'+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('singularity run -B /cvmfs:/cvmfs -B /data:/data docker://gitlab-registry.cern.ch/sft/docker/alma9-core:latest'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/install/bin/thisDualTestBeam.sh'+'\n')
    shfile.write('echo "ran thisDualTestBeam"'+'\n')
# another good direction is  "0 0.05 0.99875"  and position 0.,-7*mm,-1*cm use this for pure fiber
# DO IT BOTH PLACES!!!
    shfile.write('ddsim --compactFile=/home/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/DR'
                 +str(args.geometry)+'.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/SCEPCALsteering.py --outputFile='+outputarea+'out_'+str(args.geometry)+'_'+str(energies[i])+'GeV_e-.root --part.userParticleHandler='' -G --gun.position="'+poss+'" --gun.direction "'+direct+'" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="e-" -N '+str(args.number)+' >& '+outputarea+'sce_e_'+str(args.geometry)+'_'+str(energies[i])+'.log'+'\n')
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
    shfile.write('cd /data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/'+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('singularity run -B /cvmfs:/cvmfs -B /data:/data docker://gitlab-registry.cern.ch/sft/docker/alma9-core:latest'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/install/bin/thisDualTestBeam.sh'+'\n')
    shfile.write('echo "ran thisDualTestBeam"'+'\n')
    shfile.write('ddsim --compactFile=/home/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/DR'+str(args.geometry)+'.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/stuff4stuff/DualTestBeam/compact/SCEPCALsteering.py --outputFile='+outputarea+'out_'+str(args.geometry)+'_'+str(energies[i])+'GeV_pi-.root --part.userParticleHandler='' -G --gun.position="'+poss+'" --gun.direction "'+direct+'" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="pi-" -N '+str(args.number)+' >& '+outputarea+'sce_pi_'+str(args.geometry)+'_'+str(energies[i])+'.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")


# create the .jdl files for electrons
i=0
while (i<nenergy):
    print(i)
    jdlfile = open(hostarea+name+str(energies[i])+'-e.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+hostarea+name+str(energies[i])+"_GeV-e.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = machine == \"hepcms-henrietta.privnet\""+'\n')
#    jdlfile.write("request_memory = 15GB"+'\n')
#    jdlfile.write("RequestCpus = 4"+'\n')
    jdlfile.write("Output = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")

# create the .jdl files for pions
i=0
while (i<nenergy):
    print(i)
    jdlfile = open(hostarea+name+str(energies[i])+'-pi.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+hostarea+name+str(energies[i])+"_GeV-pi.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = machine == \"hepcms-henrietta.privnet\""+'\n')
#    jdlfile.write("request_memory = 15GB"+'\n')
#    jdlfile.write("RequestCpus = 4"+'\n')
    jdlfile.write("Output = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")



# create the submitter file
f = open("massjobs.sh",'w')
f.write('chmod 777 '+hostarea+'*'+'\n')
i=0
while (i<nenergy):
    f.write("condor_submit "+hostarea+name+str(energies[i])+'-e.jdl'+'\n')
    f.write("condor_submit "+hostarea+name+str(energies[i])+'-pi.jdl'+'\n')
    i=i+1
f.write("condor_q")
f.close()



