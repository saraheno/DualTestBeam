from array import *

hostarea="/data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/output/"



njobs=2
names=["condor-executable-sampling-20_e","condor-executable-sampling-20_pi"]

i=0
while (i<njobs):
    print i
    jdlfile = open(names[i]+'.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+names[i]+".sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Output = "+hostarea+names[i]+"_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+names[i]+"_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+names[i]+"_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print "file closed"


f = open("massjobs.sh",'w')
i=0
while (i<njobs):
    f.write("condor_submit "+names[i]+'.jdl'+'\n')
    i=i+1
f.write("condor_q")
f.close()



