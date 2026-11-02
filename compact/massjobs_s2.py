from array import *
import argparse


# python massjobs_s2.py
# -g1 DualTestBeam -g2 FSCEPonly
# -doecal 1 -dohcal 1 -ecaltype 0 -hcaltype 0 -doedge 1 -t smalldial -timecut 9999999. -gendete 3 -gendeth 3 -dodualem 2 -dodualhad 2 -ecalkappa 2 -layermode 0 -NNN 2
# -rmssc -0.0127  -rmsic 0.984  -rmsss -0.008  -rmsis 1.056
# -df 0 -of 5
# -dh 2 -oh 3
# -c /data/users/eno/CalVision/dd4hep/gitlab/DualTestBeam/compact/
# -w /data/users/eno/CalVision/dd4hep/gitlab/DualTestBeam/compact/output/
# -s /data/users/eno/CalVision/dd4hep/gitlab/DualTestBeam/compact/jobs/

argParser = argparse.ArgumentParser()
argParser.add_argument("-w", "--write", help="where to write")
argParser.add_argument("-c", "--cdarea", help="where to run (compact area)")
argParser.add_argument("-s", "--script", help="where to put scripts")
argParser.add_argument("-g1", "--geometry1", help="main geometry ")
argParser.add_argument("-g2", "--geometry2", help="hcal calibration geometry ")
argParser.add_argument("-t", "--steering", help="choose steering file ")
argParser.add_argument("-doecal","--doecal", help="0 no 1 yes")
argParser.add_argument("-dohcal","--dohcal", help="0 no 1 yes")
argParser.add_argument("-hcaltype","--hcaltype", help="0 fiber 1 sampling")
argParser.add_argument("-ecaltype","--ecaltype", help="0 2 layer crystal 1 crilin")
argParser.add_argument("-doedge","--doedge", help="0 no 1 yes")
argParser.add_argument("-timecut","--timecut", help="neg is no time cut unit ns")
argParser.add_argument("-dodualem","--dodualem", help="0 scale set by electrons 1 set by pions 2 dual correction")
argParser.add_argument("-dodualhad","--dodualhad", help="0 scale set by electrons 1 set by pions 2 dualcorrection")
argParser.add_argument("-NNN","--NNN", help="number to process")
argParser.add_argument("-gendete","--gendete", help="ecal 1 active media photons 2 photodetector  3 deposited energies 5 digit")
argParser.add_argument("-gendeth","--gendeth", help="hcal 1 active media photons 2 photodetector  3 deposited energies")
argParser.add_argument("-ecalkappa","--ecalkappa", help="ecal kappa mode see code")
argParser.add_argument("-df", "--directionf", help="0 is straight, 1 is 0.05, 2 is 0.0125 full detector")
argParser.add_argument("-of", "--originf", help="0 is 0,0,-1  1 is 0,-7,-1  2 is 0,0,-80 3 is 0,-1.5,-80 full detector")
argParser.add_argument("-dh", "--directionh", help="0 is straight, 1 is 0.05, 2 is 0.0125 hcal only")
argParser.add_argument("-oh", "--originh", help="0 is 0,0,-1  1 is 0,-7,-1  2 is 0,0,-80 3 is 0,-1.5,-80 hcal only")
argParser.add_argument("-layermode","--layermode", help="0 all 1 last")
argParser.add_argument("-rmssc","--rmssc", help="rms corrected dual parm")
argParser.add_argument("-rmsic","--rmsic", help="rms corrected dual parm")
argParser.add_argument("-rmsss","--rmsss", help="rms corrected dual parm")
argParser.add_argument("-rmsis","--rmsis", help="rms corrected dual parm")


args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry1)
print("args.name=%s" % args.geometry2)
print("args.name=%s" % args.steering)


outputarea=args.write
hostarea=args.cdarea
scriptarea=args.script

print("output area is %s" % outputarea )
print("host area is %s" % hostarea )
print("script area is %s" % scriptarea )




#nenergy=9
#energies=[10,15,20,25,30,35,40,45,50]
nenergy=1
energies=[20]
name=str(args.geometry1)+"_"+str(args.geometry2)+'_steering_'+str(args.steering)+'_timecut_'+str(args.timecut)+'_g'+str(args.gendete)+"-"+str(args.gendeth)+'_dual_'+str(args.dodualem)+"-"+str(args.dodualhad)+'_kappa_'+str(args.ecalkappa)+"_layermode_"+str(args.layermode)+"-directionf-"+str(args.directionf)+"-originf-"+str(args.originf)+"-directionh-"+str(args.directionh)+"-originh-"+str(args.originh)
name2="s2-condor-executable-"+name




# create the .sh files 
i=0
while (i<nenergy):
    print(i)
    aname1 = str(args.steering)+'_direction_'+str(args.directionf)+'_origin_'+str(args.originf)+'_energy_'+str(energies[i])+'GeV'
    aname2 = str(args.steering)+'_direction_'+str(args.directionh)+'_origin_'+str(args.originh)+'_energy_'+str(energies[i])+'GeV'
    aaaname = str(args.steering)+'_timecut_'+str(args.timecut)+'_directionf_'+str(args.directionf)+'_originf_'+str(args.originf)+'_directionh_'+str(args.directionh)+'_originh_'+str(args.originh)+'_dualem_'+str(args.dodualem)+'_dualhad_'+str(args.dodualhad)+'_kappaecal_'+str(args.ecalkappa)+'_layermode_'+str(args.layermode)+'_gendete_'+str(args.gendete)+'_gendeth_'+str(args.gendeth)+'_energy_'+str(energies[i])+'GeV'
    shfile = open(scriptarea+name2+str(energies[i])+'_GeV.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd '+hostarea+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('root -b -l -q \'Resolution.C('+str(args.NNN)+',"./output/out_'+str(args.geometry1)+'_'+str(aname1)+'_e-.root","./output/out_'+str(args.geometry1)+'_'+str(aname1)+'_pi-.root","./output/out_'+str(args.geometry2)+'_'+str(aname2)+'_e-.root","./output/out_'+str(args.geometry2)+'_'+str(aname2)+'_pi-.root",'+str(energies[i])+','+str(args.doecal)+','+str(args.dohcal)+','+str(args.ecaltype)+','+str(args.hcaltype)+','+str(args.doedge)+','+str(args.timecut)+','+str(args.gendete)+','+str(args.gendeth)+',"./output/hists_'+str(args.geometry1)+'_'+str(args.geometry2)+"_"+aaaname+'.root",0,'+str(args.dodualem)+','+str(args.dodualhad)+','+str(args.ecalkappa)+','+str(args.layermode)+','+str(args.rmssc)+','+str(args.rmsic)+','+str(args.rmsss)+','+str(args.rmsis)+')\' >& ./Figures/s2_'+str(args.geometry1)+'_'+str(args.geometry2)+'_'+aaaname+'.log \n' );
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
    jdlfile = open(scriptarea+name2+str(energies[i])+'_GeV.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+scriptarea+name2+str(energies[i])+"_GeV.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write('Requirements = machine == "hepcms-henrietta.privnet" '+'\n')
    #    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"r720\""+'\n')
    jdlfile.write("Output = "+outputarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+outputarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+outputarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")


# create the submitter file
f = open("massjobs_s2.sh",'w')
f.write('chmod 777 '+scriptarea+'*'+'\n')
i=0
while (i<nenergy):
    f.write("condor_submit "+scriptarea+name2+str(energies[i])+'_GeV.jdl'+'\n')
    i=i+1
f.write("condor_q")
f.close()



