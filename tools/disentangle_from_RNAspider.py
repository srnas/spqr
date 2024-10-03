#script for disentangling RNA structures annotated with RNAspider using SPQR simulations and all-atom refinement.

import argparse
import os
import subprocess

NATNT=5
SPQR_PATH="SPQRPATH"
ENERGPATH=SPQR_PATH+"/interactions/intrac.btb"
BINPATH=SPQR_PATH+"/bin/"
CLASHSTEPS=5
def run_commands(command, where):
    commout=[]
    for ii in command:
        pp=subprocess.run(ii.split(),capture_output=True,text=True,cwd=where)
        commout.append(pp.stdout)
    return commout
def write_params(paramd, where):
    PARAMFILE=open(where+"/params.spqr", "w")
    for ii in paramd: PARAMFILE.write(ii+" "+paramd[ii]+"\n")
    PARAMFILE.close()
def rename_files(wdir,localdir,ind):
    ENDPDB,ENDMC="final.p"+ind+".pdb","configs/chk.last.p"+ind+".mc"
    copy_file(ENDPDB, "refined_"+wdir+"_"+ind+".pdb",localdir)
    copy_file(ENDMC, "refined_"+wdir+"_"+ind+".mc",localdir)
    return ENDPDB,ENDMC
def check_clashes(wdir,rseed):
    comm=[]
    ind=str(rseed).zfill(2)
    comm.append("cp "+BINPATH+"SPQR_ENERG .")
    comm.append("./SPQR_ENERG -i final.p"+ind+".pdb -t")
    clashoutput=[]
    for ii in comm :
        pp=subprocess.run(ii.split(),capture_output=True,text=True,cwd=wdir)
        clashoutput.append(pp.stdout)
    clflag=True
    for ll in clashoutput[1].split('\n'):
        if ll.strip()=="No clashes found.": clflag=False
    return clflag
def run_simulation(wdir,exe,rseed=1):
    comm=[]
    comm.append("./"+exe+" -i "+str(rseed))
    ind=str(rseed).zfill(2)
    for ii in comm : subprocess.run(ii.split(),stdout=subprocess.PIPE,text=True,cwd=wdir)
def rerun_simulation(wdir, exe,rseed=1):
    comm=[]
    ind=str(rseed).zfill(2)
    comm.append("cp configs/chk.last.p"+ind+".mc spqr_inits/init.mc")
    for ii in comm : subprocess.run(ii.split(),stdout=subprocess.PIPE,text=True,cwd=wdir)
    run_simulation(wdir, exe,rseed)
def create_simulation(where, wdir, binary, initfile, paramdict, rseed, add_files=[]):
    comm=[]
    comm.append("mkdir " + wdir)
    comm.append("mkdir " + wdir + "/spqr_inits")
    comm.append("cp "+BINPATH+binary+" "+wdir)
    ind=str(rseed).zfill(2)
    destfile=""
    if initfile[-4:]=='.pdb': destfile=wdir+"/spqr_inits/init.p"+ind+".pdb"
    if initfile[-3:]=='.mc': destfile=wdir+"/spqr_inits/init.p"+ind+".mc"
    if destfile=="":
        print("Invalid initfile format!")
        exit(1)
    comm.append("cp "+initfile+" "+destfile)
    for ii in comm : subprocess.run(ii.split(),stdout=subprocess.PIPE,text=True,cwd=where)
    PARAMFILE=open(wdir+"/params.spqr", "w")
    for ii in paramdict : PARAMFILE.write(ii+" "+paramdict[ii]+"\n")
    PARAMFILE.close()
    comm=[]
    if add_files!=[]:
        tline="cp " 
        for ff in add_files:
            tline+=ff+" "
        tline+=wdir
        for ii in comm : subprocess.run(ii.split(),stdout=subprocess.PIPE,text=True)
def copy_file(orig, dest, wdir):
    comm=[]
    comm.append("cp "+orig+" "+dest)
    for ii in comm : subprocess.run(ii.split(),stdout=subprocess.PIPE,text=True,cwd=wdir)    
def create_template(template, ermsddir):
    blockfile=[]
    for ii in open(template).readlines():
        if len(ii)>4:
            if ii[0:4]=='ATOM': blockfile.append(ii)
    NNT=len(blockfile)
    if NNT%NATNT!=0:
        print("Invalid template")
        exit(1)
    NNT=NNT//NATNT
    TEMPLFILE=open(ermsddir+"/ermsd_frags.lst", "w")
    wline="REMARK ERMSD PARAMS 1 100"
    TEMPLFILE.write(wline+"\n")
    wline="REMARK ERMSD GROUP 50 "
    for ii in range(NNT): wline+=" "+str(ii)
    TEMPLFILE.write(wline+"\n")
    for wline in blockfile : TEMPLFILE.write(wline)
    TEMPLFILE.close()


    
parser=argparse.ArgumentParser()
parser.add_argument("-i", "--initfile", type=str, default="init.pdb", help="PDB file")
parser.add_argument("-o", "--output", type=str, default="refined", help="Output directory name")    
parser.add_argument("-s", "--spider", type=str, default="init.csv", help="RNAspider file in csv format")
parser.add_argument("-t", "--sec_struct", type=str, default="init.ss", help="Secondary structure in Vienna format")
parser.add_argument("-r", "--seed", type=int, default=1, help="Random seed")
args=parser.parse_args()   
WDIR=args.output
RSEED=args.seed
IND=str(RSEED).zfill(2)
SPIDERFILE=args.spider
INIFILE=args.initfile
SSFILE=args.sec_struct
genparams=[["TEMPERATURE", "3"],["PDB_OUTPUT", "0"],["RG_COUPL","0  0"],["MC_PH_XYZ","2"],["MC_NT_XYZ","2"],["MC_STEPS",  " 1000"],["MC_TRAJ_STEPS",   "1000000"],["MC_CHKP_STEPS","10000000"],["RANDOM_SEED","1"],["MC_NT_ANGLE","0.5   0.5"],["ENERGS_PATH",ENERGPATH],["SA_TINI","9"],["SA_TMIN","1"],["SA_TFAC","0.75"],["SA_STEP","0"],["SA_NT","50"],["SA_PREENERG", "0"],["SA_SFAC","0.75"],["SA_RTIMES","0"]]

#initial setup
comm=[]
DIRNAME=WDIR
if os.path.exists(WDIR):
    print("Please choose a different output name as " + WDIR)
    exit(1)
comm.append("mkdir "+WDIR)
out00=run_commands(comm,".")


#energy minimization
#./SPQR_REFINE
print("Removing clashes from entangled structure")
DIR01=WDIR+"/01_energy_mini/"
#CLASH REMOVAL
comm=[]
iname="unref.pdb"
pdbstruct01=DIR01+iname
comm.append("mkdir "+DIR01)
comm.append("cp "+SPQR_PATH+"/tools/python3/pdb2spqr.py "+ DIR01 )
comm.append("cp "+INIFILE+" "+ DIR01 )
out01=run_commands(comm,".")
comm=[]
comm.append("python pdb2spqr.py " + INIFILE)
out012=run_commands(comm,DIR01)
INICONF=open(pdbstruct01, "w")
for ii in out012: INICONF.write(ii)
INICONF.close()
clashparams=dict(genparams)
clashparams['SA_NT']='10'
clashparams['MC_STEPS']='3000'
clashparams['SA_TINI']='5'
CLASHDIR="01_CLASHRM"
create_simulation(".",DIR01+CLASHDIR,"SPQR_wSA",pdbstruct01,clashparams,RSEED)
clashflag=check_clashes(DIR01+CLASHDIR,RSEED)
if clashflag:
    print("Energy minimization required for removing clashes. Running...")
    run_simulation(DIR01+CLASHDIR,"SPQR_wSA",RSEED)
    clashflag=check_clashes(DIR01+CLASHDIR,RSEED)
for clash in range(1,CLASHSTEPS):
    if not clashflag: break
    else:
        print("Step ",clash)
        rerun_simulation(DIR01+CLASHDIR,"SPQR_wSA",RSEED)
        clashflag=check_clashes(DIR01+CLASHDIR,RSEED)
if clashflag:
    print("Maximum number of clash simulations reached. Try increasing CLASHSTEPS")
    exit(1)
else:
    print("Minimization successful. Structure is free of clashes.")
    pdbfile,mcfile=rename_files(WDIR,DIR01+CLASHDIR,IND)


pdbstruct01=DIR01+CLASHDIR+"/final.p"+IND+".pdb"
struct01=DIR01+CLASHDIR+"/configs/chk.last.p"+IND+".mc"
TEMPLATEFILE=pdbstruct01


#disentanglement
print("Starting disentanglement simulation.")
DIR02=WDIR+"/02_disentanglement"
comm=[]
comm.append("mkdir "+DIR02)
comm.append("mkdir "+DIR02+"/spqr_inits")
comm.append("cp "+SPQR_PATH+"/bin/SPQR_cMC "+ DIR02)
comm.append("cp "+SPQR_PATH+"/tools/python3/spider2spqr.py "+ DIR02)
comm.append("cp "+ SPIDERFILE+ " "+ DIR02)
comm.append("cp "+ INIFILE+ " "+ DIR02)
comm.append("cp "+ SSFILE+ " "+ DIR02)
comm.append("cp "+ struct01 + " " + DIR02+"/spqr_inits/init.mc")
out021=run_commands(comm,".")
disentparams=dict(genparams)
disentparams['MC_STEPS']='6000'
disentparams['TEMPERATURE']='9'
disentparams['MC_NT_ANGLE']='0.5 0.5'
disentparams['MC_NT_XYZ']='2'

write_params(disentparams,DIR02)
#preprocessing
comm=[]
comm.append("python spider2spqr.py -i "+SPIDERFILE + " -p " + INIFILE + " -t " + SSFILE)
comm.append("./SPQR_cMC")
out022=run_commands(comm,DIR02)
struct02=DIR02+"/configs/chk.last.p00.mc"
pdbstruct02=DIR02+"/final.p00.pdb"
if os.path.isfile(struct02):
    print("Disentanglement simulation finished.")
else:
    print("Disentanglement simulation failed.")
    exit(1)

#energy minimization
print("Removing clashes from disentangled structure")
DIR03=WDIR+"/03_energy_mini/"
#CLASH REMOVAL
comm=[]
comm.append("mkdir "+DIR03)
out03=run_commands(comm,".")
clashparams=dict(genparams)
clashparams['SA_NT']='5'
CLASHDIR="01_CLASHRM"
create_simulation(".",DIR03+CLASHDIR,"SPQR_wSA",pdbstruct02,clashparams,RSEED)
clashflag=check_clashes(DIR03+CLASHDIR,RSEED)
if clashflag:
    print("Energy minimization required for removing clashes. Running...")
    run_simulation(DIR03+CLASHDIR,"SPQR_wSA",RSEED)
    clashflag=check_clashes(DIR03+CLASHDIR,RSEED)
for clash in range(1,CLASHSTEPS):
    if not clashflag: break
    else:
        print("Step ",clash)
        rerun_simulation(DIR03+CLASHDIR,"SPQR_wSA",RSEED)
        clashflag=check_clashes(DIR03+CLASHDIR,RSEED)
if clashflag:
    print("Maximum number of clash simulations reached. Try increasing CLASHSTEPS")
    exit(1)
else:
    print("Minimization successful. No clashes found")
    pdbfile,mcfile=rename_files(WDIR,DIR03+CLASHDIR,IND)

#ERMSD MINIMIZATION
print("ERMSD minimization")
ermsdparams=dict(genparams)
ermsdparams['TEMPERATURE']='9'
ermsdparams['MC_STEPS']='500'
ermsdparams['SA_TFAC']='0.5'
ermsdparams['SA_NT']='5'
ERMSDDIR="/02_ERMSD/"
INIFILEMC=DIR03+CLASHDIR+"/refined_"+WDIR+"_"+IND+".mc"
create_simulation(".",DIR03+ERMSDDIR,"SPQR_mSA",INIFILEMC,ermsdparams,RSEED)
create_template(TEMPLATEFILE,DIR03+ERMSDDIR)
run_simulation(DIR03+ERMSDDIR, "SPQR_mSA",RSEED)
pdbfile,mcfile=rename_files(WDIR,DIR03+ERMSDDIR,IND)
print("Procedure finished.")
print("Backmapping structure...")
comm=[]
comm.append("cp "+SPQR_PATH+"/tools/python3/SPQR_BBACKMAP.py "+ WDIR )
outbm1=run_commands(comm,".")
comm=[]
BMAPPEDNAME="AA_"+WDIR+".pdb"
print( DIR03+ERMSDDIR+pdbfile)
comm.append("python SPQR_BBACKMAP.py -i "+ "03_energy_mini"+ERMSDDIR+pdbfile + " -o "+BMAPPEDNAME)
outbm1=run_commands(comm,WDIR)
print("All-atom structure saved in ", WDIR+"/"+BMAPPEDNAME)
