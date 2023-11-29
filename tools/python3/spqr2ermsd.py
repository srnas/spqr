from math import sqrt
import argparse
import sys
NATNT=5
ALLFILE=[]

KERMSD=50
KSSERMSD=50
RERMSD=100

parser=argparse.ArgumentParser()
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-o","--output", help="Output file",type=str,default="out_ermsd_frags.lst")
parser.add_argument("-t","--sstruct", help="Secondary structure",type=str,default="")
parser.add_argument("-r","--onlystems", help="Only stems",action="store_true")
parser.add_argument("-e","--ensurestems", help="Ensure integrity of stems",action="store_true")
parser.add_argument("-k1","--K_ERMSD", help="K ERMSD",type=float,default=KERMSD)
parser.add_argument("-k2","--K_SS_ERMSD", help="K ERMSD secondary structure",type=float,default=KSSERMSD)
parser.add_argument("-kr","--R_ERMSD", help="R ERMSD ",type=float,default=RERMSD)



args=parser.parse_args()
INPUTFILE=args.input
if INPUTFILE=="":
    print("ERROR : input file needed.")
    parser.print_help()
    exit(1)

WITHSS=False
if args.sstruct != "" :
    ssfile=open(args.sstruct,"r")
    WITHSS=True
    ssfile.readline()
    SEQ=ssfile.readline().strip()
    rawSSTR=list(ssfile.readline().strip())
ONLYSS=args.onlystems
ENSUSS=args.ensurestems
KERMSD=args.K_ERMSD
KSSERMSD=args.K_SS_ERMSD
RERMSD=args.R_ERMSD


if WITHSS:
    for nt in range(0,len(rawSSTR)):
        if rawSSTR[nt]!="." and rawSSTR[nt]!=")" and rawSSTR[nt]!="(":
            rawSSTR[nt]="."
    #SSTR="".join(rawSSTR)

    SSDICT = {
        "(":")"
    }
    bpairs=[]
    for br in SSDICT:
        pnt=-1
        restart=True
        while(restart):
            restart=False
            for nt in range(0,len(rawSSTR)):
                if(rawSSTR[nt]==br):
                    pnt=nt
                if(rawSSTR[nt]==SSDICT[br] and pnt>=0):
                    bpairs.append([pnt,nt])
                    rawSSTR[nt]="."
                    rawSSTR[pnt]="."
                    restart=True
                    break
    ssnts=[]
    bpairs=sorted(bpairs)
    for nt in bpairs:
        ssnts.append(nt[0])
        ssnts.append(nt[1])
    
    #stems
    preSTEMS,STEMS=[],[]
    CSTACK=[bpairs[0]]
    for bp in range(1,len(bpairs)):
        if(bpairs[bp][0]==bpairs[bp-1][0]+1 and bpairs[bp][1]==bpairs[bp-1][1]-1):
            CSTACK.append(bpairs[bp])
        else:
            if(len(CSTACK)>0):
                preSTEMS.append(CSTACK)
            CSTACK=[bpairs[bp]]
    if(len(CSTACK)>0):
        preSTEMS.append(CSTACK)
    for st in preSTEMS:
        unf=[]
        for pa in st:
            unf.append(pa[0])
        for pa in reversed(st):
            unf.append(pa[1])
        STEMS.append(unf)
##########

ERMSDFILE=args.output
for line in open(INPUTFILE):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)


NATS=len(ALLFILE)
NNT=NATS//NATNT

#order is base, pucker, glyc - as usual: puck 3,2=0,1; glyc A,H,S=0,1,2
NTLIST=list(range(0,len(ALLFILE)//NATNT))
if WITHSS:
    if(ONLYSS):
        NTLIST=sorted(ssnts)

COUNTAT=0
OUTFILE=open(ERMSDFILE,"w")
orig_stdout=sys.stdout
sys.stdout=OUTFILE
NGROUPS=1
if WITHSS and ONLYSS:
    NGROUPS=len(STEMS) 
print("REMARK ERMSD PARAMS "+ str(NGROUPS)+" " + str(RERMSD))

if WITHSS and ENSUSS:
    stp="REMARK ERMSD SSPAIRS " + str(KSSERMSD)
    for pair in bpairs:
        stp=stp+" " +str(pair[0])+" "+str(pair[1])
    print(stp)

if WITHSS and ONLYSS:
    for st in STEMS:
        stp="REMARK ERMSD GROUP " + str(KERMSD)
        for nt in st:
            stp=stp+" "+str(nt)
        print(stp)
            
else:
    stp="REMARK ERMSD GROUP " + str(KERMSD)
    for nt in NTLIST:
        stp=stp+" "+str(nt)
    print(stp)

for nt in NTLIST:
    for at in range(0,NATNT):
        print(ALLFILE[nt*NATNT+at].strip())

sys.stdout=orig_stdout
OUTFILE.close()

