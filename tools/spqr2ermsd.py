from math import sqrt
import argparse
import sys
NATNT=5
ALLFILE=[]

KERMSD=10
KSSERMSD=50
RERMSD=15

parser=argparse.ArgumentParser()
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-o","--output", help="Output file",type=str,default="bbm_output.pdb")
parser.add_argument("-t","--sstruct", help="Secondary structure",type=str,default="")
args=parser.parse_args()
ssfile=open(args.sstruct,"r")
ssfile.readline()
SEQ=ssfile.readline().strip()
rawSSTR=list(ssfile.readline().strip())

for nt in xrange(0,len(rawSSTR)):
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
        for nt in xrange(0,len(rawSSTR)):
            if(rawSSTR[nt]==br):
                pnt=nt
            if(rawSSTR[nt]==SSDICT[br] and pnt>=0):
                bpairs.append([pnt,nt])
                rawSSTR[nt]="."
                rawSSTR[pnt]="."
                restart=True
                break



INPUTFILE=args.input
if INPUTFILE=="":
    print "ERROR : input file needed."
    parser.print_help()
    exit(1)
ERMSDFILE=args.output
for line in open(INPUTFILE):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)


NATS=len(ALLFILE)
NNT=NATS/NATNT

#order is base, pucker, glyc - as usual: puck 3,2=0,1; glyc A,H,S=0,1,2

COUNTAT=0
OUTFILE=open(ERMSDFILE,"w")
orig_stdout=sys.stdout
sys.stdout=OUTFILE
print "REMARK ERMSD PARAMS 1 " + str(RERMSD)
stp="REMARK ERMSD SSPAIRS " + str(KSSERMSD)
for pair in bpairs:
    stp=stp+" " +str(pair[0])+" "+str(pair[1])
print stp
stp="REMARK ERMSD GROUP " + str(KERMSD)
for nt in xrange(0,len(ALLFILE)/NATNT):
    stp=stp+" "+str(nt)
print stp
for at in xrange(0,len(ALLFILE)):
    print ALLFILE[at].strip()
sys.stdout=orig_stdout
OUTFILE.close()
#print "ENDMDL"
