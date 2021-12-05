from math import sqrt
import argparse
import sys
NATNT=5
ALLFILE=[]
KERMSD=50
RERMSD=100
parser=argparse.ArgumentParser()
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-f","--fasta", help="Secondary structure file",type=str,default="input.fa")
parser.add_argument("-o","--output", help="Output name",type=str,default="ermsd_frags_ss.lst")
args=parser.parse_args()
OUTNAME=args.output
FASTA=args.fasta
SSDICT = {
    "(":")"
    #,
    #"[":"]",
    #"{":"}",
    #"<":">"
}
INPUTFILE=args.input
if INPUTFILE=="":
    print("ERROR : input file needed.")
    parser.print_help()
    exit(1)


for line in open(INPUTFILE):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)

if(FASTA!=""):
    seq=[]
    strandseq=[]
    fullss=[]
    rawseq=[]
    fastafile=open(FASTA,"r")
    fastablock=fastafile.readlines()
    fastaseq=fastablock[1].strip()
    fastasstruct=fastablock[2].strip()
    if(fastasstruct!=""):
        ssflag=True
    fastafile.close()
    for i in range(0,len(fastaseq)):
        if(fastaseq[i]!="&"):
            strandseq.append(fastaseq[i])
            rawseq.append(fastaseq[i])
            if(ssflag):
                fullss.append(fastasstruct[i])
        if(fastaseq[i]=="&"):
            seq.append(strandseq)
            strandseq=[]
    seq.append(strandseq)

bpairs=[]
for br in SSDICT:
    pnt=-1
    restart=True
    while(restart):
        restart=False
        for nt in range(0,len(fullss)):
            if(fullss[nt]==br):
                pnt=nt
            if(fullss[nt]==SSDICT[br] and pnt>=0):
                bpairs.append([pnt,nt])
                fullss[nt]="."
                fullss[pnt]="."
                restart=True
                break
opairs=bpairs.sort()
SSSTACKS=[]
CSTACK=[bpairs[0]]
for bp in range(1,len(bpairs)):
    if(bpairs[bp][0]==bpairs[bp-1][0]+1 and bpairs[bp][1]==bpairs[bp-1][1]-1):
        CSTACK.append(bpairs[bp])
    else:
        if(len(CSTACK)>0):
            SSSTACKS.append(CSTACK)
        CSTACK=[bpairs[bp]]
if(len(CSTACK)>0):
    SSSTACKS.append(CSTACK)

STACKSEQS=[]
stra=[]
for ST in range(0,len(SSSTACKS)):
    stra=[]
    for nt in range(0,len(SSSTACKS[ST])):
        stra.append(rawseq[SSSTACKS[ST][nt][0]])
    STACKSEQS.append(stra)


frags=open(OUTNAME,"w")
frags.write( "REMARK ERMSD PARAMS " + str(len(STACKSEQS)) + " "+str(RERMSD)+"\n")
for ST in range(0,len(SSSTACKS)):
    frags.write( "REMARK ERMSD GROUP "+str(KERMSD)+" ")
    for ist in SSSTACKS[ST]:
        frags.write(str(ist[0])+" ")
    for ist in reversed(SSSTACKS[ST]):
        frags.write(str(ist[1])+" ")
    frags.write("\n")

for ST in range(0,len(SSSTACKS)):
    STLEN=len(SSSTACKS[ST])
    #print SSSTACKS[ST]
    for PA in range(0,STLEN):
        for AT in range(0,NATNT):
            #print PA*NATNT+AT
            #print SSSTACKS[ST][PA][0]*NATNT+AT
            frags.write(ALLFILE[SSSTACKS[ST][PA][0]*NATNT+AT])
    for PA in range(0,STLEN):
        for AT in range(0,NATNT):
            #print (STLEN-PA-1)*NATNT+AT
            #print SSSTACKS[ST][(STLEN-PA-1)][1]*NATNT+AT
            frags.write(ALLFILE[SSSTACKS[ST][(STLEN-PA-1)][1]*NATNT+AT])
frags.close()


