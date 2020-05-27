from math import sqrt,pi
from scipy import integrate
import numpy as np
import argparse
import sys

def diloops(l1, l2):
    c1,c2=np.array([0,0,0]),np.array([0,0,0])
    for a1 in l1:
        c1=c1+a1
    c1=c1/len(l1)
    for a2 in l2:
        c2=c2+a2
    c2=c2/len(l2)
    max1,max2=0,0
    for a1 in l1:
        max1=max(max1,np.linalg.norm(a1-c1))
    for a2 in l2:
        max2=max(max2,np.linalg.norm(a2-c2))
    dcl1l2=np.linalg.norm(c1-c2)
    if(dcl1l2>max1+max2):
        ret=False
    else:
        ret=True
    return ret
    
def difloops(l1,l2):
    fl=True
    for a1 in l1:
        for a2 in l2:
            if (a1[0]==a2[0] and a1[1]==a2[1] and a1[2]==a2[2]):
                fl=False
                break
    return fl

parser=argparse.ArgumentParser()
parser.add_argument("-t","--sstruct", help="Secondary structure file. First line must be sequence. Second, secondary structure in Vienna format.",type=str,default="")
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-o","--output", help="Output file",type=str,default="linked_loops_output.lst")
args=parser.parse_args()
SSFILE=args.sstruct

#fullss=list(args.sstruct)
INPUTFILE=args.input
LINKFILE=args.output
if INPUTFILE=="" or SSFILE=="":
    print "ERROR : input or secondary structure files missing."
    parser.print_help()
    exit(1)

ALLSSFILE=[]
for line in open(SSFILE):
    if line[0]!=">":
        ALLSSFILE.append(line.strip())
fullss=list(ALLSSFILE[1])

SSDICT = {
    "(":")",
    "[":"]",
    "{":"}",
    "<":">"
}
KL=50
##FIND HAIRPINS##
HAIRPINS=[]
loops=[]
tloop=[]
op=0
for nt in xrange(0,len(fullss)):
    if fullss[nt] == ".":
        tloop.append(nt)
        op=1
    if(op==1 and fullss[nt]!="."):
        op=0
        loops.append(tloop)
        tloop=[]
if op==1:
    loops.append(tloop)
for lo in loops:
    fl=0
    linit,lend=lo[0],lo[len(lo)-1]
    if(linit > 0 and lend<len(fullss)-1):
        for br in SSDICT:
            if(fullss[linit-1]==br and fullss[lend+1]==SSDICT[br]):
                fl=1
    if(fl==1):
        HAIRPINS.append(lo)

##SORT BASE PAIRS##
bpairs=[]
for br in SSDICT:
    pnt=-1
    restart=True
    while(restart):
        restart=False
        for nt in xrange(0,len(fullss)):
            if(fullss[nt]==br):
                pnt=nt
            if(fullss[nt]==SSDICT[br] and pnt>=0):
                bpairs.append([pnt,nt])
                fullss[nt]="."
                fullss[pnt]="."
                restart=True
                break

##FIND STACKS##
opairs=bpairs.sort()
preSTEMS,STEMS=[],[]
CSTACK=[bpairs[0]]
for bp in xrange(1,len(bpairs)):
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

##GET COORDINATES##
N_PARTS_PER_NT=5
ALLFILE=[]
for line in open(INPUTFILE):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)
NATS=len(ALLFILE)
NNT=NATS/N_PARTS_PER_NT
COORDS=[]
for i in xrange(0,NNT):
    cline=ALLFILE[i*N_PARTS_PER_NT]
    sline=ALLFILE[i*N_PARTS_PER_NT+3]
    pline=ALLFILE[i*N_PARTS_PER_NT+4]
    basetyp=cline[17:20].strip()
    resind=int(cline[22:26].strip())
    chain=cline[21].strip()
    BASE=[float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())]
    SUGR=[float(sline[30:38].strip()), float(sline[38:46].strip()), float(sline[46:54].strip())]
    PHOS=[float(pline[30:38].strip()), float(pline[38:46].strip()), float(pline[46:54].strip())]
    COORDS.append([BASE,SUGR,PHOS])

##COORDINATES IN LOOPS##
HPCOORDS=[]
STCOORDS=[]
#rint HAIRPINS
#rint STEMS
for hp in HAIRPINS:
    hpc=[]
    hpc.append(np.array(COORDS[hp[0]-1][0]))
    hpc.append(np.array(COORDS[hp[0]-1][1]))
    for nt in xrange(0,len(hp)):
        hpc.append(np.array(COORDS[hp[nt]][2]))
        hpc.append(np.array(COORDS[hp[nt]][1]))
    hpc.append(np.array(COORDS[hp[len(hp)-1]+1][2]))
    hpc.append(np.array(COORDS[hp[len(hp)-1]+1][1]))
    hpc.append(np.array(COORDS[hp[len(hp)-1]+1][0]))
    hpc.append(np.array(COORDS[hp[0]-1][0]))
    HPCOORDS.append(hpc)

for st in STEMS:
    stc=[]
    stc.append(np.array(COORDS[st[0]][0]))
    stc.append(np.array(COORDS[st[0]][1]))
    for nt in xrange(1,len(st)/2):
        stc.append(np.array(COORDS[st[nt]][2]))
        stc.append(np.array(COORDS[st[nt]][1]))
    stc.append(np.array(COORDS[st[len(st)/2-1]][0]))
    stc.append(np.array(COORDS[st[len(st)/2]][0]))
    stc.append(np.array(COORDS[st[len(st)/2]][1]))
    for nt in xrange(len(st)/2+1,len(st)):
        stc.append(np.array(COORDS[st[nt]][2]))
        stc.append(np.array(COORDS[st[nt]][1]))
    stc.append(np.array(COORDS[st[len(st)-1]][0]))
    stc.append(np.array(COORDS[st[0]][0]))
    STCOORDS.append(stc)

LOOPCOORDS=[]
for hp in xrange(0,len(HPCOORDS)):
    LOOPCOORDS.append([HPCOORDS[hp],["hairpin ",hp]])
for st in xrange(0,len(STCOORDS)):
    LOOPCOORDS.append([STCOORDS[st],["stem ", st]])
SEGF=lambda t1,t2 : np.dot(((p1- q1)* t1+ q1- ((p2- q2)*t2+ q2)),(np.cross((p1- q1), (p2- q2))))/np.linalg.norm((p1- q1)*t1+ q1-((p2- q2)*t2+ q2))**3

ALLINDEXES=[]
for iloop1 in xrange(0,len(LOOPCOORDS)):
    for iloop2 in xrange(iloop1+1,len(LOOPCOORDS)):
        GAUSSINT=0
        dist,sameflag=diloops(LOOPCOORDS[iloop1][0],LOOPCOORDS[iloop2][0]),difloops(LOOPCOORDS[iloop1][0],LOOPCOORDS[iloop2][0])
        if(dist and sameflag):
            for pa1 in xrange(0,len(LOOPCOORDS[iloop1][0])-1):
                for pa2 in xrange(0,len(LOOPCOORDS[iloop2][0])-1):
                    p1,q1,p2,q2=LOOPCOORDS[iloop1][0][pa1+1],LOOPCOORDS[iloop1][0][pa1],LOOPCOORDS[iloop2][0][pa2+1],LOOPCOORDS[iloop2][0][pa2]
                    #print p1, q1, p2, q2
                    GAUSSINT=GAUSSINT+integrate.dblquad(SEGF,0,1,lambda x: 0, lambda y: 1)[0]
            GAUSSINT=GAUSSINT/(4.0*pi)
            if(GAUSSINT>0.5):
                type1,looi1=LOOPCOORDS[iloop1][1][0],LOOPCOORDS[iloop1][1][1]
                type2,looi2=LOOPCOORDS[iloop2][1][0],LOOPCOORDS[iloop2][1][1]
                print "LINK DETECTED: loops "+str(iloop1)+" ( "+type1+")"+" and "+str( iloop2)+" ( "+type2+")"+". Gauss integral = "+str(GAUSSINT)
                indexes=[]
                if(type1=="hairpin "):
                    indexes.append([HAIRPINS[looi1],"hp"])
                elif(type1=="stem "):
                    indexes.append([STEMS[looi1],"st"])
                if(type2=="hairpin "):
                    indexes.append([HAIRPINS[looi2],"hp"])
                elif(type2=="stem "):
                    indexes.append([STEMS[looi2],"st"])
                ALLINDEXES.append(indexes)
NLINKS=len(ALLINDEXES)
if NLINKS>0:
#ELEMS="ABCDEFGHIJKL"
#for hp in HPCOORDS:
#    for icoo in xrange(0,len(hp)):
#        print ELEMS[cnt],hp[icoo][0],hp[icoo][1],hp[icoo][2]
#    cnt=cnt+1
#
#for st in STCOORDS:
#    for icoo in xrange(0,len(st)):
#        print ELEMS[cnt],st[icoo][0],st[icoo][1],st[icoo][2]
#    cnt=cnt+1
    OUTFILE=open(LINKFILE,"w")
    orig_stdout=sys.stdout
    sys.stdout=OUTFILE
    print "REMARK LNKDLPS "+str(KL)+" "+str(NLINKS)
    for li in ALLINDEXES:
        typs=li[0][1]+li[1][1]
        loop1,loop2=[],[]
        if li[0][1]=='hp':
            loop1.append(li[0][0][0]-1)
            loop1.append(li[0][0][len(li[0][0])-1]+1)
        if li[0][1]=='st':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])/2-1])
            loop1.append(li[0][0][len(li[0][0])/2])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[1][1]=='hp':
            loop2.append(li[1][0][0]-1)
            loop2.append(li[1][0][len(li[1][0])-1]+1)
        if li[1][1]=='st':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])/2-1])
            loop2.append(li[1][0][len(li[1][0])/2])
            loop2.append(li[1][0][len(li[1][0])-1])
        
        print typs+ " "+' '.join(str(x) for x in loop1) + " "+' '.join(str(x) for x in loop2)
        
    OUTFILE.close()
    sys.stdout=orig_stdout
