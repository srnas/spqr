from math import sqrt,pi
from scipy import integrate
import numpy as np
import argparse
import sys

parser=argparse.ArgumentParser()
parser.add_argument("-t","--sstruct", help="Secondary structure file. First line must be sequence. Second, secondary structure in Vienna format.",type=str,default="")
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-o","--output", help="Output file",type=str,default="linked_loops_output.lst")
parser.add_argument("-p","--calc_type_pierce", help="Perform piercing calculation",action="store_true")
parser.add_argument("-v","--verbose", help="Prints loop indexes",action="store_true")
args=parser.parse_args()
SSFILE=args.sstruct

calctype="g"
if args.calc_type_pierce:
    calctype="p"
verboseflag=False
if args.verbose:
    verboseflag=True
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
NNT=len(fullss)
SSDICT = {
    "(":")"
}
#   "[":"]",
#   "{":"}",
#   "<":">"
#

KL=50

def make_triangles(vertlist):
    #remember that in the current nomenclature, the last vertex is repeated
    triangles=[]
    NL=len(vertlist)
    CM=(np.add.reduce(vertlist)-vertlist[NL-1])/float(NL-1)
    for tr in xrange(0,NL-2):
        triangles.append([vertlist[tr],vertlist[tr+1],CM])
    triangles.append([vertlist[NL-2],vertlist[0],CM])
    return triangles

def check_pierce(tri, vlist):
    npie=0
    V1,V2=tri[1]-tri[0],tri[2]-tri[0]
    NORM=np.cross(V1,V2)
    NORM=NORM/sqrt(np.dot(NORM,NORM))
    for segind in xrange(0,len(vlist)-1):
        seg=vlist[segind+1]-vlist[segind]
        tinter=np.dot(NORM,(tri[0]-vlist[segind]))/np.dot(NORM,seg)
        if(tinter>0 and tinter<1):
            ivec=tinter*seg+vlist[segind]
            check=np.dot(ivec-tri[0],NORM)
            if(check>0.001):
                print "no!", check
                exit(1)
            if(np.dot(np.cross(ivec-tri[0],tri[1]-tri[0]),np.cross(tri[2]-tri[0],tri[1]-tri[0]))>0 and
               np.dot(np.cross(ivec-tri[1],tri[2]-tri[1]),np.cross(tri[0]-tri[1],tri[2]-tri[1]))>0 and
               np.dot(np.cross(ivec-tri[2],tri[0]-tri[2]),np.cross(tri[1]-tri[2],tri[0]-tri[2]))>0):
                npie=npie+1
    return npie

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

def inddifloops(l1,l2):
    fl=True
    for a1 in l1:
        for a2 in l2:
            if (a1==a2):
                fl=False
                break
    return fl


#FIND UNSTRUCTURED LOOPS
UNSTLS=[]
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

#f calctype=="g":
for lo in loops:
    fl=1
    linit,lend=lo[0],lo[len(lo)-1]
    if(linit > 0 and lend<len(fullss)-1):
        fl=1
        for br in SSDICT:
            #if(fullss[linit-1]!=br or fullss[lend+1]!=SSDICT[br]):
            #    fl=1
            if(fullss[linit-1]==br and fullss[lend+1]==SSDICT[br]):
                fl=0
    if(fl==1 and len(lo)>2):
        if(linit>1):
            linit=linit-1
        if(lend<NNT-1):
            lend=lend+1
        UNSTLS.append(range(linit,lend+1))

#f calctype=="p":
#   for lo in loops:
#       fl=0
#       linit,lend=lo[0],lo[len(lo)-1]
#       if(linit>1):
#           linit=linit-1
#       if(lend<NNT-1):
#           lend=lend+1
#       
#       if len(lo)>2:
#           UNSTLS.append(range(linit,lend+1))


##FIND HAIRPINS##
HAIRPINS=[]
loops=[]
tloop=[]
op=0
#f calctype=="g":
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
        if(linit>1):
            linit=linit-1
        if(lend<NNT-1):
            lend=lend+1
        HAIRPINS.append(range(linit,lend+1))

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

if verboseflag:
    print "UNSTLS" , UNSTLS
    print "HAIRPINS", HAIRPINS
    print "STEMS" , STEMS

##GET COORDINATES##
N_PARTS_PER_NT=5
ALLFILE=[]
for line in open(INPUTFILE):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)
NATS=len(ALLFILE)
#NT=NATS/N_PARTS_PER_NT
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
ULCOORDS=[]

for ul in UNSTLS:
    ulc=[]
    ulc.append(np.array(COORDS[ul[0]][0]))
    ulc.append(np.array(COORDS[ul[0]][1]))
    for nt in xrange(1,len(ul)):
        ulc.append(np.array(COORDS[ul[nt]][2]))
        ulc.append(np.array(COORDS[ul[nt]][1]))
    ulc.append(np.array(COORDS[ul[len(ul)-1]][0]))
    ulc.append(np.array(COORDS[ul[0]][0]))
    ULCOORDS.append(ulc)


for hp in HAIRPINS:
    hpc=[]
    hpc.append(np.array(COORDS[hp[0]][0]))
    hpc.append(np.array(COORDS[hp[0]][1]))
    for nt in xrange(1,len(hp)):
        hpc.append(np.array(COORDS[hp[nt]][2]))
        hpc.append(np.array(COORDS[hp[nt]][1]))
    hpc.append(np.array(COORDS[hp[len(hp)-1]+1][0]))
    hpc.append(np.array(COORDS[hp[0]][0]))
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
#print STCOORDS
ALLINDEXES=[]
if calctype=="g":
    LOOPCOORDS=[]
    for hp in xrange(0,len(HPCOORDS)):
        LOOPCOORDS.append([HPCOORDS[hp],["hairpin ",hp]])
    for st in xrange(0,len(STCOORDS)):
        LOOPCOORDS.append([STCOORDS[st],["stem ", st]])
    for ul in xrange(0,len(ULCOORDS)):
        LOOPCOORDS.append([ULCOORDS[ul],["unst ", ul]])
    SEGF=lambda t1,t2 : np.dot(((p1- q1)* t1+ q1- ((p2- q2)*t2+ q2)),(np.cross((p1- q1), (p2- q2))))/np.linalg.norm((p1- q1)*t1+ q1-((p2- q2)*t2+ q2))**3

    for iloop1 in xrange(0,len(LOOPCOORDS)):
        for iloop2 in xrange(iloop1+1,len(LOOPCOORDS)):
    #for iloop1 in xrange(0,0):
    #    for iloop2 in xrange(iloop1+1,0):
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
                    print "LINK DETECTED: loops "+str(iloop1)+" ( "+type1+" )"+" and "+str( iloop2)+" ( "+type2+")"+". Gauss integral = "+str(GAUSSINT)
                    indexes=[]
                    if(type1=="hairpin "):
                        indexes.append([HAIRPINS[looi1],"hp"])
                    elif(type1=="stem "):
                        indexes.append([STEMS[looi1],"st"])
                    elif(type1=="unstl "):
                        indexes.append([UNSTLS[looi1],"ul"])
                    if(type2=="hairpin "):
                        indexes.append([HAIRPINS[looi2],"hp"])
                    elif(type2=="stem "):
                        indexes.append([STEMS[looi2],"st"])
                    elif(type2=="unstl "):
                        indexes.append([UNSTLS[looi2],"ul"])
                    ALLINDEXES.append(indexes)
    print "Link analysis finished"

#print "Calculating pierces"
#check piercings
PIERCEINDEXES=[]

if calctype=="p":
    LOOPCOORDS=[]
    LOOPINDEXS=[]
    for hp in xrange(0,len(HPCOORDS)):
        LOOPCOORDS.append([HPCOORDS[hp],["hairpin ",hp]])
        LOOPINDEXS.append(HAIRPINS[hp])
    for st in xrange(0,len(STCOORDS)):
        LOOPCOORDS.append([STCOORDS[st],["stem ", st]])
        LOOPINDEXS.append(STEMS[st])
    for ul in xrange(0,len(ULCOORDS)):
        LOOPCOORDS.append([ULCOORDS[ul],["unst ", ul]])
        LOOPINDEXS.append(UNSTLS[ul])
#    print LOOPINDEXS
    for iloop1 in xrange(0,len(LOOPCOORDS)):
        trianglist=make_triangles(LOOPCOORDS[iloop1][0])
        for iloop2 in xrange(0,len(LOOPCOORDS)):
            #sameflag=inddifloops(LOOPCOORDS[iloop1][0],LOOPCOORDS[iloop2][0])
            sameflag=inddifloops(LOOPINDEXS[iloop1],LOOPINDEXS[iloop2])
            #print sameflag, iloop1,iloop2
            #if(sameflag):
             #   print sameflag,iloop1,iloop2
            if(sameflag):
 #               print sameflag
                totpierce=0
                for triang in trianglist:
                    npierce=check_pierce(triang,LOOPCOORDS[iloop2][0])
                    totpierce=totpierce+npierce
                if totpierce>0:
                    type1,looi1=LOOPCOORDS[iloop1][1][0],LOOPCOORDS[iloop1][1][1]
                    type2,looi2=LOOPCOORDS[iloop2][1][0],LOOPCOORDS[iloop2][1][1]
                    print "LINK DETECTED: loops "+str(iloop1)+" ( "+type1+" )"+" and "+str( iloop2)+" ( "+type2+")"+" pierced " + str(totpierce) +" times."
                    indexes=[]
                    if(type1=="hairpin "):
                        indexes.append([HAIRPINS[looi1],"hp"])
                    elif(type1=="stem "):
                        indexes.append([STEMS[looi1],"st"])
                    elif(type1=="unst "):
                        indexes.append([UNSTLS[looi1],"ul"])
                    if(type2=="hairpin "):
                        indexes.append([HAIRPINS[looi2],"hp"])
                    elif(type2=="stem "):
                        indexes.append([STEMS[looi2],"st"])
                    elif(type2=="unst "):
                        indexes.append([UNSTLS[looi2],"ul"])
                    if(iloop1<iloop2):
                        PIERCEINDEXES.append(indexes)

#    for hp in xrange(0,len(ULCOORDS)):
#        trianglist=make_triangles(ULCOORDS[hp])
#        for ul in xrange(hp+1,len(ULCOORDS)):
#            totpierce=0
#            for triang in trianglist:
#                npierce=check_pierce(triang,ULCOORDS[ul])
#                totpierce=totpierce+npierce
#            if totpierce>0:
#                print "Unstructured loop ",ul, "pierces unstructured loop (",ul,")", totpierce, "times."
#                indexes=[]
#                type1,type2="ul","ul"
#                indexes.append([UNSTLS[ul],type1])
#                indexes.append([UNSTLS[hp],type2])
#                PIERCEINDEXES.append(indexes)
#for st in xrange(0,len(STCOORDS)):
#    trianglist=make_triangles(STCOORDS[st])
#    for ul in xrange(0,len(ULCOORDS)):
#        totpierce=0
#        for triang in trianglist:
#            npierce=check_pierce(triang,ULCOORDS[ul])
#            totpierce=totpierce+npierce
#        if totpierce>0:
#            print "Unstructured loop ",ul, "pierces stack (",st,")", totpierce, "times."
#            indexes=[]
#            type1,type2="ul","st"
#            indexes.append([UNSTLS[ul],type1])
#            indexes.append([STEMS[st],type2])
#            PIERCEINDEXES.append(indexes)

NLINKS,NPIERCES=len(ALLINDEXES),len(PIERCEINDEXES)

if NLINKS>0 or NPIERCES>0:
    OUTFILE=open(LINKFILE,"w")
    orig_stdout=sys.stdout
    sys.stdout=OUTFILE
    print "REMARK LNKDLPS "+str(KL)+" "+str(NLINKS+NPIERCES)
    for li in ALLINDEXES:
        typs=li[0][1]+li[1][1]
        loop1,loop2=[],[]
        if li[0][1]=='hp':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[0][1]=='ul':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[0][1]=='st':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])/2-1])
            loop1.append(li[0][0][len(li[0][0])/2])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[1][1]=='hp':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])-1])
        if li[1][1]=='ul':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])-1])
        if li[1][1]=='st':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])/2-1])
            loop2.append(li[1][0][len(li[1][0])/2])
            loop2.append(li[1][0][len(li[1][0])-1])
        print typs+ " "+' '.join(str(x) for x in loop1) + " "+' '.join(str(x) for x in loop2)
    for li in PIERCEINDEXES:
        typs=li[0][1]+li[1][1]
        loop1,loop2=[],[]
        if li[0][1]=='ul':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[0][1]=='st':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])/2-1])
            loop1.append(li[0][0][len(li[0][0])/2])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[0][1]=='hp':
            loop1.append(li[0][0][0])
            loop1.append(li[0][0][len(li[0][0])-1])
        if li[1][1]=='ul':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])-1])
        if li[1][1]=='st':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])/2-1])
            loop2.append(li[1][0][len(li[1][0])/2])
            loop2.append(li[1][0][len(li[1][0])-1])
        if li[1][1]=='hp':
            loop2.append(li[1][0][0])
            loop2.append(li[1][0][len(li[1][0])-1])
        print typs+ " "+' '.join(str(x) for x in loop1) + " "+' '.join(str(x) for x in loop2)

    OUTFILE.close()
    sys.stdout=orig_stdout
