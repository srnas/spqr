from math import sqrt,pi,fabs
from scipy import integrate
import numpy as np
import argparse
import sys

parser=argparse.ArgumentParser()
parser.add_argument("-t","--sstruct", help="Secondary structure file. First line must be sequence. Second, secondary structure in Vienna format.",type=str,default="")
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-o","--output", help="Output file",type=str,default="output")
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
OUTPUT=args.output
LINKFILE="linked_loops_"+OUTPUT+".lst"
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

KL=500

def make_triangles(vertlist):
    #remember that in the current nomenclature, the last vertex is repeated
    triangles=[]
    NL=len(vertlist)
    CM=(np.add.reduce(vertlist)-vertlist[NL-1])/float(NL-1)
    for tr in xrange(0,NL-1):
        triangles.append([vertlist[tr],vertlist[tr+1],CM])
        #xtriangles.append([vertlist[NL-2],vertlist[0],CM])
    return triangles

def check_pierce(tri, vlist):
    npie=0
    V1,V2=tri[1]-tri[0],tri[2]-tri[0]
    NORM=np.cross(V1,V2)
    NORM=NORM/sqrt(np.dot(NORM,NORM))
    for segind in xrange(0,len(vlist)-1):
        intersc=False
        for coo in tri:
            intersc=(intersc or np.array_equal(coo,vlist[segind+1]) or np.array_equal(coo,vlist[segind]))
        if(intersc==False):
            #if(np.dot(tri[0]-vlist[segind+1],tri[0]-vlist[segind+1])==0 or np.dot(tri[1]-vlist[segind+1],tri[1]-vlist[segind+1])==0 or np.dot(tri[2]-vlist[segind+1],tri[2]-vlist[segind+1])==0 or np.dot(tri[0]-vlist[segind],tri[0]-vlist[segind])==0 or np.dot(tri[1]-vlist[segind],tri[1]-vlist[segind])==0 or np.dot(tri[2]-vlist[segind],tri[2]-vlist[segind])==0):
             #   print "AIEEEEEEEE!"
             #   exit(1)
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
                #else:
            #print "someone is repeated:" , tri, vlist[segind+1], vlist[segind]
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



##FIND HAIRPINS##
HAIRPINS=[]
loops=[]
tloop=[]
op=0
#f calctype=="g":
for nt in xrange(0,NNT):
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
    if(linit > 0 and lend<NNT-1):
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
        for nt in xrange(0,NNT):
            if(fullss[nt]==br):
                pnt=nt
            if(fullss[nt]==SSDICT[br] and pnt>=0):
                bpairs.append([pnt,nt])
                fullss[nt]="."
                fullss[pnt]="."
                restart=True
                break

PAIRLIST=[]
for nt in xrange(0,NNT):
    cpair=-1
    for pair in bpairs:
        if(nt==pair[0]):
            cpair=pair[1]
        if(nt==pair[1]):
            cpair=pair[0]
    PAIRLIST.append(cpair)

##FIND STEMS##
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

#FIND UNSTRUCTURED LOOPS
UNSTLS=[]
loops=[]
tloop=[]
op=0
prevst=0

for nt in xrange(0,NNT-1):
    ist,ihp,iul,nist=-1,-1,-1,-1
    for st in xrange(0,len(STEMS)):
        if(nt in STEMS[st]):
            ist=st
        if(nt+1 in STEMS[st]):
            nist=st
    for hp in xrange(0,len(HAIRPINS)):
        if(nt in HAIRPINS[hp]):
            ihp=hp
            
    if(((ist!=-1 and nist==-1 and ihp==-1) or (ist!=-1 and nist!=-1 and nist!=ist) or (ist==-1 and nt==0)) ):
        #start internal loop
        lstart,lnt=nt,nt
        tloop.append(nt)
        lendfl=False
        lnt=lnt+1
        vpair=-1
        while lendfl == False:
            if(lnt==NNT-1):
                UNSTLS.append(tloop)
                tloop=[]
                lendfl=True
                break
            tloop.append(lnt)
            if(PAIRLIST[lnt]==-1):
                lnt=lnt+1
            else:
                UNSTLS.append(tloop)
                tloop=[]
                lendfl=True

#now we join the loops
#USE THE FULL UNSTRUCTURED LOOPS FOR GAUSSIAN INTEGRALS
FULLUNSTLS=[]
for iul in xrange(0,len(UNSTLS)):
    uloop=UNSTLS[iul]
    tloop=uloop
    lstart=uloop[0]
    cend=uloop[len(uloop)-1]
    lendfl=False
    while lendfl == False:
        #find next
        hasnext=False
        for nloop in UNSTLS:
            nstart=nloop[0]
            
            if(cend==PAIRLIST[nstart]):
                floop=nloop
                hasnext=True
                break
        if hasnext is False:
            tloop=[]
            lendfl=True
            break
        else:
            tloop=tloop+floop
        cend=floop[len(floop)-1]
        if(cend==PAIRLIST[lstart]):
            if(sorted(tloop) not in FULLUNSTLS):
                FULLUNSTLS.append(sorted(tloop))
            tloop=[]
            lendfl=True
            break
if verboseflag:
    FILENAME=OUTPUT+"_LOOPS.dat"
    tfile=open(FILENAME,"w")
    tfile.write( "UNSTRUCTURED STRANDS\n")
    for lo in UNSTLS:
        for nt in lo:
            tfile.write(str(nt)+" ")
        tfile.write("\n")
    tfile.write( "UNSTRUCTURED LOOPS\n")
    for lo in FULLUNSTLS:
        for nt in lo:
            tfile.write(str(nt)+" ")
        tfile.write("\n")
    tfile.write( "HAIRPINS\n")
    for lo in HAIRPINS:
        for nt in lo:
            tfile.write(str(nt)+" ")
        tfile.write("\n")
    tfile.write("STEMS\n")
    for lo in STEMS:
        for nt in lo:
            tfile.write(str(nt)+" ")
        tfile.write("\n")


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
FULLULCOORDS=[]

for ful in FULLUNSTLS:
    #here we exclude the first and last unstructured loops if they are free strands
#    if (0 not in ful and NNT-1 not in ful):
    #the loops here are all closed. no worries
    ulc=[]
    fflag=True
    for nt in xrange(0,len(ful)-1):
        if(fflag==True):
            if(nt>0):
                ulc.append(np.array(COORDS[ful[nt-1]][0]))
            ulc.append(np.array(COORDS[ful[nt]][0]))
            ulc.append(np.array(COORDS[ful[nt]][1]))


        else:
            ulc.append(np.array(COORDS[ful[nt]][2]))
            ulc.append(np.array(COORDS[ful[nt]][1]))
        if(ful[nt+1]-ful[nt]==1):
            fflag=False
        else:
            fflag=True
    ulc.append(np.array(COORDS[ful[len(ful)-1]][2]))
    ulc.append(np.array(COORDS[ful[len(ful)-1]][1]))
    ulc.append(np.array(COORDS[ful[len(ful)-1]][0]))
    ulc.append(np.array(COORDS[ful[0]][0]))
    FULLULCOORDS.append(ulc)


for ul in UNSTLS:
    ulc=[]
    ulc.append(np.array(COORDS[ul[0]][0]))
    ulc.append(np.array(COORDS[ul[0]][1]))
    for nt in xrange(1,len(ul)):
        ulc.append(np.array(COORDS[ul[nt]][2]))
        ulc.append(np.array(COORDS[ul[nt]][1]))
    ulc.append(np.array(COORDS[ul[len(ul)-1]][0]))
    #ulc.append(np.array(COORDS[ul[0]][0]))
    ULCOORDS.append(ulc)

for hp in HAIRPINS:
    hpc=[]
    hpc.append(np.array(COORDS[hp[0]][0]))
    hpc.append(np.array(COORDS[hp[0]][1]))
    for nt in xrange(1,len(hp)):
        hpc.append(np.array(COORDS[hp[nt]][2]))
        hpc.append(np.array(COORDS[hp[nt]][1]))
    hpc.append(np.array(COORDS[hp[len(hp)-1]][0]))
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

if verboseflag:
    lcnt=0
    cnt=0
    for loo in FULLULCOORDS:
        FILENAME=OUTPUT+"_LOOP_JL"+str(cnt)+".xyz"
        tfile=open(FILENAME,"w")
        tfile.write(str( len(loo))+"\n")
        tfile.write("Junction"+" "+str( cnt)+"\n")
        for nt in loo:
            line="C "+str( nt[0])+" "+str( nt[1])+" "+str( nt[2])+"\n"
            tfile.write(line)

        cnt,lcnt=cnt+1,lcnt+1
        tfile.close()

    cnt=0
    for loo in ULCOORDS:
        FILENAME=OUTPUT+"_LOOP_UL"+str(cnt)+".xyz"
        tfile=open(FILENAME,"w")
        tfile.write(str( len(loo))+"\n")
        tfile.write("Unstructured"+" "+str( cnt)+"\n")
        for nt in loo:
            line="C "+str( nt[0])+" "+str( nt[1])+" "+str( nt[2])+"\n"
            tfile.write(line)

        cnt,lcnt=cnt+1,lcnt+1
        tfile.close()

    cnt=0
    for loo in HPCOORDS:
        FILENAME=OUTPUT+"_LOOP_HP"+str(cnt)+".xyz"
        tfile=open(FILENAME,"w")
        tfile.write(str( len(loo))+"\n")
        tfile.write("Hairpin"+" "+str( cnt)+"\n")
        for nt in loo:
             line="C "+str( nt[0])+" "+str( nt[1])+" "+str( nt[2])+"\n"
             tfile.write(line)
        cnt,lcnt=cnt+1,lcnt+1
    cnt=0
    for loo in STCOORDS:
        FILENAME=OUTPUT+"_LOOP_ST"+str(cnt)+".xyz"
        tfile=open(FILENAME,"w")
        tfile.write(str( len(loo))+"\n")
        tfile.write("Stem"+" "+str( cnt)+"\n")
        for nt in loo:
             line="C "+str( nt[0])+" "+str( nt[1])+" "+str( nt[2])+"\n"
             tfile.write(line)
        cnt,lcnt=cnt+1,lcnt+1
    
    exit(1)
ALLINDEXES=[]
if calctype=="g":
    print "Calculating with Gaussian integrals"
    LOOPCOORDS=[]
    for hp in xrange(0,len(HPCOORDS)):
        LOOPCOORDS.append([HPCOORDS[hp],["hairpin ",hp]])
    for st in xrange(0,len(STCOORDS)):
        LOOPCOORDS.append([STCOORDS[st],["stem ", st]])
    for ul in xrange(0,len(FULLULCOORDS)):
        LOOPCOORDS.append([FULLULCOORDS[ul],["unstl ", ul]])
    SEGF=lambda t1,t2 : np.dot(((p1- q1)* t1+ q1- ((p2- q2)*t2+ q2)),(np.cross((p1- q1), (p2- q2))))/np.linalg.norm((p1- q1)*t1+ q1-((p2- q2)*t2+ q2))**3

    for iloop1 in xrange(0,len(LOOPCOORDS)):
        for iloop2 in xrange(iloop1+1,len(LOOPCOORDS)):
    #for iloop1 in xrange(0,0):
    #    for iloop2 in xrange(iloop1+1,0):
            #print "Loop " + str(iloop1)+" and " +str(iloop2)
            GAUSSINT=0
            dist,sameflag=diloops(LOOPCOORDS[iloop1][0],LOOPCOORDS[iloop2][0]),difloops(LOOPCOORDS[iloop1][0],LOOPCOORDS[iloop2][0])
            
            if(dist and sameflag):
                #print "Loop " + str(iloop1)+","+str(LOOPCOORDS[iloop1][1])+" and " +str(iloop2)+" + "+str(LOOPCOORDS[iloop2][1])+" with "+str(dist)+" "+str(sameflag)
                for pa1 in xrange(0,len(LOOPCOORDS[iloop1][0])-1):
                    for pa2 in xrange(0,len(LOOPCOORDS[iloop2][0])-1):
                        p1,q1,p2,q2=LOOPCOORDS[iloop1][0][pa1+1],LOOPCOORDS[iloop1][0][pa1],LOOPCOORDS[iloop2][0][pa2+1],LOOPCOORDS[iloop2][0][pa2]
                        #print p1, q1, p2, q2
                        GAUSSINT=GAUSSINT+integrate.dblquad(SEGF,0,1,lambda x: 0, lambda y: 1)[0]
                GAUSSINT=GAUSSINT/(4.0*pi)
                #print GAUSSINT
                if(fabs(GAUSSINT)>0.5):
                    type1,looi1=LOOPCOORDS[iloop1][1][0],LOOPCOORDS[iloop1][1][1]
                    type2,looi2=LOOPCOORDS[iloop2][1][0],LOOPCOORDS[iloop2][1][1]
                    indexes=[]
                    if(type1=="hairpin "):
                        indexes.append([HAIRPINS[looi1],"hp"])
                    elif(type1=="stem "):
                        indexes.append([STEMS[looi1],"st"])
                    elif(type1=="unstl "):
                        indexes.append([FULLUNSTLS[looi1],"il"])
                    if(type2=="hairpin "):
                        indexes.append([HAIRPINS[looi2],"hp"])
                    elif(type2=="stem "):
                        indexes.append([STEMS[looi2],"st"])
                    elif(type2=="unstl "):
                        indexes.append([FULLUNSTLS[looi2],"il"])
                    ALLINDEXES.append(indexes)
                    floop1=indexes[0]
                    floop2=indexes[1]
                    print "LINK DETECTED: loops "+str(floop1)+" ( "+type1+" )"+" and "+str( floop2)+" ( "+type2+")"+". Gauss integral = "+str(GAUSSINT)

    print "Link analysis finished"

#print "Calculating pierces"
#check piercings
PIERCEINDEXES=[]
REDUNDANTINDEXES=[]

if calctype=="p":
    LOOPCOORDS1,LOOPINDEXS1=[],[]
    LOOPCOORDS2,LOOPINDEXS2=[],[]

    for hp in xrange(0,len(HPCOORDS)):
        LOOPCOORDS1.append([HPCOORDS[hp],["hairpin ",hp]])
        LOOPINDEXS1.append(HAIRPINS[hp])
        LOOPCOORDS2.append([HPCOORDS[hp],["hairpin ",hp]])
        LOOPINDEXS2.append(HAIRPINS[hp])
    for st in xrange(0,len(STCOORDS)):
        LOOPCOORDS1.append([STCOORDS[st],["stem ", st]])
        LOOPINDEXS1.append(STEMS[st])
        LOOPCOORDS2.append([STCOORDS[st],["stem ", st]])
        LOOPINDEXS2.append(STEMS[st])
    for ul in xrange(0,len(ULCOORDS)):
        LOOPCOORDS2.append([ULCOORDS[ul],["unst ", ul]])
        LOOPINDEXS2.append(UNSTLS[ul])

    for iloop1 in xrange(0,len(LOOPCOORDS1)):
        trianglist=make_triangles(LOOPCOORDS1[iloop1][0])
        for iloop2 in xrange(0,len(LOOPCOORDS2)):
            if(iloop1!= iloop2):
                totpierce=0
                for triang in trianglist:
                    npierce=check_pierce(triang,LOOPCOORDS2[iloop2][0])
                    totpierce=totpierce+npierce
                if totpierce>0:
                    type1,looi1=LOOPCOORDS1[iloop1][1][0],LOOPCOORDS1[iloop1][1][1]
                    type2,looi2=LOOPCOORDS2[iloop2][1][0],LOOPCOORDS2[iloop2][1][1]
                    indexes=[]
                    if(type1=="hairpin "):
                        indexes.append([HAIRPINS[looi1],"hp"])
                    elif(type1=="stem "):
                        indexes.append([STEMS[looi1],"st"])
                    elif(type1=="unst "):
                        indexes.append([UNSTLS[looi1],"il"])
                    if(type2=="hairpin "):
                        indexes.append([HAIRPINS[looi2],"hp"])
                    elif(type2=="stem "):
                        indexes.append([STEMS[looi2],"st"])
                    elif(type2=="unst "):
                        indexes.append([UNSTLS[looi2],"il"])
                    #if(iloop1<iloop2):
                    #PIERCEINDEXES.append(indexes)
                    REDUNDANTINDEXES.append(indexes)
                            
for pair in REDUNDANTINDEXES:
    if(pair not in PIERCEINDEXES and [pair[1],pair[0]] not in PIERCEINDEXES):
        PIERCEINDEXES.append(pair)
        print "LINK DETECTED: ("+str(pair[0][1])+") : "+str(pair[0][0])+ "  and  ("+str(pair[1][1])+") : "+str(pair[1][0])
NLINKS,NPIERCES=len(ALLINDEXES),len(PIERCEINDEXES)
LINKOUTPUT=[]
if  NPIERCES>0 and calctype=="p":
    for li in PIERCEINDEXES:
        #typs=li[0][1]+li[1][1]
        loop1,loop2=[],[]
        if li[0][1]=='il':
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
        if li[1][1]=='il':
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
        #print typs+ " "+' '.join(str(x) for x in loop1) + " "+' '.join(str(x) for x in loop2)
        LINKOUTPUT.append([[li[0][1],' '.join(str(x) for x in loop1)],[li[1][1],' '.join(str(x) for x in loop2)]])

if(calctype=="p"):
    OUTFILE=open(LINKFILE,"w")
    orig_stdout=sys.stdout
    sys.stdout=OUTFILE
    LOOPLIST=[]
    cnt=0
    for li in LINKOUTPUT:
        for loo in li:
            if loo not in LOOPLIST:
                LOOPLIST.append(loo)
                cnt=cnt+1
    #print LINKOUTPUT
    print str(NPIERCES)+" "+str(len(LOOPLIST))+" "+str(KL)+" "+str(KL)
    for loo in LOOPLIST:
        print loo[0], loo[1]
    for li in LINKOUTPUT:
        print LOOPLIST.index(li[0]),LOOPLIST.index(li[1])
    OUTFILE.close()
    sys.stdout=orig_stdout

