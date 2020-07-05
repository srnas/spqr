from math import sqrt,pi,cos,sin
from copy import deepcopy
import numpy as np
import argparse
import sys
import random as rd
rd.seed()
#DUPLEX=False
KERMSD=50
RERMSD=100
parser=argparse.ArgumentParser()
parser.add_argument("-s","--sequences", help="Sequences",type=str,default="")
parser.add_argument("-t","--sstruct", help="Secondary structure",type=str,default="")
parser.add_argument("-c","--centers", help="Strand centers",type=str,default="")
parser.add_argument("-g","--forgi", help="Forgi vectors: coords and twists",type=str,default="")
parser.add_argument("-o","--output", help="Output name",type=str,default="init")
parser.add_argument("-p","--pairs", help="No stacks, only pairs",action='store_true',default=False)
parser.add_argument("-d","--duplex", help="Makes duplex with complementary structure",action='store_true',default=False)

args=parser.parse_args()
ssflag=False
PFLAG=args.pairs
DUPLEX=True
FORGIFLAG=False
FORGICG=args.forgi
PRDUPLEX=args.duplex
splseq=args.sequences.split("&")
splsst=args.sstruct.split("&")
OUTNAME=args.output

if(args.sequences==""):
    print "ERROR: a sequence must be entered."
    parser.print_help()
    exit(1)
NSTRANDS=len(splseq)
if NSTRANDS > 1:
    print "ERROR: More than one strand not supported yet. Please try splitting your job."
    exit(1)
SSNSTRAN=len(splsst)
CENTERS,ZCENTERS=[],[]
for st in xrange(0,NSTRANDS):
    ZCENTERS.append([0,0,0])

if(args.centers==""):
    for st in xrange(0,NSTRANDS):
        CENTERS.append([0,0,0])
else:
    CENCOORDS=args.centers.split()
    if(len(CENCOORDS)!=3*NSTRANDS):
        print "Number of centers "+str(len(CENCOORDS))+ " inconsistent with number of coordinates of strand centers "+str(3*NSTRANDS)+ " !"
        exit(1)
    else:
        for st in xrange(0,NSTRANDS):
            CENTERS.append([float(CENCOORDS[st*3]),float(CENCOORDS[st*3+1]),float(CENCOORDS[st*3+2])])


FGCOORDS,FGTWISTS=[np.array([0,0,0]),np.array([0,0,0])],[np.array([0,0,0]),np.array([0,0,0])]
fgline=FORGICG.split()
if FORGICG!="":
    if(len(fgline)!=12):
        print "Incorrect number of coordinates for forgi coordinates (must be 12)."
        exit(1)
    else:
        FORGIFLAG=True
        FGCOORDS[0]=np.array([float(fgline[0]),float(fgline[1]),float(fgline[2])])
        FGCOORDS[1]=np.array([float(fgline[3]),float(fgline[4]),float(fgline[5])])
        FGTWISTS[0]=np.array([float(fgline[6]),float(fgline[7]),float(fgline[8])])
        FGTWISTS[1]=np.array([float(fgline[9]),float(fgline[10]),float(fgline[11])])
        #FGCM=(FGCOORDS[0]+FGCOORDS[1])/2.0
        FGCM=FGCOORDS[0]
        FGAXIS=(FGCOORDS[1]-FGCOORDS[0])
        FGAXIS=FGAXIS/sqrt(np.dot(FGAXIS,FGAXIS))
        FORGIVECS=[FGCM,FGAXIS,FGTWISTS[0],FGTWISTS[1]]

if(args.sstruct!=""):
    ssflag=True
errfl=False
if(SSNSTRAN!=NSTRANDS):
    errfl=True
if(not errfl):
    for ii in xrange(0,NSTRANDS):
        if(len(splseq[ii])!=len(splsst[ii])):
            errfl=True
if(errfl and len(args.sstruct)>0):
    print "ERROR: Number of nucleotides must be consistent between secondary structure and sequence!"
    parser.print_help()
    exit(1)
if PRDUPLEX and ssflag :
    print "ERROR:It is not recommended to use options -d and -t simultaneously."
    exit(1)

seq=[]
strandseq=[]
fullss=[]
rawseq=[]
for i in xrange(0,len(args.sequences)):
    if(args.sequences[i]!="&"):
        strandseq.append(args.sequences[i])
        rawseq.append(args.sequences[i])
        if(ssflag):
            fullss.append(args.sstruct[i])
    if(args.sequences[i]=="&"):
        seq.append(strandseq)
        strandseq=[]
seq.append(strandseq)

DIM=3
NAT=5
NBA=4
STTOT=2
if (DUPLEX):
    STTOT=4
SSDICT = {
    "(":")",
    "[":"]",
    "{":"}",
    "<":">"
}

#sequence is arg1 - we read all the stacks
##########FUNCTIONS #################
def get_forgi_matrix(leftstrand, rightstrand, forgivecs,typ):
    #first we calculate the cm
    CM=[0,0,0]
    for nt in leftstrand:
        #print nt
        for at in nt[0]:
            CM[2]=CM[2]+at[2]
    for nt in rightstrand:
        for at in nt[0]:
            CM[2]=CM[2]+at[2]
    for d in xrange(0,DIM):
        CM[d]=CM[d]/(2.0*len(leftstrand)*NAT)
#        FORGIVECS=[FGCM,FGAXIS,FGTWISTS[0],FGTWISTS[1]]
    t0,   CM,b0=np.array([0,0,1]),np.array(CM),twistdict[typ]
    tgtCM,t1,b1=forgivecs[0],forgivecs[1],forgivecs[2]
    c0,c1=np.cross(t0,b0),np.cross(t1,b1)
    B1cart= np.array([[t0[0],b0[0],c0[0]],
                      [t0[1],b0[1],c0[1]],
                      [t0[2],b0[2],c0[2]]])
    rotB1 = np.array([[np.dot(t0,t1),np.dot(t0,b1),np.dot(t0,c1)],
                      [np.dot(b0,t1),np.dot(b0,b1),np.dot(b0,c1)],
                      [np.dot(c0,t1),np.dot(c0,b1),np.dot(c0,c1)]])
    #otB1=np.array([[1,0,0],[0,1,0],[0,0,1]])
    cartB1= np.array([[t0[0],t0[1],t0[2]],
                      [b0[0],b0[1],b0[2]],
                      [c0[0],c0[1],c0[2]]])
    #print np.dot(np.dot(B1cart,rotB1),cartB1)
    #return [np.dot(np.dot(B1cart,rotB1),cartB1),tgtCM-CM]
    return [np.dot(np.dot(B1cart,rotB1),cartB1),tgtCM]


def forgi_arrange(fgparams, ntc):
    nnt=[]
    for at in ntc:
        #print at
        cat=np.array(at)
        nnt.append(np.dot(fgparams[0],cat)+fgparams[1])
    #exit(1)
    return nnt

def get_bname(nt):
    ret="X"
    if(nt==0):
        ret="a"
    elif(nt==1):
        ret="u"
    elif(nt==2):
        ret="g"
    elif(nt==3):
        ret="c"
    else:
        print "Base type not recognized"
    return ret

def base_index(base):
    ret=-1
    if(base=='A' or base=='a'):
        ret=0
    elif(base=='U' or base=='u'):
        ret=1
    elif(base=='G' or base=='g'):
        ret=2
    elif(base=='C' or base=='c'):
        ret=3
    else:
        print "Base not recognized!"
    return ret

def get_complb(bas):
    ret="X"
    if(bas=="A"):
        ret="U"
    elif(bas=="U"):
        ret="A"
    elif(bas=="G"):
        ret="C"
    elif(bas=="C"):
        ret="G"
    else:
        print "Base not recognized!"
    return ret

def normalize(vec):
    norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
    return [vec[0]/norm, vec[1]/norm, vec[2]/norm]

def pdbprint(nt, resind, bas, chain,center):
    resname='{:3}'.format(bas)
    #chindex='{:1}'.format(chain)
    chindex="0"
    resindex='{:4}'.format(resind)
    record="ATOM  "
    xc,yc,zc=center[0],center[1],center[2]
    atindex='{:5}'.format(resind*NAT)
    atname="BASE"
    xtemp, ytemp, ztemp=nt[0][0],nt[0][1],nt[0][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp+xc),'{:8.3f}'.format(ytemp+yc),'{:8.3f}'.format(ztemp+zc)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*NAT+1)
    atname="XVEC"
    xtemp, ytemp, ztemp=nt[1][0],nt[1][1],nt[1][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp+xc),'{:8.3f}'.format(ytemp+yc),'{:8.3f}'.format(ztemp+zc)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*NAT+2)
    atname="YVEC"
    xtemp, ytemp, ztemp=nt[2][0],nt[2][1],nt[2][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp+xc),'{:8.3f}'.format(ytemp+yc),'{:8.3f}'.format(ztemp+zc)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*NAT+3)
    atname="SUGR"
    xtemp, ytemp, ztemp=nt[3][0],nt[3][1],nt[3][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp+xc),'{:8.3f}'.format(ytemp+yc),'{:8.3f}'.format(ztemp+zc)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*NAT+4)
    atname="PHOS"
    xtemp, ytemp, ztemp=nt[4][0],nt[4][1],nt[4][2]
    nidx, nidy, nidz=len(str(int(xtemp))),len(str(int(ytemp))),len(str(int(ztemp)))
    xpos,ypos,zpos='{:8.3f}'.format(xtemp+xc),'{:8.3f}'.format(ytemp+yc),'{:8.3f}'.format(ztemp+zc)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    return

def cross_prod(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

def get_basis(ntcoord):
    ori=[ntcoord[0][0],ntcoord[0][1],ntcoord[0][2]]
    xat=[ntcoord[1][0],ntcoord[1][1],ntcoord[1][2]]
    yat=[ntcoord[2][0],ntcoord[2][1],ntcoord[2][2]]
    xvec=[xat[0]-ori[0],xat[1]-ori[1],xat[2]-ori[2]]
    yvec=[yat[0]-ori[0],yat[1]-ori[1],yat[2]-ori[2]]
    xvec,yvec=normalize(xvec),normalize(yvec)
    zvec=cross_prod(xvec,yvec)
    zvec=normalize(zvec)
    return [ori,xvec,yvec,zvec]

def get_stacked_basis(refcoords, current_base, ntdo, ntup,info,wchain,wchain2,DUPLEX):
    new_base=deepcopy(get_basis(refcoords[ntdo][ntup][0]))
    
    #in chain 0
    new_stack=deepcopy(refcoords[ntdo][ntup][1])
    #translation to origin
    for a in xrange(0,NAT):
        for d in xrange(0,DIM):
            new_stack[a][d]=new_stack[a][d]-new_base[0][d]
    #rotation
    nnew_stack=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    for a in xrange(0,NAT):
        for i in xrange(0,DIM):
            for j in xrange(0,DIM):
                for k in xrange(0,DIM):
                    nnew_stack[a][i]=nnew_stack[a][i]+current_base[j+1][i]*new_base[j+1][k]*new_stack[a][k]
    #new translation
    for a in xrange(0,NAT):
        for d in xrange(0,DIM):
            nnew_stack[a][d]=nnew_stack[a][d]+current_base[0][d]

    #new rotation, for growth
    tppos=np.array(nnew_stack[4])
    tempstack=[np.array(nnew_stack[0])-tppos,np.array(nnew_stack[1])-tppos,np.array(nnew_stack[2])-tppos,np.array(nnew_stack[3])-tppos,np.array(nnew_stack[4])-tppos]
    randflag=0
    while randflag==0:
        randflag=1
        rand_axis,rand_angle=np.array([2*rd.random()-1,2*rd.random()-1,2*rd.random()-1]),2*pi*rd.random()
        rand_axis=rand_axis/sqrt(np.dot(rand_axis,rand_axis))
        ct,st,ux,uy,uz=cos(rand_angle),sin(rand_angle),rand_axis[0],rand_axis[1],rand_axis[2]
        rand_matrix=np.array([[ ct+ux*ux*(1-ct)   , ux*uy*(1-ct)-uz*st , ux*uz*(1-ct)+uy*st], [ uy*ux*(1-ct)+uz*st , ct+uy*uy*(1-ct)    , uy*uz*(1-ct)-ux*st], [ uz*ux*(1-ct)-uy*st , uz*uy*(1-ct)+ux*st , ct+uz*uz*(1-ct)   ]] )
        rrotnt=[np.dot(rand_matrix,tempstack[0])+tppos,np.dot(rand_matrix,tempstack[1])+tppos,np.dot(rand_matrix,tempstack[2])+tppos,np.dot(rand_matrix,tempstack[3])+tppos,np.dot(rand_matrix,tempstack[4])+tppos]
        for bl in wchain:
            prevnt=bl[0]
            co1,co2,co3=np.array(prevnt[0]),np.array(prevnt[3]),np.array(prevnt[4])
            cm1,cm2=(co1+co2+co3)/3.0,(rrotnt[0]+rrotnt[3]+rrotnt[4])/3.0
            if(np.dot(cm1-cm2,cm1-cm2)<25):
                randflag=0
            
    nnew_stack[0]=rrotnt[0].tolist()
    nnew_stack[1]=rrotnt[1].tolist()
    nnew_stack[2]=rrotnt[2].tolist()
    nnew_stack[3]=rrotnt[3].tolist()
    nnew_stack[4]=rrotnt[4].tolist()
    
    wchain.append([nnew_stack, info[0]])
    ret=get_basis(nnew_stack)
    
    ##DUPLEX
    if(DUPLEX):
    #we do the same but at the other chain, that is, index 2
        new_stack=deepcopy(refcoords[ntdo][ntup][2])
        #translation to origin
        for a in xrange(0,NAT):
            for d in xrange(0,DIM):
                new_stack[a][d]=new_stack[a][d]-new_base[0][d]
        #rotation
        nnew_stack=[[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
        for a in xrange(0,NAT):
            for i in xrange(0,DIM):
                for j in xrange(0,DIM):
                    for k in xrange(0,DIM):
                        nnew_stack[a][i]=nnew_stack[a][i]+current_base[j+1][i]*new_base[j+1][k]*new_stack[a][k]
        #new translation
        for a in xrange(0,NAT):
            for d in xrange(0,DIM):
                nnew_stack[a][d]=nnew_stack[a][d]+current_base[0][d]
        wchain2.append([nnew_stack, info[1]])
    return ret

#READ THE SEQUENCE AND PRINT THE COORDINATES
def get_chain(sequ,initnt,refcoords,strand):
    chain1,chain2=[],[]
    nnt=len(sequ)
    nst=nnt-1
    info=[]
    #print first nt
    nt=base_index(sequ[0])
    current_base=get_basis(refcoords[nt][0][0])
    info=[[initnt, sequ[0].upper(),strand], [2*nnt-1-initnt,get_complb(sequ[0].upper()),strand+1]]
    chain1.append([refcoords[nt][0][0], info[0]])
    if(DUPLEX):
        #info.append([nnt, get_complb(sequ[0].upper()),strand])
        chain2.append([refcoords[nt][0][3], info[1]])
    initnt=initnt+1
    #relabel duplex chain index
    for st in xrange(0,nst):
        print "stack", st
        stack=sequ[st:st+2]
        info=[[initnt, stack[1].upper(), strand],[2*nnt-1-initnt,get_complb(stack[1].upper()),strand+1]]
        current_base=get_stacked_basis(refcoords,current_base, base_index(stack[0]), base_index(stack[1]), info, chain1, chain2, DUPLEX)
        initnt=initnt+1
    
    return [chain1,initnt,chain2]

############################STACK COORDINATES#####################################            

FILEAA=["ATOM      0 BASE A   0   0       5.454  -1.762   0.630",
        "ATOM      1 XVEC A   0   0       6.246  -1.156   0.559",
        "ATOM      2 YVEC A   0   0       4.856  -1.014   0.341",
        "ATOM      3 SUGR A   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS A   0   0       2.847  -8.829   3.227",
        "ATOM      5 BASE A   0   1       5.541   1.463   3.440",
        "ATOM      6 XVEC A   0   1       5.881   2.401   3.369",
        "ATOM      7 YVEC A   0   1       4.635   1.770   3.151",
        "ATOM      8 SUGR A   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS A   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE U   1   2       1.974   5.444   1.783",
        "ATOM     11 XVEC U   1   2       2.955   5.639   1.774",
        "ATOM     12 YVEC U   1   2       2.162   4.508   2.080",
        "ATOM     13 SUGR U   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS U   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE U   1   3       4.602   3.515  -1.027",
        "ATOM     16 XVEC U   1   3       5.533   3.149  -1.036",
        "ATOM     17 YVEC U   1   3       4.256   2.626  -0.730",
        "ATOM     18 SUGR U   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS U   1   3       2.907   8.210  -3.750"]

FILEAC=["ATOM      0 BASE A   0   0       5.454  -1.762   0.630",
        "ATOM      1 XVEC A   0   0       6.246  -1.156   0.559",
        "ATOM      2 YVEC A   0   0       4.856  -1.014   0.341",
        "ATOM      3 SUGR A   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS A   0   0       2.847  -8.829   3.227",
        "ATOM      5 BASE C   0   1       5.791  -0.441   3.840",
        "ATOM      6 XVEC C   0   1       6.365   0.378   3.847",
        "ATOM      7 YVEC C   0   1       5.010   0.109   3.541",
        "ATOM      8 SUGR C   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS C   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE G   1   2       3.674   4.432   2.180",
        "ATOM     11 XVEC G   1   2       4.671   4.403   2.239",
        "ATOM     12 YVEC G   1   2       3.629   3.477   2.472",
        "ATOM     13 SUGR G   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS G   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE U   1   3       4.602   3.515  -1.027",
        "ATOM     16 XVEC U   1   3       5.533   3.149  -1.036",
        "ATOM     17 YVEC U   1   3       4.256   2.626  -0.730",
        "ATOM     18 SUGR U   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS U   1   3       2.907   8.210  -3.750"]

FILEAG=["ATOM      0 BASE A   0   0       5.454  -1.762   0.630",
        "ATOM      1 XVEC A   0   0       6.246  -1.156   0.559",
        "ATOM      2 YVEC A   0   0       4.856  -1.014   0.341",
        "ATOM      3 SUGR A   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS A   0   0       2.847  -8.829   3.227",
        "ATOM      5 BASE G   0   1       5.559   1.495   3.440",
        "ATOM      6 XVEC G   0   1       5.948   2.414   3.381",
        "ATOM      7 YVEC G   0   1       4.672   1.852   3.148",
        "ATOM      8 SUGR G   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS G   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE C   1   2       2.010   5.449   1.780",
        "ATOM     11 XVEC C   1   2       2.994   5.629   1.773",
        "ATOM     12 YVEC C   1   2       2.184   4.511   2.079",
        "ATOM     13 SUGR C   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS C   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE U   1   3       4.602   3.515  -1.027",
        "ATOM     16 XVEC U   1   3       5.533   3.149  -1.036",
        "ATOM     17 YVEC U   1   3       4.256   2.626  -0.730",
        "ATOM     18 SUGR U   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS U   1   3       2.907   8.210  -3.750"]

FILEAU=["ATOM      0 BASE A   0   0       5.454  -1.762   0.630",
        "ATOM      1 XVEC A   0   0       6.246  -1.156   0.559",
        "ATOM      2 YVEC A   0   0       4.856  -1.014   0.341",
        "ATOM      3 SUGR A   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS A   0   0       2.847  -8.829   3.227",
        "ATOM      5 BASE U   0   1       5.772  -0.472   3.837",
        "ATOM      6 XVEC U   0   1       6.357   0.339   3.846",
        "ATOM      7 YVEC U   0   1       4.999   0.089   3.540",
        "ATOM      8 SUGR U   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS U   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE A   1   2       3.638   4.429   2.180",
        "ATOM     11 XVEC A   1   2       4.632   4.347   2.251",
        "ATOM     12 YVEC A   1   2       3.539   3.477   2.469",
        "ATOM     13 SUGR A   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS A   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE U   1   3       4.602   3.515  -1.027",
        "ATOM     16 XVEC U   1   3       5.533   3.149  -1.036",
        "ATOM     17 YVEC U   1   3       4.256   2.626  -0.730",
        "ATOM     18 SUGR U   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS U   1   3       2.907   8.210  -3.750"]

FILECA=["ATOM      0 BASE C   0   0       4.635  -3.499   1.030",
        "ATOM      1 XVEC C   0   0       5.560  -3.119   1.037",
        "ATOM      2 YVEC C   0   0       4.275  -2.616   0.731",
        "ATOM      3 SUGR C   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS C   0   0       2.756  -8.806   2.933",
        "ATOM      5 BASE A   0   1       5.541   1.463   3.440",
        "ATOM      6 XVEC A   0   1       5.881   2.401   3.369",
        "ATOM      7 YVEC A   0   1       4.635   1.770   3.151",
        "ATOM      8 SUGR A   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS A   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE U   1   2       1.974   5.444   1.783",
        "ATOM     11 XVEC U   1   2       2.955   5.639   1.774",
        "ATOM     12 YVEC U   1   2       2.162   4.508   2.080",
        "ATOM     13 SUGR U   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS U   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE G   1   3       5.486   1.745  -0.630",
        "ATOM     16 XVEC G   1   3       6.310   1.182  -0.571",
        "ATOM     17 YVEC G   1   3       4.932   0.966  -0.338",
        "ATOM     18 SUGR G   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS G   1   3       2.907   8.210  -3.750"]

FILECC=["ATOM      0 BASE C   0   0       4.635  -3.499   1.030",
        "ATOM      1 XVEC C   0   0       5.560  -3.119   1.037",
        "ATOM      2 YVEC C   0   0       4.275  -2.616   0.731",
        "ATOM      3 SUGR C   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS C   0   0       2.756  -8.806   2.933",
        "ATOM      5 BASE C   0   1       5.791  -0.441   3.840",
        "ATOM      6 XVEC C   0   1       6.365   0.378   3.847",
        "ATOM      7 YVEC C   0   1       5.010   0.109   3.541",
        "ATOM      8 SUGR C   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS C   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE G   1   2       3.674   4.432   2.180",
        "ATOM     11 XVEC G   1   2       4.671   4.403   2.239",
        "ATOM     12 YVEC G   1   2       3.629   3.477   2.472",
        "ATOM     13 SUGR G   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS G   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE G   1   3       5.486   1.745  -0.630",
        "ATOM     16 XVEC G   1   3       6.310   1.182  -0.571",
        "ATOM     17 YVEC G   1   3       4.932   0.966  -0.338",
        "ATOM     18 SUGR G   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS G   1   3       2.907   8.210  -3.750"]


FILECG=["ATOM      0 BASE C   0   0       4.635  -3.499   1.030",
        "ATOM      1 XVEC C   0   0       5.560  -3.119   1.037",
        "ATOM      2 YVEC C   0   0       4.275  -2.616   0.731",
        "ATOM      3 SUGR C   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS C   0   0       2.756  -8.806   2.933",
        "ATOM      5 BASE G   0   1       5.559   1.495   3.440",
        "ATOM      6 XVEC G   0   1       5.948   2.414   3.381",
        "ATOM      7 YVEC G   0   1       4.672   1.852   3.148",
        "ATOM      8 SUGR G   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS G   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE C   1   2       2.010   5.449   1.780",
        "ATOM     11 XVEC C   1   2       2.994   5.629   1.773",
        "ATOM     12 YVEC C   1   2       2.184   4.511   2.079",
        "ATOM     13 SUGR C   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS C   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE G   1   3       5.486   1.745  -0.630",
        "ATOM     16 XVEC G   1   3       6.310   1.182  -0.571",
        "ATOM     17 YVEC G   1   3       4.932   0.966  -0.338",
        "ATOM     18 SUGR G   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS G   1   3       2.907   8.210  -3.750"]

FILECU=["ATOM      0 BASE C   0   0       4.635  -3.499   1.030",
        "ATOM      1 XVEC C   0   0       5.560  -3.119   1.037",
        "ATOM      2 YVEC C   0   0       4.275  -2.616   0.731",
        "ATOM      3 SUGR C   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS C   0   0       2.756  -8.806   2.933",
        "ATOM      5 BASE U   0   1       5.772  -0.472   3.837",
        "ATOM      6 XVEC U   0   1       6.357   0.339   3.846",
        "ATOM      7 YVEC U   0   1       4.999   0.089   3.540",
        "ATOM      8 SUGR U   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS U   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE A   1   2       3.638   4.429   2.180",
        "ATOM     11 XVEC A   1   2       4.632   4.347   2.251",
        "ATOM     12 YVEC A   1   2       3.539   3.477   2.469",
        "ATOM     13 SUGR A   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS A   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE G   1   3       5.486   1.745  -0.630",
        "ATOM     16 XVEC G   1   3       6.310   1.182  -0.571",
        "ATOM     17 YVEC G   1   3       4.932   0.966  -0.338",
        "ATOM     18 SUGR G   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS G   1   3       2.907   8.210  -3.750"]

FILEGA=["ATOM      0 BASE G   0   0       5.486  -1.745   0.630",
        "ATOM      1 XVEC G   0   0       6.310  -1.182   0.571",
        "ATOM      2 YVEC G   0   0       4.932  -0.966   0.338",
        "ATOM      3 SUGR G   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS G   0   0       2.481  -8.676   3.165",
        "ATOM      5 BASE A   0   1       5.541   1.463   3.440",
        "ATOM      6 XVEC A   0   1       5.881   2.401   3.369",
        "ATOM      7 YVEC A   0   1       4.635   1.770   3.151",
        "ATOM      8 SUGR A   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS A   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE U   1   2       1.974   5.444   1.783",
        "ATOM     11 XVEC U   1   2       2.955   5.639   1.774",
        "ATOM     12 YVEC U   1   2       2.162   4.508   2.080",
        "ATOM     13 SUGR U   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS U   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE C   1   3       4.635   3.499  -1.030",
        "ATOM     16 XVEC C   1   3       5.560   3.119  -1.037",
        "ATOM     17 YVEC C   1   3       4.275   2.616  -0.731",
        "ATOM     18 SUGR C   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS C   1   3       2.907   8.210  -3.750"]

FILEGC=["ATOM      0 BASE G   0   0       5.486  -1.745   0.630",
        "ATOM      1 XVEC G   0   0       6.310  -1.182   0.571",
        "ATOM      2 YVEC G   0   0       4.932  -0.966   0.338",
        "ATOM      3 SUGR G   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS G   0   0       2.481  -8.676   3.165",
        "ATOM      5 BASE C   0   1       5.791  -0.441   3.840",
        "ATOM      6 XVEC C   0   1       6.365   0.378   3.847",
        "ATOM      7 YVEC C   0   1       5.010   0.109   3.541",
        "ATOM      8 SUGR C   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS C   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE G   1   2       3.674   4.432   2.180",
        "ATOM     11 XVEC G   1   2       4.671   4.403   2.239",
        "ATOM     12 YVEC G   1   2       3.629   3.477   2.472",
        "ATOM     13 SUGR G   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS G   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE C   1   3       4.635   3.499  -1.030",
        "ATOM     16 XVEC C   1   3       5.560   3.119  -1.037",
        "ATOM     17 YVEC C   1   3       4.275   2.616  -0.731",
        "ATOM     18 SUGR C   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS C   1   3       2.907   8.210  -3.750"]

FILEGG=["ATOM      0 BASE G   0   0       5.486  -1.745   0.630",
        "ATOM      1 XVEC G   0   0       6.310  -1.182   0.571",
        "ATOM      2 YVEC G   0   0       4.932  -0.966   0.338",
        "ATOM      3 SUGR G   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS G   0   0       2.481  -8.676   3.165",
        "ATOM      5 BASE G   0   1       5.559   1.495   3.440",
        "ATOM      6 XVEC G   0   1       5.948   2.414   3.381",
        "ATOM      7 YVEC G   0   1       4.672   1.852   3.148",
        "ATOM      8 SUGR G   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS G   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE C   1   2       2.010   5.449   1.780",
        "ATOM     11 XVEC C   1   2       2.994   5.629   1.773",
        "ATOM     12 YVEC C   1   2       2.184   4.511   2.079",
        "ATOM     13 SUGR C   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS C   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE C   1   3       4.635   3.499  -1.030",
        "ATOM     16 XVEC C   1   3       5.560   3.119  -1.037",
        "ATOM     17 YVEC C   1   3       4.275   2.616  -0.731",
        "ATOM     18 SUGR C   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS C   1   3       2.907   8.210  -3.750"]

FILEGU=["ATOM      0 BASE G   0   0       5.486  -1.745   0.630",
        "ATOM      1 XVEC G   0   0       6.310  -1.182   0.571",
        "ATOM      2 YVEC G   0   0       4.932  -0.966   0.338",
        "ATOM      3 SUGR G   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS G   0   0       2.481  -8.676   3.165",
        "ATOM      5 BASE U   0   1       5.772  -0.472   3.837",
        "ATOM      6 XVEC U   0   1       6.357   0.339   3.846",
        "ATOM      7 YVEC U   0   1       4.999   0.089   3.540",
        "ATOM      8 SUGR U   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS U   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE A   1   2       3.638   4.429   2.180",
        "ATOM     11 XVEC A   1   2       4.632   4.347   2.251",
        "ATOM     12 YVEC A   1   2       3.539   3.477   2.469",
        "ATOM     13 SUGR A   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS A   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE C   1   3       4.635   3.499  -1.030",
        "ATOM     16 XVEC C   1   3       5.560   3.119  -1.037",
        "ATOM     17 YVEC C   1   3       4.275   2.616  -0.731",
        "ATOM     18 SUGR C   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS C   1   3       2.907   8.210  -3.750"]

FILEUA=["ATOM      0 BASE U   0   0       4.602  -3.515   1.027",
        "ATOM      1 XVEC U   0   0       5.533  -3.149   1.036",
        "ATOM      2 YVEC U   0   0       4.256  -2.626   0.730",
        "ATOM      3 SUGR U   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS U   0   0       2.641  -8.798   2.914",
        "ATOM      5 BASE A   0   1       5.541   1.463   3.440",
        "ATOM      6 XVEC A   0   1       5.881   2.401   3.369",
        "ATOM      7 YVEC A   0   1       4.635   1.770   3.151",
        "ATOM      8 SUGR A   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS A   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE U   1   2       1.974   5.444   1.783",
        "ATOM     11 XVEC U   1   2       2.955   5.639   1.774",
        "ATOM     12 YVEC U   1   2       2.162   4.508   2.080",
        "ATOM     13 SUGR U   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS U   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE A   1   3       5.454   1.762  -0.630",
        "ATOM     16 XVEC A   1   3       6.246   1.156  -0.559",
        "ATOM     17 YVEC A   1   3       4.856   1.014  -0.341",
        "ATOM     18 SUGR A   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS A   1   3       2.907   8.210  -3.750"]

FILEUC=["ATOM      0 BASE U   0   0       4.602  -3.515   1.027",
        "ATOM      1 XVEC U   0   0       5.533  -3.149   1.036",
        "ATOM      2 YVEC U   0   0       4.256  -2.626   0.730",
        "ATOM      3 SUGR U   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS U   0   0       2.641  -8.798   2.914",
        "ATOM      5 BASE C   0   1       5.791  -0.441   3.840",
        "ATOM      6 XVEC C   0   1       6.365   0.378   3.847",
        "ATOM      7 YVEC C   0   1       5.010   0.109   3.541",
        "ATOM      8 SUGR C   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS C   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE G   1   2       3.674   4.432   2.180",
        "ATOM     11 XVEC G   1   2       4.671   4.403   2.239",
        "ATOM     12 YVEC G   1   2       3.629   3.477   2.472",
        "ATOM     13 SUGR G   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS G   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE A   1   3       5.454   1.762  -0.630",
        "ATOM     16 XVEC A   1   3       6.246   1.156  -0.559",
        "ATOM     17 YVEC A   1   3       4.856   1.014  -0.341",
        "ATOM     18 SUGR A   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS A   1   3       2.907   8.210  -3.750"]

FILEUG=["ATOM      0 BASE U   0   0       4.602  -3.515   1.027",
        "ATOM      1 XVEC U   0   0       5.533  -3.149   1.036",
        "ATOM      2 YVEC U   0   0       4.256  -2.626   0.730",
        "ATOM      3 SUGR U   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS U   0   0       2.641  -8.798   2.914",
        "ATOM      5 BASE G   0   1       5.559   1.495   3.440",
        "ATOM      6 XVEC G   0   1       5.948   2.414   3.381",
        "ATOM      7 YVEC G   0   1       4.672   1.852   3.148",
        "ATOM      8 SUGR G   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS G   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE C   1   2       2.010   5.449   1.780",
        "ATOM     11 XVEC C   1   2       2.994   5.629   1.773",
        "ATOM     12 YVEC C   1   2       2.184   4.511   2.079",
        "ATOM     13 SUGR C   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS C   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE A   1   3       5.454   1.762  -0.630",
        "ATOM     16 XVEC A   1   3       6.246   1.156  -0.559",
        "ATOM     17 YVEC A   1   3       4.856   1.014  -0.341",
        "ATOM     18 SUGR A   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS A   1   3       2.907   8.210  -3.750"]

FILEUU=["ATOM      0 BASE U   0   0       4.602  -3.515   1.027",
        "ATOM      1 XVEC U   0   0       5.533  -3.149   1.036",
        "ATOM      2 YVEC U   0   0       4.256  -2.626   0.730",
        "ATOM      3 SUGR U   0   0       6.783  -5.879   2.766",
        "ATOM      4 PHOS U   0   0       2.641  -8.798   2.914",
        "ATOM      5 BASE U   0   1       5.772  -0.472   3.837",
        "ATOM      6 XVEC U   0   1       6.357   0.339   3.846",
        "ATOM      7 YVEC U   0   1       4.999   0.089   3.540",
        "ATOM      8 SUGR U   0   1       8.884  -1.283   5.576",
        "ATOM      9 PHOS U   0   1       6.882  -5.338   6.560",
        "ATOM     10 BASE A   1   2       3.638   4.429   2.180",
        "ATOM     11 XVEC A   1   2       4.632   4.347   2.251",
        "ATOM     12 YVEC A   1   2       3.539   3.477   2.469",
        "ATOM     13 SUGR A   1   2       2.531   8.611   0.044",
        "ATOM     14 PHOS A   1   2      -1.989   8.480  -0.940",
        "ATOM     15 BASE A   1   3       5.454   1.762  -0.630",
        "ATOM     16 XVEC A   1   3       6.246   1.156  -0.559",
        "ATOM     17 YVEC A   1   3       4.856   1.014  -0.341",
        "ATOM     18 SUGR A   1   3       6.783   5.879  -2.766",
        "ATOM     19 PHOS A   1   3       2.907   8.210  -3.750"]

twistdict={
    "A" : np.array([0.99999, -0.00442147, 0]),
    "U" : np.array([0.99999, 0.00442147, 0]),
    "G" : np.array([0.999993, 0.00370642, 0]),
    "C" : np.array([0.999993, -0.00370642, 0])
}

ALLCFILES=[[FILEAA, FILEAU, FILEAG, FILEAC],[FILEUA, FILEUU, FILEUG, FILEUC],[FILEGA, FILEGU, FILEGG, FILEGC],[FILECA,FILECU, FILECG, FILECC]]

############################HERE STARTS THE CODE##################################
refcoord=[]
for nt1 in xrange(0,NBA):
    wr2=[]
    for nt2 in xrange(0,NBA):
        ALLFILE=[]
        for line in ALLCFILES[nt1][nt2]:
            nam=line[:4]
            if(nam.strip()=="ATOM"):
                ALLFILE.append(line)
        stack=[]
        for nnt in xrange(0,STTOT):
            ntcoords=[]
            for at in xrange(0,NAT):
                x=float(ALLFILE[nnt*NAT+at][30:38].strip())
                y=float(ALLFILE[nnt*NAT+at][38:46].strip())
                z=float(ALLFILE[nnt*NAT+at][46:54].strip())
                ntcoords.append(np.array([x,y,z]))
            stack.append(ntcoords)
        wr2.append(stack)
    refcoord.append(wr2)

#get the chains
currnt=0
ALLCHAINS=[]
for strand in xrange(0,len(seq)):
    lchain,currnt,rchain=get_chain(seq[strand],currnt,refcoord,strand)
    ALLCHAINS.append([lchain,rchain])

pdbfile=open(OUTNAME+".pdb","w")
orig_stdout=sys.stdout
sys.stdout=pdbfile
for strand in xrange(0,len(seq)):
    lnnt=len(seq[strand])
    #print ALLCHAINS[strand][0],ALLCHAINS[strand][1], FORGIVECS, ALLCHAINS[strand][0][0][1][1]
    if(FORGIFLAG):
        forgiparams=get_forgi_matrix(ALLCHAINS[strand][0],ALLCHAINS[strand][1], FORGIVECS, ALLCHAINS[strand][0][0][1][1])
    for nt in xrange(0,lnnt):
        prnt=ALLCHAINS[strand][0][nt][0]
        if(FORGIFLAG):
            prnt=forgi_arrange(forgiparams,prnt)
        pdbprint(prnt, ALLCHAINS[strand][0][nt][1][0],ALLCHAINS[strand][0][nt][1][1],ALLCHAINS[strand][0][nt][1][2],CENTERS[strand])
    if(PRDUPLEX):
        for nt in xrange(0,len(seq[strand])):
            prnt=ALLCHAINS[strand][1][lnnt-1-nt][0]
            if(FORGIFLAG):
                prnt=forgi_arrange(forgiparams,prnt)
            pdbprint(prnt, ALLCHAINS[strand][1][lnnt-1-nt][1][0],ALLCHAINS[strand][1][lnnt-1-nt][1][1],ALLCHAINS[strand][1][lnnt-1-nt][1][2],CENTERS[strand])
sys.stdout=orig_stdout
pdbfile.close()

###############SECONDARY STRUCTURE FILES####################
#we look for the fragments to be built in the ermsd fashion#
#ssf=open("ss.temp","w")
if(ssflag):
    frags=open("ermsd_frags_"+OUTNAME+".lst","w")
    ERMSDFRAGS=[]
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
    
    opairs=bpairs.sort()
    SSSTACKS=[]
    CSTACK=[bpairs[0]]
    for bp in xrange(1,len(bpairs)):
        if(bpairs[bp][0]==bpairs[bp-1][0]+1 and bpairs[bp][1]==bpairs[bp-1][1]-1):
            CSTACK.append(bpairs[bp])
        else:
            if(len(CSTACK)>0):
                SSSTACKS.append(CSTACK)
            CSTACK=[bpairs[bp]]
    if(len(CSTACK)>0):
        SSSTACKS.append(CSTACK)
    
    if PFLAG:
        SSSTACKS=[]
        for pair in bpairs:
            SSSTACKS.append([pair])
    
    
    STACKSEQS=[]
    stra=[]
    for ST in xrange(0,len(SSSTACKS)):
        stra=[]
        for nt in xrange(0,len(SSSTACKS[ST])):
            stra.append(rawseq[SSSTACKS[ST][nt][0]])
        STACKSEQS.append(stra)  
    #for ST in xrange(0,len(STACKSEQS)):
    #    ssf.write(''.join(STACKSEQS[ST]))
    #    ssf.write("\n")

    frags.write( "REMARK ERMSD PARAMS " + str(len(STACKSEQS)) + " "+str(RERMSD)+"\n")
    for ST in xrange(0,len(SSSTACKS)):
        frags.write( "REMARK ERMSD GROUP "+str(KERMSD)+" ")
        #for ind in xrange(0,len(SSSTACKS[ST])):
        #    frags.write(str(SSSTACKS[ST][ind][0])+" ")
        #for ind in xrange(0,len(SSSTACKS[ST])):
        #    frags.write(str(SSSTACKS[ST][ind][1])+" ")
        for ist in SSSTACKS[ST]:
            frags.write(str(ist[0])+" ")
        for ist in reversed(SSSTACKS[ST]):
            frags.write(str(ist[1])+" ")
        frags.write("\n")

    orig_stdout=sys.stdout
    sys.stdout=frags 
    for ST in xrange(0,len(STACKSEQS)):
        lchain,currnt,rchain=get_chain(STACKSEQS[ST],0,refcoord,ST)
        lnnt=len(STACKSEQS[ST])
        for nt in xrange(0,lnnt):
            if not FORGIFLAG:
                tpcoords=lchain[nt][0]
            else:
                tpcoords=lchain[nt][0]
            pdbprint(tpcoords, lchain[nt][1][0],lchain[nt][1][1],lchain[nt][1][2],ZCENTERS[0])

        if(DUPLEX):
            for nt in xrange(0,lnnt):
                if not FORGIFLAG:
                    tpcoords=rchain[lnnt-1-nt][0]
                else:
                    tpcoords=rchain[lnnt-1-nt][0]
                pdbprint(tpcoords, rchain[lnnt-1-nt][1][0], get_complb(lchain[lnnt-1-nt][1][1]), rchain[lnnt-1-nt][1][2],ZCENTERS[0])
                #pdbprint(rchain[lnnt-1-nt][0], rchain[lnnt-1-nt][1][0], get_complb(lchain[lnnt-1-nt][1][1]), rchain[lnnt-1-nt][1][2],ZCENTERS[0])
   #ssf.close()
    frags.close()
    sys_stdout=orig_stdout
