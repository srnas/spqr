from math import sqrt
from copy import deepcopy
import numpy as np
import argparse
import sys

KERMSD=50
RERMSD=100
parser=argparse.ArgumentParser()
parser.add_argument("-g","--forgi", help="Forgi file",type=str,default="")
parser.add_argument("-o","--output", help="Output name",type=str,default="forgi_ermsd.lst")
parser.add_argument("-q","--split", help="Stems will be restricted as separate fragments",action='store_true',default=False)
parser.add_argument("-i","--initnt", help="First nucleotide index",type=int,default=0)

args=parser.parse_args()
ssflag=False
DUPLEX=True
FORGIFLAG=True
FIRSTNT=args.initnt
SPLITFLAG=args.split
FORGIFILE=args.forgi
OUTNAME=args.output
#FORGICG=args.forgi

if FORGIFILE=="":
    print "ERROR: Input file name needed!"
    parser.print_help()
    exit(1)

seq,definitions,ids,stems,coords,twists=[],[],[],[],[],[]
    
for line in open(FORGIFILE):
    splline=line.split()
    if splline==[] or line[0]=="=":
        break
    
    comm=splline[0].strip()
    if(comm=="seq"):
        seq=splline[1].strip()
    elif(comm=="define"):
        if(splline[1].strip()[0]=="s"):
            definitions.append((splline[1].strip(),[int(splline[2].strip()),int(splline[3].strip()),int(splline[4].strip()),int(splline[5].strip())]))
    elif(comm=="coord"):
        if(splline[1].strip()[0]=="s"):
            coords.append((splline[1].strip(),[float(splline[2].strip()),float(splline[3].strip()),float(splline[4].strip()),float(splline[5].strip()),float(splline[6].strip()),float(splline[7].strip()) ]))
    elif(comm=="twist"):
        if(splline[1].strip()[0]=="s"):
            twists.append((splline[1].strip(),[float(splline[2].strip()),float(splline[3].strip()),float(splline[4].strip()),float(splline[5].strip()),float(splline[6].strip()),float(splline[7].strip()) ]))
    elif(comm=="seq_ids"):
        ids=splline[1:]

stemDict=dict(definitions)
coordDict=dict(coords)
twistDict=dict(twists)

#make list
STEMLIST=[]
for st in definitions:
    STEMLIST.append(st[0])

#make sequences
SEQS=[]
for st in STEMLIST:
    stseq=seq[stemDict[st][0]-1:stemDict[st][1]]
    compst=seq[stemDict[st][2]-1:stemDict[st][3]]
    SEQS.append((st,stseq))
seqDict=dict(SEQS)
NSTEMS=len(SEQS)

FORGIPARAMS=[]
for st in STEMLIST:
    fgcoor1=np.array([coordDict[st][0],coordDict[st][1],coordDict[st][2]])
    fgcoor2=np.array([coordDict[st][3],coordDict[st][4],coordDict[st][5]])
    fgtwis1=np.array([twistDict[st][0],twistDict[st][1],twistDict[st][2]])
    fgtwis2=np.array([twistDict[st][3],twistDict[st][4],twistDict[st][5]])
    FGCM=0.5*(fgcoor1+fgcoor2)
    FGAXIS=fgcoor2-fgcoor1
    FGAXIS=FGAXIS/sqrt(np.dot(FGAXIS,FGAXIS))
    #forgipar=[FGCM,FGAXIS,fgtwis1,fgtwis2]
    forgipar=[fgcoor1,FGAXIS,fgtwis1,fgtwis2]
    FORGIPARAMS.append((st,forgipar))
fparamsDict=dict(FORGIPARAMS)

DIM=3
NAT=5
NBA=4
STTOT=4

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
    t0,   CM,b0=np.array([0,0,1]),np.array(CM),templtwistdict[typ]
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

def pdbprint(nt, resind, bas, chain):
    resname='{:3}'.format(bas)
    chindex='{:1}'.format(chain)
    #chindex="0"
    resindex='{:4}'.format(resind)
    record="ATOM  "
    xc,yc,zc=0,0,0
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
    chain2.append([refcoords[nt][0][3], info[1]])
    initnt=initnt+1
    #relabel duplex chain index
    for st in xrange(0,nst):
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

templtwistdict={
    "A" : np.array([0.99999, -0.00442147, 0]),
    "U" : np.array([0.99999, 0.00442147, 0]),
    "G" : np.array([0.999993, 0.00370642, 0]),
    "C" : np.array([0.999993, -0.00370642, 0])
}

ALLCFILES=[[FILEAA, FILEAU, FILEAG, FILEAC],[FILEUA, FILEUU, FILEUG, FILEUC],[FILEGA, FILEGU, FILEGG, FILEGC],[FILECA,FILECU, FILECG, FILECC]]

############################HERE STARTS THE CODE##################################
#get reference coordinates
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

pdbfile=open(OUTNAME,"w")
orig_stdout=sys.stdout
sys.stdout=pdbfile

#print header
#print "REMARK ERMSD PARAMS " + str(NSTEMS) + " "+str(RERMSD)
print "REMARK ERMSD PARAMS 1 " + str(RERMSD)

indexes=""
for st in STEMLIST:
    if SPLITFLAG:
        indexes=""
    for ii in xrange(stemDict[st][0],stemDict[st][1]+1):
        indexes=indexes+str(ii-1+FIRSTNT)+" "
    for ii in xrange(stemDict[st][2],stemDict[st][3]+1):
        indexes=indexes+str(ii-1+FIRSTNT)+" "
    if SPLITFLAG:
        print "REMARK ERMSD GROUP "+str(KERMSD)+" "+indexes
if not SPLITFLAG:
    print "REMARK ERMSD GROUP "+str(KERMSD)+" "+indexes

#print the chains
currnt=0
ALLCHAINS=[]
for st in STEMLIST:
    lchain,currnt,rchain=get_chain(seqDict[st],currnt,refcoord,0)
    ALLCHAINS.append((st,[lchain,rchain]))
chainDict=dict(ALLCHAINS)

for st in STEMLIST:
    lnnt=len(seqDict[st])
    forgiparams=get_forgi_matrix(chainDict[st][0],chainDict[st][1], fparamsDict[st], chainDict[st][0][0][1][1])
    ntcnt=0
    for nt in xrange(0,lnnt):
        prnt=forgi_arrange(forgiparams,chainDict[st][0][nt][0])
        pdbprint(prnt, ntcnt,chainDict[st][0][nt][1][1],chainDict[st][0][nt][1][2])
        ntcnt=ntcnt+1
    for nt in xrange(0,lnnt):
        prnt=forgi_arrange(forgiparams,chainDict[st][1][lnnt-1-nt][0])
        #print chainDict[st][1][lnnt-1-nt][1][0],chainDict[st][1][lnnt-1-nt][1][1],chainDict[st][1][lnnt-1-nt][1][2]
        pdbprint(prnt, ntcnt,chainDict[st][1][lnnt-1-nt][1][1],chainDict[st][1][lnnt-1-nt][1][2])
        ntcnt=ntcnt+1
sys.stdout=orig_stdout
pdbfile.close()
