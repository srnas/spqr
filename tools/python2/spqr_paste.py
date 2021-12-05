from math import sqrt
import argparse
import sys
import numpy as np
N_PARTS_PER_NT=5
ALLFILE1=[]
ALLFILE2=[]

def pdbprint(nt, resind, bas, chain ):
    resname='{:3}'.format(bas)
    chindex='{:1}'.format(chain)
    resindex='{:4}'.format(resind)
    record="ATOM  "
    atindex='{:5}'.format(resind*N_PARTS_PER_NT)
    atname="BASE"
    xtemp, ytemp, ztemp=nt[0][0],nt[0][1],nt[0][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*N_PARTS_PER_NT+1)
    atname="XVEC"
    xtemp, ytemp, ztemp=nt[1][0],nt[1][1],nt[1][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*N_PARTS_PER_NT+2)
    atname="YVEC"
    xtemp, ytemp, ztemp=nt[2][0],nt[2][1],nt[2][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*N_PARTS_PER_NT+3)
    atname="SUGR"
    xtemp, ytemp, ztemp=nt[3][0],nt[3][1],nt[3][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    atindex='{:5}'.format(resind*N_PARTS_PER_NT+4)
    atname="PHOS"
    xtemp, ytemp, ztemp=nt[4][0],nt[4][1],nt[4][2]
    nidx, nidy, nidz=len(str(int(xtemp))),len(str(int(ytemp))),len(str(int(ztemp)))
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    return

parser=argparse.ArgumentParser()
parser.add_argument("-f","--first", help="First input file",type=str,default="")
parser.add_argument("-s","--second", help="Second inpu tfile",type=str,default="")
args=parser.parse_args()

INPUTFILE1=args.first
INPUTFILE2=args.second
if INPUTFILE1=="" or INPUTFILE2=="":
    print "ERROR : input file needed."
    parser.print_help()
    exit(1)

for line in open(INPUTFILE1):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE1.append(line)
for line in open(INPUTFILE2):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE2.append(line)

NATS1=len(ALLFILE1)
NNT1=NATS1/N_PARTS_PER_NT

NATS2=len(ALLFILE2)
NNT2=NATS2/N_PARTS_PER_NT

SELNT=NATS1-N_PARTS_PER_NT
cline=ALLFILE1[SELNT]
xline=ALLFILE1[SELNT+1]
yline=ALLFILE1[SELNT+2]

CM=np.array([float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())])
rXV=np.array([float(xline[30:38].strip())-CM[0], float(xline[38:46].strip())-CM[1], float(xline[46:54].strip())-CM[2]])
rYV=np.array([float(yline[30:38].strip())-CM[0], float(yline[38:46].strip())-CM[1], float(yline[46:54].strip())-CM[2]])
mXV=sqrt(np.dot(rXV,rXV))
mYV=sqrt(np.dot(rYV,rYV))
XV=rXV/float(mXV)
YV=rYV/float(mYV)
rZV=np.cross(XV,YV)
mZV=sqrt(np.dot(rZV,rZV))
ZV=rZV/float(mZV)


SELNT2=0
cline=ALLFILE2[SELNT2]
xline=ALLFILE2[SELNT2+1]
yline=ALLFILE2[SELNT2+2]
oCM=np.array([float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())])
rXV=np.array([float(xline[30:38].strip())-oCM[0], float(xline[38:46].strip())-oCM[1], float(xline[46:54].strip())-oCM[2]])
rYV=np.array([float(yline[30:38].strip())-oCM[0], float(yline[38:46].strip())-oCM[1], float(yline[46:54].strip())-oCM[2]])
mXV=sqrt(np.dot(rXV,rXV))
mYV=sqrt(np.dot(rYV,rYV))
oXV=rXV/float(mXV)
oYV=rYV/float(mYV)
rZV=np.cross(oXV,oYV)
mZV=sqrt(np.dot(rZV,rZV))
oZV=rZV/float(mZV)

#rotation matrix
#ROTMAT=np.array([[np.dot(XV,oXV),np.dot(XV,oYV),np.dot(XV,oZV)],[np.dot(YV,oXV),np.dot(YV,oYV),np.dot(YV,oZV)],[np.dot(ZV,oXV),np.dot(ZV,oYV),np.dot(ZV,oZV)]])
#ROTMAT=np.array([[1,0,0],[0,1,0],[0,0,1]])
#ROTMAT=np.transpose(oROTMAT)



for i in xrange(0,NNT2):
    cline=ALLFILE2[i*N_PARTS_PER_NT]
    xline=ALLFILE2[i*N_PARTS_PER_NT+1]
    yline=ALLFILE2[i*N_PARTS_PER_NT+2]
    sline=ALLFILE2[i*N_PARTS_PER_NT+3]
    pline=ALLFILE2[i*N_PARTS_PER_NT+4]
    basetyp=cline[17:20].strip()
    resind=int(cline[22:26].strip())
    chain=cline[21].strip()
    puck="3"
    glyc="A"
    temppuck="X"
    tempglyc="X"
    if(len(cline)>=57):
        temppuck=cline[55]
        tempglyc=cline[56]
    if((temppuck=="3" or temppuck=="2") and (tempglyc=="A" or tempglyc=="H" or tempglyc=="S")):
        puck=temppuck
        glyc=tempglyc
    #for resline in open(coordfile):
    ccoo=np.array([float(cline[30:38].strip())-oCM[0], float(cline[38:46].strip())-oCM[1], float(cline[46:54].strip())-oCM[2]])
    xcoo=np.array([float(xline[30:38].strip())-oCM[0], float(xline[38:46].strip())-oCM[1], float(xline[46:54].strip())-oCM[2]])
    ycoo=np.array([float(yline[30:38].strip())-oCM[0], float(yline[38:46].strip())-oCM[1], float(yline[46:54].strip())-oCM[2]])
    scoo=np.array([float(sline[30:38].strip())-oCM[0], float(sline[38:46].strip())-oCM[1], float(sline[46:54].strip())-oCM[2]])
    pcoo=np.array([float(pline[30:38].strip())-oCM[0], float(pline[38:46].strip())-oCM[1], float(pline[46:54].strip())-oCM[2]])

    rp=np.array([np.dot(ccoo,oXV),np.dot(ccoo,oYV),np.dot(ccoo,oZV)])
    nccoo=(rp[0]*XV+rp[1]*YV+rp[2]*ZV)+CM

    rp=np.array([np.dot(xcoo,oXV),np.dot(xcoo,oYV),np.dot(xcoo,oZV)])
    nxcoo=(rp[0]*XV+rp[1]*YV+rp[2]*ZV)+CM
    
    rp=np.array([np.dot(ycoo,oXV),np.dot(ycoo,oYV),np.dot(ycoo,oZV)])
    nycoo=(rp[0]*XV+rp[1]*YV+rp[2]*ZV)+CM

    rp=np.array([np.dot(scoo,oXV),np.dot(scoo,oYV),np.dot(scoo,oZV)])
    nscoo=(rp[0]*XV+rp[1]*YV+rp[2]*ZV)+CM

    rp=np.array([np.dot(pcoo,oXV),np.dot(pcoo,oYV),np.dot(pcoo,oZV)])
    npcoo=(rp[0]*XV+rp[1]*YV+rp[2]*ZV)+CM
    newnt=[nccoo,nxcoo,nycoo,nscoo,npcoo]
    pdbprint(newnt,i+NNT1,basetyp,0)
