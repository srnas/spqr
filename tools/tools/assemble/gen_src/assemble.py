from math import sqrt
from copy import deepcopy
import sys
DIM=3
NAT=5
NBA=4
STTOT=2
DUPLEX=0
#sequence is arg1 - we read all the stacks
if(sys.argv[1]=="-d"):
    DUPLEX=1
    STTOT=4
chain1=[]
chain2=[]
##########FUNCTIONS #################
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
    resindex='{:4}'.format(resind)
    record="ATOM  "
    
    atindex='{:5}'.format(resind*NAT)
    atname="BASE"
    xtemp, ytemp, ztemp=nt[0][0],nt[0][1],nt[0][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    atindex='{:5}'.format(resind*NAT+1)
    atname="XVEC"
    xtemp, ytemp, ztemp=nt[1][0],nt[1][1],nt[1][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    atindex='{:5}'.format(resind*NAT+2)
    atname="YVEC"
    xtemp, ytemp, ztemp=nt[2][0],nt[2][1],nt[2][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    atindex='{:5}'.format(resind*NAT+3)
    atname="SUGR"
    xtemp, ytemp, ztemp=nt[3][0],nt[3][1],nt[3][2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    atindex='{:5}'.format(resind*NAT+4)
    atname="PHOS"
    xtemp, ytemp, ztemp=nt[4][0],nt[4][1],nt[4][2]
    nidx, nidy, nidz=len(str(int(xtemp))),len(str(int(ytemp))),len(str(int(ztemp)))
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
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

def get_stacked_basis(current_base, ntdo, ntup,info):
    new_base=deepcopy(get_basis(refcoord[ntdo][ntup][0]))
    
    #in chain 0
    new_stack=deepcopy(refcoord[ntdo][ntup][1])
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
    chain1.append([nnew_stack, info])
    ret=get_basis(nnew_stack)
    
    if(DUPLEX):
        #we do the same but at the other chain, that is, index 2
        new_stack=deepcopy(refcoord[ntdo][ntup][2])
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
        chain2.append([nnew_stack, info])
    
    return ret

############################HERE STARTS THE CODE##################################            
#READ THE TEMPLATES
refcoord=[]
for nt1 in xrange(0,NBA):
    wr2=[]
    for nt2 in xrange(0,NBA):
        ant1, ant2=get_bname(nt1), get_bname(nt2)
        nam="gen_src/"+ant1+ant2+"_spqr.pdb"
        ALLFILE=[]
        for line in open(nam):
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
                ntcoords.append([x,y,z])
            stack.append(ntcoords)
        wr2.append(stack)
    refcoord.append(wr2)

#READ THE SEQUENCE AND PRINT THE COORDINATES
if(DUPLEX):
    seq=sys.argv[2]
else:
    seq=sys.argv[1]
nnt=len(seq)
nst=nnt-1
#print first nt
nt=base_index(seq[0])
current_base=get_basis(refcoord[nt][0][0])
info=[0, seq[0].upper(),0]
chain1.append([refcoord[nt][0][0], info])
if(DUPLEX):
    info=[nnt, get_complb(seq[0].upper()),1]
    chain2.append([refcoord[nt][0][3], info])

for st in xrange(0,nst):
    stack=seq[st:st+2]
    info=[st+1, stack[1].upper(), 0]
    current_base=get_stacked_basis(current_base, base_index(stack[0]), base_index(stack[1]), info)

for st in xrange(0,nnt):
    pdbprint(chain1[st][0], chain1[st][1][0],chain1[st][1][1],chain1[st][1][2])
if(DUPLEX):
    #for st in xrange(0,nnt):
    #   info2=[st+nnt, get_complb(chain1[st][1][1]), 1]
    #  pdbprint(chain2[st][0], info2[0], info2[1], info2[2])
    for st in xrange(0,nnt):
        info2=[nnt+st, get_complb(chain1[nnt-1-st][1][1]), 1]
        #resindex must be counted in the opposite direction, following the 5'-3' convention
        pdbprint(chain2[nnt-1-st][0], info2[0], info2[1], info2[2])
