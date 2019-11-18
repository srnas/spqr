import sys
from math import sqrt
N_PARTS_PER_NT=5
ALLFILE=[]
for line in open(sys.argv[1]):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)

NATS=len(ALLFILE)
NNT=NATS/N_PARTS_PER_NT

#print "MODEL"
COUNTAT=0
for i in xrange(0,NNT):
    cline=ALLFILE[i*N_PARTS_PER_NT]
    xline=ALLFILE[i*N_PARTS_PER_NT+1]
    yline=ALLFILE[i*N_PARTS_PER_NT+2]
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

    CM=[float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())]
    rXV=[float(xline[30:38].strip())-CM[0], float(xline[38:46].strip())-CM[1], float(xline[46:54].strip())-CM[2]]
    rYV=[float(yline[30:38].strip())-CM[0], float(yline[38:46].strip())-CM[1], float(yline[46:54].strip())-CM[2]]
    mXV=sqrt(rXV[0]*rXV[0]+rXV[1]*rXV[1]+rXV[2]*rXV[2])
    mYV=sqrt(rYV[0]*rYV[0]+rYV[1]*rYV[1]+rYV[2]*rYV[2])
    XV=[rXV[0]/mXV,rXV[1]/mXV,rXV[2]/mXV]
    YV=[rYV[0]/mYV,rYV[1]/mYV,rYV[2]/mYV]
    rZV=[XV[1]*YV[2]-XV[2]*YV[1], XV[2]*YV[0]-XV[0]*YV[2], XV[0]*YV[1]-XV[1]*YV[0]]
    mZV=sqrt(rZV[0]*rZV[0]+rZV[1]*rZV[1]+rZV[2]*rZV[2])
    ZV=[rZV[0]/mZV,rZV[1]/mZV,rZV[2]/mZV]
    
    coordfile="bbm_templates/"+basetyp+"_"+glyc+puck+".rco"
    for resline in open(coordfile):
        reslist=resline.split()
        attype=reslist[0].strip()
        atx=float(reslist[1].strip())
        aty=float(reslist[2].strip())
        atz=float(reslist[3].strip())
        newx=atx*XV[0]+aty*YV[0]+atz*ZV[0]+CM[0]
        newy=atx*XV[1]+aty*YV[1]+atz*ZV[1]+CM[1]
        newz=atx*XV[2]+aty*YV[2]+atz*ZV[2]+CM[2]
        strX=("%0.3f" % newx).rjust(7)
        strY=("%0.3f" % newy).rjust(7)
        strZ=("%0.3f" % newz).rjust(7)
        if(len(strX)>6):
            strX=("%0.2f" % newx).rjust(7)
        if(len(strY)>6):
            strY=("%0.2f" % newy).rjust(7)
        if(len(strZ)>6):
            strZ=("%0.2f" % newz).rjust(7)
        print "ATOM ", str(COUNTAT ).rjust(5), attype.rjust(4),  str(basetyp).rjust(3), chain, str(resind).rjust(3),"   ",  strX, strY, strZ," 1.00  0.00"
        COUNTAT=COUNTAT+1
        
#print "ENDMDL"
