from math import sin , cos , pi
import argparse


N_PARTS_PER_NT=5
parser=argparse.ArgumentParser()
parser.add_argument("pdbfile", help="The spqr-pdb file name")
parser.add_argument('-t', dest='pdbtemplate', help="The template name, for indexing atoms", default='nofile')
args=parser.parse_args()
FILENAME=args.pdbfile
TEMPLATE=args.pdbtemplate

NTOT=[22,20,23,20]
C2IND=[16,15,16,15]
C4IND=[18,12,19,12]
C6IND=[13,10,13,10]
PHIND=[1,1,1,1]

R2X=1.29
R4X=-1.31*cos(pi*63.82/180)
R4Y=-1.31*sin(pi*63.82/180)
R6X=-1.38*sin(pi*31.18/180)
R6Y= 1.38*cos(pi*31.18/180)

Y2X=1.4
Y4X=-1.42*sin(pi*31.2/180)
Y4Y=1.42*cos(pi*31.2/180)
Y6X=-1.3*cos(pi*61.31/180)
Y6Y=-1.3*sin(pi*61.31/180)


ALLFILE=[]
for line in open(FILENAME):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)

NATS=len(ALLFILE)
NNT=NATS/N_PARTS_PER_NT
curr=0
for i in xrange(0,NNT):
    cline=ALLFILE[i*N_PARTS_PER_NT]
    xline=ALLFILE[i*N_PARTS_PER_NT+1]
    yline=ALLFILE[i*N_PARTS_PER_NT+2]
    pline=ALLFILE[i*N_PARTS_PER_NT+4]
    base=cline[17:20].strip()
    resind=int(cline[22:26].strip())
    chain=cline[21].strip()
    #atname=cline[12:16].strip()
    cpos=[float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())]
    xpos=[float(xline[30:38].strip()), float(xline[38:46].strip()), float(xline[46:54].strip())]
    ypos=[float(yline[30:38].strip()), float(yline[38:46].strip()), float(yline[46:54].strip())]
    PHPOS=[float(pline[30:38].strip()), float(pline[38:46].strip()), float(pline[46:54].strip())]
    
    Xv=[xpos[0]-cpos[0], xpos[1]-cpos[1],xpos[2]-cpos[2]]
    Yv=[ypos[0]-cpos[0], ypos[1]-cpos[1],ypos[2]-cpos[2]]
    
    if(base=='A' or base=='G'):
        C2POS=[Xv[0]*R2X+cpos[0], Xv[1]*R2X+cpos[1], Xv[2]*R2X+cpos[2] ]
        C4POS=[Xv[0]*R4X+Yv[0]*R4Y+cpos[0], Xv[1]*R4X+Yv[1]*R4Y+cpos[1], Xv[2]*R4X+Yv[2]*R4Y+cpos[2]]
        C6POS=[Xv[0]*R6X+Yv[0]*R6Y+cpos[0], Xv[1]*R6X+Yv[1]*R6Y+cpos[1], Xv[2]*R6X+Yv[2]*R6Y+cpos[2]]
    elif(base=='U' or base=='C'):
        C2POS=[Xv[0]*Y2X+cpos[0], Xv[1]*Y2X+cpos[1], Xv[2]*Y2X+cpos[2] ]
        C4POS=[Xv[0]*Y4X+Yv[0]*Y4Y+cpos[0], Xv[1]*Y4X+Yv[1]*Y4Y+cpos[1], Xv[2]*Y4X+Yv[2]*Y4Y+cpos[2]]
        C6POS=[Xv[0]*Y6X+Yv[0]*Y6Y+cpos[0], Xv[1]*Y6X+Yv[1]*Y6Y+cpos[1], Xv[2]*Y6X+Yv[2]*Y6Y+cpos[2]]
        
    C2I,C4I,C6I=0,0,0
    if(TEMPLATE != 'nofile'):
        for teline in open(TEMPLATE):
            TEATT=teline[0:5].strip()
            if(TEATT=="ATOM"):
                TEIND=teline[6:11].strip()
                TENAM=teline[12:16].strip()
                TERES=int(teline[22:26].strip())
                
                if(TERES==i+1 and TENAM=="C2"):
                    C2I=int(TEIND)
                if(TERES==i+1 and TENAM=="C4"):
                    C4I=int(TEIND)
                if(TERES==i+1 and TENAM=="C6"):
                    C6I=int(TEIND)
    else:
        NTTYP="X"
        if(base=='A'):
            C2I,C4I,C6I,PHI=C2IND[0]+curr,C4IND[0]+curr,C6IND[0]+curr,PHIND[0]+curr
        if(base=='U'):
            C2I,C4I,C6I,PHI=C2IND[1]+curr,C4IND[1]+curr,C6IND[1]+curr,PHIND[0]+curr
        if(base=='G'):
            C2I,C4I,C6I,PHI=C2IND[2]+curr,C4IND[2]+curr,C6IND[2]+curr,PHIND[0]+curr
        if(base=='C'):
            C2I,C4I,C6I,PHI=C2IND[3]+curr,C4IND[3]+curr,C6IND[3]+curr,PHIND[0]+curr

        #this is the nucleotide is the first of its chain and therefore it lacks the phosphate group
        if(curr==0):
            C2I,C4I,C6I,curr,PHI=C2I-3,C4I-3,C6I-3,curr-3,-1
    #############
    
    NTTYP="X"
    if(base=='A'):
        NTTYP="A"
        curr=curr+NTOT[0]
    if(base=='U'):
        NTTYP="U"
        curr=curr+NTOT[1]
    if(base=='G'):
        NTTYP="G"
        curr=curr+NTOT[2]
    if(base=='C'):
        NTTYP="C"
        curr=curr+NTOT[3]
    if(base=='C' or base=='U'):
        if(PHI>=0):
            print "ATOM  "+str(PHI).rjust(5)+"  P   "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(PHPOS[0])+"{0:8.3f}".format(PHPOS[1])+"{0:8.3f}".format(PHPOS[2])+"  1.00                 P"
        print "ATOM  "+str(C6I).rjust(5)+"  C6  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C6POS[0])+"{0:8.3f}".format(C6POS[1])+"{0:8.3f}".format(C6POS[2])+"  1.00                 C"
        print "ATOM  "+str(C4I).rjust(5)+"  C4  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C4POS[0])+"{0:8.3f}".format(C4POS[1])+"{0:8.3f}".format(C4POS[2])+"  1.00                 C"
        print "ATOM  "+str(C2I).rjust(5)+"  C2  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C2POS[0])+"{0:8.3f}".format(C2POS[1])+"{0:8.3f}".format(C2POS[2])+"  1.00                 C"
    else:
        if(PHI>=0):
            print "ATOM  "+str(PHI).rjust(5)+"  P   "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(PHPOS[0])+"{0:8.3f}".format(PHPOS[1])+"{0:8.3f}".format(PHPOS[2])+"  1.00                 P"
        print "ATOM  "+str(C2I).rjust(5)+"  C2  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C2POS[0])+"{0:8.3f}".format(C2POS[1])+"{0:8.3f}".format(C2POS[2])+"  1.00                 C"
        print "ATOM  "+str(C4I).rjust(5)+"  C4  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C4POS[0])+"{0:8.3f}".format(C4POS[1])+"{0:8.3f}".format(C4POS[2])+"  1.00                 C"
        
