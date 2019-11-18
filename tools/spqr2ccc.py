from math import sin , cos , pi
import sys
#FILENAME=str(sys.argv[1])
N_PARTS_PER_NT=5

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
for line in open(sys.argv[1]):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        ALLFILE.append(line)

NATS=len(ALLFILE)
NNT=NATS/N_PARTS_PER_NT

for i in xrange(0,NNT):
    cline=ALLFILE[i*N_PARTS_PER_NT]
    xline=ALLFILE[i*N_PARTS_PER_NT+1]
    yline=ALLFILE[i*N_PARTS_PER_NT+2]
    base=cline[17:20].strip()
    resind=int(cline[22:26].strip())
    chain=cline[21].strip()
    #atname=cline[12:16].strip()
    cpos=[float(cline[30:38].strip()), float(cline[38:46].strip()), float(cline[46:54].strip())]
    xpos=[float(xline[30:38].strip()), float(xline[38:46].strip()), float(xline[46:54].strip())]
    ypos=[float(yline[30:38].strip()), float(yline[38:46].strip()), float(yline[46:54].strip())]
        
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
    for teline in open("TEMPLATE.pdb"):
        TEATT=teline[0:5].strip()
        if(TEATT=="ATOM"):
            TEIND=teline[6:11].strip()
            TENAM=teline[12:16].strip()
            TERES=int(teline[22:26].strip())
       
            if(TEATT=="ATOM" and TERES==i+1 and TENAM=="C2"):
                C2I=TEIND
            if(TEATT=="ATOM" and TERES==i+1 and TENAM=="C4"):
                C4I=TEIND
            if(TEATT=="ATOM" and TERES==i+1 and TENAM=="C6"):
                C6I=TEIND
                #print TENAM, TERES, TEIND
            
    NTTYP="X"
    if(base=='A'):
        NTTYP="A"
    if(base=='U'):
        NTTYP="U"
    if(base=='G'):
        NTTYP="G"
    if(base=='C'):
        NTTYP="C"

    print "ATOM  "+str(C2I).rjust(5)+"  C2  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C2POS[0])+"{0:8.3f}".format(C2POS[1])+"{0:8.3f}".format(C2POS[2])+"  1.00                 C"
    print "ATOM  "+str(C4I).rjust(5)+"  C4  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C4POS[0])+"{0:8.3f}".format(C4POS[1])+"{0:8.3f}".format(C4POS[2])+"  1.00                 C"
    print "ATOM  "+str(C6I).rjust(5)+"  C6  "+NTTYP.rjust(3)+" A"+str(i).rjust(4)+"    "+"{0:8.3f}".format(C6POS[0])+"{0:8.3f}".format(C6POS[1])+"{0:8.3f}".format(C6POS[2])+"  1.00                 C"
