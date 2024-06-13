import numpy as np
import argparse
import random as rd
from collections import Counter

NATNT=5
IPHOS=4
MAXBP=1.7
#["N9",	"C1'",   "N*", "CT"     ,1.3830,   354803.2,0,0],
#["C1'",	"C2'",	 "CT", "CT"	,1.5260,   259408.0,0,0],
#       ["C2'",	"C3'",	 "CT", "CT"	,1.5260,   259408.0,0,0],
#       ["C2'",	"O2'",	 "CT", "OH"	,1.4100,   267776.0,0,0],
#       ["C3'",	"O3'",	 "CT", "OS"	,1.4100,   267776.0,0,0],
#       ["C4'",	"C3'",	 "CT", "CT"	,1.5260,   259408.0,0,0],
#       ["C4'",	"O4'",	 "CT", "OS"	,1.4100,   267776.0,0,0],
#       ["C1'",	"O4'",	 "CT", "OS"	,1.4100,   267776.0,0,0],
#       ["P",	"O3'",	 "P",  "OS"	,1.6100,   192464.0,1,0],
#       ["O5'",	"C5'",	 "OH", "CI"	,1.4100,   267776.0,1,1],
#       ["P",	"O5'",	 "P",  "OS"	,1.6100,   192464.0,1,1],
#       ["OP1",	"P",	 "O2", "P"	,1.4800,   439320.0,1,1],
#       ["OP2",	"P",	 "O2", "P"	,1.4800,   439320.0,1,1],
#       ["C4'",	"C5'",	 "CT", "CI"	,1.5260,   259408.0,1,1]]

BONDS=[["C1'",	"C2'",	 "CT", "CT"	,1.5260,   400000.0,0,0],
       ["C2'",	"C3'",	 "CT", "CT"	,1.5260,   400000.0,0,0],
       ["C2'",	"O2'",	 "CT", "OH"	,1.4100,   400000.0,0,0],
       ["C3'",	"O3'",	 "CT", "OS"	,1.4100,   400000.0,0,0],
       ["C4'",	"C3'",	 "CT", "CT"	,1.5260,   400000.0,0,0],
       ["C4'",	"O4'",	 "CT", "OS"	,1.4100,   400000.0,0,0],
       ["C1'",	"O4'",	 "CT", "OS"	,1.4100,   400000.0,0,0],
       ["C4'",	"C5'",	 "CT", "CI"	,1.5260,   400000.0,0,0],
       ["P",	"O3'",	 "P",  "OS"	,1.6100,   400000.0,1,0],
       ["O5'",	"C5'",	 "OH", "CI"	,1.4100,   400000.0,1,1],
       ["P",	"O5'",	 "P",  "OS"	,1.6100,   400000.0,1,1],
       ["OP1",	"P",	 "O2", "P"	,1.4800,   400000.0,1,1],
       ["OP2",	"P",	 "O2", "P"	,1.4800,   400000.0,1,1],
       ["C4'",	"C5'",	 "CT", "CI"	,1.5260,   400000.0,1,1]]

ANGLES=[["C1'",	"C2'",	"O2'",	"CT","CT","OH",	109.500*np.pi/180,    10000,0,0,0],
        ["C1'",	"C2'",	"C3'",	"CT","CT","CT",	109.500*np.pi/180,    10000,0,0,0],
        ["C2'",	"C3'",	"O3'",	"CT","CT","OS",	109.500*np.pi/180,    10000,0,0,0],
        ["C2'",	"C3'"	,"C4'",	"CT","CT","CT",	109.500*np.pi/180,    10000,0,0,0],
        ["C3'",	"O3'",	"P",	"CT","OS","P",	120.500*np.pi/180,    10000,0,0,0],
        ["C3'",	"C4'",	"C5'",	"CT","CT","OS",	109.500*np.pi/180,    10000,0,0,0],
        ["C3'",	"C4'",	"O4'",	"CT","CT","OS",	109.500*np.pi/180,    10000,0,0,0],
	["C4'",	"O4'",	"C1'",	"CT","OS","CT",	109.500*np.pi/180,    10000,0,0,0],
        ["O4'",	"C1'",	"C2'",	"OS","CT","CT",	109.500*np.pi/180,    10000,0,0,0],
	["O2'",	"C2'",	"C3'",	"OH","CT","CT",	109.500*np.pi/180,    10000,0,0,0],
        ["O3'",	"C3'",	"C4'",	"OS","CT","CT",	109.500*np.pi/180,    10000,0,0,0],
        
	["O3'",	"P",	"OP1",	"OS","P","O2",	108.230*np.pi/180,    10000,0,1,1],
        ["O3'",	"P",	"OP2",	"OS","P","O2",	108.230*np.pi/180,    10000,0,1,1],
        ["O3'",	"P",	"O5'",	"OS","P","OS",	102.600*np.pi/180,    10000,0,1,1],
        ["P",	"O5'",	"C5'",	"P","OS","CI",	120.500*np.pi/180,    10000,1,1,1],
        ["O5'",	"C5'",	"C4'",	"OS","CI","CT",	109.500*np.pi/180,    10000,1,1,1],
        ["OP1",	"P",	"O5'",	"OS","P","OS",	102.600*np.pi/180,    10000,1,1,1],
        ["OP2",	"P",	"O5'",	"OS","P","OS",	102.600*np.pi/180,    10000,1,1,1],
        ["OP1",	"P",	"OP2",	"OS","P","OS",	102.600*np.pi/180,    10000,1,1,1]]


#first 5 are from pucker c3' endo, last 5 are from pucker c2' endo
DIHEDRALS=[["C1'",	"C2'",	"C3'","C4'",	"CT","CT","CT",	"CT", 39.21*np.pi/180,    10000,0,0,0,0,3],
           ["C2'",	"C3'",	"C4'","O4'",	"CT","CT","CT",	"OS",-39.23*np.pi/180,    10000,0,0,0,0,3],
           ["C3'",	"C4'",	"O4'","C1'",	"CT","CT","OS",	"CT", 23.65*np.pi/180,    10000,0,0,0,0,3],
           ["C4'",	"O4'",	"C1'","C2'",	"CT","OS","CT",	"CT",  1.86*np.pi/180,    10000,0,0,0,0,3],
           ["O4'",	"C1'",	"C2'","C3'",	"OS","CT","CT",	"CT",-26.32*np.pi/180,    10000,0,0,0,0,3],
           ["C1'",	"C2'",	"C3'","C4'",	"CT","CT","CT",	"CT",-37.03*np.pi/180,    10000,0,0,0,0,2],
           ["C2'",	"C3'",	"C4'","O4'",	"CT","CT","CT",	"OS", 24.54*np.pi/180,    10000,0,0,0,0,2],
           ["C3'",	"C4'",	"O4'","C1'",	"CT","CT","OS",	"CT", -0.81*np.pi/180,    10000,0,0,0,0,2],
           ["C4'",	"O4'",	"C1'","C2'",	"CT","OS","CT",	"CT",-23.38*np.pi/180,    10000,0,0,0,0,2],
           ["O4'",	"C1'",	"C2'","C3'",	"OS","CT","CT",	"CT", 37.61*np.pi/180,    10000,0,0,0,0,2],
           ["C1'",	"C2'",	"C3'","C4'",	"CT","CT","CT",	"CT", 39.21*np.pi/180,    10000,1,1,1,1,3],
           ["C2'",	"C3'",	"C4'","O4'",	"CT","CT","CT",	"OS",-39.23*np.pi/180,    10000,1,1,1,1,3],
           ["C3'",	"C4'",	"O4'","C1'",	"CT","CT","OS",	"CT", 23.65*np.pi/180,    10000,1,1,1,1,3],
           ["C4'",	"O4'",	"C1'","C2'",	"CT","OS","CT",	"CT",  1.86*np.pi/180,    10000,1,1,1,1,3],
           ["O4'",	"C1'",	"C2'","C3'",	"OS","CT","CT",	"CT",-26.32*np.pi/180,    10000,1,1,1,1,3],
           ["C1'",	"C2'",	"C3'","C4'",	"CT","CT","CT",	"CT",-37.03*np.pi/180,    10000,1,1,1,1,2],
           ["C2'",	"C3'",	"C4'","O4'",	"CT","CT","CT",	"OS", 24.54*np.pi/180,    10000,1,1,1,1,2],
           ["C3'",	"C4'",	"O4'","C1'",	"CT","CT","OS",	"CT", -0.81*np.pi/180,    10000,1,1,1,1,2],
           ["C4'",	"O4'",	"C1'","C2'",	"CT","OS","CT",	"CT",-23.38*np.pi/180,    10000,1,1,1,1,2],
           ["O4'",	"C1'",	"C2'","C3'",	"OS","CT","CT",	"CT", 37.61*np.pi/180,    10000,1,1,1,1,2]
]


def get_dist(vec1,vec2):
    d2=np.dot(vec1-vec2,vec1-vec2)
    return np.sqrt(d2)

#OPEN BACKMAPPED
def read_pdb(filename):
    alldata=[]
    coords,pdbdata,mobile,glpdata=[],[],[],[]
    if filename=="":
        print("File not found!")
        exit(1)
    for line in open(filename):
        nam=line[:4]
        if(nam.strip()=="ATOM"):
            alldata.append(line)
    currnt=int(alldata[0][22:26].strip())
    ntcoord,ntdata=[],[]
    for atom in alldata:
        thisnt=int(atom[22:26].strip())
        name=atom[12:16].strip()
        resname=atom[17:20].strip()
        puck=3
        if(len(atom)>=60):
            puck=int(atom[57].strip())
        if thisnt!=currnt:
            currnt=thisnt
            coords.append(ntcoord)
            pdbdata.append(ntdata)
            ntcoord,ntdata=[],[]
        ntcoord.append(np.array([float(atom[30:38].strip()), float(atom[38:46].strip()), float(atom[46:54].strip())]) )
        ntdata.append([thisnt,name,resname])
        if atom[12:16].strip()=="BASE":
            glpdata.append(puck)
    coords.append(ntcoord)
    pdbdata.append(ntdata)
    return coords,pdbdata,len(coords),glpdata

#CREATE BONDLISTS AND ANGLELISTS
def create_interaction_lists(coords,atdata,glpdata):
    DMAX=4
    firstnt=atdata[0][0]
    bondlists, anglelists, dihedrallists, evlists=[],[],[],[]
    for ii in range(len(coords)):
        ibonds,abonds,dbonds,evneighs=[],[],[],[]
        for jj in range(len(coords)):
            for bb in range(len(BONDS)):
                if ( ((BONDS[bb][0]==atdata[ii][1] and BONDS[bb][6]==atdata[ii][0]-firstnt) and (BONDS[bb][1]==atdata[jj][1] and BONDS[bb][7]==atdata[jj][0]-firstnt))
                     or
                     ((BONDS[bb][1]==atdata[ii][1] and BONDS[bb][7]==atdata[ii][0]-firstnt) and (BONDS[bb][0]==atdata[jj][1] and BONDS[bb][6]==atdata[jj][0]-firstnt)) ):
                    ibonds.append([jj,bb])
        
        for aa in range(len(ANGLES)):
            tri=[[ANGLES[aa][0],ANGLES[aa][8]],[ANGLES[aa][1],ANGLES[aa][9]],[ANGLES[aa][2],ANGLES[aa][10]]]
            #if atdata[ii][1]=="O3'" and atdata[ii][0]-firstnt==0:
            #    print([atdata[ii][1],atdata[ii][0]-firstnt])
            #    print(tri)
            if [atdata[ii][1],atdata[ii][0]-firstnt] in tri:
                #print("yess")
                indta=[]
                for tt in tri:
                    for kk in range(len(coords)):
                        if atdata[kk][1]==tt[0] and atdata[kk][0]-firstnt==tt[1]:
                            indta.append(kk)
                            break
                if len(indta)==3:
                    abonds.append([indta,aa])
        for dd in range(len(DIHEDRALS)):
            tet=[[DIHEDRALS[dd][0],DIHEDRALS[dd][10],DIHEDRALS[dd][14]],
                 [DIHEDRALS[dd][1],DIHEDRALS[dd][11],DIHEDRALS[dd][14]],
                 [DIHEDRALS[dd][2],DIHEDRALS[dd][12],DIHEDRALS[dd][14]],
                 [DIHEDRALS[dd][3],DIHEDRALS[dd][13],DIHEDRALS[dd][14]]]
            if [atdata[ii][1],atdata[ii][0]-firstnt,glpdata[atdata[ii][0]-firstnt]] in tet:
                indtd=[]
                for tt in tet:
                    for kk in range(len(coords)):
                        if atdata[kk][1]==tt[0] and atdata[kk][0]-firstnt==tt[1] and glpdata[atdata[kk][0]-firstnt]==tt[2]:
                            indtd.append(kk)
                            break
                if len(indtd)==4:
                    dbonds.append([indtd,dd])
            
        for jj in range(len(coords)):
            bflag=True
            if ii==jj: bflag=False
            for bb in ibonds:
                if bb[0]==jj:
                    bflag=False
                    break
            for aa in abonds:
                if jj in aa[0]:
                    bflag=False
                    break
            if bflag:
                dis=get_dist(coords[ii],coords[jj])
                if dis<DMAX: evneighs.append([jj,0])
        #print("before", firstnt)
        #if atdata[ii][1]=="O3'" and atdata[ii][0]-firstnt==0:
        #if atdata[ii][1]=="P" and atdata[ii][0]-firstnt==1:
        #    print("here")
        #    for jj in evneighs:
        #        print(atdata[jj[0]])
        #    print("bonds")
        #    for jj in ibonds:
        #        print(atdata[jj[0]])
        #    print("angles")
        #    for jj in abonds:
        #        print(atdata[jj[0][0]],atdata[jj[0][1]],atdata[jj[0][2]])
        #    exit(1)
        bondlists.append(ibonds)
        anglelists.append(abonds)
        dihedrallists.append(dbonds)
        evlists.append(evneighs)
        
    return bondlists, anglelists, dihedrallists, evlists

def get_local_energy(nt, coords, atdata, bondlists, anglelists, dihedrallists, nblists, refcoords, refdata):
    energ=0
    #bonded loop
    for bb in bondlists[nt]:
        btype,ntneigh=bb[1],bb[0]
        benerg=0.5*BONDS[btype][5]*(get_dist(coords[nt],coords[ntneigh])-BONDS[btype][4])**2
        #print(atdata[nt], atdata[ntneigh], get_dist(coords[nt],coords[ntneigh]),benerg)
        energ+=benerg
    for bb in anglelists[nt]:
        btype, nt1,nt2,nt3 = bb[1],bb[0][0],bb[0][1],bb[0][2]
        vec1,vec2=coords[nt1]-coords[nt2], coords[nt3]-coords[nt2]
        ang=np.arccos(np.dot(vec1,vec2)/np.sqrt(np.dot(vec1,vec1)*np.dot(vec2,vec2)))
        aenerg=0.5*ANGLES[btype][7]*(ang-ANGLES[btype][6])**2
        energ+=aenerg
    for bb in dihedrallists[nt]:
        btype,nt1,nt2,nt3,nt4=bb[1],bb[0][0],bb[0][1],bb[0][2],bb[0][3]
        vec1,vec2,vec3=coords[nt2]-coords[nt1],coords[nt3]-coords[nt2],coords[nt4]-coords[nt3]
        vcr1,vcr2=np.cross(vec1,vec2),np.cross(vec2,vec3)
        #norm1,norm2,vn2=np.sqrt(np.dot(vcr1,vcr1)),np.sqrt(np.dot(vcr2,vcr2)),np.sqrt(np.dot(vec2,vec2))
        vn2=np.sqrt(np.dot(vec2,vec2))
        cphi=np.dot(vcr1,vcr2)
        sphi=np.dot(vec2,np.cross(vcr1,vcr2))/(vn2)
        dih=np.arctan2(sphi,cphi)
        denerg=0.5*DIHEDRALS[btype][9]*(dih-DIHEDRALS[btype][8])**2
        #print(nt, bb,dih*180.0/np.pi)
        energ+=denerg
        
    #non-bonded loop
    for neigh in nblists[nt]:
        etype,ntneigh=neigh[1],neigh[0]
        evenerg=get_EV_energy(get_dist(coords[nt],coords[ntneigh]),etype)
        energ+=evenerg
    #push phosphate to reference
    firstnt=atdata[0][0]
    if atdata[nt][1]=="P" and atdata[nt][0]==firstnt+1:
        lnt=0
        REFPHOS=refcoords[1*NATNT+IPHOS]
        dp=get_dist(REFPHOS, coords[nt])
        energ+=0.5*1000000*dp**2
    return energ

def get_EV_energy(dis,etype):
    DCUT=2.5
    EMAX=1000
    #EMAX=0
    evenerg=0
    if dis<DCUT:
        evenerg=EMAX*(1-dis/DCUT)
    return evenerg

def write_coords(coords,atdata,datafile):
    datafile.write("MODEL 1\n")
    for ii in range(len(coords)):
        record="ATOM  "
        atindex='{:5}'.format(ii)
        atname='{:4}'.format(atdata[ii][1])
        resname='{:3}'.format(atdata[ii][2])
        chindex="0"
        resindex='{:4}'.format(atdata[ii][0])
        xpos,ypos,zpos='{:8.3f}'.format(coords[ii][0]),'{:8.3f}'.format(coords[ii][1]),'{:8.3f}'.format(coords[ii][2])
        line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos+"  1.00  0.00\n"
        datafile.write(line)
    datafile.write("ENDMDL\n")

def write_full_coords(coords,atdata,datafile):
    datafile.write("MODEL 1\n")
    cnt=0
    for ii in range(len(coords)):
        for jj in range(len(coords[ii])):
            record="ATOM  "
            atindex='{:5}'.format(cnt)
            atname='{:4}'.format(atdata[ii][jj][1])
            resname='{:3}'.format(atdata[ii][jj][2])
            chindex="0"
            resindex='{:4}'.format(atdata[ii][jj][0])
            xpos,ypos,zpos='{:8.3f}'.format(coords[ii][jj][0]),'{:8.3f}'.format(coords[ii][jj][1]),'{:8.3f}'.format(coords[ii][jj][2])
            line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos+"  1.00  0.00\n"
            cnt=cnt+1
            datafile.write(line)
    datafile.write("ENDMDL\n")
    
def mc_move_atom(pos):
    amp=0.1
    delta=np.array([rd.random()-0.5,rd.random()-0.5, rd.random()-0.5])
    return pos+delta

def mc_sweep(kt, coords, atdata, bondlists, anglelists, dihedrallists, evlists, mobile, refcoords, refdata):
    for ii in mobile:
        ene_init=get_local_energy(ii,coords,atdata,bondlists,anglelists, dihedrallists,evlists, refcoords, refdata)
        tempcoord=coords[ii]
        coords[ii]=mc_move_atom(coords[ii])
        ene_trial=get_local_energy(ii,coords,atdata,bondlists,anglelists, dihedrallists,evlists, refcoords, refdata)
        accflag=False
        if kt==0:
            if ene_trial<ene_init: accflag=True
        else:
            if rd.random() < np.exp(-(ene_trial-ene_init)/kt): accflag=True
        if accflag==False:
            coords[ii]=tempcoord

            
def need_phos_relax(coords,atdata,refcoords,refdata):
    phos,c3p,c5p=[],[],[]
    firstnt=atdata[0][0]
    #fetch p, o3' and o5'
    for ii in range(len(atdata)):
        if atdata[ii][1]=="P" and atdata[ii][0]==firstnt+1: phos=coords[ii]
        if atdata[ii][1]=="O5'" and atdata[ii][0]==firstnt+1: c5p=coords[ii]
        if atdata[ii][1]=="O3'" and atdata[ii][0]==firstnt: c3p=coords[ii]
    refphos=refcoords[IPHOS]
    res=True
    if get_dist(phos,c3p)<MAXBP and get_dist(phos,c5p)<MAXBP:
        res=False
    return res

    
def relax_dinucleotide(coords,atdata,mobile,refcoords,refdata,glpdata):
    KT=0.0
    BONDLISTS, ANGLELISTS, DIHEDRALLISTS, EVLISTS=create_interaction_lists(coords,atdata,glpdata)
    for ii in range(10):
        if need_phos_relax(coords,atdata,refcoords,refdata):
            for jj in range(100):
                mc_sweep(KT,coords,atdata,BONDLISTS,ANGLELISTS,DIHEDRALLISTS,EVLISTS,mobile,refcoords,refdata)
        else: break
    return coords

def get_dinucleotide(dnt,fullcoords,fulldata,refcoords,refdata,glpdata):
    ntcoords,ntdata,ntmobile,ntrefcoords,ntrefdata,ntglpdata=[],[],[],[],[],[]
    for ii in fullcoords[dnt]: ntcoords.append(ii)
    for ii in fullcoords[dnt+1]: ntcoords.append(ii)
    for ii in fulldata[dnt]: ntdata.append(ii)
    for ii in fulldata[dnt+1]: ntdata.append(ii)
    for ii in refcoords[dnt]: ntrefcoords.append(ii)
    for ii in refcoords[dnt+1]: ntrefcoords.append(ii)
    for ii in refdata[dnt]: ntrefdata.append(ii)
    for ii in refdata[dnt+1]: ntrefdata.append(ii)
    ntglpdata=[glpdata[dnt],glpdata[dnt+1]]
    cnt=0
    firstnt=ntdata[0][0]
    for ii in ntdata:
        name=ii[1]
        ntid=ii[0]
        #if name in ["OP2","P","OP1","O5'","C5'","C4'","O4'","C1'","C2'","O3'","O2'","C3'"] : ntmobile.append(cnt)
        if name in ["OP2","P","OP1","O5'","C5'"] and ntid==firstnt+1 : ntmobile.append(cnt)
        if name in ["C3'","O3'","C4'","O4'","C2'","O2'"] and ntid==firstnt : ntmobile.append(cnt)
        cnt=cnt+1
    
    return ntcoords,ntdata,ntmobile,ntrefcoords,ntrefdata,ntglpdata
    
def update_coords(fullcoords,dnt,coords,atdata):
    nt1,nt2=[],[]
    firstnt=atdata[0][0]
    for ii in range(len(atdata)):
        if atdata[ii][0]==firstnt:
            nt1.append(coords[ii])
        else :
            nt2.append(coords[ii])
    fullcoords[dnt]=nt1
    fullcoords[dnt+1]=nt2
    
    
    return fullcoords
    
###MAIN###
parser=argparse.ArgumentParser()
parser.add_argument("-i","--input", help="Input file",type=str,default="")
parser.add_argument("-r","--reference", help="Reference spqr pdb file",type=str,default="")
parser.add_argument("-o", "--output", help="Output file", type=str, default="mc_mini.pdb")
args=parser.parse_args()
if args.input=="" or args.reference=="":
    parser.print_help()
    exit(1)
FULLCOORDS,FULLDATA,NNT, TEMPDATA=read_pdb(args.input)
REFCOORDS,REFDATA,REFNNT,GLPDATA=read_pdb(args.reference)

OUTFILE=open(args.output, "w")
write_full_coords(FULLCOORDS,FULLDATA,OUTFILE)
#big loop
for dnt in range(NNT-1):
    COORDS,ATDATA,MOBILE,REFDNTCOORDS,REFDNTDATA,CGGLP=get_dinucleotide(dnt,FULLCOORDS,FULLDATA,REFCOORDS,REFDATA,GLPDATA)
    print("DNT ", dnt)
    relax_dinucleotide(COORDS,ATDATA,MOBILE,REFDNTCOORDS,REFDNTDATA,CGGLP)
    FULLCOORDS=update_coords(FULLCOORDS,dnt,COORDS,ATDATA)
    write_full_coords(FULLCOORDS,FULLDATA,OUTFILE)
#boltzmann constant is 0.0083144621kJ mol−1K−1 . 300K makes kT=2.49433863
OUTFILE.close()

OUTFILE=open("last"+args.output, "w")
write_full_coords(FULLCOORDS,FULLDATA,OUTFILE)
OUTFILE.close()
