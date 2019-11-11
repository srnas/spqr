from math import sqrt,pi,atan2,cos
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-r","--renum",help="Renumbering of atoms", action="store_true")
parser.add_argument("filename", help="The pdb file name")
parser.add_argument("-n","--nframes", help="Number of frames",type=int,default=1)
args=parser.parse_args()

if(args.renum):
    RENUM=1
else:
    RENUM=0
NFRAMES=args.nframes
FILENAME=args.filename



visited = {}
ln=0
c2pos={}
c4pos={}
c6pos={}
ppos={}
prevncl=1
fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
fch3,fch4=0,0
fd1,fd2,fd4=0,0,0
N_PARTS_PER_NT=5
RINDS=[]

#the following are in the order: A3, H3, S3, A2, H2, S2
RSUG_A=[-1.56584, -4.8467,0.709821]
RSUG_H=[-1.56584, -4.8467,0.709821]
RSUG_S=[-0.0632599, -4.345263,0.0595195]
YSUG_A=[1.05136, -3.42329,0.652454]
YSUG_H=[1.05136, -3.42329,0.652454]

RPHO_A3=[ -6.10965 , -4.45512, 0.776946 ]
RPHO_H3=[ -4.25369 , -4.3032 , -3.52153 ]
RPHO_S3=[ 2.37143  , -3.66344, 4.03475  ]
RPHO_A2=[ -5.62207 , -4.84927, 0.322851 ]
RPHO_H2=[ -4.72923 , -5.16856, -2.15611 ]
RPHO_S2=[ 3.9117, -3.80641, 1.08371  ]

YPHO_A3=[ -3.36428, -4.5741 , 0.573644 ]
YPHO_H3=[ -1.47508, -4.3344 , -3.66775 ]
YPHO_A2=[ -2.66998, -5.03771, 0.298707 ]
YPHO_H2=[ -2.2033 ,-5.17071 , -1.47934 ]

PHO_ARR_R=[[RPHO_A3,RPHO_A2],[RPHO_H3,RPHO_H2],[RPHO_S3,RPHO_S2]]
PHO_ARR_Y=[[YPHO_A3,YPHO_A2],[YPHO_H3,YPHO_H2]]
SUG_ARR_R=[RSUG_A,RSUG_H,RSUG_S]
SUG_ARR_Y=[YSUG_A,YSUG_H]

DEFAULT_CHI=179.0*pi/180.0
DEFAULT_DEL=80.0*pi/180.0
DEL_DIV=100.0*pi/180.0

def mod3(x,y,z):
    nx=float(x)
    ny=float(y)
    nz=float(z)
    return sqrt(nx*nx+ny*ny+nz*nz)
def mod(vec):
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
def dot_pro(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]        
def cro_pro(xv,yv):
    return [xv[1]*yv[2]-xv[2]*yv[1],xv[2]*yv[0]-xv[0]*yv[2],xv[0]*yv[1]-xv[1]*yv[0]]
def sqdist(v1, v2):
    return (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) + (v1[2]-v2[2])*(v1[2]-v2[2])
def get_dih(a,b,c,d):
    b1=[b[0]-a[0],b[1]-a[1],b[2]-a[2]]
    b2=[c[0]-b[0],c[1]-b[1],c[2]-b[2]]
    b3=[d[0]-c[0],d[1]-c[1],d[2]-c[2]]
    le=dot_pro(b2,cro_pro(cro_pro(b1,b2),cro_pro(b2,b3)))/mod(b2)
    ri=dot_pro(cro_pro(b1,b2),cro_pro(b2,b3))
    return atan2(le, ri)

def class_glp(chi, delta, res):
    #gl: 0=A, 1=H, 2=S
    #pk: 0=C3'endo, 1=C2'endo
    gl=0
    pk=0
    dchi=chi*180.0/pi
    ddel=delta*180.0/pi
    if((dchi>=-180 and dchi<=-120) or (dchi>155 and dchi<=180)):
        gl=0
    elif(dchi>=-120 and dchi<=-10):
        gl=1
    elif(dchi>=35 and dchi <= 145 and (res=="A" or res=="G")):
        gl=2
    if(delta<DEL_DIV):
        pk=0
    else:
        pk=1

    return [gl,pk]

def get_cg_beads(a,b,c,sug, ph,bas, resind, chain, cg_glp):
    cmx=(float(a[0])+float(b[0])+float(c[0]))/3.0
    cmy=(float(a[1])+float(b[1])+float(c[1]))/3.0
    cmz=(float(a[2])+float(b[2])+float(c[2]))/3.0
    cm=[cmx,cmy,cmz]
    Xx=float(a[0])-cmx
    Xy=float(a[1])-cmy
    Xz=float(a[2])-cmz
    modX=mod3(Xx,Xy,Xz)
    Xx=Xx/modX
    Xy=Xy/modX
    Xz=Xz/modX
    
    #if(sug[0]==0 and sug[1]==0 and sug[2]==0):
    #    sug=[cmx, cmy, cmz]
    
    if bas == 'A' or bas == 'G' :
        Yx=float(c[0])-cmx
        Yy=float(c[1])-cmy
        Yz=float(c[2])-cmz
    elif bas == 'C' or bas == 'U':
        Yx=float(b[0])-cmx
        Yy=float(b[1])-cmy
        Yz=float(b[2])-cmz
    else:
        print "Base ", bas, " not recognized!"
        exit
    xdoty=Xx*Yx+Xy*Yy+Xz*Yz
    Yx=Yx-xdoty*Xx
    Yy=Yy-xdoty*Xy
    Yz=Yz-xdoty*Xz
    
    modY=mod3(Yx,Yy,Yz)
    Yx=Yx/modY
    Yy=Yy/modY
    Yz=Yz/modY
    
    Zx=Xy*Yz-Xz*Yy
    Zy=Xz*Yx-Xx*Yz
    Zz=Xx*Yy-Xy*Yx
    
    CV=[cmx, cmy,cmz]
    XV=[Xx,Xy,Xz]
    YV=[Yx,Yy,Yz]
    ZV=[Zx,Zy,Zz]
    if bas == 'A':
        typ="0"
    if bas == 'U':
        typ="1"
    if bas == 'G':
        typ="2"
    if bas == 'C':
        typ="3"

    
    sugx, sugy, sugz=sug[0], sug[1], sug[2]
    Xvecx,Xvecy,Xvecz=Xx+cmx,Xy+cmy,Xz+cmz
    Yvecx,Yvecy,Yvecz=Yx+cmx,Yy+cmy,Yz+cmz
    
    
    #phosphate
    if(ph[0]==0 and ph[1]==0 and ph[2]==0):
        if bas == 'A' or bas == 'G' :
            pho_pos=PHO_ARR_R[glp[0]][glp[1]]
        else:
            pho_pos=PHO_ARR_Y[glp[0]][glp[1]]
        phox=pho_pos[0]*Xx + pho_pos[1]*Yx + pho_pos[2]*Zx + cmx
        phoy=pho_pos[0]*Xy + pho_pos[1]*Yy + pho_pos[2]*Zy + cmy
        phoz=pho_pos[0]*Xz + pho_pos[1]*Yz + pho_pos[2]*Zz + cmz
    else:
        phox=ph[0]
        phoy=ph[1]
        phoz=ph[2]

    if bas == 'A' or bas == 'G' :
        sug_pos=SUG_ARR_R[glp[0]]
    else:
        sug_pos=SUG_ARR_Y[glp[0]]
    sugx=sug_pos[0]*Xx + sug_pos[1]*Yx + sug_pos[2]*Zx + cmx
    sugy=sug_pos[0]*Xy + sug_pos[1]*Yy + sug_pos[2]*Zy + cmy
    sugz=sug_pos[0]*Xz + sug_pos[1]*Yz + sug_pos[2]*Zz + cmz
    
    coords=[CV,[Xvecx,Xvecy,Xvecz],[Yvecx,Yvecy,Yvecz],[sugx,sugy,sugz],[phox,phoy,phoz]]
    return[coords,bas,chain,resind]
    
def print_cg_res(cgres, cur_at, cg_glp):
 #build the pdbline
    coords=cgres[0]
    bas=cgres[1]
    chain=cgres[2]
    resind=int(cgres[3])
    cm=coords[0]
    xvec=coords[1]
    yvec=coords[2]
    sug=coords[3]
    pho=coords[4]
    
    resname='{:3}'.format(bas)
    if(RENUM==1):
        chindex='{:1}'.format(chain)
        resindex='{:4}'.format(resind)
    else:
        chindex='{:1}'.format(RINDS[resind][1])
        resindex='{:4}'.format(RINDS[resind][0])

    gl,pk='A','3'
    if(cg_glp[0]==0):
        gl="A"
    elif(cg_glp[0]==1):
        gl="H"
    else:
        gl="S"
    if(cg_glp[1]==0):
        pk="3"
    else:
        pk="2"
    glp=" "+gl+pk+"NA"
    record="ATOM  "
    
    #atindex='{:5}'.format(resind*N_PARTS_PER_NT)
    atindex='{:5}'.format(cur_at)
    atname="BASE"
    xtemp, ytemp, ztemp=cm[0], cm[1], cm[2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos+" "+glp
    print line
    
    #atindex='{:5}'.format(resind*N_PARTS_PER_NT+1)
    atindex='{:5}'.format(cur_at+1)
    atname="XVEC"
    xtemp, ytemp, ztemp=xvec[0],xvec[1],xvec[2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    #atindex='{:5}'.format(resind*N_PARTS_PER_NT+2)
    atindex='{:5}'.format(cur_at+2)
    atname="YVEC"
    xtemp, ytemp, ztemp=yvec[0],yvec[1],yvec[2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    #atindex='{:5}'.format(resind*N_PARTS_PER_NT+3)
    atindex='{:5}'.format(cur_at+3)
    atname="SUGR"
    xtemp, ytemp, ztemp=sug[0],sug[1],sug[2]
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
    
    #atindex='{:5}'.format(resind*N_PARTS_PER_NT+4)
    atindex='{:5}'.format(cur_at+4)
    atname="PHOS"
    xtemp, ytemp, ztemp=pho[0],pho[1],pho[2]
    nidx, nidy, nidz=len(str(int(xtemp))),len(str(int(ytemp))),len(str(int(ztemp)))
    xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
    line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
    print line
        
    return


ALLFILE=[]
eff_ats=0
modflag=0

#EMAP_SUG=0
#ODIF_NTS=0

#for line in open(sys.argv[1]):
for line in open(FILENAME):
    nam=line[:4]
    if(nam.strip()=="ATOM"):
        base=line[17:20].strip()
        if(base=="A" or base=="G" or base=="U" or base=="C" or base=="RGN"):
            ALLFILE.append(line)
            if(modflag==0):
                eff_ats=eff_ats+1
    elif(nam.strip()=="ENDM"):
        modflag=1

            
INDSEQS=[]                
INDS=[]
SEQ=[]
n9fl,o6fl, o4fl, n4fl=0,0,0,0
prevchain=0
prevind=0
first=0
nind=0
at_renum=0
res_renum=0
res_list=[]
single_res=[]
model=0
RINDS=[]
endchk=0
filecnt=0

for line in ALLFILE:
    #if(line.strip().split()[0]=="MODEL"):
    if(filecnt==0):
        #opening="REMARK FRAME "+line.split()[1]
        #reset stuff
        INDSEQS=[]                
        INDS=[]
        RINDS=[]
        SEQ=[]
        n9fl,o6fl, o4fl, n4fl=0,0,0,0
        prevchain=0
        prevind=0
        first=0
        nind=0
        at_renum=0
        res_renum=0
        res_list=[]
        single_res=[]
        #opening="MODEL FRAME "+line.split()[1]
        model=model+1
        #print opening
    if(line.strip().split()[0]=="ATOM"):
        filecnt=filecnt+1
        #atind=line[6:11].strip()
        base=line[17:20].strip()
        resind=int(line[22:26].strip())
        chain=line[21].strip()
        if(first==0):
            first=1
            prevind=resind
            prevchain=chain
            SEQ.append(base)
            INDS.append(nind)
            RINDS.append([resind,chain])
        atname=line[12:16].strip()
        if(resind!=prevind or prevchain!=chain):
            nind=nind+1
            SEQ.append(base)
            INDSEQS.append([n9fl,o6fl, o4fl, n4fl])
            n9fl,o6fl, o4fl, n4fl=0,0,0,0
            prevind=resind
            prevchain=chain
            INDS.append(nind)
            RINDS.append([resind,chain])
            res_list.append(single_res)
            single_res=[]
        single_res.append(at_renum)
        at_renum=at_renum+1
        if atname == "N9":
            n9fl=1
        if atname == "O6":
            o6fl=1
        if atname == "O4":
            o4fl=1
        if atname == "N4":
            n4fl=1
    if(filecnt==eff_ats):
        filecnt=0
        #break
INDSEQS.append([n9fl,o6fl, o4fl, n4fl])
res_list.append(single_res)

NATS=0
for re in xrange(0,len(INDS)):
    NATS=NATS+len(res_list[re])

curr_at=0
for frame in xrange(0,NFRAMES):
    prev_chain=0
    chain=0
    new_resind=0
    opening="MODEL FRAME "+str(frame+1)
    #print opening
    for re in xrange(0,len(INDS)):
        for at in res_list[re]:
            line=ALLFILE[at+frame*(2+NATS)]
            if(at==res_list[re][0]):
                this_chain=line[21:22]
                if(re==0):
                    prev_chain=this_chain
                if(prev_chain!=this_chain):
                    chain=chain+1
            resind=int(line[22:26].strip())
            atname=line[12:16].strip()
            base=line[17:20].strip()
            x=float(line[30:38].strip())
            y=float(line[38:46].strip())
            z=float(line[46:54].strip())
            if atname == "C2":
                if(fc2<1):
                    c2pos = [x,y,z]
                    fc2=1
#                    if(base=="A" or base=="G"):
                    if(base=="C" or base=="U"):
                        chi4pos=[x,y,z]
                        fch4=1
            elif atname == "C4":
                if(fc4<1):
                    c4pos = [x,y,z]
                    fc4=1
                    if(base=="A" or base=="G"):
                    #if(base=="C" or base=="U"):
                        chi4pos=[x,y,z]
                        fch4=1
            elif atname == "C6":
                if(fc6<1):
                    c6pos = [x,y,z]
                    fc6=1
            elif atname == "P":
                if(fp<1):
                    ppos = [x,y,z]
                    fp=1
            elif(atname == "C1\'" or atname == "C1*"):
                if(fs1<1):
                    s1pos = [x,y,z]
                    fs1=1
            elif(atname == "C2\'" or atname == "C2*"):
                if(fs2<1):
                    s2pos = [x,y,z]
                    fs2=1
            elif(atname == "C3\'" or atname == "C3*"):
                if(fs3<1):
                    s3pos = [x,y,z]
                    fs3=1
            elif(atname == "C4\'" or atname == "C4*"):
                if(fs4<1):
                    s4pos = [x,y,z]
                    d2pos = [x,y,z]
                    fs4=1
                    fd2=1
            elif(atname == "O4\'" or atname == "O4*"):
                if(fs5<1):
                    s5pos = [x,y,z]
                    fs5=1
            elif(atname == "C5\'" or atname == "C5*"):
                if(fd1<1):
                    d1pos = [x,y,z]
                    fd1=1
#            elif(atname == "C4\'" or atname == "C4*"):
#                if(fd2<1):
#                    d2pos = [x,y,z]
#                    fd2=1
            elif(atname == "O3\'" or atname == "O3*"):
                if(fd4<1):
                    d4pos = [x,y,z]
                    fd4=1
            if(base=="A" or base=="G"):
                if(atname == "N9"):
                    if(fch3<1):
                        chi3pos = [x,y,z]
                        fch3=1
            elif(base=="U" or base=="C"):
                if(atname == "N1"):                
                    if(fch3<1):
                        chi3pos = [x,y,z]
                        fch3=1
        prev_chain=this_chain
        #if(fc2==1 and fc4==1 and fc6==1 and sp1==1 and sp2==1 and sp3==1 and sp4==1 and sp5==1):
        if(fc2==1 and fc4==1 and fc6==1):
            if(fp==0):
                ppos=[0,0,0]
            if(fs1==0 or fs2==0 or fs3==0 or fs4==0 or fs5==0):
                sugpos=[0,0,0]
            else:
                sugposx=(s1pos[0]+s2pos[0]+s3pos[0]+s4pos[0]+s5pos[0])/5.0
                sugposy=(s1pos[1]+s2pos[1]+s3pos[1]+s4pos[1]+s5pos[1])/5.0
                sugposz=(s1pos[2]+s2pos[2]+s3pos[2]+s4pos[2]+s5pos[2])/5.0
                sugpos=[sugposx,sugposy,sugposz]
            newbas='X'
            if(INDSEQS[new_resind][0]==1):
                if(INDSEQS[new_resind][1]==1):
                    newbas='G'
                else:
                    newbas='A'
            else:
                if(INDSEQS[new_resind][2]==1):
                    newbas='U'
                else:
                    newbas='C'

            #calculate dihedral angles
            if(fs5==1 and fs1==1 and fch3==1 and fch4==1):
                chi=get_dih(s5pos, s1pos, chi3pos, chi4pos)
            else:
                chi=DEFAULT_CHI
            fch3,fch4=0,0
            if(fd1==1 and fd2==1 and fs3==1 and fd4==1):
                delta=get_dih(d1pos, d2pos, s3pos, d4pos)
            else:
                delta=DEFAULT_DEL
            fd1,fd2,fd4=0,0,0
            glp=class_glp(chi,delta,newbas)
            #print newbas, new_resind
            cg_res=get_cg_beads(c2pos,c4pos,c6pos,sugpos,ppos,newbas,new_resind,chain,glp)
            print_cg_res(cg_res,curr_at,glp)
            curr_at=curr_at+N_PARTS_PER_NT
            fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
            new_resind=new_resind+1
        else:
            #print "Residue " + str(re) + " incomplete!"
            fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
            new_resind=new_resind+1
    ending="ENDMDL"
    #print ending
