from math import sqrt,pi,atan2,cos

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

GTEMPLATE=["P","OP1","OP2","O5\'","C5\'","C4\'","O4\'","C3\'","O3\'","C2\'","O2\'","C1\'","N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"]
#CTEMPLATE=["P","O1P","O2P","O5\'","C5\'","C4\'","O4\'","C1\'","N1","C6","C5","C4","N4","N3","C2","O2","C3\'","C2\'","O2\'","O3\'"]
#UTEMPLATE=["P","O1P","O2P","O5\'","C5\'","C4\'","O4\'","C1\'","N1","C6","C5","C4","O4","N3","C2","O2","C3\'","C2\'","O2\'","O3\'"]
#ATEMPLATE=["P","O1P","O2P","O5\'","C5\'","C4\'","O4\'","C1\'","N9","C8","N7","C5","C6","N6","N1","C2","N3","C4","C3\'","C2\'","O2\'","O3\'"]
CTEMPLATE=["P","OP1","OP2","O5\'","C5\'","C4\'","O4\'","C1\'","N1","C6","C5","C4","N4","N3","C2","O2","C3\'","C2\'","O2\'","O3\'"]
UTEMPLATE=["P","OP1","OP2","O5\'","C5\'","C4\'","O4\'","C1\'","N1","C6","C5","C4","O4","N3","C2","O2","C3\'","C2\'","O2\'","O3\'"]
ATEMPLATE=["P","OP1","OP2","O5\'","C5\'","C4\'","O4\'","C1\'","N9","C8","N7","C5","C6","N6","N1","C2","N3","C4","C3\'","C2\'","O2\'","O3\'"]
DEF_ORDER=[ATEMPLATE,UTEMPLATE,GTEMPLATE,CTEMPLATE]
#DEF_OPTIONAL=[[],[],["OP1","OP2"],[]]

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

def res_as_int(base):
    #A=0,U=1,G=2,C=3
    ret=-1
    if(base=="A"):
        ret=0
    if(base=="U"):
        ret=1
    if(base=="G"):
        ret=2
    if(base=="C"):
        ret=3
    return ret

def trans_glp(glp):
    pu,gl=-1,-1
    if(glp[0]=="A"):
        gl=0
    elif(glp[0]=="H"):
        gl=1
    elif(glp[0]=="S"):
        gl=2
    if(glp[1]=="3"):
        pu=0
    if(glp[1]=="2"):
        pu=1
    return [gl,pu]
    
def class_glp(chi, delta, res):
    #gl: 0=A, 1=H, 2=S
    #pk: 0=C3'endo, 1=C2'endo
    gl="X"
    pk="X"
    dchi=chi*180.0/pi
    ddel=delta*180.0/pi
    if((dchi>=-180 and dchi<=-120) or (dchi>155 and dchi<=180)):
        gl="A"
    elif(dchi>=-120 and dchi<=-10):
        gl="H"
    elif(dchi>=35 and dchi <= 145 and (res=="A" or res=="G")):
        gl="S"
    if(delta<DEL_DIV):
        pk="3"
    else:
        pk="2"
    #print res,dchi, gl, ddel, pk
    return [gl,pk,chi*180/pi,delta*180/pi]

def get_cg_beads(a,b,c,sug, ph,bas, resind, chain, glp):
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
    cg_glp=trans_glp(glp)
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
            pho_pos=PHO_ARR_R[cg_glp[0]][cg_glp[1]]
        else:
            pho_pos=PHO_ARR_Y[cg_glp[0]][cg_glp[1]]
        phox=pho_pos[0]*Xx + pho_pos[1]*Yx + pho_pos[2]*Zx + cmx
        phoy=pho_pos[0]*Xy + pho_pos[1]*Yy + pho_pos[2]*Zy + cmy
        phoz=pho_pos[0]*Xz + pho_pos[1]*Yz + pho_pos[2]*Zz + cmz
    else:
        phox=ph[0]
        phoy=ph[1]
        phoz=ph[2]

    if bas == 'A' or bas == 'G' :
        sug_pos=SUG_ARR_R[cg_glp[0]]
    else:
        sug_pos=SUG_ARR_Y[cg_glp[0]]
    sugx=sug_pos[0]*Xx + sug_pos[1]*Yx + sug_pos[2]*Zx + cmx
    sugy=sug_pos[0]*Xy + sug_pos[1]*Yy + sug_pos[2]*Zy + cmy
    sugz=sug_pos[0]*Xz + sug_pos[1]*Yz + sug_pos[2]*Zz + cmz
    
    coords=[CV,[Xvecx,Xvecy,Xvecz],[Yvecx,Yvecy,Yvecz],[sugx,sugy,sugz],[phox,phoy,phoz]]
    return[coords,bas,chain,resind]
    
def read_pdb_atoms(ALLFILE,cgflag=0):
    fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
    fch3,fch4=0,0
    fd1,fd2,fd4=0,0,0
    ##we check the number of nucleotides
    NTDAT=[]
    BLOCKRES=[]
    for lin in xrange(0,len(ALLFILE)):
        base=ALLFILE[lin][17:20].strip()
        resind=int(ALLFILE[lin][22:26].strip())
        chain=ALLFILE[lin][21].strip()
        icode=ALLFILE[lin][26].strip()
        dat=[base, resind, chain,icode]
        if dat not in NTDAT:
            NTDAT.append(dat)
    ATLIST=[]
    ATLABEL=[]
    for at in xrange(0,len(ALLFILE)):
        ATLIST.append([])
        ATLABEL.append([])
    for lin in xrange(0,len(ALLFILE)):
        resind=int(ALLFILE[lin][22:26].strip())
        atname=ALLFILE[lin][12:16].strip()
        for res in xrange(0,len(NTDAT)):
            if(resind==NTDAT[res][1]):
                ATLIST[res].append(lin)
                ATLABEL[res].append([lin,atname])
    for re in xrange(0,len(NTDAT)):
        NTCOORDS=[]
        for at in ATLIST[re]:
            line=ALLFILE[at]
            resind=int(line[22:26].strip())
            atname=line[12:16].strip()
            base=line[17:20].strip()
            x=float(line[30:38].strip())
            y=float(line[38:46].strip())
            z=float(line[46:54].strip())
            NTCOORDS.append([x,y,z])
            if atname == "C2":
                if(fc2<1):
                    c2pos = [x,y,z]
                    fc2=1
                    if(base=="C" or base=="U"):
                        chi4pos=[x,y,z]
                        fch4=1
            elif atname == "C4":
                if(fc4<1):
                    c4pos = [x,y,z]
                    fc4=1
                    #if(base=="C" or base=="U"):
                    if(base=="A" or base=="G"):
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
        ###TO GET THE SPQR COORDINATES###
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
            glp=class_glp(chi,delta,base)
            if(cgflag):
                cg_res=get_cg_beads(c2pos,c4pos,c6pos,sugpos,ppos,base,resind,chain,glp)
            fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
        else:
            print "Residue " + str(re) + " incomplete!"
            quit(1)
            fc2,fc4,fc6,fp,fs1,fs2,fs3,fs4,fs5=0,0,0,0,0,0,0,0,0
            new_resind=new_resind+1
        #######PACK EVERYTHING#########
        if(len(ATLIST[re])!=len(NTCOORDS)):
           print "ERROR: missing atoms in residue "+str(re)
           quit(2)
           
        NEWATLABEL=[]
        NEWNTCOORDS=[]
        intres=res_as_int(NTDAT[re][0])
        for ia in xrange(0,len(DEF_ORDER[intres])):
            for newi in xrange(0,len(ATLABEL[re])):
                #if(ATLABEL[re][newi][1] in DEF_OPTIONAL[intres]):
                #    NEWATLABEL.append(ATLABEL[re][newi])
                #    NEWNTCOORDS.append(NTCOORDS[newi])
                #    break
                if(ATLABEL[re][newi][1]==DEF_ORDER[intres][ia]):
                    NEWATLABEL.append(ATLABEL[re][newi])
                    NEWNTCOORDS.append(NTCOORDS[newi])
                    break
                
        WRESIDUE=[NTDAT[re],NEWATLABEL,NEWNTCOORDS,glp]
        if(cgflag):
           WRESIDUE.append(cg_res)
        BLOCKRES.append(WRESIDUE)
    return BLOCKRES

def mimic_nt(res,template):
    NEWLABEL=[]
    NEWNTCOORDS=[]
    for da in xrange(0,len(template[1])):
        for ea in xrange(0,len(res[1])):
            if(template[1][da][1]==res[1][ea][1]):
                NEWLABEL.append(res[1][ea])
                NEWNTCOORDS.append(res[2][ea])
                break
    RET=[res[0],NEWLABEL,NEWNTCOORDS,res[3]]
    if(len(res)>4):
        RET.append(res[4])
    return RET

def print_nt(WRESIDUE):
    for at in xrange(0,len(WRESIDUE[1])):
        record="ATOM  "
        atindex='{:5}'.format(int(WRESIDUE[1][at][0]))
        atname='{:4}'.format(WRESIDUE[1][at][1])
        resname='{:3}'.format(WRESIDUE[0][0])
        chindex=WRESIDUE[0][2]
        resindex='{:4}'.format(WRESIDUE[0][1])
        xtemp, ytemp, ztemp=WRESIDUE[2][at][0],WRESIDUE[2][at][1],WRESIDUE[2][at][2]
        xpos,ypos,zpos='{:8.3f}'.format(xtemp),'{:8.3f}'.format(ytemp),'{:8.3f}'.format(ztemp)
        line=record+atindex+" "+atname+" "+resname+" "+chindex+resindex+"    "+xpos+ypos+zpos
        print line
