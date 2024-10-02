import numpy as np
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-i","--infile", help="Input RNAspider file",type=str,default="")
parser.add_argument("-o","--output", help="Output SPQR file",type=str,default="linked_loops.lst")
parser.add_argument("-p","--pdbfile", help="Input PDB file",type=str,default="")
parser.add_argument("-t","--sstruct", help="Secondary structure file", type=str, default="")
parser.add_argument("-v","--verbose", help="Verbose option",action='store_true')
parser.add_argument("-K1","--Kcmcm", help="Repulsion constant cm-cm",type=float,default=500)
parser.add_argument("-K2","--Kcmlp", help="Repulsion constant cm-loop",type=float,default=500)
ssflag,pdbflag=False,False

args=parser.parse_args()


Kcmcm,Kcmlp=args.Kcmcm,args.Kcmlp
HPCUTOFF,HPMAXSIZE=4,10
currresind=-1

REDUCIBLE_LINKS=['L&L', 'D&L', 'L(S)', 'D(S)']

###open 3d structure
pdbcoords,ntcoords,list3d2d=[],[],[]
cnt2d=0
if args.pdbfile!="":
    pdbflag=True
    for ii in open(args.pdbfile).readlines():
        ll=ii.split()
        ntc=[]
        if ll[0]=='ATOM':
            resind=int(ii[22:26].strip())
            chain=ii[21:22].strip()
            if cnt2d==0: currchain=chain
            if resind!=currresind and currresind>-1:
                pdbcoords.append([ntcoords,currresind])
                ntcoords=[]
                list3d2d.append(((currresind,currchain),cnt2d))
                cnt2d+=1
            atname=ii[12:16].strip()
            if atname=='P' or atname=="C4'":
                x,y,z=float(ii[30:38].strip()),float(ii[38:46].strip()),float(ii[46:54].strip())
                ntcoords.append([np.array([x,y,z]),atname])
            currresind=resind
            currchain=chain
    pdbcoords.append([ntcoords,currresind])
    list3d2d.append(((currresind,chain),cnt2d))
    cnt2d+=1
#for uu in list3d2d:
#    print(uu)

dict3d2d=dict(list3d2d)

if args.infile=="":
    print("ERROR: invalid file name")
    parser.print_help()
    exit(1)

###open secondary structure
if args.sstruct!="": ssflag=True
if ssflag:
    ALLSSFILE=[]
    for line in open(args.sstruct):
        if line[0]!=">" and line[0]!="#":
            ALLSSFILE.append(line.strip())
    fullss=list(ALLSSFILE[1])
    origss=list(ALLSSFILE[1])
    sequence=list(ALLSSFILE[0])
    NNT=len(fullss)
    SSDICT = {
        "(":")",
       "[":"]",
       "{":"}",
       "<":">"
    }
    bpairs=[]
    for br in SSDICT:
        pnt=-1
        restart=True
        while(restart):
            restart=False
            for nt in range(0,NNT):
                if(fullss[nt]==br):
                    pnt=nt
                if(fullss[nt]==SSDICT[br] and pnt>=0):
                    bpairs.append([pnt,nt])
                    fullss[nt]="."
                    fullss[pnt]="."
                    restart=True
                    break
    PAIRLIST=[]
    for nt in range(0,NNT):
        cpair=-1
        for pair in bpairs:
            if(nt==pair[0]):
                cpair=pair[1]
            if(nt==pair[1]):
                cpair=pair[0]
        PAIRLIST.append(cpair)
    ##FIND STEMS##
    opairs=bpairs.sort()
    preSTEMS,STEMS=[],[]
    CSTACK=[bpairs[0]]
    for bp in range(1,len(bpairs)):
        if(bpairs[bp][0]==bpairs[bp-1][0]+1 and bpairs[bp][1]==bpairs[bp-1][1]-1):
            CSTACK.append(bpairs[bp])
        else:
            if(len(CSTACK)>0):
                preSTEMS.append(CSTACK)
            CSTACK=[bpairs[bp]]
    if(len(CSTACK)>0):
        preSTEMS.append(CSTACK)
    for st in preSTEMS:
        unf=[]
        for pa in st:
            unf.append(pa[0])
        for pa in reversed(st):
            unf.append(pa[1])
        STEMS.append(unf)
###########################

def get_distsq(vec1,vec2):
    return np.dot(vec1-vec2,vec1-vec2)

def in_loop(nt, ntlist):
    ret = False
    for strand in range(len(ntlist)):
        if in_strand(nt,ntlist[strand]):
            ret=True
            break
    return ret

def in_strand(nt, strand):
    ret = False
    if nt<=strand[1] and nt >= strand[0]:
        ret = True
    return ret

def look_for_unpaired(ntmin,ntmax,ntl,red=True):
    segl=[]
    #here we look for unpaired nts
    for ii in range(ntmin, ntmax+1):
        if origss[ii]==".": segl.append(ii)
    #if there are not, we look for hinges
    if segl==[]:
        nta,ntb=ntmin-1,ntmax+1
        if nta<0:nta=ntmin
        if ntb>=NNT:ntb=ntmax
        for ipi in range(nta,ntb+1):
            if not belong_to_same_stem(ipi,ipi+1):
                if ipi not in segl: segl.append(ipi)
                if ipi+1 not in segl:  segl.append(ipi+1)
    #if empty, loop cannot be reduced
    if segl==[]: return segl
    #in the second case, if ip belongs to loop, we return it in the form [a,b]
    if red:
        mn,mx=min(segl),max(segl)
        segl=[mn-1,mx+1]
    return segl

def get_min_dist(ps1,ps2):
    distssq=[]
    for ii in range(2):
        for jj in range(2):
            distssq.append(get_distsq(ps1[ii][0],ps2[jj][0]))
    return np.sqrt(min(distssq))

def expand_list(li):
    ret=[]
    for ii in range(li[0],li[len(li)-1]+1):ret.append(ii)
    return ret

def find_stem(nt):
    ret=-1
    for ss in range(len(STEMS)):
        if nt in STEMS[ss]: ret=ss
    return ret

def belong_to_same_stem(nt1,nt2):
    ret=False
    st1,st2=find_stem(nt1),find_stem(nt2)
    if st1==st2: ret=True
    return ret

def reduce_loops(spqrtype, ntlist, phosugcoords,ipoints,fullntlist,modflags,ord_loops):
    #exclude cases
    if spqrtype=='st' or ((spqrtype=='il' or spqrtype=='hp') and  (ntlist[0][1]+1-ntlist[0][0]) < HPMAXSIZE ) : return spqrtype,ntlist,False
    newloop=ntlist
    newtype=spqrtype

    #we have to make sure that the ips belong to the same loop and classify them according to the strand
    inner_ips=[]
    for ip in ipoints:
        if in_loop(ip,ntlist):
            inner_ips.append(ip)
    ip_list=[]
    for stra in range(len(ntlist)):
        ipstr=[]
        for ip in inner_ips:
            if in_strand(ip,ntlist[stra]):
                ipstr.append(ip)
        ip_list.append(ipstr)
    IPINLOOP=False
    for stra in range(len(ntlist)):
        if len(ip_list[stra])>0 :IPINLOOP=True

    #LOOP DOES NOT CONTAIN INTERSECTION POINTS
    #add "fake" intersection points to strands which do not contain them
    #this requires that the other loop, which contains the IPs, has been already reduced
    
    if not IPINLOOP :
        modloop=[]
        tloop=ntlist
        unpaired_nts=[]
        for cstrand in ntlist:
            elist=look_for_unpaired(cstrand[0],cstrand[len(cstrand)-1],ntlist,red=False)
            unpaired_nts=unpaired_nts+elist

        if unpaired_nts!=[]:
            #we assume that the ips belong to the same loop and strand, which should hold for D(S) and L(S) types of link
            ipcentroid=np.zeros(3)
            set_for_centroid=ipoints
            if modflags[ord_loops[0]]==True and len(fullntlist[ord_loops[0]])==1: set_for_centroid=expand_list(fullntlist[ord_loops[0]][0])
            #this guarantees that the centroid will be the reduced loop which contains the ips, and not just the ips
            nips=0
            for ntip in set_for_centroid:
                if ntip< len(phosugcoords) and ntip>=0:
                    ipcentroid=ipcentroid+phosugcoords[ntip][0][0][0]+phosugcoords[ntip][0][1][0]
                    nips=nips+2
                ipcentroid=ipcentroid/nips
            
            selnt=unpaired_nts[0]
            mindist=get_distsq(phosugcoords[unpaired_nts[0]][0][0][0],ipcentroid)
            for ntunp in unpaired_nts:
                for at in range(2):
                    rdist=get_distsq(phosugcoords[ntunp][0][at][0],ipcentroid)
                    if rdist<mindist:
                        mindist=rdist
                        selnt=ntunp
                        tloop=[selnt-1,selnt+1]
            if tloop==ntlist:
                modloop.append(tloop)
            else:
                modloop.append(tloop)
        chflag=False
        if modloop!=ntlist and modloop!=[] : chflag=True
        if chflag :
            newloop=modloop
            if spqrtype=='hp': newtype='il'
            if spqrtype=='ju' and len(newloop)==1:
                newtype='il'
        
    if IPINLOOP:
    #CONTAINS INTERSECTION POINTS
        #we reduce the list of intersection points
        redloop=[]
        chflag=True
        for sips in ip_list:
            if len(sips)>0:
                ipmin,ipmax=sips[0],sips[0]
                if len(sips)>1:
                    for iip in sips:
                        if iip<ipmin: ipmin=iip
                        if iip>ipmax: ipmax=iip
                #now we look for unpaired regions close to the intersection point
                segment=look_for_unpaired(ipmin,ipmax,ntlist)
                if segment==[]:
                    print("Loop cannot be reduced")
                    chflag=False
                redloop.append(segment)
        if redloop==[]:chflag=False
        if chflag:
            newloop=redloop
            if spqrtype=='hp': newtype='il'
            if spqrtype=='ju' and len(newloop)==1:
                newtype='il'
    return newtype,newloop,chflag

def increase_stem(nts):
    a,b,c,d=nts[0][0],nts[0][1],nts[1][0],nts[1][1]
    if a>0 and d<NNT-1  : a,d=a-1,d+1
    if b<NNT-1 and c>0: b,c=b+1,c-1
    return [[a,b],[c,d]]

def find_loops_with_ips(ntlist, IPS):
    ret=[False,False]
    for loo in range(2):
        fl=False
        for stra in ntlist[loo]:
            for ip in IPS:
                if ip >= stra[0] and ip <=stra[len(stra)-1]: fl=True
        ret[loo]=fl

    li=[0,1]
    if ret[0]==False: li=[1,0]
    return li 


#def merge_links(orig_links):
#    new_links=orig_links
#    for ii in orig_links:
#        print(ii)
#    exit(1)
#    return new_links

def unpack_links(rawdat):
    allnts, allips, allnams,alletypes=[],[],[],[]
    for rr in rawdat:
        ll=rr[0].split(";")
        etype=ll[1]
        verbetype=ll[2]
        elem1=ll[3]+ll[4]

        ll=rr[3].split(";")
        elem2=ll[3]+ll[4]
        elem1=elem1.replace('\"','')
        elem2=elem2.replace('\"','')
        nam1=elem1.split("(")[0].split()
        nam2=elem2.split("(")[0].split()

        rnts1=elem1.split("2D:")[1].split("3D")[0]
        rnts2=elem2.split("2D:")[1].split("3D")[0]
        nts1,nts2=[],[]

        for ii in rnts1.split(","):
            jj=ii.split("-")
            nta,ntb=int(jj[0].strip()),int(jj[1].strip())
            nts1.append([nta,ntb])
        for ii in rnts2.split(","):
            jj=ii.split("-")
            nta,ntb=int(jj[0].strip()),int(jj[1].strip())
            nts2.append([nta,ntb])

        nts=[nts1,nts2]
        nams=[nam1,nam2]

        #intersection points
        IPS=[]
        CHAIN_IPS=[]
        for ii in range(6):
            ll=rr[ii].split("IP")
            if len(ll)>1:
                ipll=ll[1].split(":")
                print(ipll)

                #input numbering is in 3D style
                ipnt1,ipch1=int(ipll[1].split(",")[0].strip()[1:]),ipll[0][-1]
                ipnt1_2d=dict3d2d[(ipnt1,ipch1)]
                #ipnt1=int(ipll[1].split(",")[0].strip()[1:])-1 
                if ipnt1_2d not in IPS :
                    IPS.append(ipnt1_2d)
                    CHAIN_IPS.append(ipll[0][-1])
                ipnt2,ipch2=int(ipll[2].split(",")[0].strip()[1:]),ipll[1][-1]
                ipnt2_2d=dict3d2d[(ipnt2,ipch2)]
                #ipnt2=int(ipll[2].split(",")[0].strip()[1:])-1

                if ipnt2_2d not in IPS :
                    IPS.append(ipnt2_2d)
                    CHAIN_IPS.append(ipll[1][-1])
                print(ipnt1_2d,ipnt2_2d)
        if args.verbose:
            if  (nams[0][0]=='Loop' and len(nts1)==1) or nams[0][0]=="Single" : length1=nts1[0][1]- nts1[0][0]+1
            else : length1=nts1
            if  (nams[1][0]=='Loop' and len(nts2)==1) or nams[1][0]=="Single" : length2=nts2[0][1]- nts2[0][0]+1
            else : length2=nts2
        allnts.append(nts)
        allnams.append(nams)
        allips.append(IPS)
        alletypes.append(etype)
    return allnts,allips,allnams,alletypes

def read_links(infile):
    temp,rawdat=[],[]
    for ii in open(infile).readlines():
        fc=ii[0]
        if fc=="E": temp=[]
        temp.append(ii)
        if ii[:4]==";;;;":
            rawdat.append(temp)
    return rawdat


def order_and_reduce_loops(allnts_,allIPS_,allnams_,alletypes_):
    links,rings=[],[]
    for li in range(len(allnts_)):
        nts,IPS,nams,etype=allnts_[li],allIPS_[li],allnams_[li],alletypes_[li]
        spqrtypes=["X","X"]
        #loops with ips are considered first
        ordered_loops=find_loops_with_ips(nts, IPS)
        modified_flags=[False,False]

        for ii in ordered_loops:
            #identify loop types according to spqr
            if nams[ii][0]=="Loop":
                if len(nts[ii])==1: #hairpin
                    spqrtypes[ii]='hp'
                    #elif len(nts2)==2: #internal loop
                    #    spqrtyp1="il"
                else:
                    spqrtypes[ii]='ju'

            elif nams[ii][0]=="Dinucleotide": spqrtypes[ii]='st'
            elif nams[ii][0]=="Single": spqrtypes[ii]='il'
            else:
                print("Unrecognized loop type ", nams[ii])
                exit(1)


            if args.verbose: print("before",[spqrtypes[ii],nts[ii]])

            ##### REDUCE LINKS
            if pdbflag and ssflag and (etype in REDUCIBLE_LINKS):
                #here reduce first the loop which contains the intersection points
                spqrtypes[ii],nts[ii],modified_flags[ii]=reduce_loops(spqrtypes[ii], nts[ii], pdbcoords,IPS,nts,modified_flags,ordered_loops)
            if args.verbose: print("Contacts: ",IPS)
            if spqrtypes[ii]=='st' and ssflag:nts[ii]=increase_stem(nts[ii])

            ring_to_append=[spqrtypes[ii],nts[ii]]
            if args.verbose:print("after",ring_to_append)
            if ring_to_append not in rings: rings.append(ring_to_append)

        #correct once more
        if etype=="L(L)":
            if spqrtypes[0]=='il' and spqrtypes[1]=='il': spqrtypes=['ju','ju']
        #restore original order
        firstloop,secondloop=ordered_loops[0],ordered_loops[1]
        links.append([spqrtypes,nts,etype,[rings.index([spqrtypes[firstloop],nts[firstloop]]),rings.index([spqrtypes[secondloop],nts[secondloop]])]])
    return links, rings

def write_lst_file(links_,rings_,outname_):
    outfile=open(outname_,"w")
    outfile.write(str(len(links_))+" "+str(len(rings_))+" "+str(Kcmcm)+" "+str(Kcmlp)+"\n")
    for ii in rings_:
        lin=ii[0]
        if ii[0]=='ju':
            lin+=" "+str(len(ii[1]))
        for jj in ii[1]:
            for kk in jj:
                lin+=" "+str(kk)
        outfile.write(str(lin)+"\n")

    for ii in links_:
        ll=ii[3]
        lin=str(ll[0])+" "+str(ll[1])
        outfile.write(str(lin)+"\n")
    outfile.close()

def loops_are_adjacent(loop1, loop2, name1, name2):
    #the returning value is the number of strands that match
    #we do not include junctions so far
    ret=[]
    nstr1,nstr2=len(loop1),len(loop2)
    if nstr1==1 and nstr2==1:
        f1,e1=loop1[0][0],loop1[0][len(loop1[0])-1]
        f2,e2=loop2[0][0],loop2[0][len(loop2[0])-1]
        if e1==f2: ret=[[f1,e2]]
        if e2==f1: ret=[[f2,e1]]
    elif nstr1==1 and nstr2==2:
        f1,e1=loop1[0][0],loop1[0][len(loop1[0])-1]
        f2a,e2a=loop2[0][0],loop2[0][len(loop2[0])-1]
        f2b,e2b=loop2[1][0],loop2[1][len(loop2[1])-1]
        if e1==f2a: ret=[[f1,e2a],[f2b,e2b]]
        if e1==f2b: ret=[[f2a,e2a],[f1,e2b]]
        if f1==e2a: ret=[[f2a,e1],[f2b,e2b]]
        if f1==e2b: ret=[[f2a,e2a],[f2b,e1]]
    elif nstr1==2 and nstr2==1:
        f2,e2=loop2[0][0],loop2[0][len(loop2[0])-1]
        f1a,e1a=loop1[0][0],loop1[0][len(loop1[0])-1]
        f1b,e1b=loop1[1][0],loop1[1][len(loop1[1])-1]
        if e2==f1a: ret=[[f2,e1a],[f1b,e1b]]
        if e2==f1b: ret=[[f1a,e1a],[f2,e1b]]
        if f2==e1a: ret=[[f1a,e2],[f1b,e1b]]
        if f2==e1b: ret=[[f1a,e1a],[f1b,e2]]
    elif nstr1==2 and nstr2==2:
        f1a,e1a=loop1[0][0],loop1[0][len(loop1[0])-1]
        f1b,e1b=loop1[1][0],loop1[1][len(loop1[1])-1]
        f2a,e2a=loop2[0][0],loop2[0][len(loop2[0])-1]
        f2b,e2b=loop2[1][0],loop2[1][len(loop2[1])-1]
        if e1a==f2a and f1b==e2b: ret=[[f1a,e2a],[f2b,e1b]]
        if f1a==e2a and e1b==f2b: ret=[[f2a,e1a],[f1b,e2b]]
        if e1a==f2b and f1b==e2a: ret=[[f1a,e2b],[f2a,e1b]]
        if f1a==e2b and e1b==f2a: ret=[[f2b,e1a],[f1b,e2a]]
    if len(ret)==2:
        if ret[0][1]==ret[1][0]: ret=[[ret[0][0],ret[1][1]]]
    fnam=''
    if {name1,name2}=={'Dinucleotide','Dinucleotide'}: fnam='Dinucleotide'
    if {name1,name2}=={'Dinucleotide','Loop'}: fnam='Loop'
    if {name1,name2}=={'Dinucleotide','Single'}: fnam='Loop'
    if {name1,name2}=={'Single','Single'}: fnam='Single'
    if {name1,name2}=={'Single','Loop'}:
        if len(ret)==1: fnam='Single'
        else: fnam='Loop'
    if {name1,name2}=={'Loop','Loop'}: fnam='Loop'
    
    return ret,fnam
    


def opposite(num):
    if num==1: return 0
    if num==0: return 1

def loops_are_contained(loop1,loop2,name1,name2):
    cflag=False
    container=[]
    contname=''
    if loop1==loop2: return loop1,True,name1
    if len(loop1)==1 and len(loop2)==1:
        if loop1[0][0]<=loop2[0][0] and loop1[0][-1]>=loop2[0][-1]: container,cflag,contname=loop1,True,name1
        if loop1[0][0]>=loop2[0][0] and loop1[0][-1]<=loop2[0][-1]: container,cflag,contname=loop2,True,name2
    if len(loop1)==1 and len(loop2)==2:
        #loop1 can not contain loop2
        if loop1[0][0]>=loop2[0][0] and loop1[0][-1]<=loop2[0][-1]: container,cflag,contname=loop2,True,name2
        if loop1[0][0]>=loop2[1][0] and loop1[0][-1]<=loop2[1][-1]: container,cflag,contname=loop2,True,name2
    if len(loop2)==1 and len(loop1)==2:
        #loop2 can not contain loop1
        if loop2[0][0]>=loop1[0][0] and loop2[0][-1]<=loop1[0][-1]: container,cflag,contname=loop1,True,name1
        if loop2[0][0]>=loop1[1][0] and loop2[0][-1]<=loop1[1][-1]: container,cflag,contname=loop1,True,name1
    if len(loop2)==2 and len(loop1)==2:
        if loop2[0][0]>=loop1[0][0] and loop2[0][-1]<=loop1[0][-1] and loop2[1][0]>=loop1[1][0] and loop2[1][-1]<=loop1[1][-1] :container,cflag,contname=loop1,True,name1
        if loop1[0][0]>=loop2[0][0] and loop1[0][-1]<=loop2[0][-1] and loop1[1][0]>=loop2[1][0] and loop1[1][-1]<=loop2[1][-1] :container,cflag,contname=loop2,True,name2
    return container,cflag,contname

def merge_link_etype(typ1,typ2):
    ret=typ1
    if typ1 in REDUCIBLE_LINKS: ret=typ1
    if typ2 in REDUCIBLE_LINKS: ret=typ2
    return ret

def merge_links_if_possible(thislink,otherlink,thisinfo,otherinfo):
    #links are merged if a) they have a common loop and the other two are adjacent or contained
    newlink,newinfo,retval=thislink,thisinfo,-1

    thisname,othername=thisinfo[2][0][0],otherinfo[2][0][0]
    common_loop,common_flag,common_name=loops_are_contained(thislink[0],otherlink[0],thisname,othername)
    if common_flag: loop1tomerge,loop2tomerge,l1tmname,l2tmname=thislink[1],otherlink[1],thisinfo[2][1][0],otherinfo[2][1][0]
    else:
        thisname,othername=thisinfo[2][1][0],otherinfo[2][0][0]
        common_loop,common_flag,common_name=loops_are_contained(thislink[1],otherlink[0],thisname,othername)
        if common_flag: loop1tomerge,loop2tomerge,l1tmname,l2tmname=thislink[0],otherlink[1],thisinfo[2][0][0],otherinfo[2][1][0]
        else:
            thisname,othername=thisinfo[2][0][0],otherinfo[2][1][0]
            common_loop,common_flag,common_name=loops_are_contained(thislink[0],otherlink[1],thisname,othername)
            if common_flag: loop1tomerge,loop2tomerge,l1tmname,l2tmname=thislink[1],otherlink[0],thisinfo[2][1][0],otherinfo[2][0][0]
            else:
                thisname,othername=thisinfo[2][1][0],otherinfo[2][1][0]
                common_loop,common_flag,common_name=loops_are_contained(thislink[1],otherlink[1],thisname,othername)
                if common_flag: loop1tomerge,loop2tomerge,l1tmname,l2tmname=thislink[0],otherlink[0],thisinfo[2][0][0],otherinfo[2][0][0]
    if not common_flag : return newlink,thisinfo,retval
    merged_loop,merged_name=loops_are_adjacent(loop1tomerge,loop2tomerge,l1tmname,l2tmname)
    if merged_loop==[]: merged_loop,merged_flag,merged_name=loops_are_contained(loop1tomerge,loop2tomerge,l1tmname,l2tmname)
    if merged_loop!=[]:
        newlink=[common_loop,merged_loop]
        #here we merge the information
        thisIPS,otherIPS=thisinfo[0],otherinfo[0]
        
        #get names of common loop and loops to merge
        thisetype,otheretype=thisinfo[1],otherinfo[1]
        thisnames,othernames=thisinfo[2],otherinfo[2]
        newIPS=thisIPS+[data for data in otherIPS if data not in thisIPS]
        newIPS.sort()
        newnames,newetype=[[common_name,'',''],[merged_name,'','']],merge_link_etype(thisinfo[1],otherinfo[1])
        newinfo=[newIPS,newetype,newnames]
        retval=1
    return newlink,newinfo,retval
                            
def do_merge_links(allnts, allIPS,allnams,alletypes):
    newlist,newIPS,newnams,newetypes=[],[],[],[]
    newlinks,looplist,linklist,loopnames=[],[],[],[]
    #we have a list with all the loops
    for li in range(len(allnts)):
        for iili in range(2):
            if allnts[li][iili] not in looplist:
                looplist.append(allnts[li][iili])
                loopnames.append(allnams[li][iili])
    for li in range(len(allnts)):
        #search indexes of loops
        loop1,loop2=-1,-1
        for iloo in range(len(looplist)):
            if allnts[li][0]==looplist[iloo]:
                loop1=iloo
            if allnts[li][1]==looplist[iloo]:
                loop2=iloo
        linklist.append([[loop1,loop2],allIPS[li],alletypes[li],allnams[li]])

    templist=linklist
    newlist=[]

    while len(templist)>0:
        currlink=[looplist[templist[0][0][0]],looplist[templist[0][0][1]]]
        currinfo=[templist[0][1],templist[0][2],templist[0][3]]
        mergetimes=1
        
        while mergetimes>0:
            mergetimes=0
            towipe=[]
            for li in range(1,len(templist)):
                exttemplink=[looplist[templist[li][0][0]],looplist[templist[li][0][1]]]
                extinfo=[templist[li][1],templist[li][2],templist[li][3]]
                currlink,currinfo,merged=merge_links_if_possible(currlink,exttemplink,currinfo,extinfo)
                if merged!=-1:
                    mergetimes+=1
                    towipe.append(li)
            tlist=[]
            for rl in range(1,len(templist)):
                if rl not in towipe: tlist.append(templist[rl])
            
            templist=[currlink]+tlist
            #if tlist==[]:
            #    mergetimes=0
            #    templist=[currlink]+tlist

        newlist.append(currlink)
        newIPS.append(currinfo[0])
        newetypes.append(currinfo[1])
        newnams.append(currinfo[2])
        if templist!=[]:
            templist.pop(0)

    return newlist,newIPS,newnams,newetypes

############## MAIN ##############
#READ LINKS
rawdata=read_links(args.infile)
allnts, allIPS, allnams,alletypes=unpack_links(rawdata)

NEWnts,NEWIPS,NEWnams,NEWetypes=do_merge_links(allnts, allIPS,allnams,alletypes)
#print(NEWnts)



#print (merge_links(allnts, allIPS,allnams,alletypes))

links, rings=order_and_reduce_loops(NEWnts,NEWIPS,NEWnams,NEWetypes)

write_lst_file(links, rings, args.output)

######################
