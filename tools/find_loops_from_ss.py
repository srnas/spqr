import sys
KPHOS=500

SSTR=sys.argv[1].strip()
#print SSTR

OUT="REMARK PHPULL "+ str(KPHOS)
#out="REMARK GROUP "
#ntl=[]
#ntl.append(out)

LOOPS=[]
THISLOOP=[]
for c in xrange(0,len(SSTR)):
    if(SSTR[c]=="."):
        THISLOOP.append(c)
    if((SSTR[c]=="(" or SSTR[c]==")")):
        
        if(( 0 not in THISLOOP) and len(THISLOOP)>=5 and (SSTR[len(SSTR)-1] not in THISLOOP)):
            LOOPS.append(THISLOOP)
        THISLOOP=[]
        #    if(SSTR[c]=="(" or SSTR[c]==")"):
#        ntl.append(" "+str(c)+" ")
#print "".join(ntl)

for elem in LOOPS:
    OUT=OUT+" "+str(elem[0])
    OUT=OUT+" "+str(elem[len(elem)-1])

print OUT
