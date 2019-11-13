set -e
VER=0.9.5.1
SUFF=e
SPQR_PATH=/home/spoblete/source_codes/SPQR_$VER
ENERGY_PATH=$SPQR_PATH/interactions
SRC_PATH=$SPQR_PATH/src_$VER$SUFF
SRC_PP_PATH=$SPQR_PATH/src_${VER}${SUFF}_nl_fp
TOOLS_PATH=$SPQR_PATH/tools
NPROC=3

PPULL_NMC=20000
PPULL_NSA=10
PPULL_TSA=1

CLASH_NMC=20000
CLASH_NSA=10
CLASH_TSA=1

ERMSD_NMC=5000
ERMSD_NSA=20
ERMSD_TSA=5
ERMSD_NTO=$((ERMSD_NMC*ERMSD_NSA))

ENERG_NMC=1000
ENERG_NSA=5
ENERG_TSA=1

BINCLASH=MPI_SPQR_wrmSA_$VER$SUFF
BINPPCLASH=MPI_SPQR_wrmPHSA_$VER${SUFF}NL
BINANNOT=SPQR_ENERG_$VER$SUFF
BINERMSD=MPI_SPQR_ermsdSA_$VER$SUFF
BINENERG=MPI_SPQR_SA_$VER$SUFF
ERMSD_PREF=10000
ERMSD_CUTOFF=10

######DEFAULTS######
PDBSTRUCT=""
FSPQRCONV=""
OVERWRITE=""
SECSTRNAME=""
TARGETSTRCT=""
NAME="j`printf %06d $RANDOM`"
USAGE="Usage : ./do_mini.bash -i <input.pdb> -o <jobname> [-t <target.pdb>] [-d <dssrfile>] [-k]"
CRINIT=init.pdb
MPROC=`echo $NPROC | awk '{print $1-1}'`
######LABELS########
LKR=KR
LCR=CR
LEP=EP
LEM=EM
##########READ ARGUMENTS##############
while (( "$#" )); do
  case "$1" in
    -i)
	PDBSTRUCT=$2
	shift 2
	;;
    -o)
	NAME=$2
	shift 2
	;;
    -c)
	FSPQRCONv=1
	shift 2
	;;
    -t)
	TARGETSTRCT=$2
	shift 2
	;;
    -d)
	SECSTRNAME=$2
	shift 2
	;;
    -k)
	OVERWRITE=1 
	shift 1
	;;
    -h)
	echo $USAGE
	shift 2
	;;
    --) # end argument parsing
	shift
	break
	;;
    -*|--*=) # unsupported flags
    echo "Error: Unsupported flag $1" >&2
    exit 1
    ;;
  esac
done
WORKDIR=opt_$NAME
if [ !  $PDBSTRUCT ] || [ ! -f $PDBSTRUCT ] ; then echo "ERROR: Input file not found!" >&2 ; echo $USAGE ; exit 1 ; fi
if [ -d $WORKDIR ] && [ ! $OVERWRITE  ] ; then echo "ERROR: There is already a minimization with the name $NAME. Rename that directory or overwrite with the -k option." >&2 ; echo $USAGE ; exit 1 ; fi
if [ ! $TARGETSTRCT ] ; then TARGETSTRCT=$PDBSTRUCT ; fi
if [ ! -f $TARGETSTRCT ] ; then echo "ERROR: File $TARGETSTRCT not found!" >&2 ; echo $USAGE ; exit 1 ; fi
if [ $SECSTRNAME ] && [ ! -f $SECSTRNAME ] ; then echo "ERROR: File $SECSTRNAME not found!" >&2; echo $USAGE ; exit 1 ; fi
if [ $OVERWRITE ] ; then echo "WARNING: Overwriting previous result : Deleting directory $WORKDIR!" ; fi

######################################

if [ $OVERWRITE ] ; then rm -r $WORKDIR ; fi
mkdir $WORKDIR
if [ $FSPQRCONV ] ; then
    echo "Creating SPQR pdb format structure..."
    cp $TOOLS_PATH/pdb2spqr.py .
    python pdb2spqr.py $PDBSTRUCT > $WORKDIR/init.pdb
    rm pdb2spqr.py
else cp $PDBSTRUCT $WORKDIR/init.pdb
fi

############ STEP 0 : KNOT REMOVAL #########
if [ $SECSTRNAME ] ; then
    SIND=$LKR
    echo "Starting with knot removal. Good luck."
    echo "Creating files..."
    mkdir $WORKDIR/$SIND
    cp files/phpull.spqr $WORKDIR/$SIND/input.spqr
    mkdir $WORKDIR/$SIND/configs $WORKDIR/$SIND/pdb_inits
    cp $WORKDIR/init.pdb $WORKDIR/$SIND/pdb_inits
    cp $SRC_PP_PATH/$BINPPCLASH $WORKDIR/$SIND
    ss=`cat $SECSTRNAME | head -1`
    cp $TOOLS_PATH/find_loops_from_ss.py .
    python find_loops_from_ss.py "$ss" > $WORKDIR/$SIND/phpull_loops.lst
    rm find_loops_from_ss.py
    
    echo "Setting up the parameters for knot removal..."
    sed -i "s/kr_steps/$PPULL_NMC/g" $WORKDIR/$SIND/input.spqr
    sed -i "s/kr_nt/$PPULL_NSA/g" $WORKDIR/$SIND/input.spqr
    sed -i "s/kr_ti/$PPULL_TSA/g" $WORKDIR/$SIND/input.spqr
    sed -i "s|ENERPTH|$ENERGY_PATH|g" $WORKDIR/$SIND/input.spqr
    
    echo "Running knot removal annealing..."
    cd $WORKDIR/$SIND
    mpirun -np $NPROC ./$BINPPCLASH input.spqr > ppclash_$NAME.out
    echo "Annealing finished!"
    touch temppp ; rm temppp
    for iconf in `seq 0 $MPROC`;do
	IND=`printf %02d $iconf`
	ENER=`head -1 final.p$IND.pdb| awk '{print $3}'`
	echo $IND $ENER >> temppp
    done
    sort -g -k 2 temppp > temp3;
    SELECTED=`head -1 temp3 | awk '{print $1}'` ;
    rm temppp temp3
    cp final.p$SELECTED.pdb ../final.$LKR.pdb
    CRINIT=final.$LKR.pdb
    cd ../..
    echo "Knot removal finished!"
fi

########### STEP 1 : CLASH REMOVAL ################
echo "Creating files for clash removal..."
SIND=$LCR
mkdir $WORKDIR/$SIND
cp files/opt.spqr $WORKDIR/$SIND/input.spqr
mkdir $WORKDIR/$SIND/configs $WORKDIR/$SIND/pdb_inits

cp $WORKDIR/$CRINIT $WORKDIR/$SIND/pdb_inits/init.pdb

cp $SRC_PATH/$BINCLASH $WORKDIR/$SIND
cp $SRC_PATH/$BINANNOT $WORKDIR/$SIND

echo "Setting up the parameters for clashing removal..."
sed -i "s/mc_steps/$CLASH_NMC/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_nt/$CLASH_NSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_ti/$CLASH_TSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s|ENERPTH|$ENERGY_PATH|g" $WORKDIR/$SIND/input.spqr
echo "Running clash removal annealing..."
cd $WORKDIR/$SIND
mpirun -np $NPROC ./$BINCLASH input.spqr > clash_$NAME.out
touch temp1; rm temp1 ; touch temp2 ; rm temp2
for iconf in `seq 0 $MPROC`;do
    IND=`printf %02d $iconf`
    mcconf=`ls configs/chk*p$IND.mc | tail -1`
    ./$BINANNOT $mcconf -t > temp2
    fc=`grep "No clashes found" temp2 | wc | awk '{print $1}'`
    en=`grep "TOTAL ENERGY" temp2 | awk '{print $4}'`
    if [ "$fc" -eq "1" ] ; then echo $mcconf $en >>  temp1 ; fi
done
if [ ! -f temp1 ] ; then echo "There are still clashes in the initial condition. Try increasing the CLASH_NMC variable in the script." ; exit ; fi
NSEL=`wc temp1 | awk '{print $1}'`
sort -g -k 2 temp1 > temp3;
CRSEL=`head -1 temp3 | awk '{print $1}'` ;
rm temp1 temp2 temp3
cp $CRSEL sel.$LCR.mc
cd ../..
echo "Clashes removed successfuly!"


########### STEP 2 : ERMSD PULLING ###################
echo "Structure refinement: minimization of ERMSD with respect to original structure"
SIND=$LEP
mkdir $WORKDIR/$SIND ; cp files/opt.spqr $WORKDIR/$SIND/input.spqr ; mkdir $WORKDIR/$SIND/configs $WORKDIR/$SIND/pdb_inits
sed -i "s/mc_steps/$ERMSD_NMC/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_nt/$ERMSD_NSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_ti/$ERMSD_TSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s|ENERPTH|$ENERGY_PATH|g" $WORKDIR/$SIND/input.spqr
cp $WORKDIR/$LCR/$CRSEL $WORKDIR/$SIND/pdb_inits/init.mc
echo  "REMARK ERMSD PARAMS 1 $ERMSD_PREF $ERMSD_CUTOFF" > ermsd_frags.lst
NATS=`wc $TARGETSTRCT | awk '{print $1}'`; NGROUPS=`echo $NATS | awk '{print int($1/5)-1}'` 
echo "REMARK ERMSD GROUP `seq 0 $NGROUPS | tr "\n" " " `" >> ermsd_frags.lst
grep -v REMARK $TARGETSTRCT >> ermsd_frags.lst ; mv ermsd_frags.lst $WORKDIR/$SIND
cp $SRC_PATH/$BINERMSD $WORKDIR/$SIND
cd $WORKDIR/$SIND

mpirun -np $NPROC ./$BINERMSD input.spqr > ermsd_$NAME.out

cd ../..

echo "ERMSD minimization succesful!"


########### STEP 3 : ENERGY MINIMIZATION #############
echo "Structure refinement: minimization of SPQR energy"
SIND=$LEM
mkdir $WORKDIR/$SIND ; cp files/eopt.spqr $WORKDIR/$SIND/input.spqr ; mkdir $WORKDIR/$SIND/configs $WORKDIR/$SIND/pdb_inits
sed -i "s/mc_steps/$ENERG_NMC/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_nt/$ENERG_NSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s/sa_ti/$ENERG_TSA/g" $WORKDIR/$SIND/input.spqr
sed -i "s|ENERPTH|$ENERGY_PATH|g" $WORKDIR/$SIND/input.spqr
PMAX=$((NPROC - 1))
for p in `seq 0 $PMAX` ; do
    pind=`printf %02d $p`
    SEL_ERMSD_FILE=`ls -t $WORKDIR/$LEP/configs/chk*.p$pind.mc | head -1`
    cp $SEL_ERMSD_FILE $WORKDIR/$SIND/pdb_inits/init.p$pind.mc
done
cp $SRC_PATH/$BINENERG $WORKDIR/$SIND

cd $WORKDIR/$SIND
mpirun -np $NPROC ./$BINENERG input.spqr > energ_$NAME.out

cp final.p* ..
cd ../..
echo "Energy minimization finished!"

########################################################
echo "Optimization ended successfully."


