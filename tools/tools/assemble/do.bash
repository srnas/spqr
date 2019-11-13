SEQ=$1
SST=$2
python gen_src/fullassemble.py -s $SEQ -t $SST > init.pdb
while read p ; do
python gen_src/assemble.py -d $p >> ermsd_frags.lst    
done < gen_src/ss.temp
