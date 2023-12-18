from __future__ import print_function
import forgi
import forgi.threedee.classification.aminor as ftca
import sys
rnas = forgi.load_rna(sys.argv[1])
for rna in rnas:
    #print ("DSSR", rna.name, list(rna.dssr.aminor_interactions()))
    print ("forgi", rna.name, list(ftca.all_interactions(rna)))
#    print ("DSSR", rna.dssr._dssr)
