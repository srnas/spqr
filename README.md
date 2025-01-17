# SPQR - SPlit and conQueR for RNA structure prediction and refinement

SPQR is a coarse-grained representation and energy function for the prediction of the three-dimensional structure of RNA from the knowledge of its sequence [1].

A nucleotide is represented as two particles: an anisotropic one which stands for the nucleoside and a point particle for the phosphate group.
<p align="center">
<img src="https://github.com/srnas/spqr/blob/master/doc/src/rnacg.png" width="200">
<p>

## Features

### Interactions
Our energy function makes a geometrically detailed description of canonical and non-canonical base pairs, stacking interactions and even base-phosphate interactions! They are constructed as knowledge-based potentials inspired by the ESCORE function [2].

### Base degrees of freedom
Also, SPQR allows the definition of different glycosidic bond angle conformations: ANTI, HIGH-ANTI and SYN, and it is able to distinguish between sugar pucker conformations: C3' endo and C2' endo.

### Backmapping procedure
Due to the anisotropic representation of the nucleobases, SPQR allows the reintroduction of atomistic details in the predicted structures through steered-Molecular Dynamics [3].

SPQR has successfully been tested on several structures, including duplexes, hairpins and junctions [1].

### The code
SPQR is implemented in a versatile and clear, bug free C code which is available to download for free. The simulations can run constant temperature Monte Carlo simulations, and use Simulated Annealing and Simulated Tempering algorithms for finding the minimum energy conformations. It also includes a set of analysis tools for making the whole process easier.

We include here a small tutorial for the folding of a small tetraloop hairpin, and a user manual.

# Getting started

## Installing

After downloading the source, go to the `spqr/src` directory and run
```
autoconf
./configure
make
make install
```

The binary files will be stored in the `spqr/bin` directory and the tools, in the `spqr/tools` file. A more detailed manual is stored in `spqr/doc`.

## Quick energy minimization and scoring
### For optimizing a structure, remove its topological artifacts or generate a suitable initial configuration for a SPQR simulation

A simple minimization can be run on a pdb file containing a single-strand of RNA. Copy the `SPQR_REFINE` script into any directory and run
```
./SPQR_REFINE -i <PDBFILE> [ -o <OUTPUTNAME> -t <SECONDARYSTRUCTURE> ]
```
This script runs internally several simulations which minimize and remove artifacts in order to provide a suitable structure for simulations.
The script will create and store the results in the folder refSPQR<PDBFILE>. Optionally, the name can be specified with the `-o <OUTPUTNAME>` option.
The folder will contain a pdb file in SPQR format (`refSPQR<PDBFILE>.pdb`), an all-atom reconstruction (`at_refSPQR<PDBFILE>.pdb`) and a .mc file (`refSPQR<PDBFILE>.mc`), which can be used for further SPQR simulations with full precision (recommended). The all-atom reconstruction is suitable for MD simulations as well.
The minimization of energy is a good starting point for evaluating the SPQR energy, which is contained in the `.ene` file in the working directory. In addition, it removes clashes and fixes broken bonds, which might be of help when dealing with structures generated by fragment assembly methods, such as the ones generated by Ernwin [4,5].

## Disentanglement
RNA structures generated by a variety of experimental and computational procedures can present artificial entanglements between their secondary structure elements. RNAspider is a tool which allows to quickly identify, visualize and classify such artifacts [6], with a webserver available at [this link] (https://rnaspider.cs.put.poznan.pl). 
With an entangled pdb structure, its secondary structure (.ss file) and the output file of RNAspider (in csv format), SPQR is able to remove the topological artifacts by carefully displacing the precise nucleotides in the right direction. Starting from an all-atom pdb structure, one must run
```
python disentangle_from_RNAspider.py -i $SPIDERFILE.csv -o refined -p $PDBFILE.pdb -t $SEC_STRUCT.ss
```
The output structures will be in the directory `refined`. Its all-atom version, in `refined/AA_refined.pdb` can be improved a bit before an energy minization with Molecular Dynamics, which is highly recommended. By running
```
python SPQR_BMAP_RELAX.py -i AA_refined.pdb -r SPQR_refined.pdb -o AA_better.pdb
```
with the tools provided in the `tools/` directory.
An example of an entangled structure, together with the csv and secondary structure, is found in the directory `tutorials/RNAspider` , ready to be disentangled!

## Visualization
SPQR pdb files can be visualized in VMD. In order to load the topology, a tcl script must be loaded from the console by running
```
source at.tcl
```
which can be found in the `tools/` directory.

## A single run
Three files are required for running a SPQR simulation : an initial condition, a parameters file and the binary to execute. An initial condition can be generated using the `SPQR_ASSEMBLE.py` script in the `tools/` directory or by preparing an initial condition through the minimization procedure aforementioned. They must be stored in a directory named `spqr_inits/` as `init.pdb` or `init.mc`.
The generation of a single strand can be done by running
```
python SPQR_ASSEMBLE.py -s <SEQUENCE>
```
which generates an `init.pdb` file.

Note that the pdb structure must be in SPQR format, which can also be generated using the script `spqr2pdb.py` found in `tools/` as well.

The parameters file must be named `params.spqr`, and specifies the length, saving frequency, energy function and other restraints such as surface or radius of gyration potentials (see documentation) [6]. An example of this file is contained in the `tutorials/single_run/` directory. For a single run, copy the `SPQR_MC` binary from the `bin/` directory and run
```
./SPQR_MC [ -i <INDEX> ]
```
The INDEX provided is a positive integer number, optional, which serves as an identifier for running different simulations with the same initial condition but different random seeds. More initial conditions can be specified with the index by naming them properly in the `spqr_inits/` directory (see the user's guide in the `doc/` directory).

The output files are contained in the `configs/` directory, while a copy in pdb format is created in the working directory under the name `final.p<INDEX>.pdb`.

## Simulated annealing
By properly specifying the parameters in the `params.spqr`, one can run a simulated annealing simulation. The binary file is `SPQR_SA`, and must be run by invoking
```
./SPQR_SA [ -i <INDEX> ]
```

## Secondary structure and three-dimensional restraints
SPQR supports ERMSD restraints [3,7] for imposing secondary structure restraints or imposing three-dimensional templates on the simulation. These tools are extremely powerful, and demand the additional file `ermsd_frags.dat`. Such a file can be created from the knowledge of a secondary structure or a template all-atom pdb. For details, see the documentation.






## References
[1] S. Poblete, S. Bottaro and G. Bussi, A nucleobase-centric coarse-grained model for structure prediction of RNA motifs. _Nucleic Acids Res._ 46(4), 1674-1683 (2018).  [Link] https://doi.org/10.1093/nar/gkx1269 

[2] S. Bottaro, F. Di Palma and G. Bussi, The role of nucleobase interactions in RNA structure and dynamics. _Nucleic Acids Res._ 42(21), 13306-13314 (2014).

[3] S. Poblete, S. Bottaro and G. Bussi, Effects and limitations of a nucleobase-driven backmapping procedure for nucleic acids using steered molecular dynamics. _Biophys. Biochem. Res. Comm._ 498 (2), 352-358 (2018).

[4] P. Kerpedjiev, C. Höner zu Siederdissen and I. L. Hofacker, Predicting RNA 3D structure using a coarse-grain helix-centered model. _RNA_ 21(6):1110-1121, (2015)

[5] B. Thiel, G. Bussi, S. Poblete and I. L. Hofacker. _Nucleic Acids Res._ 52(16), e73 (2024) [Link] https://doi.org/10.1093/nar/gkae602

[6] S. Poblete, A. Bozic, M. Kanduc, R. Podgornik and H. V. Guzman, RNA Secondary Structures Regulate Adsorption of Fragments onto Flat Substrates. _ACS Omega_ 6(48):32823-32831 (2021).

[7] S. Poblete and H. V. Guzman, Structural 3D Domain Reconstruction of the RNA Genome from Viruses with Secondary Structure Models. _Viruses_ 13(8):1555 (2021).
