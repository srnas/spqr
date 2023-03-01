# SPQR - SPlit and conQueR for RNA structure prediction

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

## Quick energy minimization

A simple minimization can be run on a pdb file containing a single-strand of RNA. Copy the `SPQR_REFINE` script into any directory and run
```
./SPQR_REFINE -i <PDBFILE> [ -o <OUTPUTNAME> -t <SECONDARYSTRUCTURE> ]
```
This script runs internally several simulations which minimize and remove artifacts in order to provide a suitable structure for simulations.
The script will create and store the results in the folder refSPQR<PDBFILE>. Optionally, the name can be specified with the `-o <OUTPUTNAME>` option.
The folder will contain a pdb file in SPQR format, an all-atom reconstruction and a .mc file, which can be used for further SPQR simulations with full precision (recommended). The all-atom reconstruction is suitable for MD simulations as well.
The minimization of energy is a good starting point for evaluating the SPQR energy, which is contained in the `.ene` file in the working directory. In addition, it removes clashes and fixes broken bonds, which might be of help when dealing with structures generated by fragment assembly methods, such as the ones generated by Ernwin [4].
If the secondary structure elements such as hairpins, stems and internal loops are entangled, these artifacts can also be removed automatically [5]. In this case, the secondary structure must be provided in a fasta file containing three lines : a comment line, the sequence and the secondary structure in Vienna format, without pseudoknots. The option `-t <SECONDARYSTRUCTURE>` takes care of this.
Other options can be found in the user's guide and the help message displayed with the `-h` flag.

## A single run
Three files are required for running a SPQR simulation : an initial condition, a parameters file and the binary to execute. An initial condition can be generated using the `SPQR_ASSEMBLE.py` script in the `tools/` directory or by preparing an initial condition through the minimization procedure aforementioned. They must be stored in a directory named `spqr_inits/` as `init.pdb` or `init.mc`. Note that the pdb structure must be in SPQR format, which can also be generated using the script `spqr2pdb.py` found in `tools/`. The parameters file must be named `params.spqr`, and specifies the length, saving frequency, energy function and other restraints such as surface or radius of gyration potentials (see documentation) [6]. An example of this file is contained in the `examples/` directory. For a single run, copy the `SPQR_MC` binary from the `bin/` directory and run
```
./SPQR_MC [ -i <INDEX> ]
```
The INDEX provided is a positive integer number, optional, which serves as an identifier for running different simulations with the same initial condition but different random seeds.
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

[5] B. Thiel, G. Bussi, S. Poblete and I. L. Hofacker, https://www.biorxiv.org/content/10.1101/2022.07.02.498583v1

[6] S. Poblete, A. Bozic, M. Kanduc, R. Podgornik and H. V. Guzman, RNA Secondary Structures Regulate Adsorption of Fragments onto Flat Substrates. _ACS Omega_ 6(48):32823-32831 (2021).

[7] S. Poblete and H. V. Guzman, Structural 3D Domain Reconstruction of the RNA Genome from Viruses with Secondary Structure Models. _Viruses_ 13(8):1555 (2021).
