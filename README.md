## SPQR - SPlit and conQueR for RNA structure prediction!

SPQR is a coarse-grained representation and energy function for the prediction of the three-dimensional structure of RNA from the knowledge of its sequence [1].

A nucleotide is represented as two particles: an anisotropic one which stands for the nucleoside and a point particle for the phosphate group.

<img src="http://people.sissa.it/~spoblete/webpage/multiscale.png" width="500">

# Features

## Interactions
Our energy function makes a geometrically detailed description of canonical and non-canonical base pairs, stacking interactions and even base-phosphate interactions! They are constructed as knowledge-based potentials inspired by the ESCORE function [2].

## Base degrees of freedom
Also, SPQR allows the definition of different glycosidic bond angle conformations: ANTI, HIGH-ANTI and SYN, and it is able to distinguish between sugar pucker conformations: C3' endo and C2' endo.

## Backmapping procedure
Due to the anisotropic representation of the nucleobases, SPQR allows the reintroduction of atomistic details in the predicted structures through steered-Molecular Dynamics [3].

SPQR has successfully been tested on several structures, including duplexes, hairpins and junctions [1].

# The code
SPQR is implemented in a versatile and clear, bug free C code which is available to download for free. The simulations can run constant temperature Monte Carlo simulations, and use Simulated Annealing and Simulated Tempering algorithms for finding the minimum energy conformations. It also includes a set of analysis tools for making the whole process easier.

We include here a small tutorial for the folding of a small tetraloop hairpin, and a user manual.






## References
[1] S. Poblete, S. Bottaro and G. Bussi, A nucleobase-centric coarse-grained model for structure prediction of RNA motifs. _Nucleic Acids Res._ 46(4), 1674-1683 (2018).  [Link] https://doi.org/10.1093/nar/gkx1269 

[2] S. Bottaro, F. Di Palma and G. Bussi, The role of nucleobase interactions in RNA structure and dynamics. _Nucleic Acids Res._ 42(21), 13306-13314 (2014).

[3] S. Poblete, S. Bottaro and G. Bussi, Effects and limitations of a nucleobase-driven backmapping procedure for nucleic acids using steered molecular dynamics. _Biophys. Biochem. Res. Comm._ 498 (2), 352-358 (2018).
