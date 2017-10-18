## SPQR - SPlit and conQueR for RNA structure prediction!

SPQR is a coarse-grained representation and energy function for the prediction of the three-dimensional structure of RNA from the knowledge of its sequence.

A nucleotide is represented as two particles: an anisotropic one which stands for the nucleoside and a point particle for the phosphate group.

<img src="http://people.sissa.it/~spoblete/webpage/multiscale.png" width="500">

# Features

## Interactions
Our energy function makes a geometrically detailed description of canonical and non-canonical base pairs, stacking interactions and even base-phosphate interactions! They are constructed as knowledge-based potentials inspired by the ESCORE function [1].

## Base degrees of freedom
Also, SPQR allows the definition of different glycosidic bond angle conformations: ANTI, HIGH-ANTI and SYN, and it is able to distinguish between sugar pucker conformations: C3' endo and C2' endo.

## Backmapping procedure
Due to the anisotropic representation of the nucleobases, SPQR allows the reintroduction of atomistic details in the predicted structures through steered-Molecular Dynamics [2].

SPQR has successfully been tested on several structures, including duplexes, hairpins and junctions [3]

# The code
SPQR is implemented in a versatile and clear, bug free C code which is available to download for free. The simulations can run constant temperature Monte Carlo simulations, and use Simulated Annealing and Simulated Tempering algorithms for finding the minimum energy conformations. It also includes a set of analysis tools for making the whole process easier.

We include here a small tutorial for the folding of a small tetraloop hairpin, and a brief user manual.






## References
[1] S. Bottaro, F. Di Palma and G. Bussi, _Nucleic Acids Res._ 42(21), 13306-13314 (2014).

[2] S. Poblete, S. Bottaro and G. Bussi, [Link] (https://arxiv.org/abs/1709.08691).

[3] S. Poblete, S. Bottaro and G. Bussi, _submitted_.
You can use the [editor on GitHub](https://github.com/srnas/spqr/edit/master/README.md) to maintain and preview the content for your website in Markdown files.
