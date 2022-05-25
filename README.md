# phyBWT

phyBWT is a new alignment-, assembly-, and reference-free method that builds a partition tree without relying on the pairwise comparison of sequences, thus avoiding to use a distance matrix to infer phylogeny.

It applies the properties of the Extended Burrows-Wheeler Transform (EBWT) to the idea of decomposition for phylogenetic inference. 
In particular, it hinges the combinatorial properties of the *EBWT positional clustering* framework recently introduced, overcoming the limitations of employing *k*-mers with a priori fixed size. 
Finally, phyBWT infers the tree structure by comparing all the sequences simultaneously, instead of performing their pairwise comparisons.
