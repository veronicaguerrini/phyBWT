# phyBWT

phyBWT is a new alignment-, assembly-, and reference-free method that builds a partition tree without relying on the pairwise comparison of sequences, thus avoiding to use a distance matrix to infer phylogeny.

It applies the properties of the Extended Burrows-Wheeler Transform (EBWT) to the idea of decomposition for phylogenetic inference. 
In particular, it hinges the combinatorial properties of the *EBWT positional clustering* framework recently introduced , overcoming the limitations of employing *k*-mers with a priori fixed size. 
Finally, phyBWT infers the tree structure by comparing all the sequences simultaneously, instead of performing their pairwise comparisons.

### Install

```sh
git clone https://github.com/veronicaguerrini/phyBWT
cd phyBWT
```

### Compile
phyBWT can handle datasets of different types (short reads, long reads, contigs, or entire genomes). 

For short reads, phyBWT must be compiled by using

```sh
make
```
while for sequences longer than 250 bp, set SHORT=0

```sh
make SHORT=0
```

### Quick test

After compiling with parameter SHORT=0

```sh
./phyBWT prasinovirus/prasinoviruses.fasta prasinovirus/prasinoviruses.txt prasino.out 12 0.5 6
```

## References

Guerrini V., Conte A., Grossi R., Liti G., Rosone G., Tattini L., phyBWT: Alignment-Free Phylogeny via eBWT Positional Clustering. *Submitted paper*.
