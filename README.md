# phyBWT

phyBWT is a new alignment-, assembly-, and reference-free method that builds a partition tree without relying on the pairwise comparison of sequences, thus avoiding to use a distance matrix to infer phylogeny.

It applies the properties of the Extended Burrows-Wheeler Transform (EBWT) to the idea of decomposition for phylogenetic inference. 
In particular, it hinges the combinatorial properties of the *EBWT positional clustering* framework recently introduced , overcoming the limitations of employing *k*-mers with a priori fixed size. 
Finally, phyBWT infers the tree structure by comparing all the sequences simultaneously, instead of performing their pairwise comparisons.

Let *S={S<sub>1</sub>,...,S<sub>n</sub>}* be the input collection of sequences, where each *S<sub>i</sub>* is a multiset of strings representing an organism (e.g. sequencing reads, contigs, genome). The tool phyBWT takes as input the following data structures:
- the extended Burrows–Wheeler transform (ebwt), or multi-string BWT, of collection *S*;
- the longest common prefix array (lcp) of collection *S*;
- the color document array (cda) of collection *S*.

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

### Preprocessing steps

The required data structures eBWT, LCP and DA can be built independently from phyBWT. 
This is a good feature that allows the user to choose the most appropriate tool according to the resources available and the dataset composition (short reads or longer sequences).

For instance, to build .ebwt, .lcp, and .da files from scratch from a single fasta file, one could use BCR [https://github.com/giovannarosone/BCR_LCP_GSA] for short reads and for longer sequences, gsufsort [https://github.com/felipelouza/gsufsort]. Note that gsufsort tool returns the output files with slightly different filename extensions.

To install and compile BCR and gsufsort for the preprocessing, as well as phyBWT one could run

```sh
Install.sh
```
Note that by default the above script compile phyBWt for short reads. To correctly compile phyBWT for longer sequences, the parameter SHORT inside the script must be set to 0.

To obtain the color document array (CDA) from the DA file (fastaFile.da), one could use

```sh
./create_cda fastaFile fileInfo
```
where fileInfo describes the number of sequences in each multiset *S<sub>i</sub>* of the collection *S*.


### Quick test

After compiling with parameter SHORT=0

```sh
./phyBWT prasinovirus/prasinoviruses.fasta prasinovirus/prasinoviruses.txt prasino.out 13 0.5 7
```

## References

Guerrini V., Conte A., Grossi R., Liti G., Rosone G., Tattini L., phyBWT: Alignment-Free Phylogeny via eBWT Positional Clustering. *Submitted paper*.
