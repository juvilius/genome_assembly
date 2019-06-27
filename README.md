# De Bruijn Genome Assembly

Creates random reads from a genome sequence and tries to put it back together.

## Getting Started

```
git clone https://github.com/molnxx/genome_assembly.git
```

### Prerequisites

You need numpy and a preferably circular DNA sequence.

## Usage

###Create reads

Create random reads with a coverage **n** and variable length of around **l**:

```
python de_bruijn_assembly.py create -i genome.fasta -n 5 -l 200 -o outfile
```

*outfile* determines the fasta file with all the reads.

### Assemble genome

Assemble the genome sequence with de Bruijn graphs:

```
python de_bruijn_assembly.py assemble -i reads.fasta -k 43 -p orig_genome.fasta -o result_file
```

*k* is the length of the kmers. *orig_genome.fasta* is the genome sequence to proof the results (which are stored in *result_file*).

##Examples

Included is the genome of pUC19, a well-known vector. Also the example from the paper referenced below.

## License

Public Domain

## Acknowledgments

* Phillip E C Compeau et al. for *How to apply de Bruijn graphs to genome assembly* [doi](https://doi.org/10.1038/nbt.2023)
* Atom Community (https://github.com/atom/atom)
* NumPy devs (https://www.numpy.org/)
