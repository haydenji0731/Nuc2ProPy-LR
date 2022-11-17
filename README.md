<<<<<<< HEAD
# Protein Search Project

We aim to convert DNA sequence reads into amino acid sequences and then search against large-scale protein databases.
=======
# Long DNA Read to Protein Aligner

We aim to convert DNA sequence reads into amino acid sequences and then search against large-scale protein databases.

**11/16/22** - Please note that the commit history might look a little bit odd. [haydenji0731](https://github.com/haydenji0731) had previously committed a substantial number of times with the wrong username and email, and thus these commits were fixed and pushed again.

### How to Run ###

We recommend running the program with a local python interpreter. Create a python venv and install the following packages. 

```
pyfastx, numpy, pyinstaller
```

Then, to build the program into an executable, run:

```
$ pyinstaller main.py --onefile -n lonblast
$ cd dist
```

Finally, in your command line prompt, run:

```
$ ./lonblast -h
```
Proceed by reading the help message. Additionally, there are some benchmark datasets available in the data directory. An examplar run command is:

```
$ ./lonblast --out-dir ../results --out-file test_flanked.aln --aa-file refseq_grch38_unflanked_cds_head_AA.n16.aa.fa ../data/refseq_grch38_flanked_cds_head.n16.fa ../data/GCF_000001405.40_GRCh38.p14_protein.fasta
```


>>>>>>> 3a291efe2cc427a6e05183c86709d95762995408
