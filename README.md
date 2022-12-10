# Long DNA Read to Protein Aligner

We aim to convert DNA sequence reads into amino acid sequences and then search against large-scale protein databases.

**11/16/22** - Please note that the commit history might look a little bit odd. [haydenji0731](https://github.com/haydenji0731) had previously committed a substantial number of times with the wrong username and email, and thus these commits were fixed and pushed again.

### How to Run ###

We recommend running the program with a local python interpreter. These are the required packatges to run this tool. 

```
pyfastx, numpy, pyinstaller
```
You can install them manually or simply run:
```
pip install -r requirements.txt
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
$ ./lonblast --out-dir ../results --out-file simulated_reads_1_2_combined.aln --aa-file simulated_reads_1_2_combined.aa.fa ../ontSimulation/reads/simulated_reads_1_2_combined.fasta ../data/MANE.GRCh38.v1.0.refseq_protein.faa --kmer-size 9 --gap-open 11 --gap-extend 1
```
