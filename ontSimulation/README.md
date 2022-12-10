# ONT Read Simulation

We used [Nanosim](https://github.com/bcgsc/NanoSim) to generate in silico ONT reads.
The reads we generated are intended to simulated a typical genomic read for H. sapiens. 
Please follow instruction in the Nanosim repo to replicate the process through which we generated our reads. In total, we generated 123 in silico ONT reads, consiting of 103 aligned reads and 20 unaligned reads. We used the reference fasta that only contained CDS regions of the human genome. The reference fasta can be obtained [here](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/110/GCF_000001405.40_GRCh38.p14/).
We used the pre-trained **human_NA12878_DNA_FAB49712_guppy** error profile. The command that was used is as follows:

```
$ simulator.py genome -dna_type linear -rg ./test/GCF_000001405.40_GRCh38.p14_cds_from_genomic_altered.fna -c ./profile/human_NA12878_DNA_FAB49712_guppy/training -b guppy -max 140000 -min 10000 -n 100
```
Additionally, for benchmarking speed of the aligner with a single read, we extract the very first read in our set of reads by running:
```
$ head -n 2 ./reads/simulated_reads_1_2_combined.fasta
```
We also annotated a true label to each aligned read and the unaligned reads got assigned with '.' instead. 
This ground truth alignment results file can be found here:
```
./reads/simulated_reads_1_2_combined.aln
```
This directory also contains helper python scripts that were used to annotate each read with its true label. 

1. alterFasta.py --> alter the reference fasta for the identifier starts with a gene name so it can later be used to identify which protein the read should align to. By default, Nanosim takes a prefix of the identifier in the reference fasta so this step was necessary
2. findProtByGene.py --> find the protein that correspond to a gene
3. processUnaligned.py --> assign '.' to each unaligned read
