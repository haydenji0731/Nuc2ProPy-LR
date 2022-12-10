# ONT Read Simulation

We used [Nanosim](https://github.com/bcgsc/NanoSim) to generate in silico ONT reads.
The reads we generated are intended to simulated a typical genomic read for H. sapiens. 
Please follow instruction in the Nanosim repo to replicate the process through which we generated our reads. 
We used the pre-trained **human_NA12878_DNA_FAB49712_guppy** error profile. The command that was used is as follows:

```
simulator.py genome -dna_type linear -rg ./test/GCF_000001405.40_GRCh38.p14_cds_from_genomic_altered.fna -c ./profile/human_NA12878_DNA_FAB49712_guppy/training -b guppy -max 140000 -min 10000 --fastq -n 100
```
