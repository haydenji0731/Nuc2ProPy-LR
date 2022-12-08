# Benchmark Results Explained

This directory contains the lonblast outputted alignment results for **nanosim**-generated reads. For this project, we generated 123 reads. 
Of these reads, 20 aren't "unaligned," indicating they shouldn't find a match in the protein database and the rest (103) are "aligned."
<br/><br/>
To test on a local machine, the read file was truncated into 4 chunks (30, 30, 30, 33). 
The alignment results that correspond to each of these chunks are named "simulated_reads_1_2_combined.nXX_X.aln" such that X is some number.
Hence, for evaluation, you will need to concat all these alignment results **in order** (e.g., 1 before 2, 2 before 3, etc.)
<br/><br/>
Under this directory, the **blastx** directory contains the alignment results generated using the same set of **nanosim**-generated reads. 
For evaluation purposes, please use:
```
simulated_reads_1_2_combined_custom_processed.aln
```
This alignment file is a further processed version of the default blastx output that comes in the two-column format (query, target). 
This way, for evaluation purposes, the blastx results can be compared to the ground truth alignment results comes in this format. 
<br/><br/>
The ground truth alignment results is located here:
```
cs647_proj_group62/ontSimulation/reads/simulated_reads_1_2_combined.aln
```
