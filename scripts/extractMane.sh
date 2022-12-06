# output gtf only contains start_codon entry (one per gene)
awk '$3 == "start_codon" { print $0 }' GCF_000001405.40_GRCh38.p14_genomic.gtf | grep "tag \"MANE Select\"" > GCF_000001405.40_GRCh38.p14_genomic_mane_select.gtf


