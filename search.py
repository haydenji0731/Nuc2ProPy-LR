# import packages
import findOrf
import buildIndex
import extractAA
import align
import time
import pyfastx
import translationTable
import sys
# supports BLOSUM 45, 50, 62, 80, & 90
import blosum as bl

# define global variables
aln_dict = {}
seq_dict = {}


def main(query_file, target_file, out_dir, out_file, aa_file, extract_orf, kmer_size, orf_file=None):
    print("successfully entered search main!")
    start = time.time()
    reads = pyfastx.Fastx(query_file)
    if extract_orf:
        findOrf.find_orfs(reads, translationTable.trans_table)
        print(len(findOrf.orf_dict.keys()))
        findOrf.write_output_fasta(out_dir, out_file)
        duration = time.time() - start
        print("Processed %.0f sequences in %.4fs. Of all, %.0f sequences did not contain a valid orf."
              % (findOrf.num_seq, duration, findOrf.invalid_seq))
    else:
        extractAA.extract_aa(reads, translationTable.trans_table)
        extractAA.write_output_fasta(out_dir, aa_file)
        duration = time.time() - start
        print("Processed %.0f sequences in %.4fs." % (extractAA.num_seq, duration))
    proteins = pyfastx.Fastx(target_file)
    prot_db = buildIndex.db_index(kmer_size, proteins)
    sub_mat = bl.BLOSUM(62)
    global aln_dict
    global seq_dict
    for read_name in extractAA.aa_dict.keys():
        aa_seqs = extractAA.aa_dict[read_name]
        max_rf = 0
        global_max = -1
        global_aln = ()
        global_mat = None
        # align all 6 reading frames
        for query in aa_seqs:
            rf = query[1]
            seq = query[0]
            if len(seq) < kmer_size:
                print("Query amino acid sequence is too short compared to the kmer size. No significant sequence "
                      "similarity is expected.")
            else:
                query_pref = extractAA.kmer_prefilter(seq, prot_db, kmer_size)
                if len(query_pref) > 0:
                    local_max = -1
                    local_aln = ()
                    local_mat = None
                    for target_idx in query_pref:
                        sw_mat, aln_score = extractAA.sw_align(seq, prot_db.seq_index[target_idx][0], sub_mat)
                        if aln_score > local_max:
                            local_max = aln_score
                            local_mat = sw_mat
                            # target name, score
                            local_aln = (prot_db.seq_index[target_idx][1], aln_score)
                    if local_max > global_max:
                        global_max = local_max
                        global_mat = local_mat
                        global_aln = local_aln
                        max_rf = rf
                else:
                    print("No matching kmer was found between query and target in reading frame %.0f" % rf)
        if global_max != -1:
            aln_dict[read_name] = global_aln
            print("Optimal alignment for read " + read_name + " was found in reading frame %.0f" % max_rf)
        else:
            print("No alignment was found for read " + read_name)















