# import packages
import findOrf
import buildIndex
import extractAA
import align
import time
import pyfastx
import translationTable
import sys
import os
# supports BLOSUM 45, 50, 62, 80, & 90
import blosum as bl

# define global variables
aln_dict = {}
# seq_dict = {}
num_seq = 0
no_aln_seq = 0


def write_aln(out_dir, out_file):
    global aln_dict
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fa_path = os.path.join(out_dir, out_file)
    with open(out_fa_path, 'w') as fh:
        for query in aln_dict.keys():
            aln = aln_dict[query]
            fh.write(query + "\t" + aln[0] + "\t" + aln[1] + "\t" + str(aln[2]) + "\n")


def main(query_file, target_file, out_dir, out_file, aa_file, kmer_size,
         gap_open, gap_extend, extract_orf, orf_file=None):
    print("successfully entered search main!")
    start = time.time()
    reads = pyfastx.Fastx(query_file)
    if extract_orf:
        findOrf.find_orfs(reads, translationTable.trans_table)
        print(len(findOrf.orf_dict.keys()))
        findOrf.write_output_fasta(out_dir, orf_file)
        duration = time.time() - start
        print("Extracted orfs from %.0f reads in %.4fs. Of all, %.0f reads did not contain a valid orf."
              % (findOrf.num_seq, duration, findOrf.invalid_seq))
    else:
        extractAA.extract_aa(reads, translationTable.trans_table)
        extractAA.write_output_fasta(out_dir, aa_file)
        duration = time.time() - start
        print("Translated %.0f reads to amino acid sequences in %.4fs." % (extractAA.num_seq, duration))
        start = time.time()
    proteins = pyfastx.Fastx(target_file)
    prot_db = buildIndex.db_index(kmer_size, proteins)
    duration = time.time() - start
    print("Indexed target db containing %.0f protein sequences in %.4fs." % (prot_db.num_seq, duration))
    start = time.time()
    sub_mat = bl.BLOSUM(62)
    global aln_dict
    global num_seq
    global no_aln_seq
    # global seq_dict
    # set up output file
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fa_path = os.path.join(out_dir, out_file)
    fh = open(out_fa_path, 'w')
    for read_name in extractAA.aa_dict.keys():
        num_seq += 1
        aa_seqs = extractAA.aa_dict[read_name]
        max_rf = 0
        global_max = -1
        global_aln = ()
        # report top 5 scoring alignments per read
        global_max_arr = []
        # global_mat = None
        # align all 6 reading frames
        for query in aa_seqs:
            rf = query[1]
            print("Processing reading frame %.0f" % rf)
            seq = query[0]
            # print(seq)
            if len(seq) < kmer_size:
                print("Query amino acid sequence is too short compared to the kmer size. No significant sequence "
                      "similarity is expected.")
            else:
                query_pref = align.kmer_prefilter(seq, prot_db, kmer_size)
                if len(query_pref) > 0:
                    # local_max = -1
                    # local_aln = ()
                    # local_mat = None
                    for target_idx in query_pref:
                        # sw_mat, aln_score = align.sw_align(seq, prot_db.seq_index[target_idx][1], sub_mat)
                        sw_mat, aln_score = align.sw_align_affine(seq, prot_db.seq_index[target_idx][1],
                                                                  gap_open, gap_extend, sub_mat)
                        aln = (aln_score, prot_db.seq_index[target_idx][0], prot_db.seq_index[target_idx][2], rf)
                        if len(global_max_arr) < 5:
                            global_max_arr.append(aln)
                            global_max_arr.sort(reverse=True)
                        else:
                            # print(global_max_arr[-1:][0][0])
                            if global_max_arr[-1:][0][0] < aln_score:
                                global_max_arr.append(aln)
                                global_max_arr.sort(reverse=True)
                                global_max_arr = global_max_arr[:5]
                        # print(global_max_arr)
                        # default
                        # if aln_score > local_max:
                        #     local_max = aln_score
                        #     # local_mat = sw_mat
                        #     # protein ID, protein name, aln score
                        #     local_aln = (prot_db.seq_index[target_idx][0], prot_db.seq_index[target_idx][2], aln_score)

                    # default
                    # if local_max > global_max:
                    #     global_max = local_max
                    #     # global_mat = local_mat
                    #     global_aln = local_aln
                    #     max_rf = rf
                else:
                    print("No matching kmer was found between query and target in reading frame %.0f" % rf)
        if len(global_max_arr) != 0:
            print("Optimal alignment:")
            print("\tProtein ID is: " + global_max_arr[0][1])
            print("\tProtein product is: " + global_max_arr[0][2])
            print("\tAlignment score is: " + str(global_max_arr[0][0]))
            print("Optimal alignment for read " + read_name + " was found in reading frame %.0f" % global_max_arr[0][3])
            for i in range(len(global_max_arr)):
                fh.write(read_name + "\t" + global_max_arr[i][1] + "\t" + global_max_arr[i][2] + "\t" +
                         str(global_max_arr[i][0]) + "\t" + str(global_max_arr[i][3]) + "\n")
        else:
            no_aln_seq += 1
            print("No alignment was found for read " + read_name)

    #     if global_max != -1:
    #         aln_dict[read_name] = global_aln
    #         print("Optimal alignment:")
    #         print("\tProtein ID is: " + global_aln[0])
    #         print("\tProtein product is: " + global_aln[1])
    #         print("\tAlignment score is: " + str(global_aln[2]))
    #         print("Optimal alignment for read " + read_name + " was found in reading frame %.0f" % max_rf)
    #     else:
    #         no_aln_seq += 1
    #         print("No alignment was found for read " + read_name)
    # write_aln(out_dir, out_file)
    fh.close()
    duration = time.time() - start
    print("Aligned %.0f reads / sequences in %.4fs. Of all, %.0f reads / sequences did not find an alignment."
          % (num_seq, duration, no_aln_seq))














