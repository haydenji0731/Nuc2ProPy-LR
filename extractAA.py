# import packages
import sys
import translationTable
import pyfastx
import time
import os

num_seq = 0
aa_dict = {}
templ = "ACGT"
rev_compl = "TGCA"


def extract_aa(reads, trans_table):
    global aa_dict
    for read in reads:
        global num_seq
        num_seq += 1
        seq = read[1].upper()
        rev_trans = str.maketrans(templ, rev_compl)
        seq_rev = seq.translate(rev_trans)[::-1]
        aa_seqs = []
        for i in range(0, 3):
            aa_seq = ""
            aa_seq_rev = ""
            for j in range(0, len(seq), 3):
                if i + j + 2 <= len(seq) - 1:
                    codon = seq[i + j:i + j + 3]
                    codon_rev = seq_rev[i + j:i + j + 3]
                    # check if stop codon
                    aa_seq += trans_table[0][codon][0]
                    aa_seq_rev += trans_table[0][codon_rev][0]
                else:
                    break
            aa_seqs.append((aa_seq, i + 1))
            aa_seqs.append((aa_seq_rev, -(i + 1)))
        aa_dict[read[0]] = aa_seqs
    # return aa_seqs


def write_output_fasta(out_dir, fa_name):
    global aa_dict
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fa_path = os.path.join(out_dir, fa_name)
    # print(out_fa_path)
    with open(out_fa_path, 'w') as fh:
        for read_name in aa_dict.keys():
            aa_seqs = aa_dict[read_name]
            for seq in aa_seqs:
                fh.write(">" + read_name + "\n")
                fh.write(seq[0] + "\n")


def main(input_file, out_dir, out_file):
    start = time.time()
    reads = pyfastx.Fastx(input_file)
    extract_aa(reads, translationTable.trans_table)
    write_output_fasta(out_dir, out_file)
    duration = time.time() - start
    print("Processed %.0f sequences in %.4fs" % (num_seq, duration))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
