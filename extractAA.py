# import packages
import sys
import translationTable
import pyfastx
import time

num_seq = 0
templ = "ACGT"
rev_compl = "TGCA"


def extract_aa(reads, trans_table):

    aa_dict = {}
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
            aa_seqs.append((aa_seq, i+1))
            aa_seqs.append((aa_seq_rev, -(i+1)))
        aa_dict[read[0]] = aa_seqs
    return aa_seqs


def main(input_file):
    start = time.time()
    reads = pyfastx.Fastx(input_file)
    aa_seqs = extract_aa(reads, translationTable.trans_table)
    # print(aa_seqs)
    duration = time.time() - start
    print("Processed %.0f sequences in %.4fs" % (num_seq, duration))



if __name__ == '__main__':
    main(sys.argv[1])