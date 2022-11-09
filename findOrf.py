# import packages
import pyfastx
import sys
import os
import translationTable
# TODO: report duration of orf search
import time

num_seq = 0
templ = "ACGT"
rev_compl = "TGCA"


def insert_start(codon_dict, rf, start_idx):
    if rf in codon_dict.keys():
        codon_dict[rf]["start"].append(start_idx)
    else:
        codon_dict[rf] = {"start": [start_idx],
                          "stop": []}


def insert_stop(codon_dict, rf, stop_idx):
    if rf in codon_dict.keys():
        codon_dict[rf]["stop"].append(stop_idx)
    else:
        codon_dict[rf] = {"start": [],
                          "stop": [stop_idx]}


def find_orfs(reads, trans_table, orf_dict):
    for read in reads:
        global num_seq
        num_seq += 1
        codon_dict = {}
        seq = read[1]
        rev_compl_trans = str.maketrans(templ, rev_compl)
        seq_rev_compl = seq.translate(rev_compl_trans)[::-1]
        for i in range(len(seq) - 2):
            codon = seq[i:i+3]
            codon_rev_compl = seq_rev_compl[i:i+3]
            rf = i % 3 + 1
            rf_rev_compl = -(i % 3 + 1)
            if trans_table[0][codon][1] == 1:
                insert_start(codon_dict, rf, i)
            if trans_table[0][codon_rev_compl][1] == 1:
                insert_start(codon_dict, rf_rev_compl, i)
            if trans_table[0][codon][1] == 2:
                insert_stop(codon_dict, rf, i)
            if trans_table[0][codon_rev_compl][1] == 2:
                insert_stop(codon_dict, rf, i)
        longest_orf = find_longest_orf(codon_dict)
        if longest_orf["len"] == 0:
            # allow partial 3' as well?
            print(read[0] + " doesn't contain a valid orf!")
        else:
            aa_seq = extract_aa_seq(longest_orf, seq, seq_rev_compl, trans_table)
            orf_dict[read[0]] = {"start": longest_orf["start"],
                                 "stop": longest_orf["stop"],
                                 "len": longest_orf["len"],
                                 "seq": aa_seq}


def find_longest_orf(codon_dict):
    max_orf = {"start": -1, "stop": -1, "len": 0, "frame": 0}
    for i in range(-3, 4):
        if i != 0 and i in codon_dict.keys():
            start_l = codon_dict[i]["start"]
            stop_l = codon_dict[i]["stop"]
            if len(start_l) != 0 and len(stop_l) != 0:
                for start in start_l:
                    for stop in stop_l:
                        if stop > start:
                            orf_len = stop - start
                            if orf_len > max_orf["len"]:
                                max_orf["start"] = start
                                max_orf["stop"] = stop
                                max_orf["len"] = orf_len
                                max_orf["frame"] = i
                            break
    return max_orf


def extract_aa_seq(orf, seq, seq_rev_compl, trans_table):
    aa_seq = ""
    if orf["frame"] >= 0:
        for i in range(orf["start"], orf["stop"] - 2, 3):
            codon = seq[i:i+3]
            aa_seq += trans_table[0][codon][0]
    else:
        for i in range(orf["start"], orf["stop"] - 2, 3):
            codon = seq_rev_compl[i:i+3]
            aa_seq += trans_table[0][codon][0]
    return aa_seq


def write_output_fasta(out_dir, fa_name, orf_dict):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fa_path = os.path.join(out_dir, fa_name)
    with open(out_fa_path, 'w') as fh:
        for read_name in orf_dict.keys():
            orf = orf_dict[read_name]
            fh.write(">" + read_name + " " + str(orf["start"]) + " " + str(orf["stop"]) + " " + str(orf["len"]) + "\n")
            fh.write(orf["seq"] + "\n")


def main(input_file, out_dir, out_file):
    start = time.time()
    reads = pyfastx.Fastx(input_file)
    orf_dict = {}
    find_orfs(reads, translationTable.trans_table, orf_dict)
    write_output_fasta(out_dir, out_file, orf_dict)
    duration = time.time() - start
    print("Processed %.0f sequences in %.4fs" % (num_seq, duration))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
