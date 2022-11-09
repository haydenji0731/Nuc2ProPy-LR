# import packages
import pyfastx
import sys
import time

num_seq = 0

class db_index:
    def __init__(self, kmer_size):
        self.kmer_table = dict()
        self.kmer_size = kmer_size

    def insert(self, protein):
        name = protein[0]
        seq = protein[1]
        for i in range(0, len(seq) - self.kmer_size + 1):
            kmer = seq[i:i + 5]
            if kmer in self.kmer_table.keys():
                self.kmer_table[kmer].add(name)
            else:
                name_set = {name}
                self.kmer_table[kmer] = name_set


def main(input_file, kmer_size, out_dir, out_file):
    start = time.time()
    proteins = pyfastx.Fastx(input_file)
    prot_db = db_index(kmer_size)
    global num_seq
    for protein in proteins:
        num_seq += 1
        prot_db.insert(protein)
    duration = time.time() - start
    print("Processed %.0f sequences in %.2fs" % (num_seq, duration))


if __name__ == '__main__':
    # TODO: replace this with flags
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
