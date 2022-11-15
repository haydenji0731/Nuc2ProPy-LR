# import packages
import pyfastx
import sys
import time


class db_index:
    # overload the constructor
    def __init__(self, kmer_size, proteins=None):
        self.seq_index = dict()
        self.kmer_table = dict()
        self.kmer_size = kmer_size
        self.num_seq = 0
        if proteins:
            idx = 0
            for protein in proteins:
                # seq, seq info
                self.seq_index[idx] = (protein[1], protein[2])
                self.insert(protein[1], idx)
                idx += 1

    def insert(self, seq, idx):
        for i in range(0, len(seq) - self.kmer_size + 1):
            kmer = seq[i:i + self.kmer_size]
            if self.contains(kmer):
                self.kmer_table[kmer].add((idx, i))
            else:
                name_set = {(idx, i)}
                self.kmer_table[kmer] = name_set
        self.num_seq += 1

    def contains(self, kmer):
        if kmer in self.kmer_table.keys():
            return True
        return False


def main(input_file, kmer_size, out_dir, out_file):
    start = time.time()
    proteins = pyfastx.Fastx(input_file)
    prot_db = db_index(kmer_size, proteins)
    duration = time.time() - start
    print("Processed %.0f sequences in %.2fs" % (prot_db.num_seq, duration))


if __name__ == '__main__':
    # TODO: replace this with flags
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
