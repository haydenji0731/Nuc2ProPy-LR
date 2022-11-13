# import packages
import pyfastx
import sys
import time


class db_index:
    # overload the constructor
    def __init__(self, kmer_size, proteins=None):
        self.kmer_table = dict()
        self.kmer_size = kmer_size
        self.num_seq = len(proteins)
        if proteins:
            for protein in proteins:
                print(protein.seq)
                self.insert(protein)

    def insert(self, protein):
        seq = protein.seq
        for i in range(0, len(protein.seq) - self.kmer_size + 1):
            kmer = protein.seq[i:i + self.kmer_size]
            if self.contains(kmer):
                # target seq identifier, pos
                self.kmer_table[kmer].add((protein.id, i))
            else:
                name_set = {(protein.id, i)}
                self.kmer_table[kmer] = name_set

    def contains(self, kmer):
        if kmer in self.kmer_table.keys():
            return True
        return False


def main(input_file, kmer_size, out_dir, out_file):
    start = time.time()
    # build_index = True by default
    proteins = pyfastx.Fasta(input_file)
    prot_db = db_index(kmer_size, proteins)
    # prot_db_alt = db_index(kmer_size)
    # print(prot_db_alt.num_seq)
    # duration = time.time() - start
    # print("Processed %.0f sequences in %.2fs" % (prot_db.num_seq, duration))


if __name__ == '__main__':
    # TODO: replace this with flags
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
