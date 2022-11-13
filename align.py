# import packages

import numpy
import buildIndex
import sys
import pyfastx

def kmer_prefilter(aa_seq, prot_db, kmer_size):
    # TODO: implement consecutive diagonal matches?
    prefilter_set = set()
    for i in range(0, len(aa_seq) - kmer_size + 1):
        kmer = aa_seq[i:i + kmer_size]
        if kmer in prot_db.contains(kmer):
            # merge 2 sets
            updated_set = prefilter_set | prot_db[kmer]

    # diagonal double match prefilter



def smith_waterman_align(x, y, cost_func):
    V = numpy.zeros(len(x) + 1, len(y) + 1, dtype=int)
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            V[i, j] = max(V[i - 1, j - 1] + cost_func(x[i - 1], y[j - 1]),  # diagonal
                          V[i - 1, j] + cost_func(x[i - 1], '-'),  # vertical
                          V[i, j - 1] + cost_func('-', y[j - 1]),  # horizontal
                          0)
    argmax = numpy.where(V == V.max())
    return V, int(V[argmax])


def compute_cost(xc, yc, go, ge, sub, match, is_open):
    if xc == yc:
        return match
    if xc == '-' or yc == '-':
        if is_open:
            return ge
        else:
            return go
    return sub


def main(input_file, kmer_size):
    proteins = pyfastx.Fastx(input_file)
    prot_db = buildIndex.db_index(kmer_size, proteins)

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
