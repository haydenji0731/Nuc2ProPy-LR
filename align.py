# import packages
import numpy
import findOrf

# def kmer_prefilter(aa_seq, prot_db):
#     prot_db.


def smith_waterman_align(x, y, s):
    V = numpy.zeros(len(x) + 1, len(y) + 1, dtype=int)
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            V[i, j] = max(V[i - 1, j - 1] + s(x[i - 1], y[j - 1]),  # diagonal
                          V[i - 1, j] + s(x[i - 1], '-'),  # vertical
                          V[i, j - 1] + s('-', y[j - 1]),  # horizontal
                          0)
    argmax = numpy.where(V == V.max())
    return V, int(V[argmax])


def compute_cost(xc, yc):
    if xc == yc:
        return 1
    # uniform gap cost (open / extend)
    if xc == '-' or yc == '-':
        return -1
    return -1


def main():
    print("implement sw alignment")


if __name__ == '__main__':
    main()
