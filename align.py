# import packages
import numpy
import buildIndex
import sys
import pyfastx
# supports BLOSUM 45, 50, 62, 80, & 90
import blosum as bl
import time

num_seqs = 0


def kmer_prefilter(aa_seq, prot_db, kmer_size):
    # TODO: implement consecutive diagonal matches?
    prefilter_set = set()
    # prefilter_diag_set = set()
    # diag_prev = dict()
    for i in range(0, len(aa_seq) - kmer_size + 1):
        kmer = aa_seq[i:i + kmer_size]
        if prot_db.contains(kmer):
            # default prefilter
            for target_idx, target_pos in prot_db.kmer_table[kmer]:
                prefilter_set.add(target_idx)

            # diagonal double match prefilter
            # target_set = prot_db.kmer_table[kmer]
            # for target_idx, target_pos in target_set:
            #     if target_idx in diag_prev.keys():
            #         if diag_prev[target_idx] == target_pos - i:
            #             prefilter_diag_set.add(target_idx)
            #         else:
            #             diag_prev[target_idx] = target_pos - i
            #     else:
            #         diag_prev[target_idx] = target_pos - i

    # return prefilter_set, prefilter_diag_set
    return prefilter_set


def sw_align(x, y, sub_mat):
    V = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            diag_c = sub_mat[x[i - 1] + y[j - 1]]
            vert_c = sub_mat[x[i - 1] + "*"]
            horz_c = sub_mat["*" + y[j - 1]]
            V[i, j] = max(V[i - 1, j - 1] + sub_mat[x[i - 1] + y[j - 1]],  # diagonal
                          V[i - 1, j] + sub_mat[x[i - 1] + "*"],  # vertical
                          V[i, j - 1] + sub_mat["*" + y[j - 1]],  # horizontal
                          0)
    # print_mat(V)
    argmax = numpy.where(V == V.max())
    if len(argmax[0]) >= 2:
        argmax = (numpy.array([argmax[0][0]]), numpy.array([argmax[1][0]]))
    return V, int(V[argmax])


# smith-waterman-gotoh algorithm
# TODO: check the score matrix
def sw_align_affine(x, y, go, ge, sub_mat):
    # ge = ge = -4, then normal sw alignment
    # end spaces are free
    H = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
    E = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
    F = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            E[i, j] = max(E[i, j - 1] - ge, H[i, j - 1] - go)
            F[i, j] = max(F[i - 1, j] - ge, H[i - 1, j] - go)
            H[i, j] = max(0, E[i, j], F[i, j], H[i - 1, j - 1] - sub_mat[x[i - 1] + y[j - 1]])
    # print_mat(H)
    argmax = numpy.where(H == H.max())
    if len(argmax[0]) >= 2:
        argmax = (numpy.array([argmax[0][0]]), numpy.array([argmax[1][0]]))
    return H, int(H[argmax])


# source: ben langmead
def trace_back(V, x, y, sub_mat):
    # get i, j for maximal cell
    i, j = numpy.unravel_index(numpy.argmax(V), V.shape)
    escript, alx, aly, alm = [], [], [], []
    while (i > 0 or j > 0) and V[i, j] != 0:
        diag, vert, horz = 0, 0, 0
        if i > 0 and j > 0:
            diag - V[i - 1, j - 1] + sub_mat[x[i - 1] + y[j - 1]]
        if i > 0:
            vert = V[i - 1, j] + sub_mat[x[i - 1] + "*"]
        if j > 0:
            horz = V[i, j - 1] + sub_mat["*" + y[j - 1]]
        if diag >= vert and diag >= horz:
            match = x[i - 1] == y[j - 1]
            escript.append('M' if match else 'S')
            alm.append('|' if match else ' ')
            alx.append(x[i - 1])
            aly.append(y[j - 1])
            i -= 1
            j -= 1
        elif vert >= horz:
            escript.append('D')
            alx.append(x[i - 1])
            aly.append('-')
            alm.append(' ')
            i -= 1
        else:
            escript.append('I')
            aly.append(y[j - 1])
            alx.append('-')
            alm.append(' ')
            j -= 1
    escript = (''.join(escript))[::-1]
    final_aln = '\n'.join(map(lambda x: ''.join(x), [alx[::-1], alm[::-1], aly[::-1]]))
    return escript, final_aln


def print_mat(V):
    num_rows, num_cols = V.shape
    for i in range(0, num_rows):
        print(V[i])


def main(query_file, target_file, kmer_size):
    start = time.time()
    querys = pyfastx.Fastx(query_file)
    proteins = pyfastx.Fastx(target_file)
    prot_db = buildIndex.db_index(kmer_size, proteins)
    sub_mat = bl.BLOSUM(62)
    aln_dict = {}
    seq_dict = {}
    for query in querys:
        global num_seqs
        num_seqs += 1
        print("Query is " + query[0])
        seq_dict[query[0]] = query[1]
        query_seq = query[1]
        if len(query_seq) < kmer_size:
            print("Query amino acid sequence is too short. No significant sequence similarity is expected.")
        else:
            query_pref = kmer_prefilter(query_seq, prot_db, kmer_size)
            max_aln_score = -1
            opt_aln = ()
            if len(query_pref) > 0:
                for target_idx in query_pref:
                    # print("Target is " + prot_db.seq_index[target_idx][1])
                    sw_mat, aln_score = sw_align(query_seq, prot_db.seq_index[target_idx][0], sub_mat)
                    # matches blastp parameter
                    # sw_mat, aln_score = sw_align_affine(query_seq, prot_db.seq_index[target_idx][0], 11, 1, sub_mat)
                    if aln_score > max_aln_score:
                        # store target name, score
                        max_aln_score = aln_score
                        opt_aln = (prot_db.seq_index[target_idx][1], aln_score)
                # store query name and its optimal aln
                aln_dict[query[0]] = opt_aln
                print("Optimal alignment:")
                print("\tProtein is: " + aln_dict[query[0]][0])
                print("\tAlignment score is: " + str(aln_dict[query[0]][1]))
            else:
                print("No matching kmer was found between query and target")
    duration = time.time() - start
    print("Processed %.0f sequences in %.2fs" % (prot_db.num_seq, duration))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
