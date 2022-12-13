import sys
import os
import pyfastx
import re

tp = 0
fp = 0
fn = 0
tn = 0


def main(read_file, aln_file, truth_file, out_dir, out_file):
    reads = pyfastx.Fastx(read_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    global tp
    global fp
    global fn
    global tn
    for read in reads:
        qname = read[0]
        found = False
        aln_fh = open(aln_file, 'r')
        for aln_l in aln_fh:
            aln_l_clean = aln_l.replace('\n', '')
            aln_qname = aln_l_clean.split("\t")[0]
            pred = aln_l_clean.split("\t")[1]
            if qname == aln_qname:
                found = True
                break
        aln_fh.close()
        truth_fh = open(truth_file, 'r')
        for truth_l in truth_fh:
            truth_l_clean = truth_l.replace('\n', '')
            truth_qname = truth_l_clean.split("\t")[0]
            label = truth_l_clean.split("\t")[1]
            if qname == truth_qname:
                break
        truth_fh.close()
        if not found:
            if label != ".":
                fn += 1
            else:
                tn += 1
        else:
            if pred == label:
                tp += 1
            # includes the case label == ".'
            else:
                fp += 1

    acc = (tp + tn) / (tp + tn + fp + fn)
    spec = tn / (tn + fp)
    sens = tp / (tp + fn)
    prec = tp / (tp + fp)
    f1 = (2 * sens * prec) / (sens + prec)
    out_fh.write(str(acc) + "\t" + str(spec) + "\t" + str(sens) +
                 "\t" + str(prec) + "\t" + str(f1))
    print("true positive is: " + str(tp))
    print("true negative is: " + str(tn))
    print("false positive is: " + str(fp))
    print("false negative is: " + str(fn))
    out_fh.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])