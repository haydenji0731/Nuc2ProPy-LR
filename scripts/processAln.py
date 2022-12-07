import sys
import os


def main(aln_file, out_dir, out_file):
    aln_fh = open(aln_file, 'r')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    qid = ""
    for line in aln_fh:
        fields = line.split("\t")
        if qid != fields[0]:
            qid = fields[0]
            out_fh.write(fields[0] + "\t" + fields[1] + "\n")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])