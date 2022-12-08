import sys
import pyfastx
import os


def main(fa_file, out_dir, out_file):
    seqs = pyfastx.Fastx(fa_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    for seq in seqs:
        desc = seq[2]
        fields = desc.split(' ')
        gene = fields[0]
        gene = gene.replace('[', '')
        gene = gene.replace(']', '')
        id = gene.split('=')[1]
        new_seq_id = id + "_" + seq[0]
        out_fh.write(">" + new_seq_id + "\n" + seq[1] + "\n")
        print(new_seq_id)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])