import sys
import os
import pyfastx
import re


def main(in_file, gtf_file, kmer_size, out_dir, out_file):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    reads = pyfastx.Fastx(in_file)
    for read in reads:
        name = read[0].split(':')[0] + ":" + read[0].split(':')[1]
        gene = read[0].split('=')[1]
        if len(read[1]) < int(kmer_size) * 3:
            out_fh.write(name + "\t" + "." + "\n")
            continue
        gtf_fh = open(gtf_file, 'r')
        found = False
        for line in gtf_fh:
            if re.search(gene, line):
                found = True
                fields = line.split('\t')
                if fields[2] == "CDS":
                    desc = fields[8].split(';')
                    for f in desc:
                        f = f.strip().split(' ')
                        field_name = f[0]
                        if field_name == "protein_id":
                            protein_id = f[1].replace('"', '')
                            out_fh.write(name + "\t" + protein_id + "\n")
                            out_fh.flush()
                            break
                    break
        if not found:
            out_fh.write(name + "\t" + "." + "\n")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])