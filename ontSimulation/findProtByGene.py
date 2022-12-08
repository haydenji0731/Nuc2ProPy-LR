import sys
import os
import re
import pyfastx


def main(in_file, gtf_file, out_dir, out_file, read_file):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    read_path = os.path.join(out_dir, read_file)
    read_fh = open(read_path, 'w')
    reads = pyfastx.Fastx(in_file)
    for read in reads:
        name = read[0]
        tmp = name.split('|')[0]
        gene = tmp.split('-')[0]
        # print("gene is: " + gene)
        gtf_fh = open(gtf_file, 'r')
        for l in gtf_fh:
            if re.search(gene, l):
                annot_field = l.split('\t')
                type = annot_field[2]
                if type == "CDS":
                    desc = annot_field[8].split(';')
                    for f in desc:
                        f = f.strip().split(' ')
                        field_name = f[0]
                        if field_name == "protein_id":
                            protein_id = f[1].replace('"', '')
                            # print("protein is: " + protein_id)
                            out_fh.write(name + "\t" + protein_id + "\n")
                            read_fh.write(">" + read[0] + "\n" + read[1] + "\n" + "+\n" + read[2] + "\n")
                            break
                    break
        gtf_fh.close()
    out_fh.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])