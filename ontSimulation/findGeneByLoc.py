import sys
import os
import re


def main(in_file, gtf_file_1, gtf_file_2, out_dir, out_file):
    in_fh = open(in_file, 'r')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    for line in in_fh:
        if line[0] == '@':
            read = line[5:]
            fields = read.split(';')
            loc = fields[0].split('_')
            cds_id = loc[0].replace('-', '_')
            refseq_fh = open(gtf_file_1, 'r')
            for l in refseq_fh:
                if re.search(cds_id, l):
                    annot_fields = l.split('\t')
                    desc = annot_fields[8].split(';')
                    for d in desc:
                        d = d.strip()
                        ids = d.split(' ')
                        if ids[0] == "gene":
                            gene_id = ids[1].replace('"', '')
                            break
                    mane_fh = open(gtf_file_2, 'r')
                    for l2 in mane_fh:
                        if re.search(gene_id, l2):
                            annot_fields = l2.split('\t')
                            type = annot_fields[2]
                            if type == "CDS":
                                desc = annot_fields[8].split(';')
                                for f in desc:
                                    f = f.strip().split(' ')
                                    field_name = f[0]
                                    if field_name == "protein_id":
                                        protein_id = f[1].replace('"', '')
                                        print(protein_id)
                                        out_fh.write(line[1:].replace('\n', '') + '\t' + protein_id + '\n')
                                        break
                                break
                    mane_fh.close()
                    break
            refseq_fh.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])