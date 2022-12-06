# import packages

import sys
import os
import pyfastx


def main(annot_file, protein_file, out_dir, out_file):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fa_path = os.path.join(out_dir, out_file)
    out_fh = open(out_fa_path, 'w')
    annot_fh = open(annot_file, 'r')
    cnt = 0
    for line in annot_fh:
        desc = line.split('\t')[8]
        ids = desc.split(';')
        for id in ids:
            tmp = id.strip().split(' ')
            if tmp[0] == 'protein_id':
                protein_id = tmp[1].replace('"', '')
                found = False
                proteins = pyfastx.Fastx(protein_file)
                for protein in proteins:
                    if protein[0] == protein_id:
                        found = True
                        cnt += 1
                        out_fh.write(">" + protein[0] + "\t" + protein[2] + "\n" +
                                     protein[1] + "\n")
                        break
                if not found:
                    print("something might be wrong")
                    print(protein_id)
    annot_fh.close()
    out_fh.close()
    print(cnt)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


