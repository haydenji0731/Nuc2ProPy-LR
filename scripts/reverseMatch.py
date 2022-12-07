# import packages

import sys
import re


def main(annot_file, protein_file):
    annot_fh = open(annot_file, 'r')
    cnt = 0
    protein_dict = {}
    for l in annot_fh:
        desc = l.split('\t')[8]
        ids = desc.split(';')
        found = False
        for id in ids:
            tmp = id.strip().split(' ')
            if tmp[0] == 'protein_id':
                protein_id = tmp[1].replace('"', '')
                if protein_id in protein_dict.keys():
                    print("something might be wrong")
                    print(protein_id)
                else:
                    protein_dict[protein_id] = 1
                protein_fh = open(protein_file, 'r')
                for l2 in protein_fh:
                    if re.search(protein_id, l2):
                        cnt += 1
                        found = True
                        break
                protein_fh.close()
                if not found:
                    print(protein_id)
    annot_fh.close()
    print(cnt)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])