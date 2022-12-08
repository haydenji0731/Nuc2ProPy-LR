import sys
import pyfastx
import os


def main(in_file, out_dir, out_file):
    if not os.path.exists(out_dir):
        os.makesirs(out_dir)
    out_path = os.path.join(out_dir, out_file)
    out_fh = open(out_path, 'w')
    reads = pyfastx.Fastx(in_file)
    for read in reads:
        out_fh.write(read[0] + "\t" + "." + "\n")
    out_fh.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])