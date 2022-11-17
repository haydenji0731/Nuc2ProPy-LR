import pyfastx
import sys


def main(input_file):
    proteins = pyfastx.Fastx(input_file)
    for protein in proteins:
        print(protein)


if __name__ == '__main__':
    main(sys.argv[1])