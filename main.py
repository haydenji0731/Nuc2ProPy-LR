# import packages
import argparse
import sys
import search


def main():
    arg_parser = argparse.ArgumentParser(description='lonblast: long DNA read based protein search tool')
    usage = """
    lonblast [<options>] <query-file> <target-file>
    lonblast reports protein alignments between the query DNA reads (long) and target protein database. Input DNA file
    must either be in FASTQ or FASTA format.
    """
    # TODO: what doest store do?
    arg_parser.add_argument("--out-dir", help="Directory where the output file will be written to", action="store")
    arg_parser.add_argument("--out-file", help="Output file name that will store alignment results", action="store")
    arg_parser.add_argument("--aa-file", help="Output file name that will store all 6 reading frames amino acid "
                                              "sequence", action="store")
    arg_parser.add_argument("--extract-orf", help="Extract valid open reading frames from the input reads before "
                                                  "alignment (default = False)", default=False)
    arg_parser.add_argument("--orf-file", help="Output file name that will store the amino acid sequence after orf "
                                               "extraction")
    arg_parser.add_argument("--kmer-size", help="Kmer size that is used to index the target (default = 7)", default=7)
    arg_parser.add_argument("--gap-open", help="Gap open penalties used in scoring alignments", default=4)
    arg_parser.add_argument("--gap-extend", help="Gap extend penalties used in scoring alignments", default=4)
    arg_parser.add_argument("query", help="The input DNA read file in FASTA/FASTQ format", action="store")
    arg_parser.add_argument("target", help="Target protein database in FASTA format", action="store")
    args = arg_parser.parse_args()
    # print(args)
    # if infile not specified, exit
    query_file = args.query
    if not query_file:
        arg_parser.print_help()
        raise Exception("Must provide a query DNA read file.")
        sys.exit(1)

    # if target db not provided, exit
    target_file = args.target
    if not target_file:
        arg_parser.print_help()
        raise Exception("Must provide a target protein DB file.")
        sys.exit(1)

    # if out dir not provided, exit
    out_dir = args.out_dir
    if not out_dir:
        arg_parser.print_help()
        raise Exception("Must provide an output directory.")
        sys.exit(1)

    # if out file not provided, exit
    out_file = args.out_file
    if not out_file:
        arg_parser.print_help()
        raise Exception("Must provide an output file name.")
        sys.exit(1)

    extract_orf = args.extract_orf
    # print(extract_orf)

    # if extract orf is true but if orf file not provided, exit
    orf_file = args.orf_file
    if extract_orf and not orf_file:
        arg_parser.print_help()
        raise Exception("If extract-orf option is used, must provide an orf file name.")
        sys.exit(1)

    aa_file = args.aa_file
    if not extract_orf and not aa_file:
        arg_parser.print_help()
        raise Exception("If extract-orf option is not used, must provide an aa file name.")
        sys.exit(1)

    kmer_size = int(args.kmer_size)
    gap_open = int(args.gap_open)
    gap_extend = int(args.gap_extend)

    search.main(query_file, target_file, out_dir, out_file, aa_file, kmer_size,
                gap_open, gap_extend, extract_orf, orf_file)


if __name__ == '__main__':
    main()
