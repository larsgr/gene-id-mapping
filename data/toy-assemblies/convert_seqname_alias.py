# Convert seqname in gff or other tab separated file (from stdin) 
# using alias file. Output to stdout.
#
# alias file is tab separated table with one row per 
# sequence name, with different aliases in each column.
# e.g.:
# ssa01	NC_027300.1	CM003279.1	ssa01
# ssa02	NC_027301.1	CM003280.1	ssa02
# ssa03	NC_027302.1	CM003281.1	ssa03
#
# The first column is typically the name we want to convert to.
# Note that the same name/alias can be repeated several times 
# in a row, but never occur on multiple rows.
#
# Include options to select the column in the alias file (default to first), and 
# which column in the input file which contains the seqname (default to first)
#
# If any seqname is not found in the alias table exit with error.

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Convert seqname in tab-separated file using alias file.")
    parser.add_argument('alias_file', help='Tab-separated alias file')
    parser.add_argument('-a', '--alias-col', type=int, default=0, help='Column in alias file to convert to (0-based, default: 0)')
    parser.add_argument('-i', '--input-col', type=int, default=0, help='Column in input file with seqname (0-based, default: 0)')
    return parser.parse_args()

def build_alias_map(alias_file, alias_col):
    alias_map = {}
    with open(alias_file) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            target = cols[alias_col]
            for alias in cols:
                alias_map[alias] = target
    return alias_map

def convert_seqnames(alias_map, input_col):
    for line in sys.stdin:
        if not line.strip():
            print(line, end='')
            continue
        cols = line.rstrip('\n').split('\t')
        seqname = cols[input_col]
        if seqname not in alias_map:
            sys.stderr.write(f"Error: seqname '{seqname}' not found in alias table\n")
            sys.exit(1)
        cols[input_col] = alias_map[seqname]
        print('\t'.join(cols))

def main():
    args = parse_args()
    alias_map = build_alias_map(args.alias_file, args.alias_col)
    convert_seqnames(alias_map, args.input_col)

if __name__ == '__main__':
    main()

