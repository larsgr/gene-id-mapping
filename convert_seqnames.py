#!/usr/bin/env python3
"""
Convert sequence names in a GFF file using a mapping table.
"""

import sys
import argparse
import gzip
from pathlib import Path


def load_mapping(tsv_file, from_col, to_col):
    """Load sequence name mapping from TSV file.

    Args:
        tsv_file: Path to TSV file with header
        from_col: Column name to map from
        to_col: Column name to map to

    Returns:
        Dictionary mapping from_col values to to_col values
    """
    mapping = {}
    with open(tsv_file) as f:
        header = f.readline().strip().split('\t')

        try:
            from_idx = header.index(from_col)
            to_idx = header.index(to_col)
        except ValueError as e:
            print(f"Error: Column not found in header: {e}", file=sys.stderr)
            print(f"Available columns: {header}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(from_idx, to_idx):
                mapping[fields[from_idx]] = fields[to_idx]

    return mapping


def convert_gff(input_file, output_file, mapping):
    """Convert sequence names in GFF file.

    Args:
        input_file: Input GFF file (can be .gz)
        output_file: Output GFF file
        mapping: Dictionary mapping old names to new names
    """
    # Open input file (handle gzip)
    if input_file.endswith('.gz'):
        in_f = gzip.open(input_file, 'rt')
    else:
        in_f = open(input_file, 'r')

    # Open output file
    out_f = open(output_file, 'w')

    unmapped_count = 0
    unmapped_seqs = set()
    converted_count = 0

    for line in in_f:
        if line.startswith('#'):
            # Pass through comments
            out_f.write(line)
        else:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 9:
                seqid = fields[0]
                if seqid in mapping:
                    fields[0] = mapping[seqid]
                    converted_count += 1
                else:
                    unmapped_count += 1
                    unmapped_seqs.add(seqid)

                out_f.write('\t'.join(fields) + '\n')
            else:
                # Malformed line, pass through
                out_f.write(line)

    in_f.close()
    out_f.close()

    print(f"Converted {converted_count} lines", file=sys.stderr)
    if unmapped_count > 0:
        print(f"Warning: {unmapped_count} lines with unmapped sequence IDs", file=sys.stderr)
        print(f"Unmapped sequences: {sorted(unmapped_seqs)[:10]}...", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Convert sequence names in GFF file using a mapping table'
    )
    parser.add_argument('gff_input', help='Input GFF file (can be .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file')
    parser.add_argument('-m', '--mapping', required=True, help='TSV mapping file')
    parser.add_argument('--from-col', default='RefSeq-Accn',
                        help='Column name to map from (default: RefSeq-Accn)')
    parser.add_argument('--to-col', default='Ensembl-Name',
                        help='Column name to map to (default: Ensembl-Name)')

    args = parser.parse_args()

    # Load mapping
    print(f"Loading mapping from {args.mapping}...", file=sys.stderr)
    mapping = load_mapping(args.mapping, args.from_col, args.to_col)
    print(f"Loaded {len(mapping)} mappings", file=sys.stderr)

    # Convert GFF
    print(f"Converting {args.gff_input} -> {args.output}...", file=sys.stderr)
    convert_gff(args.gff_input, args.output, mapping)
    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
