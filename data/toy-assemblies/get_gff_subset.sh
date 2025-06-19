#!/bin/zsh

set -euo pipefail

usage() {
    echo "Usage: $0 [-p download_path] <FASTA_FILE> <GFF_URL|GFF_PATH> <OUTPUT_GFF>" >&2
    echo "  -p, --download-path   Directory to store downloaded GFF (default: temp)" >&2
    exit 1
}

# Default values
DOWNLOAD_PATH="temp"

# Parse options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--download-path)
            DOWNLOAD_PATH="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        -* )
            echo "Unknown option: $1" >&2
            usage
            ;;
        * )
            break
            ;;
    esac
done

# Check positional arguments
if [[ $# -lt 3 ]]; then
    usage
fi
FASTA_FILE="$1"
GFF_ARG="$2"
GFF_OUT="$3"

# Determine if GFF_ARG is a URL or file path
if [[ "$GFF_ARG" =~ ^https?:// ]]; then
    GFF_URL="$GFF_ARG"
    GFF_PATH="$DOWNLOAD_PATH/$(basename "$GFF_URL")"
    # Download if not present
    if [ ! -f "$GFF_PATH" ]; then
        mkdir -p "$DOWNLOAD_PATH"
        echo "Downloading $GFF_URL to $GFF_PATH"
        curl -L "$GFF_URL" -o "$GFF_PATH"
    else
        echo "$GFF_PATH already exists, skipping download."
    fi
else
    GFF_PATH="$GFF_ARG"
    if [ ! -f "$GFF_PATH" ]; then
        echo "Error: GFF file '$GFF_PATH' does not exist." >&2
        exit 1
    fi
fi

# Extract regions from fasta headers into a variable (no intermediate BED file)
get_bed_regions() {
    local fasta="$1"
    grep '^>' "$fasta" | sed 's/^>//' | \
    awk -F'[: -]' '{
        if (NF != 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2) {
            print "Error: Invalid FASTA header format (expected >chrom:start-end, with end >= start): " $0 > "/dev/stderr";
            exit 1;
        }
        chrom=$1;
        start=$2-1; # convert to 0-based
        end=$3;
        region=chrom":"(start+1)"-"end;
        print chrom"\t"start"\t"end"\t"region
    }'
}

# --- Extract GFF subset using bedtools ---
# This function extracts GFF features that are fully contained within the regions defined in the FASTA headers.
# It works as follows:
#   1. Prepends 3 BED columns (chrom, start, end) to each GFF feature.
#   2. Uses bedtools intersect (-f 1.0) to require full overlap with the regions (BED), appending 4 BED columns at the end.
#   3. Outputs the original GFF columns (columns 4-12), but replaces the seqname (column 4) with the region name (column 16 from BED).
extract_gff_subset() {
    local gff_file="$1"
    local bed_regions="$2"
    if [[ "$gff_file" == *.gz ]]; then
        gff_stream="gunzip -c \"$gff_file\""
    else
        gff_stream="cat \"$gff_file\""
    fi
    eval $gff_stream | \
        awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $0}' | \
        bedtools intersect -a stdin -b <(echo "$bed_regions") -wa -wb -f 1.0 | \
        awk 'BEGIN{FS="\t"; OFS="\t"} {$4=$16; print $4, $5, $6, $7, $8, $9, $10, $11, $12}'
}

# Convert coordinates in GFF file to be relative to region start (reads GFF from stdin)
convert_gff_coordinates() {
    local out=$1
    # Process the GFF file from stdin
    # region_start is 1-based, need to subtract the 0-based
    while IFS=$'\t' read -r seqname source type start end score strand phase attributes; do
        region_start=$(echo "$seqname" | awk -F'[:\-]' '{print $2}')
        new_start=$((start - (region_start-1)))
        new_end=$((end - (region_start-1)))
        echo -e "$seqname\t$source\t$type\t$new_start\t$new_end\t$score\t$strand\t$phase\t$attributes"
    done > "$out"
}


# Get regions in memory
BED_REGIONS="$(get_bed_regions "$FASTA_FILE")"

echo "Extracting subset of GFF..."

# Pipe extract_gff_subset into convert_gff_coordinates
extract_gff_subset "$GFF_PATH" "$BED_REGIONS" | convert_gff_coordinates "$GFF_OUT"

echo "Done"

