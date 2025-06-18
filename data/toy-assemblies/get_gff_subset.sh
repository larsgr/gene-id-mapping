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
        print chrom"\t"start"\t"end
    }'
}

# --- Extract GFF subset using bedtools ---
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
        bedtools intersect -a stdin -b <(echo "$bed_regions") -wa -f 1.0 | \
        cut -f4-
}

# Convert coordinates in GFF file to be relative to region start (reads GFF from stdin)
convert_gff_coordinates() {
    local bed_regions="$1"
    local out="$2"
    # Create an associative array with region information
    declare -A regions
    while IFS=$'\t' read -r chrom start end; do
        regions["$chrom"]="$chrom:$((start+1))-$end $start"
    done <<< "$bed_regions"
    # Process the GFF file from stdin
    while IFS=$'\t' read -r seqname source type start end score strand phase attributes; do
        if [[ -n "${regions["$seqname"]}" ]]; then
            read -r new_seqname region_start <<< "${regions["$seqname"]}"
            new_start=$((start - region_start))
            new_end=$((end - region_start))
            echo -e "$new_seqname\t$source\t$type\t$new_start\t$new_end\t$score\t$strand\t$phase\t$attributes"
        fi
    done > "$out"
}

# TODO: handle multiple regions on same chromosome
# Current


# Get regions in memory
BED_REGIONS="$(get_bed_regions "$FASTA_FILE")"

echo "Extracting subset of GFF..."

# Pipe extract_gff_subset into convert_gff_coordinates
extract_gff_subset "$GFF_PATH" "$BED_REGIONS" | convert_gff_coordinates "$BED_REGIONS" "$GFF_OUT"

echo "Done"

