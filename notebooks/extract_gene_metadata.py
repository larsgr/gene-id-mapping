#!/usr/bin/env python3
"""
Extract gene-level metadata from Ensembl and NCBI GFF3 files into compact TSVs.

Avoids loading 700-800MB GFFs in R by pre-extracting relevant fields.

Usage:
    python3 notebooks/extract_gene_metadata.py

Outputs:
    notebooks/ens_gene_metadata.tsv
    notebooks/ncbi_gene_metadata.tsv
"""

import os
import sys
import re
from urllib.parse import unquote

# Paths (relative to repo root)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

ENS_GFF = os.path.join(REPO_ROOT, "data/genomes/AtlanticSalmon/Ssal_v3.1_Ens.gff3")
NCBI_GFF = os.path.join(REPO_ROOT, "data/genomes/AtlanticSalmon/Ssal_v3.1_NCBI.gff3")

ENS_OUT = os.path.join(SCRIPT_DIR, "ens_gene_metadata.tsv")
NCBI_OUT = os.path.join(SCRIPT_DIR, "ncbi_gene_metadata.tsv")

# Ensembl gene-level feature types
ENS_GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}


def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes column into a dict, URL-decoding values."""
    attrs = {}
    for item in attr_string.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = unquote(val)
    return attrs


def extract_description_source(description):
    """Extract the source from Ensembl description like 'foo [Source:ZFIN;Acc:ZDB-...]'."""
    m = re.search(r'\[Source:([^;]+)', description)
    return m.group(1) if m else "NA"


def extract_ensembl_metadata():
    """Extract gene metadata from Ensembl GFF3."""
    print(f"Reading {ENS_GFF} ...")
    count = 0
    with open(ENS_GFF) as fin, open(ENS_OUT, "w") as fout:
        fout.write("gene_id\tfeature_type\tbiotype\tname\tdescription\t"
                   "description_source\tseqid\tstart\tend\tstrand\n")
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ENS_GENE_TYPES:
                continue

            seqid = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attrs = parse_gff_attributes(parts[8])

            # Gene ID: strip "gene:" prefix from ID attribute
            gene_id = attrs.get("ID", "")
            if gene_id.startswith("gene:"):
                gene_id = gene_id[5:]

            biotype = attrs.get("biotype", "NA")
            name = attrs.get("Name", "NA")
            description = attrs.get("description", "NA")
            desc_source = extract_description_source(description) if description != "NA" else "NA"

            fout.write(f"{gene_id}\t{feature_type}\t{biotype}\t{name}\t"
                       f"{description}\t{desc_source}\t{seqid}\t{start}\t{end}\t{strand}\n")
            count += 1

    print(f"  Wrote {count} genes to {ENS_OUT}")


def extract_ncbi_metadata():
    """Extract gene metadata from NCBI GFF3."""
    print(f"Reading {NCBI_GFF} ...")
    count = 0
    with open(NCBI_GFF) as fin, open(NCBI_OUT, "w") as fout:
        fout.write("gene_id\tbiotype\tname\tis_loc\tncbi_gene_id\t"
                   "description\tseqid\tstart\tend\tstrand\n")
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue

            seqid = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attrs = parse_gff_attributes(parts[8])

            gene_id = attrs.get("ID", "")
            biotype = attrs.get("gene_biotype", "NA")
            name = attrs.get("Name", "NA")
            is_loc = "TRUE" if name.startswith("LOC") else "FALSE"
            description = attrs.get("description", "NA")

            # Extract GeneID from Dbxref (e.g. "GeneID:123729278")
            dbxref = attrs.get("Dbxref", "")
            ncbi_gene_id = "NA"
            for ref in dbxref.split(","):
                if ref.startswith("GeneID:"):
                    ncbi_gene_id = ref.split(":", 1)[1]
                    break

            fout.write(f"{gene_id}\t{biotype}\t{name}\t{is_loc}\t{ncbi_gene_id}\t"
                       f"{description}\t{seqid}\t{start}\t{end}\t{strand}\n")
            count += 1

    print(f"  Wrote {count} genes to {NCBI_OUT}")


if __name__ == "__main__":
    extract_ensembl_metadata()
    extract_ncbi_metadata()
    print("Done.")
