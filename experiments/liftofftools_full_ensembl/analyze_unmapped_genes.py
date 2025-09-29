#!/usr/bin/env python3
"""Investigate Liftoff-unmapped genes in the ICSASG_v2 Ensembl annotation."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

ATTR_PATTERN = re.compile(r"([^=;]+)=([^;]*)")
GENE_ID_PATTERN = re.compile(r"ID=gene:([^;]+)")
TRANSCRIPT_PARENT_PATTERN = re.compile(r"Parent=gene:([^;]+)")


def parse_attributes(attr_field: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in ATTR_PATTERN.finditer(attr_field)}


def load_gene_metadata(gff_path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    with gff_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, feature_type, _, _, _, _, strand, attrs = parts
            if "gene" in feature_type:
                match = GENE_ID_PATTERN.search(attrs)
                if not match:
                    continue
                gene_id = match.group(1)
                attr_map = parse_attributes(attrs)
                metadata[gene_id] = {
                    "seqid": seqid,
                    "strand": strand,
                    "biotype": attr_map.get("biotype", ""),
                    "name": attr_map.get("Name", ""),
                    "description": attr_map.get("description", ""),
                }
            elif feature_type in {"mRNA", "transcript"}:
                match = TRANSCRIPT_PARENT_PATTERN.search(attrs)
                if not match:
                    continue
                gene_id = match.group(1)
                if gene_id not in metadata:
                    metadata[gene_id] = {
                        "seqid": seqid,
                        "strand": strand,
                        "biotype": "",
                        "name": "",
                        "description": "",
                    }
    return metadata


def load_unmapped_genes(best_effect_path: Path) -> set[str]:
    unmapped: set[str] = set()
    with best_effect_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["best_effect"] == "unmapped":
                unmapped.add(row["gene_id"])
    return unmapped


def classify_seqid(seqid: str) -> str:
    if not seqid:
        return "unknown"
    if seqid.startswith("ssa") and seqid[3:].isdigit():
        return "anchored_chromosome"
    return "unplaced_scaffold"


def summarize(metadata: dict[str, dict[str, str]], unmapped: set[str]) -> dict[str, Counter]:
    location_counter = Counter()
    biotype_counter = Counter()
    name_counter = Counter()
    missing_metadata = []

    for gene_id in unmapped:
        info = metadata.get(gene_id)
        if info is None:
            missing_metadata.append(gene_id)
            continue
        location_counter[classify_seqid(info.get("seqid", ""))] += 1
        biotype_counter[info.get("biotype", "") or "(missing)"] += 1
        name_counter["has_name" if info.get("name") else "no_name"] += 1

    return {
        "location": location_counter,
        "biotype": biotype_counter,
        "name": name_counter,
        "missing": Counter({"missing_metadata": len(missing_metadata)}),
    }


def write_detailed_table(output_path: Path, metadata: dict[str, dict[str, str]], unmapped: set[str]) -> None:
    with output_path.open("w") as handle:
        handle.write("gene_id\tseqid\tstrand\tbiotype\tname\tdescription\n")
        for gene_id in sorted(unmapped):
            info = metadata.get(gene_id, {})
            handle.write(
                "\t".join(
                    [
                        gene_id,
                        info.get("seqid", ""),
                        info.get("strand", ""),
                        info.get("biotype", ""),
                        info.get("name", ""),
                        info.get("description", ""),
                    ]
                )
                + "\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best-effects", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata = load_gene_metadata(args.gff)
    unmapped_genes = load_unmapped_genes(args.best_effects)
    summary = summarize(metadata, unmapped_genes)

    detailed_path = args.outdir / "unmapped_genes_annotation.tsv"
    write_detailed_table(detailed_path, metadata, unmapped_genes)

    summary_path = args.outdir / "unmapped_genes_summary.json"
    with summary_path.open("w") as handle:
        json.dump({k: dict(v) for k, v in summary.items()}, handle, indent=2, sort_keys=True)

    print(f"Unmapped genes analysed: {len(unmapped_genes)}")
    for category, counter in summary.items():
        print(f"\n{category.capitalize()} counts:")
        for label, count in counter.most_common():
            print(f"  {label:20s} {count:6d}")

    print(f"Detailed annotations written to {detailed_path}")
    print(f"Summary JSON written to {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
