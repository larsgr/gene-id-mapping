#!/usr/bin/env python3
"""Summarise LiftoffTools transcript-level variant effects at the gene level."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEVERITY_ORDER = [
    "identical",
    "synonymous",
    "NA",
    "inframe_insertion",
    "inframe_deletion",
    "nonsynonymous",
    "frameshift",
    "start_lost",
    "stop_gained",
    "5'_truncated",
    "3'_truncated",
    "unmapped",
]
SEVERITY_RANK = {effect: idx for idx, effect in enumerate(SEVERITY_ORDER)}


def load_transcript_gene_map(gff_path: Path) -> dict[str, str]:
    """Return transcript -> gene mapping from an Ensembl-style GFF3."""
    mapping: dict[str, str] = {}
    pattern = re.compile(r"ID=transcript:([^;]+);Parent=gene:([^;]+)")
    with gff_path.open() as handle:
        for line in handle:
            if "\t" not in line or line.startswith("#"):
                continue
            match = pattern.search(line)
            if match:
                transcript_id, gene_id = match.groups()
                mapping[transcript_id] = gene_id
    return mapping


def classify_gene(effects: list[str], *, optimistic: bool) -> str:
    """Return best (optimistic) or worst (pessimistic) effect."""
    rank_default = len(SEVERITY_RANK)
    key = (min if optimistic else max)
    return key(effects, key=lambda eff: SEVERITY_RANK.get(eff, rank_default))


def aggregate_counts(per_gene: dict[str, list[str]], *, optimistic: bool) -> Counter[str]:
    summary: Counter[str] = Counter()
    for gene, effects in per_gene.items():
        if not effects:
            continue
        summary[classify_gene(effects, optimistic=optimistic)] += 1
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--variant-effects",
        required=True,
        type=Path,
        help="Path to LiftoffTools variant_effects file",
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="Reference annotation GFF3 used for the lift (provides transcript parentage)",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for the generated TSV outputs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    transcript_map = load_transcript_gene_map(args.gff)
    if not transcript_map:
        print("[warn] No transcript→gene mappings found in GFF; check file path", file=sys.stderr)

    per_gene: defaultdict[str, list[str]] = defaultdict(list)
    missing = Counter()

    with args.variant_effects.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("transcript:"):
                transcript_id = row[0].split(":", 1)[1]
            else:
                continue
            gene_id = transcript_map.get(transcript_id)
            if gene_id is None:
                missing[transcript_id] += 1
                continue
            if len(row) >= 5:
                effect = row[4] or "NA"
            elif len(row) == 2 and row[1] == "unmapped":
                effect = "unmapped"
            else:
                effect = "NA"
            per_gene[gene_id].append(effect)

    if missing:
        print(f"[warn] {len(missing)} transcripts missing gene mapping (skipped)", file=sys.stderr)

    args.outdir.mkdir(parents=True, exist_ok=True)
    best_path = args.outdir / "variant_effects_gene_best.tsv"
    worst_path = args.outdir / "variant_effects_gene_worst.tsv"

    with best_path.open("w") as handle:
        handle.write("gene_id\tbest_effect\n")
        for gene, effects in per_gene.items():
            handle.write(f"{gene}\t{classify_gene(effects, optimistic=True)}\n")

    with worst_path.open("w") as handle:
        handle.write("gene_id\tworst_effect\n")
        for gene, effects in per_gene.items():
            handle.write(f"{gene}\t{classify_gene(effects, optimistic=False)}\n")

    best_counts = aggregate_counts(per_gene, optimistic=True)
    worst_counts = aggregate_counts(per_gene, optimistic=False)
    total_genes = len(per_gene)

    def report(label: str, counts: Counter[str]) -> None:
        print(f"\n{label} ({total_genes} genes)")
        for effect, count in counts.most_common():
            pct = 100.0 * count / total_genes if total_genes else 0.0
            print(f"  {effect:15s} {count:6d} ({pct:5.2f}%)")

    report("Optimistic aggregation", best_counts)
    report("Pessimistic aggregation", worst_counts)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
