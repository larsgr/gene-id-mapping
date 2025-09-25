#!/usr/bin/env python3
"""Compare per-gene Liftoff effects between named and unnamed genes."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import Counter
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent / ".matplotlib-cache"),
)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")  # headless backend for scripted runs
import matplotlib.pyplot as plt
import numpy as np

AGGREGATED_ORDER = [
    "identical",
    "synonymous",
    "protein_change",
    "truncation",
    "unmapped",
    "noncoding",
    "other",
]
PROTEIN_CHANGE = {
    "nonsynonymous",
    "frameshift",
    "start_lost",
    "stop_gained",
    "inframe_deletion",
    "inframe_insertion",
}
TRUNCATIONS = {"5'_truncated", "3'_truncated"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-effects", required=True, type=Path, help="TSV with columns gene_id,best_effect")
    parser.add_argument("--gff", required=True, type=Path, help="Annotation GFF3 containing gene Name attributes")
    parser.add_argument("--output-prefix", required=True, type=Path, help="Prefix for PNG/TSV outputs")
    return parser.parse_args()


def load_gene_names(gff_path: Path) -> dict[str, bool]:
    pattern = re.compile(r"ID=gene:([^;]+)")
    gene_has_name: dict[str, bool] = {}
    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#") or "\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            match = pattern.search(attrs)
            if not match:
                continue
            gene_id = match.group(1)
            name_value = None
            for attr in attrs.split(";"):
                if attr.startswith("Name="):
                    name_value = attr.split("=", 1)[1]
                    break
            gene_has_name[gene_id] = bool(name_value and name_value != gene_id)
    return gene_has_name


def aggregate(effect: str) -> str:
    if effect == "identical":
        return "identical"
    if effect == "synonymous":
        return "synonymous"
    if effect in PROTEIN_CHANGE:
        return "protein_change"
    if effect in TRUNCATIONS:
        return "truncation"
    if effect == "unmapped":
        return "unmapped"
    if effect == "NA":
        return "noncoding"
    return "other"


def proportion_z_test(a_success: int, a_total: int, b_success: int, b_total: int) -> tuple[float, float]:
    """Return z-score and two-sided p-value for difference in proportions."""
    p1 = a_success / a_total
    p2 = b_success / b_total
    pooled = (a_success + b_success) / (a_total + b_total)
    se = math.sqrt(pooled * (1.0 - pooled) * (1 / a_total + 1 / b_total))
    z = (p1 - p2) / se if se else float("inf")
    # Survival function for standard normal
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return z, p


def main() -> int:
    args = parse_args()
    output_dir = args.output_prefix.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_has_name = load_gene_names(args.gff)

    counts = {"named": Counter(), "unnamed": Counter()}
    raw_effects = {"named": Counter(), "unnamed": Counter()}
    missing = 0

    with args.gene_effects.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"]
            effect = row["best_effect"]
            has_name = gene_has_name.get(gene_id)
            if has_name is None:
                has_name = False
                missing += 1
            group = "named" if has_name else "unnamed"
            raw_effects[group][effect] += 1
            counts[group][aggregate(effect)] += 1

    totals = {group: sum(counts[group].values()) for group in counts}

    print("Gene counts by naming status:")
    for group in counts:
        print(f"  {group:7s}: {totals[group]:6d} genes")
    if missing:
        print(f"  [info] {missing} genes missing explicit Name attribute treated as unnamed")

    summary_rows = []
    percentages = {group: {} for group in counts}
    for group in counts:
        total = totals[group]
        for cat in AGGREGATED_ORDER:
            value = counts[group][cat]
            pct = 100.0 * value / total if total else 0.0
            percentages[group][cat] = pct
            summary_rows.append((group, cat, value, pct))

    # Proportion test for identical category
    z, p = proportion_z_test(
        raw_effects["named"]["identical"],
        totals["named"],
        raw_effects["unnamed"]["identical"],
        totals["unnamed"],
    )
    print(f"\nIdentical proportion: named={percentages['named']['identical']:.2f}% vs unnamed={percentages['unnamed']['identical']:.2f}%")
    print(f"Two-sided z-test: z={z:.2f}, p={p:.2e}")

    summary_tsv = args.output_prefix.parent / f"{args.output_prefix.name}_named_vs_unnamed.tsv"
    with summary_tsv.open("w") as handle:
        handle.write("group\tcategory\tgenes\tpercent\n")
        for group, cat, value, pct in summary_rows:
            handle.write(f"{group}\t{cat}\t{value}\t{pct:.4f}\n")
        handle.write(f"z_test\tidentical\t{z:.4f}\t{p:.4e}\n")

    category_labels = ["identical", "synonymous", "protein_change", "truncation", "unmapped", "noncoding"]
    idx = np.arange(len(category_labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    named_vals = [percentages["named"][cat] for cat in category_labels]
    unnamed_vals = [percentages["unnamed"][cat] for cat in category_labels]

    ax.bar(idx - width / 2, named_vals, width, label="Named genes", color="#4C72B0")
    ax.bar(idx + width / 2, unnamed_vals, width, label="No Name attribute", color="#DD8452")
    ax.set_xticks(idx)
    ax.set_xticklabels([label.replace("_", " ").title() for label in category_labels], rotation=18, ha="right")
    ax.set_ylabel("Genes (%)")
    max_height = max(named_vals + unnamed_vals) if (named_vals or unnamed_vals) else 0
    ax.set_ylim(0, max_height * 1.2 if max_height else 1)
    ax.legend(frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.set_title("Gene-level Liftoff outcomes by naming status")
    fig.tight_layout()

    chart_path = args.output_prefix.with_suffix(".png")
    fig.savefig(chart_path, dpi=180)
    print(f"Saved bar chart to {chart_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
