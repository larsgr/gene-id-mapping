#!/usr/bin/env python3
"""Sort a GFF3 file by gene blocks.

The script ensures that every gene (and all of its descendant features) is
emitted together and that gene blocks are ordered by `(seqid, gene_start)`.
Comment lines and other features that do not belong to a gene are preserved in
their original order ahead of the sorted gene blocks.

Example
-------

```
gff_block_sort.py input.gff3 -o input.sorted.gff3
```

The implementation uses only the standard library and keeps the original order
of features inside each gene block.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


GFF_COLUMNS = 9

GENE_TYPES = {
    "gene",              # Standard gene feature (NCBI & Ensembl protein-coding + some non-coding)
    "ncRNA_gene",        # Ensembl non-coding RNA genes (lncRNA, rRNA, snRNA, snoRNA, miRNA, etc.)
    "pseudogene",        # Pseudogenes (both NCBI & Ensembl)
    "V_gene_segment",    # Immunoglobulin/T-cell receptor V segments
    "D_gene_segment",    # Immunoglobulin/T-cell receptor D segments
    "J_gene_segment",    # Immunoglobulin/T-cell receptor J segments
}


@dataclass
class GeneBlock:
    """Container holding all lines that belong to a gene."""

    gene_id: str
    seqid: str
    start: int
    order: int
    lines: List[str] = field(default_factory=list)

    def add_line(self, line: str) -> None:
        self.lines.append(line)


class BlockSorterError(RuntimeError):
    """Raised when the GFF input violates expectations."""


def parse_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in raw.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, value = part, ""
        attrs[key] = value
    return attrs


def assign_gene_id(
    feature_id: Optional[str],
    parents: Sequence[str],
    feature_to_gene: Dict[str, str],
) -> Optional[str]:
    """Return the gene ID for the feature, updating the lookup if possible."""

    gene_id: Optional[str] = None
    for parent in parents:
        gene_id = feature_to_gene.get(parent)
        if gene_id:
            break

    if gene_id and feature_id:
        feature_to_gene[feature_id] = gene_id

    return gene_id


def sort_gff_blocks(input_path: str) -> Tuple[List[str], List[GeneBlock], List[str]]:
    """Return header lines, sorted gene blocks, and trailing lines."""

    header_lines: List[str] = []
    trailing_lines: List[str] = []
    blocks: Dict[str, GeneBlock] = {}
    feature_to_gene: Dict[str, str] = {}
    order_counter = 0
    seen_gene = False
    current_gene_id: Optional[str] = None
    pending_by_gene: Dict[str, List[str]] = {}

    with open(input_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            if not line or line.startswith("#"):
                if not seen_gene:
                    header_lines.append(line)
                elif current_gene_id and current_gene_id in blocks:
                    blocks[current_gene_id].add_line(line)
                else:
                    trailing_lines.append(line)
                continue

            parts = line.split("\t")
            if len(parts) != GFF_COLUMNS:
                if not seen_gene:
                    header_lines.append(line)
                else:
                    trailing_lines.append(line)
                continue

            seqid, _source, feature_type, start_s, _end_s, _score, _strand, _phase, attr_raw = parts
            attrs = parse_attributes(attr_raw)
            feature_id = attrs.get("ID")
            parents = attrs.get("Parent", "").split(",") if attrs.get("Parent") else []

            if feature_type in GENE_TYPES and feature_id:
                try:
                    start = int(start_s)
                except ValueError as exc:
                    raise BlockSorterError(f"Invalid start coordinate for {feature_type} {feature_id}: {start_s}") from exc

                block = GeneBlock(gene_id=feature_id, seqid=seqid, start=start, order=order_counter)
                block.add_line(line)
                blocks[feature_id] = block
                feature_to_gene[feature_id] = feature_id
                order_counter += 1
                seen_gene = True
                current_gene_id = feature_id
                if feature_id in pending_by_gene:
                    for pending_line in pending_by_gene.pop(feature_id):
                        block.add_line(pending_line)
                continue

            if not seen_gene:
                header_lines.append(line)
                continue

            gene_id = assign_gene_id(feature_id, parents, feature_to_gene)

            if gene_id is None:
                if parents:
                    parent_gene = parents[0]
                    pending_by_gene.setdefault(parent_gene, []).append(line)
                    if feature_id:
                        feature_to_gene[feature_id] = parent_gene
                    current_gene_id = None
                    continue
                trailing_lines.append(line)
                current_gene_id = None
                continue

            block = blocks.get(gene_id)
            if block is None:
                pending_by_gene.setdefault(gene_id, []).append(line)
                if feature_id:
                    feature_to_gene[feature_id] = gene_id
                current_gene_id = None
                continue

            block.add_line(line)
            current_gene_id = gene_id

    sorted_blocks = sorted(blocks.values(), key=lambda b: (b.seqid, b.start, b.order))
    return header_lines, sorted_blocks, trailing_lines


def write_sorted_gff(
    header_lines: Iterable[str],
    blocks: Iterable[GeneBlock],
    trailing_lines: Iterable[str],
    output_path: Optional[str],
) -> None:
    out_handle = open(output_path, "w", encoding="utf-8") if output_path else sys.stdout

    def write_lines(lines: Iterable[str]) -> None:
        for text in lines:
            out_handle.write(text)
            out_handle.write("\n")

    try:
        header_list = list(header_lines)
        block_list = list(blocks)
        trailing_list = list(trailing_lines)

        if header_list:
            write_lines(header_list)

        for block in block_list:
            write_lines(block.lines)

        if trailing_list:
            write_lines(trailing_list)
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sort GFF3 gene blocks by genomic position")
    parser.add_argument("input", help="Input GFF3 file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output path (defaults to stdout)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    try:
        header, blocks, trailing = sort_gff_blocks(args.input)
    except BlockSorterError as exc:
        sys.stderr.write(f"Error: {exc}\n")
        return 1

    write_sorted_gff(header, blocks, trailing, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
