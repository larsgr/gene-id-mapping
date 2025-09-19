#!/usr/bin/env python3
"""Within-assembly GFF3 comparison parser.

This script compares two or more sorted GFF3 annotation files from the same
assembly and reports per-gene and per-transcript overlap statistics for every
pair of overlapping genes across annotations. It has no third-party
dependencies and focuses on structural and CDS-aware similarities.
"""

from __future__ import annotations

import argparse
import collections
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

Interval = Tuple[int, int]


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


def merge_intervals(intervals: Sequence[Interval]) -> List[Interval]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: List[List[int]] = [[sorted_intervals[0][0], sorted_intervals[0][1]]]
    for start, end in sorted_intervals[1:]:
        last = merged[-1]
        if start <= last[1] + 1:
            if end > last[1]:
                last[1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def intervals_overlap(intervals_a: Sequence[Interval], intervals_b: Sequence[Interval]) -> bool:
    i = j = 0
    while i < len(intervals_a) and j < len(intervals_b):
        a_start, a_end = intervals_a[i]
        b_start, b_end = intervals_b[j]
        if a_end < b_start:
            i += 1
        elif b_end < a_start:
            j += 1
        else:
            return True
    return False


def intersection_length(intervals_a: Sequence[Interval], intervals_b: Sequence[Interval]) -> int:
    i = j = 0
    total = 0
    while i < len(intervals_a) and j < len(intervals_b):
        a_start, a_end = intervals_a[i]
        b_start, b_end = intervals_b[j]
        start = max(a_start, b_start)
        end = min(a_end, b_end)
        if start <= end:
            total += end - start + 1
        if a_end <= b_end:
            i += 1
        else:
            j += 1
    return total


def intervals_length(intervals: Sequence[Interval]) -> int:
    return sum(end - start + 1 for start, end in intervals)


# ---------------------------------------------------------------------------
# Domain objects
# ---------------------------------------------------------------------------

@dataclass
class CodingSegment:
    start: int
    end: int
    strand: str
    phase_start: int
    phase_end: int

    def phase_at(self, pos: int) -> Optional[int]:
        if pos < self.start or pos > self.end:
            return None
        if self.strand == "+":
            return (self.phase_start + (pos - self.start)) % 3
        if self.strand == "-":
            return (self.phase_end + (self.end - pos)) % 3
        return None


class Transcript:
    __slots__ = (
        "id",
        "gene_id",
        "seqid",
        "strand",
        "span_start",
        "span_end",
        "exons",
        "cds_raw",
        "merged_exons",
        "intron_boundaries",
        "exon_bp",
        "coding_segments",
        "coding_bp",
    )

    def __init__(
        self,
        transcript_id: str,
        gene_id: str,
        seqid: str,
        strand: str,
        span_start: Optional[int] = None,
        span_end: Optional[int] = None,
    ) -> None:
        self.id = transcript_id
        self.gene_id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.span_start = span_start if span_start is not None else sys.maxsize
        self.span_end = span_end if span_end is not None else -sys.maxsize
        self.exons: List[Interval] = []
        self.cds_raw: List[Interval] = []
        self.merged_exons: List[Interval] = []
        self.intron_boundaries: List[Tuple[int, int]] = []
        self.exon_bp: int = 0
        self.coding_segments: List[CodingSegment] = []
        self.coding_bp: int = 0

    def add_exon(self, start: int, end: int) -> None:
        self.exons.append((start, end))
        if start < self.span_start:
            self.span_start = start
        if end > self.span_end:
            self.span_end = end

    def add_cds(self, start: int, end: int) -> None:
        self.cds_raw.append((start, end))
        if start < self.span_start:
            self.span_start = start
        if end > self.span_end:
            self.span_end = end

    def finalize(self) -> None:
        if not self.exons:
            self.merged_exons = []
            self.intron_boundaries = []
            self.exon_bp = 0
        else:
            self.exons.sort(key=lambda x: (x[0], x[1]))
            self.merged_exons = merge_intervals(self.exons)
            self.exon_bp = intervals_length(self.exons)
            self.intron_boundaries = []
            for first, second in zip(self.exons, self.exons[1:]):
                prev_end = first[1]
                next_start = second[0]
                if prev_end < next_start:
                    self.intron_boundaries.append((prev_end, next_start))

        if not self.cds_raw:
            self.coding_segments = []
            self.coding_bp = 0
            return

        self.coding_bp = sum(end - start + 1 for start, end in self.cds_raw)
        cds_sorted = sorted(self.cds_raw, key=lambda x: (x[0], x[1]))
        segments: List[CodingSegment] = []

        if self.strand == "+":
            offset = 0
            for start, end in cds_sorted:
                length = end - start + 1
                phase_start = offset % 3
                phase_end = (phase_start + length - 1) % 3
                segments.append(CodingSegment(start, end, self.strand, phase_start, phase_end))
                offset += length
        elif self.strand == "-":
            offset = 0
            temp: List[CodingSegment] = []
            for start, end in reversed(cds_sorted):
                length = end - start + 1
                phase_end = offset % 3
                phase_start = (phase_end + length - 1) % 3
                temp.append(CodingSegment(start, end, self.strand, phase_start, phase_end))
                offset += length
            segments = list(sorted(temp, key=lambda x: (x.start, x.end)))
        else:
            # Unknown strand: still record segments but without reliable phase
            segments = [
                CodingSegment(start, end, self.strand, 0, 0)
                for start, end in cds_sorted
            ]
        self.coding_segments = segments


class Gene:
    __slots__ = (
        "id",
        "seqid",
        "strand",
        "start",
        "end",
        "transcripts",
        "merged_exons",
        "exon_bp",
    )

    def __init__(self, gene_id: str, seqid: str, strand: str, start: Optional[int] = None, end: Optional[int] = None) -> None:
        self.id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.start = start if start is not None else sys.maxsize
        self.end = end if end is not None else -sys.maxsize
        self.transcripts: Dict[str, Transcript] = {}
        self.merged_exons: List[Interval] = []
        self.exon_bp: int = 0

    def add_transcript(self, transcript: Transcript) -> None:
        self.transcripts[transcript.id] = transcript
        if transcript.span_start != sys.maxsize:
            self.start = min(self.start, transcript.span_start)
        if transcript.span_end != -sys.maxsize:
            self.end = max(self.end, transcript.span_end)

    def update_bounds(self, start: int, end: int) -> None:
        self.start = min(self.start, start)
        self.end = max(self.end, end)

    def finalize(self) -> None:
        exon_intervals: List[Interval] = []
        for transcript in self.transcripts.values():
            exon_intervals.extend(transcript.exons)
        self.merged_exons = merge_intervals(exon_intervals)
        self.exon_bp = intervals_length(self.merged_exons)
        if self.start == sys.maxsize:
            if self.merged_exons:
                self.start = self.merged_exons[0][0]
            else:
                self.start = 0
        if self.end == -sys.maxsize:
            if self.merged_exons:
                self.end = self.merged_exons[-1][1]
            else:
                self.end = 0


class Annotation:
    def __init__(self, name: str) -> None:
        self.name = name
        self.genes_by_id: Dict[str, Gene] = {}
        self.transcripts_by_id: Dict[str, Transcript] = {}
        self.genes_by_seqid: Dict[str, List[Gene]] = collections.defaultdict(list)

    def add_gene(self, gene: Gene) -> None:
        self.genes_by_id[gene.id] = gene
        self.genes_by_seqid[gene.seqid].append(gene)

    def finalize(self) -> None:
        for transcript in self.transcripts_by_id.values():
            transcript.finalize()
        for gene in self.genes_by_id.values():
            gene.finalize()
        for seqid in list(self.genes_by_seqid.keys()):
            self.genes_by_seqid[seqid] = sorted(self.genes_by_seqid[seqid], key=lambda g: g.start)


# ---------------------------------------------------------------------------
# Parsing logic
# ---------------------------------------------------------------------------

TRANSCRIPT_TYPES = {
    "mRNA",
    "transcript",
    "ncRNA",
    "lnc_RNA",
    "miRNA",
    "rRNA",
    "tRNA",
    "snRNA",
    "snoRNA",
    "primary_transcript",
    "pseudogenic_transcript",
}


def load_annotation(path: str, name: Optional[str] = None) -> Annotation:
    if name is None:
        base = os.path.basename(path)
        name = os.path.splitext(base)[0]
    annotation = Annotation(name)

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start_s, end_s, score, strand, phase, attr_raw = parts
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            attrs = parse_attributes(attr_raw)
            feature_id = attrs.get("ID")
            parents_raw = attrs.get("Parent")
            parents = parents_raw.split(",") if parents_raw else []

            if feature_type == "gene" and feature_id:
                gene = annotation.genes_by_id.get(feature_id)
                if gene is None:
                    gene = Gene(feature_id, seqid, strand, start, end)
                    annotation.add_gene(gene)
                else:
                    gene.update_bounds(start, end)
                    gene.seqid = seqid
                    gene.strand = strand if strand != "." else gene.strand
                continue

            if feature_type in TRANSCRIPT_TYPES and feature_id:
                if not parents:
                    continue
                gene_id = parents[0]
                gene = annotation.genes_by_id.get(gene_id)
                if gene is None:
                    gene = Gene(gene_id, seqid, strand, start, end)
                    annotation.add_gene(gene)
                transcript = annotation.transcripts_by_id.get(feature_id)
                if transcript is None:
                    transcript = Transcript(feature_id, gene_id, seqid, strand, start, end)
                    annotation.transcripts_by_id[feature_id] = transcript
                    gene.add_transcript(transcript)
                else:
                    gene.add_transcript(transcript)
                    transcript.gene_id = gene_id
                    transcript.seqid = seqid
                    transcript.strand = strand if strand != "." else transcript.strand
                    if start < transcript.span_start:
                        transcript.span_start = start
                    if end > transcript.span_end:
                        transcript.span_end = end
                continue

            if feature_type == "exon" and parents:
                for parent in parents:
                    transcript = annotation.transcripts_by_id.get(parent)
                    if transcript is None:
                        # Transcript feature may be missing; create stub
                        transcript = Transcript(parent, parent, seqid, strand, start, end)
                        annotation.transcripts_by_id[parent] = transcript
                    transcript.add_exon(start, end)
                continue

            if feature_type == "CDS" and parents:
                for parent in parents:
                    transcript = annotation.transcripts_by_id.get(parent)
                    if transcript is None:
                        transcript = Transcript(parent, parent, seqid, strand, start, end)
                        annotation.transcripts_by_id[parent] = transcript
                    transcript.add_cds(start, end)
                continue

    annotation.finalize()
    return annotation


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

@dataclass
class TranscriptComparison:
    transcript_a: Transcript
    transcript_b: Transcript
    stats: Dict[str, int]


@dataclass
class GeneComparison:
    gene_a: Gene
    gene_b: Gene
    transcript_pairs: List[TranscriptComparison]
    stats: Dict[str, int]


def transcript_overlap(transcript_a: Transcript, transcript_b: Transcript) -> bool:
    if transcript_a.seqid != transcript_b.seqid:
        return False
    return intervals_overlap(transcript_a.merged_exons, transcript_b.merged_exons)


def matching_exons(transcript_a: Transcript, transcript_b: Transcript) -> int:
    set_a = set(transcript_a.exons)
    set_b = set(transcript_b.exons)
    return len(set_a & set_b)


def matching_introns(transcript_a: Transcript, transcript_b: Transcript) -> int:
    set_a = set(transcript_a.intron_boundaries)
    set_b = set(transcript_b.intron_boundaries)
    return len(set_a & set_b)


def cds_overlap_same_phase(transcript_a: Transcript, transcript_b: Transcript) -> int:
    if not transcript_a.coding_segments or not transcript_b.coding_segments:
        return 0
    if transcript_a.strand != transcript_b.strand:
        return 0
    i = j = 0
    total = 0
    segs_a = transcript_a.coding_segments
    segs_b = transcript_b.coding_segments
    while i < len(segs_a) and j < len(segs_b):
        seg_a = segs_a[i]
        seg_b = segs_b[j]
        start = max(seg_a.start, seg_b.start)
        end = min(seg_a.end, seg_b.end)
        if start <= end:
            phase_a = seg_a.phase_at(start)
            phase_b = seg_b.phase_at(start)
            if phase_a is not None and phase_b is not None and phase_a == phase_b:
                total += end - start + 1
        if seg_a.end <= seg_b.end:
            i += 1
        else:
            j += 1
    return total


def compare_transcripts(transcript_a: Transcript, transcript_b: Transcript) -> TranscriptComparison:
    exons_a = len(transcript_a.exons)
    exons_b = len(transcript_b.exons)
    introns_a = len(transcript_a.intron_boundaries)
    introns_b = len(transcript_b.intron_boundaries)
    matched_exons = matching_exons(transcript_a, transcript_b)
    matched_introns = matching_introns(transcript_a, transcript_b)
    bp_a = transcript_a.exon_bp
    bp_b = transcript_b.exon_bp
    bp_overlap = intersection_length(transcript_a.merged_exons, transcript_b.merged_exons)
    cds_bp_a = transcript_a.coding_bp
    cds_bp_b = transcript_b.coding_bp
    cds_overlap = cds_overlap_same_phase(transcript_a, transcript_b)
    stats = {
        "exonsA": exons_a,
        "exonsB": exons_b,
        "match_exons": matched_exons,
        "intronsA": introns_a,
        "intronsB": introns_b,
        "match_introns": matched_introns,
        "bpA": bp_a,
        "bpB": bp_b,
        "bp_overlap": bp_overlap,
        "cds_bpA": cds_bp_a,
        "cds_bpB": cds_bp_b,
        "cds_bp_overlap_same_phase": cds_overlap,
    }
    return TranscriptComparison(transcript_a, transcript_b, stats)


def compare_gene_pair(gene_a: Gene, gene_b: Gene) -> GeneComparison:
    transcript_pairs: List[TranscriptComparison] = []
    matching_transcript_pairs = 0
    for transcript_a in gene_a.transcripts.values():
        for transcript_b in gene_b.transcripts.values():
            if not transcript_overlap(transcript_a, transcript_b):
                continue
            comparison = compare_transcripts(transcript_a, transcript_b)
            transcript_pairs.append(comparison)
            if (
                transcript_a.exons == transcript_b.exons
                and transcript_a.coding_bp == transcript_b.coding_bp
                and transcript_a.coding_bp == comparison.stats["cds_bp_overlap_same_phase"]
                and comparison.stats["bp_overlap"] == transcript_a.exon_bp == transcript_b.exon_bp
            ):
                matching_transcript_pairs += 1

    exon_overlap_bp = intersection_length(gene_a.merged_exons, gene_b.merged_exons)
    gene_overlap_bp = max(0, min(gene_a.end, gene_b.end) - max(gene_a.start, gene_b.start) + 1)
    stats = {
        "txA": len(gene_a.transcripts),
        "txB": len(gene_b.transcripts),
        "overlap_tx_pairs": len(transcript_pairs),
        "matching_tx_pairs": matching_transcript_pairs,
        "exon_bp_A": gene_a.exon_bp,
        "exon_bp_B": gene_b.exon_bp,
        "exon_bp_overlap": exon_overlap_bp,
        "gene_bp_overlap": gene_overlap_bp,
    }
    strand_match = 1 if gene_a.strand == gene_b.strand else 0
    stats["strand_match"] = strand_match
    return GeneComparison(gene_a, gene_b, transcript_pairs, stats)


def genes_overlap(gene_a: Gene, gene_b: Gene) -> bool:
    if gene_a.seqid != gene_b.seqid:
        return False
    if gene_a.end < gene_b.start or gene_b.end < gene_a.start:
        return False
    return intervals_overlap(gene_a.merged_exons, gene_b.merged_exons)


def generate_gene_pairs(annotation_a: Annotation, annotation_b: Annotation) -> Iterator[Tuple[Gene, Gene]]:
    for seqid in sorted(set(annotation_a.genes_by_seqid) & set(annotation_b.genes_by_seqid)):
        genes_a = annotation_a.genes_by_seqid[seqid]
        genes_b = annotation_b.genes_by_seqid[seqid]
        j = 0
        for gene_a in genes_a:
            while j < len(genes_b) and genes_b[j].end < gene_a.start:
                j += 1
            k = j
            while k < len(genes_b) and genes_b[k].start <= gene_a.end:
                gene_b = genes_b[k]
                if genes_overlap(gene_a, gene_b):
                    yield gene_a, gene_b
                k += 1


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def stats_to_string(stats: Dict[str, int]) -> str:
    return ";".join(f"{key}={value}" for key, value in stats.items())


def write_gene_comparisons(
    annotation_a: Annotation,
    annotation_b: Annotation,
    out_handle,
) -> None:
    for gene_a, gene_b in generate_gene_pairs(annotation_a, annotation_b):
        comparison = compare_gene_pair(gene_a, gene_b)
        if not comparison.transcript_pairs:
            continue
        gene_stats = comparison.stats.copy()
        gene_stats.setdefault("strand_match", 1 if gene_a.strand == gene_b.strand else 0)
        out_handle.write(
            "\t".join(
                [
                    "gene",
                    annotation_a.name,
                    gene_a.id,
                    ".",
                    annotation_b.name,
                    gene_b.id,
                    ".",
                    stats_to_string(gene_stats),
                ]
            )
            + "\n"
        )
        for tx_comp in comparison.transcript_pairs:
            tx_stats = tx_comp.stats
            out_handle.write(
                "\t".join(
                    [
                        "transcript",
                        annotation_a.name,
                        tx_comp.transcript_a.gene_id,
                        tx_comp.transcript_a.id,
                        annotation_b.name,
                        tx_comp.transcript_b.gene_id,
                        tx_comp.transcript_b.id,
                        stats_to_string(tx_stats),
                    ]
                )
                + "\n"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare GFF3 annotations within the same assembly")
    parser.add_argument("gff", nargs="+", help="Sorted GFF3 annotation files")
    parser.add_argument(
        "-l",
        "--labels",
        nargs="+",
        help="Annotation labels (same order as GFF files)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file (TSV). Defaults to stdout",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.labels and len(args.labels) != len(args.gff):
        sys.stderr.write("Error: number of labels must match number of GFF files\n")
        return 1

    labels = args.labels if args.labels else [None] * len(args.gff)
    annotations = [load_annotation(path, label) for path, label in zip(args.gff, labels)]

    out_handle = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    try:
        out_handle.write("feature\tannA\tgeneA\ttxA\tannB\tgeneB\ttxB\tstats\n")
        for idx_a, annotation_a in enumerate(annotations):
            for idx_b, annotation_b in enumerate(annotations):
                if idx_a == idx_b:
                    continue
                write_gene_comparisons(annotation_a, annotation_b, out_handle)
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
