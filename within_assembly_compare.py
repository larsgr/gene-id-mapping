#!/usr/bin/env python3
"""Streamed within-assembly GFF3 comparison parser.

This implementation compares sorted GFF3 annotation files for the same genome
assembly. It streams gene records from the input, computes transcript-level
metrics, rolls them up to per-gene summaries, and writes the results while
keeping only the active locus in memory. Transcript rows are optional and can
be interleaved with gene summaries on request.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple, Union


# ---------------------------------------------------------------------------
# Types and constants
# ---------------------------------------------------------------------------

Interval = Tuple[int, int]

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

GENE_TYPES = {
    "gene",              # Standard gene feature (NCBI & Ensembl protein-coding + some non-coding)
    "ncRNA_gene",        # Ensembl non-coding RNA genes (lncRNA, rRNA, snRNA, snoRNA, miRNA, etc.)
    "pseudogene",        # Pseudogenes (both NCBI & Ensembl)
    "V_gene_segment",    # Immunoglobulin/T-cell receptor V segments
    "D_gene_segment",    # Immunoglobulin/T-cell receptor D segments
    "J_gene_segment",    # Immunoglobulin/T-cell receptor J segments
}

SEVERITY_ORDER = {
    "NotMapped": 0,
    "Red": 1,
    "Yellow": 2,
    "Green": 3,
}

# Threshold constants (can be surfaced as CLI options later if needed)
GREEN_JACCARD_CDS = 0.98
GREEN_JACCARD_CDS_SECONDARY = 0.95
GREEN_JUNCTION_CDS = 0.95

GREEN_NONCODING_JACCARD = 0.90
GREEN_NONCODING_JUNCTION = 0.95
GREEN_MONOEXONIC_JACCARD = 0.95
GREEN_MONOEXONIC_RATIO_DELTA = 0.10

YELLOW_JACCARD_CDS_LOW = 0.50
YELLOW_JUNCTION_CDS_LOW = 0.50
YELLOW_JACCARD_CDS_MIN = 0.40
YELLOW_JACCARD_NONCODING_LOW = 0.50
YELLOW_JUNCTION_NONCODING_LOW = 0.50
YELLOW_MONOEXONIC_JACCARD_LOW = 0.50

RED_JACCARD_EXON_LOW = 0.10


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------


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


def merge_intervals(intervals: Iterable[Interval]) -> List[Interval]:
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    if not sorted_intervals:
        return []
    merged: List[List[int]] = [[sorted_intervals[0][0], sorted_intervals[0][1]]]
    for start, end in sorted_intervals[1:]:
        last = merged[-1]
        if start <= last[1] + 1:
            if end > last[1]:
                last[1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def intervals_length(intervals: Iterable[Interval]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def intervals_overlap(a: Sequence[Interval], b: Sequence[Interval]) -> bool:
    i = j = 0
    while i < len(a) and j < len(b):
        a_start, a_end = a[i]
        b_start, b_end = b[j]
        if a_end < b_start:
            i += 1
        elif b_end < a_start:
            j += 1
        else:
            return True
    return False


def intersection_length(a: Sequence[Interval], b: Sequence[Interval]) -> int:
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        a_start, a_end = a[i]
        b_start, b_end = b[j]
        start = max(a_start, b_start)
        end = min(a_end, b_end)
        if start <= end:
            total += end - start + 1
        if a_end <= b_end:
            i += 1
        else:
            j += 1
    return total


def union_length(a: Sequence[Interval], b: Sequence[Interval]) -> int:
    if not a and not b:
        return 0
    return intervals_length(merge_intervals(list(a) + list(b)))


def f1_score(matches: int, count_a: int, count_b: int) -> float:
    if count_a == 0 and count_b == 0:
        return 1.0
    denom = count_a + count_b
    if denom == 0:
        return 0.0
    return 2.0 * matches / denom


def float_to_str(value: float) -> str:
    return f"{value:.4f}" if not math.isnan(value) else "nan"


def stats_to_string(stats: Dict[str, Union[int, float, str]]) -> str:
    parts: List[str] = []
    for key in sorted(stats.keys()):
        value = stats[key]
        if isinstance(value, float):
            parts.append(f"{key}={float_to_str(value)}")
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts)


def severity_max(a: str, b: str) -> str:
    return a if SEVERITY_ORDER[a] >= SEVERITY_ORDER[b] else b


def severity_min(a: str, b: str) -> str:
    return a if SEVERITY_ORDER[a] <= SEVERITY_ORDER[b] else b


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
        "merged_cds",
        "intron_boundaries",
        "cds_intron_boundaries",
        "exon_bp",
        "coding_segments",
        "coding_bp",
        "monoexonic",
    )

    def __init__(self, transcript_id: str, gene_id: str, seqid: str, strand: str) -> None:
        self.id = transcript_id
        self.gene_id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.span_start = sys.maxsize
        self.span_end = -sys.maxsize
        self.exons: List[Interval] = []
        self.cds_raw: List[Interval] = []
        self.merged_exons: List[Interval] = []
        self.merged_cds: List[Interval] = []
        self.intron_boundaries: List[Tuple[int, int]] = []
        self.cds_intron_boundaries: List[Tuple[int, int]] = []
        self.exon_bp: int = 0
        self.coding_segments: List[CodingSegment] = []
        self.coding_bp: int = 0
        self.monoexonic = True

    def update_span(self, start: int, end: int) -> None:
        if start < self.span_start:
            self.span_start = start
        if end > self.span_end:
            self.span_end = end

    def add_exon(self, start: int, end: int) -> None:
        self.exons.append((start, end))
        self.update_span(start, end)

    def add_cds(self, start: int, end: int) -> None:
        self.cds_raw.append((start, end))
        self.update_span(start, end)

    def finalize(self) -> None:
        if not self.exons:
            self.merged_exons = []
            self.intron_boundaries = []
            self.exon_bp = 0
            self.monoexonic = True
        else:
            self.exons.sort(key=lambda x: (x[0], x[1]))
            self.merged_exons = merge_intervals(self.exons)
            self.exon_bp = intervals_length(self.exons)
            introns: List[Tuple[int, int]] = []
            for first, second in zip(self.exons, self.exons[1:]):
                if first[1] < second[0]:
                    introns.append((first[1], second[0]))
            self.intron_boundaries = introns
            self.monoexonic = len(introns) == 0

        if not self.cds_raw:
            self.coding_segments = []
            self.coding_bp = 0
            self.cds_intron_boundaries = []
            self.merged_cds = []
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
            segments = [
                CodingSegment(start, end, self.strand, 0, 0)
                for start, end in cds_sorted
            ]

        self.coding_segments = segments
        self.merged_cds = merge_intervals(self.cds_raw)

        cds_introns: List[Tuple[int, int]] = []
        for first, second in zip(cds_sorted, cds_sorted[1:]):
            if first[1] < second[0]:
                cds_introns.append((first[1], second[0]))
        self.cds_intron_boundaries = cds_introns

    @property
    def is_coding(self) -> bool:
        return self.coding_bp > 0


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

    def __init__(self, gene_id: str, seqid: str, strand: str, start: int, end: int) -> None:
        self.id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts: Dict[str, Transcript] = {}
        self.merged_exons: List[Interval] = []
        self.exon_bp: int = 0

    def add_transcript(self, transcript: Transcript) -> None:
        self.transcripts[transcript.id] = transcript
        if transcript.span_start < self.start:
            self.start = transcript.span_start
        if transcript.span_end > self.end:
            self.end = transcript.span_end

    def update_bounds(self, start: int, end: int) -> None:
        if start < self.start:
            self.start = start
        if end > self.end:
            self.end = end

    def finalize(self) -> None:
        for transcript in self.transcripts.values():
            transcript.finalize()
        exon_intervals: List[Interval] = []
        for transcript in self.transcripts.values():
            exon_intervals.extend(transcript.merged_exons)
        self.merged_exons = merge_intervals(exon_intervals)
        self.exon_bp = intervals_length(self.merged_exons)


# ---------------------------------------------------------------------------
# Streaming parser
# ---------------------------------------------------------------------------


class GFFSortError(RuntimeError):
    pass


def stream_genes(path: str) -> Iterator[Gene]:
    with open(path, "r", encoding="utf-8") as handle:
        current_gene: Optional[Gene] = None
        transcripts: Dict[str, Transcript] = {}
        last_gene_key: Optional[Tuple[str, int]] = None

        def finalize_current_gene() -> Optional[Gene]:
            nonlocal current_gene, transcripts
            if current_gene is None:
                return None
            for transcript in transcripts.values():
                current_gene.add_transcript(transcript)
            current_gene.finalize()
            gene_to_return = current_gene
            current_gene = None
            transcripts = {}
            return gene_to_return

        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
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

            if feature_type in GENE_TYPES and feature_id:
                if last_gene_key is not None:
                    last_seqid, last_start = last_gene_key
                    if (seqid < last_seqid) or (seqid == last_seqid and start < last_start):
                        raise GFFSortError(
                            f"GFF not sorted: {feature_type} {feature_id} at {seqid}:{start} precedes {last_seqid}:{last_start}"
                        )
                last_gene_key = (seqid, start)
                gene_done = finalize_current_gene()
                if gene_done is not None:
                    yield gene_done
                current_gene = Gene(feature_id, seqid, strand, start, end)
                transcripts = {}
                continue

            if current_gene is None:
                # Ignore features before the first gene-like feature declaration
                continue

            if feature_type in TRANSCRIPT_TYPES and feature_id:
                if not parents:
                    raise RuntimeError(f"Transcript {feature_id} missing Parent at {seqid}:{start}")
                parent_gene = parents[0]
                if parent_gene != current_gene.id:
                    raise RuntimeError(
                        f"Transcript {feature_id} refers to gene {parent_gene} while parsing {current_gene.id}"
                    )
                transcript = transcripts.get(feature_id)
                if transcript is None:
                    transcript = Transcript(feature_id, current_gene.id, seqid, strand)
                    transcripts[feature_id] = transcript
                transcript.update_span(start, end)
                current_gene.update_bounds(start, end)
                continue

            if feature_type == "exon" and parents:
                for parent in parents:
                    transcript = transcripts.get(parent)
                    if transcript is None:
                        transcript = Transcript(parent, current_gene.id, seqid, strand)
                        transcripts[parent] = transcript
                    transcript.add_exon(start, end)
                    current_gene.update_bounds(start, end)
                continue

            if feature_type == "CDS" and parents:
                for parent in parents:
                    transcript = transcripts.get(parent)
                    if transcript is None:
                        transcript = Transcript(parent, current_gene.id, seqid, strand)
                        transcripts[parent] = transcript
                    transcript.add_cds(start, end)
                    current_gene.update_bounds(start, end)
                continue

        gene_done = finalize_current_gene()
        if gene_done is not None:
            yield gene_done


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------


def transcript_pair_overlap(transcript_a: Transcript, transcript_b: Transcript) -> bool:
    if transcript_a.seqid != transcript_b.seqid:
        return False
    return intervals_overlap(transcript_a.merged_exons, transcript_b.merged_exons)


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


def cds_overlap_total(transcript_a: Transcript, transcript_b: Transcript) -> int:
    return intersection_length(transcript_a.merged_cds, transcript_b.merged_cds)


@dataclass
class TranscriptPairResult:
    transcript_a: Transcript
    transcript_b: Transcript
    stats: Dict[str, Union[int, float, str]]
    classification: str
    category: str
    notes: Dict[str, Union[int, float, str]] = field(default_factory=dict)


@dataclass
class GenePairRecord:
    label_a: str
    gene_a: Gene
    label_b: str
    gene_b: Gene
    stats: Dict[str, Union[int, float, str]]
    classification: str
    transcript_pairs: List[TranscriptPairResult]
    notes: Dict[str, Union[int, float, str]] = field(default_factory=dict)
    pending: Set[str] = field(default_factory=lambda: {"A", "B"})


@dataclass
class ComparisonOptions:
    include_transcripts: bool
    antisense_red_threshold: float


def compute_transcript_stats(transcript_a: Transcript, transcript_b: Transcript) -> Dict[str, Union[int, float, str]]:
    stats: Dict[str, Union[int, float, str]] = {}
    exons_a = len(transcript_a.exons)
    exons_b = len(transcript_b.exons)
    introns_a = len(transcript_a.intron_boundaries)
    introns_b = len(transcript_b.intron_boundaries)
    matched_exons = len(set(transcript_a.exons) & set(transcript_b.exons))
    matched_introns = len(set(transcript_a.intron_boundaries) & set(transcript_b.intron_boundaries))
    matched_cds_introns = len(set(transcript_a.cds_intron_boundaries) & set(transcript_b.cds_intron_boundaries))

    exon_overlap_bp = intersection_length(transcript_a.merged_exons, transcript_b.merged_exons)
    exon_union_bp = union_length(transcript_a.merged_exons, transcript_b.merged_exons)
    cds_overlap_bp = cds_overlap_total(transcript_a, transcript_b)
    cds_union_bp = union_length(transcript_a.merged_cds, transcript_b.merged_cds)
    cds_phase_overlap_bp = cds_overlap_same_phase(transcript_a, transcript_b)

    stats.update(
        {
            "exonsA": exons_a,
            "exonsB": exons_b,
            "intronsA": introns_a,
            "intronsB": introns_b,
            "match_exons": matched_exons,
            "match_introns": matched_introns,
            "cds_match_introns": matched_cds_introns,
            "bpA": transcript_a.exon_bp,
            "bpB": transcript_b.exon_bp,
            "bp_overlap": exon_overlap_bp,
            "bp_union": exon_union_bp,
            "cds_bpA": transcript_a.coding_bp,
            "cds_bpB": transcript_b.coding_bp,
            "cds_bp_overlap": cds_overlap_bp,
            "cds_bp_union": cds_union_bp,
            "cds_bp_overlap_same_phase": cds_phase_overlap_bp,
            "strand_agree": 1 if transcript_a.strand == transcript_b.strand else 0,
            "monoexonicA": 1 if transcript_a.monoexonic else 0,
            "monoexonicB": 1 if transcript_b.monoexonic else 0,
        }
    )

    stats["jaccard_exon"] = (exon_overlap_bp / exon_union_bp) if exon_union_bp else 0.0
    stats["jaccard_cds_phase"] = (
        cds_phase_overlap_bp / cds_union_bp if cds_union_bp else 0.0
    )
    stats["junction_f1_all"] = f1_score(matched_introns, introns_a, introns_b)
    stats["junction_f1_cds"] = f1_score(matched_cds_introns, len(transcript_a.cds_intron_boundaries), len(transcript_b.cds_intron_boundaries))
    stats["shared_introns_fraction"] = (
        matched_introns / max(introns_a, introns_b) if max(introns_a, introns_b) else 1.0
    )
    if stats["strand_agree"]:
        stats["antisense_overlap"] = 0.0
    else:
        stats["antisense_overlap"] = (exon_overlap_bp / exon_union_bp) if exon_union_bp else 0.0

    return stats


def classify_transcript_pair(
    transcript_a: Transcript,
    transcript_b: Transcript,
    stats: Dict[str, Union[int, float, str]],
) -> Tuple[str, str, Dict[str, Union[int, float, str]]]:
    notes: Dict[str, Union[int, float, str]] = {}

    strand_agree = bool(stats.get("strand_agree", 0))
    if not strand_agree:
        notes["reason"] = "antisense"
        return "NotMapped", "antisense", notes

    is_coding_a = transcript_a.is_coding
    is_coding_b = transcript_b.is_coding

    if is_coding_a and not is_coding_b or is_coding_b and not is_coding_a:
        notes["reason"] = "discordant_biotype"
        return "NotMapped", "mixed", notes

    jaccard_exon = float(stats.get("jaccard_exon", 0.0))
    jaccard_cds_phase = float(stats.get("jaccard_cds_phase", 0.0))
    junction_all = float(stats.get("junction_f1_all", 0.0))
    junction_cds = float(stats.get("junction_f1_cds", 0.0))
    mono_a = bool(stats.get("monoexonicA", 0))
    mono_b = bool(stats.get("monoexonicB", 0))
    bp_a = int(stats.get("bpA", 0))
    bp_b = int(stats.get("bpB", 0))
    shared_introns = int(stats.get("match_introns", 0))

    if is_coding_a and is_coding_b:
        category = "coding"
        if jaccard_cds_phase >= GREEN_JACCARD_CDS or (
            junction_cds >= GREEN_JUNCTION_CDS and jaccard_cds_phase >= GREEN_JACCARD_CDS_SECONDARY
        ):
            return "Green", category, notes
        if (
            YELLOW_JACCARD_CDS_LOW <= jaccard_cds_phase < GREEN_JACCARD_CDS
            or (
                junction_cds >= YELLOW_JUNCTION_CDS_LOW
                and jaccard_cds_phase >= YELLOW_JACCARD_CDS_MIN
            )
        ):
            return "Yellow", category, notes
        if jaccard_exon >= RED_JACCARD_EXON_LOW or (shared_introns > 0 and junction_all < YELLOW_JUNCTION_NONCODING_LOW):
            notes["warning"] = "weak_coding_support"
            return "Red", category, notes
        return "NotMapped", category, notes

    category = "noncoding"
    if not mono_a and not mono_b:
        if junction_all >= GREEN_NONCODING_JUNCTION or (
            jaccard_exon >= GREEN_NONCODING_JACCARD and junction_all >= 0.90
        ):
            return "Green", category, notes
    else:
        if jaccard_exon >= GREEN_MONOEXONIC_JACCARD:
            longer = max(bp_a, bp_b) if max(bp_a, bp_b) else 1
            shorter = min(bp_a, bp_b)
            ratio = shorter / longer if longer else 0.0
            if ratio >= 1.0 - GREEN_MONOEXONIC_RATIO_DELTA:
                return "Green", category, notes
    if (
        YELLOW_JACCARD_NONCODING_LOW <= jaccard_exon < GREEN_NONCODING_JACCARD
        or YELLOW_JUNCTION_NONCODING_LOW <= junction_all < GREEN_NONCODING_JUNCTION
    ):
        return "Yellow", category, notes
    if jaccard_exon >= RED_JACCARD_EXON_LOW or (shared_introns > 0 and junction_all < YELLOW_JUNCTION_NONCODING_LOW):
        notes["warning"] = "weak_noncoding_support"
        return "Red", category, notes
    return "NotMapped", category, notes


def classify_gene_pair(
    gene_a: Gene,
    gene_b: Gene,
    transcript_pairs: List[TranscriptPairResult],
    options: ComparisonOptions,
) -> Tuple[str, Dict[str, Union[int, float, str]], Dict[str, Union[int, float, str]]]:
    notes: Dict[str, Union[int, float, str]] = {}
    stats: Dict[str, Union[int, float, str]] = {}

    exon_overlap_bp = intersection_length(gene_a.merged_exons, gene_b.merged_exons)
    exon_union_bp = union_length(gene_a.merged_exons, gene_b.merged_exons)
    stats["overlap_bp"] = exon_overlap_bp
    stats["exon_union_bp"] = exon_union_bp
    stats["strand_match"] = 1 if gene_a.strand == gene_b.strand else 0

    if gene_a.strand != gene_b.strand and gene_a.strand != "." and gene_b.strand != ".":
        antisense_overlap = (exon_overlap_bp / exon_union_bp) if exon_union_bp else 0.0
        stats["antisense_overlap"] = antisense_overlap
        notes["antisense_conflict"] = 1
        if antisense_overlap >= options.antisense_red_threshold:
            stats["class_initial"] = "Red"
            return "Red", stats, notes
        stats["class_initial"] = "NotMapped"
        return "NotMapped", stats, notes

    best_coding = "NotMapped"
    best_noncoding = "NotMapped"
    best_cds_jaccard = 0.0
    best_exon_jaccard = 0.0
    discordant_biotype = False

    for pair in transcript_pairs:
        stats_local = pair.stats
        jaccard_exon = float(stats_local.get("jaccard_exon", 0.0))
        jaccard_cds = float(stats_local.get("jaccard_cds_phase", 0.0))
        best_exon_jaccard = max(best_exon_jaccard, jaccard_exon)
        best_cds_jaccard = max(best_cds_jaccard, jaccard_cds)
        if pair.category == "coding":
            best_coding = severity_max(best_coding, pair.classification)
        elif pair.category == "noncoding":
            best_noncoding = severity_max(best_noncoding, pair.classification)
        elif pair.category == "mixed":
            discordant_biotype = True

    stats["best_coding"] = best_coding
    stats["best_noncoding"] = best_noncoding
    stats["best_jaccard_exon"] = best_exon_jaccard
    stats["best_jaccard_cds_phase"] = best_cds_jaccard

    if not transcript_pairs:
        stats["class_initial"] = "NotMapped"
        return "NotMapped", stats, notes

    combined = "NotMapped"
    for candidate in (best_coding, best_noncoding):
        combined = severity_max(combined, candidate)
    stats["class_initial"] = combined

    if discordant_biotype and SEVERITY_ORDER[combined] > SEVERITY_ORDER["Yellow"]:
        combined = "Yellow"
        notes["discordant_biotype"] = 1

    return combined, stats, notes


def downgrade_class(current: str, target: str) -> str:
    return severity_min(current, target)


def analyse_gene_pair(
    label_a: str,
    gene_a: Gene,
    label_b: str,
    gene_b: Gene,
    options: ComparisonOptions,
) -> Optional[GenePairRecord]:
    transcript_pairs: List[TranscriptPairResult] = []
    for transcript_a in gene_a.transcripts.values():
        for transcript_b in gene_b.transcripts.values():
            if not transcript_pair_overlap(transcript_a, transcript_b):
                continue
            stats = compute_transcript_stats(transcript_a, transcript_b)
            classification, category, notes = classify_transcript_pair(transcript_a, transcript_b, stats)
            stats["class"] = classification
            result = TranscriptPairResult(transcript_a, transcript_b, stats, classification, category, notes)
            transcript_pairs.append(result)

    if not transcript_pairs:
        return None

    classification, stats, notes = classify_gene_pair(gene_a, gene_b, transcript_pairs, options)
    stats["class"] = classification
    return GenePairRecord(label_a, gene_a, label_b, gene_b, stats, classification, transcript_pairs, notes)


# ---------------------------------------------------------------------------
# Streaming pairwise comparison
# ---------------------------------------------------------------------------


class RecordBuffer:
    def __init__(self, out_handle, options: ComparisonOptions) -> None:
        self.out_handle = out_handle
        self.options = options
        self.records_by_gene_a: Dict[str, List[GenePairRecord]] = {}
        self.records_by_gene_b: Dict[str, List[GenePairRecord]] = {}
        self.partner_sets_a: Dict[str, Set[str]] = {}
        self.partner_sets_b: Dict[str, Set[str]] = {}

    def add_record(self, record: GenePairRecord) -> None:
        self.records_by_gene_a.setdefault(record.gene_a.id, []).append(record)
        self.records_by_gene_b.setdefault(record.gene_b.id, []).append(record)
        self.partner_sets_a.setdefault(record.gene_a.id, set()).add(record.gene_b.id)
        self.partner_sets_b.setdefault(record.gene_b.id, set()).add(record.gene_a.id)

    def finalize_gene(self, side: str, gene_id: str) -> None:
        if side == "A":
            partner_set = self.partner_sets_a.get(gene_id, set())
            records = self.records_by_gene_a.pop(gene_id, [])
        else:
            partner_set = self.partner_sets_b.get(gene_id, set())
            records = self.records_by_gene_b.pop(gene_id, [])

        partner_count = len(partner_set)
        split_flag = 1 if partner_count > 1 else 0
        for record in records:
            if side == "A":
                record.stats["partner_countA"] = partner_count
                if split_flag:
                    record.notes["split_or_merge_A"] = partner_count
                    record.classification = downgrade_class(record.classification, "Yellow")
                    record.stats["class"] = record.classification
            else:
                record.stats["partner_countB"] = partner_count
                if split_flag:
                    record.notes["split_or_merge_B"] = partner_count
                    record.classification = downgrade_class(record.classification, "Yellow")
                    record.stats["class"] = record.classification
            record.pending.discard(side)
            if not record.pending:
                self.write_record(record)

        if side == "A":
            self.partner_sets_a.pop(gene_id, None)
        else:
            self.partner_sets_b.pop(gene_id, None)

    def flush_remaining(self) -> None:
        for records in itertools.chain(self.records_by_gene_a.values(), self.records_by_gene_b.values()):
            for record in records:
                record.pending.clear()
                self.write_record(record)
        self.records_by_gene_a.clear()
        self.records_by_gene_b.clear()

    def write_record(self, record: GenePairRecord) -> None:
        stats_combined = dict(record.stats)
        for key, value in record.notes.items():
            stats_combined[f"note_{key}"] = value
        stats_str = stats_to_string(stats_combined)
        line = "\t".join(
            [
                "gene",
                record.label_a,
                record.gene_a.id,
                ".",
                record.label_b,
                record.gene_b.id,
                ".",
                stats_str,
            ]
        )
        self.out_handle.write(line + "\n")
        if self.options.include_transcripts:
            for tx_pair in record.transcript_pairs:
                tx_stats = dict(tx_pair.stats)
                for key, value in tx_pair.notes.items():
                    tx_stats[f"note_{key}"] = value
                tx_stats_str = stats_to_string(tx_stats)
                tx_line = "\t".join(
                    [
                        "transcript",
                        record.label_a,
                        tx_pair.transcript_a.gene_id,
                        tx_pair.transcript_a.id,
                        record.label_b,
                        tx_pair.transcript_b.gene_id,
                        tx_pair.transcript_b.id,
                        tx_stats_str,
                    ]
                )
                self.out_handle.write(tx_line + "\n")


def compare_pair_streaming(
    path_a: str,
    label_a: str,
    path_b: str,
    label_b: str,
    out_handle,
    options: ComparisonOptions,
) -> None:
    stream_a = stream_genes(path_a)
    stream_b = stream_genes(path_b)
    buffer = RecordBuffer(out_handle, options)

    active_b: List[Gene] = []
    next_b: Optional[Gene] = None

    try:
        next_b = next(stream_b)
    except StopIteration:
        next_b = None

    for gene_a in stream_a:
        # Remove B genes that can no longer overlap with current or future A genes
        remaining_b: List[Gene] = []
        for gene_b in active_b:
            if gene_b.seqid < gene_a.seqid or (gene_b.seqid == gene_a.seqid and gene_b.end < gene_a.start):
                buffer.finalize_gene("B", gene_b.id)
            else:
                remaining_b.append(gene_b)
        active_b = remaining_b

        # Pull in B genes that might overlap gene_a
        while next_b is not None and (
            next_b.seqid < gene_a.seqid
            or (next_b.seqid == gene_a.seqid and next_b.start <= gene_a.end)
        ):
            active_b.append(next_b)
            try:
                next_b = next(stream_b)
            except StopIteration:
                next_b = None

        # Compare gene_a against active gene_b entries
        for gene_b in active_b:
            if gene_b.seqid != gene_a.seqid:
                continue
            if gene_b.start > gene_a.end or gene_b.end < gene_a.start:
                continue
            record = analyse_gene_pair(label_a, gene_a, label_b, gene_b, options)
            if record is not None:
                buffer.add_record(record)

        buffer.finalize_gene("A", gene_a.id)

    # All genes from A processed; remaining B genes cannot overlap future genes
    for gene_b in active_b:
        buffer.finalize_gene("B", gene_b.id)

    for gene_b in stream_b:
        buffer.finalize_gene("B", gene_b.id)

    buffer.flush_remaining()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare sorted GFF3 annotations within the same assembly")
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
    parser.add_argument(
        "--include-transcripts",
        action="store_true",
        help="Include transcript-level mappings interleaved with gene summaries",
    )
    parser.add_argument(
        "--antisense-red-threshold",
        type=float,
        default=0.5,
        help="Fraction of exon union that triggers a Red warning for antisense overlaps",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if len(args.gff) < 2:
        sys.stderr.write("Error: provide at least two GFF files\n")
        return 1
    if args.labels and len(args.labels) != len(args.gff):
        sys.stderr.write("Error: number of labels must match number of GFF files\n")
        return 1

    labels = args.labels if args.labels else [None] * len(args.gff)

    out_handle = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    options = ComparisonOptions(
        include_transcripts=args.include_transcripts,
        antisense_red_threshold=args.antisense_red_threshold,
    )

    try:
        out_handle.write("feature\tannA\tgeneA\ttxA\tannB\tgeneB\ttxB\tstats\n")
        for idx_a, (path_a, label_a) in enumerate(zip(args.gff, labels)):
            label_a = label_a or f"ann{idx_a+1}"
            for idx_b in range(idx_a + 1, len(args.gff)):
                path_b = args.gff[idx_b]
                label_b = (labels[idx_b] if labels[idx_b] else f"ann{idx_b+1}")
                compare_pair_streaming(path_a, label_a, path_b, label_b, out_handle, options)
    except GFFSortError as exc:
        sys.stderr.write(f"Error: {exc}\n")
        return 1
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
