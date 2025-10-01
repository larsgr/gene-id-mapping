# Comparison Script Design

This note refines the custom comparison script that will extend the
`within_assembly_compare.py` prototype described in `docs/experiments.md` and
`docs/AI-usage.md`. The existing parser already produces transcript-level
overlaps across annotations on the same assembly. The next phase is to keep the
per-transcript detail while rolling it up into per-gene mapping summaries that
drive the Green / Yellow / Red / NotMapped labels.

## Metrics to Compute Per Transcript Pair

The transcript pass still computes the raw measures. They drive both the
transcript report and the gene aggregation logic, so persist the raw counts and
derived fractions.

- `Jaccard_exon`: overlapping exon bp ÷ union of exon bp. Here the "union"
  means the combined genomic footprint of the two transcripts being compared.
- `Jaccard_CDS_phase`: overlapping CDS bp (same strand, same inferred codon
  phase) ÷ union of CDS bp for that transcript pair; ignore segments where one
  side lacks CDS.
- `Junction_F1_all`: F1-score over all splice junctions.
- `Junction_F1_CDS`: F1-score over CDS junctions only.
- `Strand_agree`: boolean flag; false means antisense.
- `Shared_introns`: absolute count and fraction of introns that match exactly.
- `Monoexonic`: true if transcript has no introns.
- `Has_split_merge`: flag plus partner counts when overlaps span multiple
  genes.
- `Antisense_overlap`: fraction of exon union overlapped on the opposite strand.

Keep every raw count (bp, intron counts, exon counts) alongside the fractions;
thresholds below are still provisional.

## Current Per-Transcript Flow

The within-assembly parser emits one record per overlapping transcript pair with
the metrics above. That output already distinguishes coding vs. non-coding
transcripts, handles antisense overlaps, and records CDS phase matches by
recomputing reading frame from CDS features.

## Transition to Per-Gene Mapping

Rolling up to genes introduces a second pass:

1. **Group transcript pairs by gene pair** (`annA/geneA` × `annB/geneB`).
2. **Classify transcript pairs** using the same Green / Yellow / Red thresholds
   that the script applies today.
3. **Select representative evidence** for the gene: a coding gene can use any
   coding transcript pair; non-coding genes use non-coding pairs. Record which
   transcript pair(s) justified the gene score.
4. **Aggregate downgrade triggers** (splits, merges, competing partners, etc.)
   across all transcripts for the gene before finalising the gene-level label.
5. **Emit both summaries**: one gene-level record per partner pair plus the
   existing transcript-level rows.

This extra pass can reuse the in-memory structures already built for overlap
detection in `within_assembly_compare.py`: extend the gene object so it tracks
its transcript mappings and exposes helper methods for the aggregation logic.

### Gene Exon Footprint vs. Transcript Union

Per-gene statistics may still need a "gene exon footprint" for each annotation:
the union of exons across every transcript isoform for a gene. This concept is
distinct from the per-transcript "union of exon bp" used in the Jaccard metrics
above, which considers only the two transcripts being compared. When aggregating
to genes, make the distinction explicit in code and naming so downstream
consumers do not confuse the gene-level footprint with the transcript-pair
union.

## Gene-Level Classification Rules

Representative transcript rules apply separately to coding vs. non-coding genes
before combining them into a single class for the gene pair.

### Green
- **Coding** (strand agreement not required for many-to-many):
  - `Jaccard_CDS_phase ≥ 0.98`, or
  - `Junction_F1_CDS ≥ 0.95` *and* `Jaccard_CDS_phase ≥ 0.95`.
  - Ignore UTR- or non-CDS-only differences.
- **Non-coding**:
  - If both are spliced: `Junction_F1_all ≥ 0.95` or the combination of
    `Jaccard_exon ≥ 0.90` and `Junction_F1_all ≥ 0.90`.
  - If either isoform is monoexonic: `Jaccard_exon ≥ 0.95` and length ratio in
    `[0.9, 1.1]`.
- Antisense mappings never qualify as Green.

### Yellow
- **Coding** (strand must agree):
  - `0.50 ≤ Jaccard_CDS_phase < 0.98`, or
  - `Junction_F1_CDS ≥ 0.5` with `Jaccard_CDS_phase ≥ 0.40`.
- **Non-coding**:
  - `0.50 ≤ Jaccard_exon < 0.90`, or
  - `0.50 ≤ Junction_F1_all < 0.95`.
- **Monoexonic pairs**: `0.50 ≤ Jaccard_exon < 0.95`.
- Antisense overlaps covering ≥ 50% of exon union can be reported as Yellow only
  if antisense families are explicitly allowed; otherwise flag `Conflict` and
  demote to NotMapped.

### Red
- Same strand but weak evidence:
  - `0.10 ≤ Jaccard_exon < 0.50`, or
  - At least one shared intron while `Junction_F1_all < 0.50`.
- Antisense overlaps with `Jaccard_exon ≥ 0.10` go to the notes as `Conflict` by
  default. They can be promoted to Red if the workflow decides to count them.

### NotMapped
- Different chromosomes, zero overlap, or metrics below the Red thresholds.
- Any antisense case when antisense mappings are disallowed.

## Downgrade Rules (Green → Yellow)

Gene-level downgrades fire if any transcript pair triggers them; record the
trigger in the notes and demote the gene summary:

- **Split / Merge detection**: overlaps with ≥ 2 distinct partners where each
  subpair is at least Yellow.
- **Competing partner within Δ**: another partner whose key metric is within the
  delta tolerance (e.g. `ΔJaccard_CDS_phase ≤ 0.02`).
- **Frame inconsistency**: mixed in-phase and out-of-phase transcript pairs.
- **Discordant biotypes**: coding vs. non-coding disagreement within the same
  locus.
- **Problematic introns**: shared CDS junction fraction < 1.0 when
  `Jaccard_CDS_phase` barely passes Green.
- **Extensive non-CDS divergence**: `Jaccard_exon < 0.70` despite a Green CDS
  verdict.

## Secondary Flags and Notes

Carry these annotations up to the gene level; expose them in both transcript and
gene rows when present.

- `split_parts` / `merge_parts`: partner IDs and per-part metrics.
- `ambiguous_partners`: partner IDs plus metric deltas.
- `antisense_conflict`: true with overlap fraction.
- `monoexonic_caveat`: true if one side is monoexonic and the other spliced.
- `phase_disagreement_segments`: list CDS intervals excluded by phase rule.
- `utr_divergence`: exon Jaccard outside CDS.

## Decision Sketch

1. If `Strand_agree` is false:
   - If `Jaccard_exon ≥ 0.50` and antisense is allowed, return Yellow with
     `antisense_conflict`.
   - Otherwise mark NotMapped with the conflict note.
2. Evaluate coding transcript pairs:
   - If any pair meets Green coding, set `Green_coding_candidate`.
   - Else if any meets Yellow coding, set `Yellow_coding_candidate`.
3. Evaluate non-coding pairs similarly.
4. Combine per-gene verdicts:
   - Any Green candidate → provisional Green.
   - Else any Yellow candidate → Yellow.
   - Else any Red evidence → Red.
   - Otherwise NotMapped.
5. Apply downgrade rules before finalising the gene record.

## Open Design Decisions

- **Data model changes**: extend the current parser to store transcript metrics
  grouped by gene without excessive memory use. Decide whether to stream gene
  summaries immediately after processing each locus or buffer until all pairs
  are seen.
- **Output format**: confirm whether gene-level summaries should be TSV rows
  interleaved with transcript rows or split into dedicated sections/files.
- **Threshold calibration**: the provisional cut-offs above come from the
  transcript prototype; we need to validate or tune them once real assemblies
  are processed.
- **Antisense policy**: final decision on whether antisense overlaps qualify as
  Yellow/Red or get shunted to NotMapped with `Conflict` notes.
- **Representative transcript selection**: decide if the gene summary should
  cite a single best transcript pair or retain all supporting pairs.
- **Terminology for unions**: ensure the code and outputs clearly separate the
  per-transcript union used for Jaccard metrics from the gene-wide exon
  footprint aggregation.

These decisions can be iterated on once the per-gene aggregation is hooked into
`within_assembly_compare.py` and we can inspect real-world output.
