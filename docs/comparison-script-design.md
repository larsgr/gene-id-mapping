# Comparison Script Design

This note now documents the streamed `within_assembly_compare.py`
implementation. The parser keeps only the active locus in memory, computes
transcript-level overlaps, rolls them into per-gene summaries, and emits one TSV
row per overlapping gene pair (with optional transcript rows interleaved on
request).

## Streaming architecture

- **Input contract**: each GFF must be sorted by `(seqid, gene start)` with gene
  blocks kept intact. The script validates this ordering and exits with an error
  if it detects regression. Use `./gff_block_sort.py` (documented in
  `docs/somparison-script-experiments.md`) to pre-sort inputs safely.
- **Gene streaming**: the parser yields one `Gene` object at a time (with all
  transcripts resolved) and maintains only the overlapping partner genes from
  the other annotation in memory.
- **Output buffering**: per-gene records are buffered until both partners have
  no further overlaps, at which point split/merge partner counts can safely be
  applied before writing.

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

For every overlapping transcript pair the parser logs the metrics above, plus
derived values such as Jaccard scores, junction F1, monoexonic flags, and a
pair-level classification. Transcript rows are written only when
`--include-transcripts` is specified; otherwise they remain internal evidence
for the gene-level label.

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

The streaming implementation adds a lightweight buffering layer that collects
all transcript comparisons for a gene pair, records partner counts for split /​
merge detection, and writes the gene row only once the overlap window closes.

### Gene Exon Footprint vs. Transcript Union

Per-gene statistics still use a "gene exon footprint" for each annotation: the
union of exons across every transcript isoform for a gene. This concept is
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
- Antisense overlaps are no longer allowed; see "Antisense handling" below.

### Red
- Same strand but weak evidence:
  - `0.10 ≤ Jaccard_exon < 0.50`, or
  - At least one shared intron while `Junction_F1_all < 0.50`.
- Antisense overlaps with a large fractional overlap are treated as Red
  warnings (threshold defaults to `0.5`).

### NotMapped
- Different chromosomes, zero overlap, or metrics below the Red thresholds.
- Any antisense case when antisense mappings are disallowed.

### Antisense handling

- Set `note_antisense_conflict` on every opposite-strand overlap.
- If exon-overlap / exon-union ≥ antisense threshold (default 0.5) downgrade to
  Red and emit the gene row; otherwise suppress the mapping (NotMapped) but keep
  the warning in the notes.

## Downgrade Rules (Green → Yellow)

Gene-level downgrades fire if any transcript pair triggers them; record the
trigger in the notes and demote the gene summary:

- **Split / Merge detection**: overlaps with ≥ 2 distinct partners (currently
  triggers an automatic downgrade to Yellow and records `note_split_or_merge_*`).
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

- `split_or_merge_A` / `split_or_merge_B`: partner counts when we detect a split
  or merge and downgrade the class.
- `discordant_biotype`: set when coding/non-coding evidence conflicts.
- Additional flags from the earlier wish list (e.g. `ambiguous_partners`,
  `phase_disagreement_segments`) remain TODO.

## Decision Sketch

1. If `Strand_agree` is false:
   - Compute exon overlap fraction. If it exceeds the antisense threshold,
     classify as Red with an `antisense_conflict` warning; otherwise return
     NotMapped.
2. Evaluate coding transcript pairs:
   - If any pair meets Green coding, set `Green_coding_candidate`.
   - Else if any meets Yellow coding, set `Yellow_coding_candidate`.
3. Evaluate non-coding pairs similarly.
4. Combine per-gene verdicts:
   - Any Green candidate → provisional Green.
   - Else any Yellow candidate → Yellow.
   - Else any Red evidence → Red.
   - Otherwise NotMapped.
5. Apply downgrade rules before finalising the gene record (currently split /
   merge and discordant biotype are wired in; others remain on the backlog).

## Open Design Decisions

- **Input preparation**: provide or automate a safe GFF block sorter so users
  do not break gene grouping when satisfying the sorted-input requirement.
- **Threshold calibration**: the provisional cut-offs above come from the
  transcript prototype; we need to validate or tune them once real assemblies
  are processed.
- **Additional downgrade rules**: competing partner deltas, frame
  inconsistency, problematic CDS introns, and UTR divergence are still pending.
- **Representative transcript selection**: decide if the gene summary should
  cite a single best transcript pair or retain all supporting pairs (today all
  supporting pairs are kept when `--include-transcripts` is used).
- **Terminology for unions**: ensure the code and outputs clearly separate the
  per-transcript union used for Jaccard metrics from the gene-wide exon
  footprint aggregation.

These decisions can be iterated on once the per-gene aggregation is hooked into
`within_assembly_compare.py` and we can inspect real-world output.
