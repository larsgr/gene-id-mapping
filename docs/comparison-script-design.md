Metrics to compute and store

Jaccard_exon: bp overlap of exon unions divided by bp union.

Jaccard_CDS_phase: bp overlap of CDS unions counted only where both CDS are annotated and codon phase matches at each genomic position.

Junction_F1_all: F1 for exact splice junctions over all introns.

Junction_F1_CDS: F1 for junctions within CDS only.

Strand_agree: boolean.

Shared_introns: count and fraction.

Monoexonic: boolean for each gene.

Has_split_merge: flags with partner counts.

Antisense_overlap: fraction of exon union overlapped on opposite strand.

Keep all raw counts too. Thresholds are placeholders.

Representative transcript rule

Coding gene: evaluate across all transcript pairs. If any pair meets Green, the gene is Green unless a downgrade rule fires. Same for Yellow.

Noncoding gene: same approach using noncoding pairs.

Main classes
Green

Coding: Strand_agree and reciprocal at gene level is not required since many-to-many is allowed. Any transcript pair satisfies:

Jaccard_CDS_phase ≥ 0.98 or

Junction_F1_CDS ≥ 0.95 and Jaccard_CDS_phase ≥ 0.95.
UTR ignored. Non-CDS introns ignored.

Noncoding:

If spliced in both: Junction_F1_all ≥ 0.95 or (Jaccard_exon ≥ 0.90 and Junction_F1_all ≥ 0.90).

If monoexonic in either: Jaccard_exon ≥ 0.95 and length ratio within 0.9–1.1.

Antisense never Green.

Yellow

Coding: Same strand. Any transcript pair satisfies:

0.50 ≤ Jaccard_CDS_phase < 0.98 or

Junction_F1_CDS ≥ 0.5 with Jaccard_CDS_phase ≥ 0.40.

Noncoding:

0.50 ≤ Jaccard_exon < 0.90 or

0.50 ≤ Junction_F1_all < 0.95.

Monoexonic: 0.50 ≤ Jaccard_exon < 0.95.

Antisense overlap ≥ 0.50 of exon union can be Yellow only if you explicitly allow antisense families. Otherwise mark Conflict in notes and treat as NotMapped.

Red

Same strand and weak evidence:

0.10 ≤ Jaccard_exon < 0.50, or

at least one shared intron but Junction_F1_all < 0.50.

Antisense with Jaccard_exon ≥ 0.10 goes to notes as Conflict unless you want to count it as Red. Your call. Default NotMapped with Conflict note.

NotMapped

Different chromosomes, no overlap, or below Red.

Any antisense case if you choose to disallow antisense mappings.

Downgrade rules from Green to Yellow

Store triggers in notes and set class to Yellow if any fire:

Split or merge detected: one gene overlaps two or more distinct partners with nonoverlapping exon blocks where each sub-pair is at least Yellow.

Competing partner within delta: another partner within Δ of the main metric (e.g., ΔJaccard_CDS_phase ≤ 0.02 or ΔJunction_F1_CDS ≤ 0.02).

Frame inconsistency across transcripts: some transcripts align in-phase, others only out-of-phase or with zero Jaccard_CDS_phase.

Discordant biotypes: coding vs noncoding between sources for the same locus.

Problematic intron: shared CDS junction fraction < 1.0 when Jaccard_CDS_phase barely passes.

Extensive non-CDS differences: Jaccard_exon < 0.70 despite Green CDS.

Notes and secondary flags

split_parts or merge_parts with partner IDs and per-part metrics.

ambiguous_partners with metric deltas.

antisense_conflict: true with overlap fraction.

monoexonic_caveat: true when one side is monoexonic and the other is spliced.

phase_disagreement_segments: list of CDS intervals excluded by phase rule.

utr_divergence: exon Jaccard outside CDS.

Decision sketch

If Strand_agree is false:

If Jaccard_exon ≥ 0.50 and you allow antisense, set Yellow with antisense_conflict. Else NotMapped with note.

Else compute class for coding transcript pairs:

If any pair meets Green coding, set Green_coding_candidate.

Else if any meets Yellow coding, set Yellow_coding_candidate.

Compute class for noncoding pairs similarly.

Combine:

If any Green candidate exists, class = Green initially.

Else if any Yellow candidate exists, class = Yellow.

Else if any Red condition exists, class = Red.

Else NotMapped.

Apply downgrade rules to Green. If any trigger, class → Yellow.