# Comparison Script Experiments

## 2025-10-01 Block-sorted toy annotations (Ensembl vs NCBI)

Goal: sanity-check the streaming comparison on the toy ICSASG_v2 annotation pair.

### Preparation

The parser insists on sorted input and on keeping gene blocks intact. I wrote a
quick Python helper (`python - <<'PY' ...`) that buckets every gene with its
nested features and rewrites the file ordered by `(seqid, gene_start)`.
Sorted copies now live under `experiments/sorted/` with the suffix
`.blocksorted.gff`.

### Run

```
./within_assembly_compare.py \
  experiments/sorted/ICSASG_v2_hoxca_Ens.blocksorted.gff \
  experiments/sorted/ICSASG_v2_hoxca_NCBI.blocksorted.gff \
  -l ens ncbi \
  -o experiments/comparison_test.tsv
```

Re-running with `--include-transcripts` confirmed that transcript rows are only
emitted when explicitly requested.

### Findings

* Streaming works: resident memory stays near a single locus, yet the split-
  merge downgrades still fire (see `note_split_or_merge_*` in the output).
* Most genes classify `Green`; loci with multiple partners correctly downgrade
  to `Yellow`.
* No antisense warnings surfaced in this dataset; will need a dedicated test
  case later.

The TSV artifact is stored as `experiments/comparison_test.tsv`
(and `comparison_test_with_tx.tsv` for the transcript-inclusive variant).
