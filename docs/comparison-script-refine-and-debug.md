# Refining and debugging the comparison script

We have implemented a working script that can compare gff files within an assembly. This has only been tested with a toy-dataset so we need to test more broadly to find edge-cases. E.g. compare full gff files, e.g. between the old ICSASG_v2 salmon assembly and the current Ssal_v3.1, and between NCBI and Ensembl annotations.

Experimentation will probably reveal issues that needs to fixed or features that needs to added or refined

### Cross-Assembly Liftoff Comparisons

Multiple Atlantic salmon assemblies have been lifted to Ssal_v3.1 and compared with the native Ensembl annotation:

- **ICSASG_v2 → Ssal_v3.1**: Full comparison documented in `experiments/comparison_runs/ens_lift_vs_native.tsv`
- **Ssal_Brian_v1.0 → Ssal_v3.1**: Comparison output in `experiments/comparison_runs/brian_lift_vs_native.tsv`

A comprehensive Rmarkdown analysis comparing both liftoffs is available at `notebooks/compare_liftoff_assemblies.Rmd`, which evaluates:
- Classification distributions (Green/Yellow/Red/NotMapped)
- Quality metrics (CDS-phase and exon Jaccard indices)
- Stable Ensembl ID preservation
- Split/merge events

This analysis helps determine which source assembly provides better annotation continuity for Salmobase integration.

### Within-Assembly Annotation Comparison (NCBI vs Ensembl)

A direct comparison of NCBI and Ensembl annotations for the Ssal_v3.1 assembly has been completed:

- **Location**: `experiments/comparison_runs/ncbi_vs_ensembl_Ssal_v3.1/`
- **Report**: `comparison_report.html` (interactive HTML report with visualizations)
- **Data**: 42,041 NCBI genes compared with 40,589 Ensembl genes
- **Key findings**:
  - 31.2% Green (high-quality matches)
  - 37.9% Yellow (significant differences)
  - 7.6% Red (weak evidence)
  - 23.1% NotMapped (no sufficient overlap)

This comparison required converting NCBI RefSeq chromosome names to Ensembl format using the Salmobase chromosome mapping table. The analysis includes detailed metrics on CDS overlap, junction agreement, and split/merge events, providing valuable insights for cross-annotation linking in Salmobase.

## Fix comparison script to handle all gene-like features (2026-02-03)

**What was done:**
* Added `GENE_TYPES` constant to `within_assembly_compare.py` and `gff_block_sort.py`
  - Includes: `gene`, `ncRNA_gene`, `pseudogene`, `V_gene_segment`, `D_gene_segment`, `J_gene_segment`
* Updated gene detection logic from `if feature_type == "gene"` to `if feature_type in GENE_TYPES`
* Modified error messages to include actual feature type
* Re-sorted NCBI and Ensembl GFF files: `Ssal_v3.1_NCBI.sorted.gff3` and `Ssal_v3.1_Ens.sorted.gff3`
* Re-ran comparison generating `experiments/comparison_runs/ncbi_vs_ensembl_Ssal_v3.1/comparison_v2.tsv`
* Re-rendered comparison report: `experiments/comparison_runs/ncbi_vs_ensembl_Ssal_v3.1/comparison_report.html`

**Findings:**
The fix dramatically expanded the scope of the comparison:
- **NCBI genes compared**: 54,016 (previously ~42k)
- **Ensembl genes compared**: 54,268 (previously ~40k)
- **Total comparison rows**: 67,546

**Coding status distribution** (transcript-level):
- Both coding: 452,094
- **Both non-coding: 15,170** (NEW - was entirely missing before!)
- NCBI coding, Ensembl non-coding: 10,220
- NCBI non-coding, Ensembl coding: 28,912

**Source GFF files contain**:
- NCBI: 65,475 gene-like features
- Ensembl: 69,551 gene-like features
  - Including 21,576 `ncRNA_gene` features (lncRNAs, rRNAs, snRNAs, etc.)
  - Including 608 `pseudogene` features

The comparison now properly processes all gene-like features from both annotations, providing a complete and fair comparison. The appearance of the "Both non-coding" category with over 15,000 entries confirms that the script was previously missing a substantial portion of the annotations.

