# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project develops a bioinformatics workflow for mapping gene annotations across and within genome assemblies. The primary goal is to convert gene IDs from one annotation to another (e.g., from ICSASG_v2 to Ssal_v3.1) with accompanying metadata on mapping confidence and annotation similarity. The workflow will ultimately generate mapping tables for [Salmobase](https://salmobase.org).

## Development Environment

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate idmap

# Note: Some tools (ParsEval, Liftoff/LiftoffTools) are not available on osx-arm64 via conda
# Use Docker for those tools if needed on Mac ARM
```

The conda environment includes:
- Python 3.10
- bedtools, agat, gffcompare, gffread
- JupyterLab for analysis notebooks
- pandas, matplotlib, seaborn for data analysis

## Architecture

The workflow follows a two-step approach:

### 1. Cross-Assembly Mapping (Liftoff/LiftoffTools)
- Lifts gene annotations from source to target assembly
- Generates variant and synteny information
- Outputs stored in `experiments/liftoff_*/` and `experiments/liftofftools_*/` directories

### 2. Within-Assembly Comparison (Custom Script)
- Compares lifted annotations with native target annotations on the same assembly
- Uses streaming architecture to handle large genome-scale files
- Classifies gene mappings as Green/Yellow/Red based on overlap quality
- Main implementation: `within_assembly_compare.py`

## Key Scripts

### `gff_block_sort.py`
Sorts GFF3 files by `(seqid, gene_start)` while preserving gene blocks (all features belonging to a gene stay together).

**Usage:**
```bash
./gff_block_sort.py input.gff3 -o output.sorted.gff3
```

**Critical requirement:** All GFF3 files MUST be sorted this way before running `within_assembly_compare.py`. The comparison script validates this ordering and will error if it detects unsorted input.

### `within_assembly_compare.py`
Streams sorted GFF3 annotation files and computes gene/transcript-level overlap metrics.

**Usage:**
```bash
# Basic gene-level comparison
./within_assembly_compare.py annA.sorted.gff3 annB.sorted.gff3 -o comparison_output.tsv

# Include transcript-level details
./within_assembly_compare.py annA.sorted.gff3 annB.sorted.gff3 -o comparison_output.tsv --include-transcripts
```

**design document**
`docs/comparison-script-design.md` contains a detailed description. ALWAYS READ AND UPDATE THIS IF CHANGES ARE MADE!

**Architecture:**
- Streaming parser: only keeps the active locus in memory
- Processes one gene at a time with overlapping partner genes
- Computes transcript-level metrics, then rolls up to gene-level summaries
- Emits one TSV row per overlapping gene pair

## GFF3 Processing Conventions

### Gene ID Format
Gene IDs typically follow the pattern: `ID=gene:ENSSASG00000000001`

### Sorting Requirements
- All GFF files must be sorted by `(seqid, gene start)` with gene blocks intact
- Use `gff_block_sort.py` before any comparison operations
- The script validates ordering and will error on unsorted input

### Transcript Types
The comparison script recognizes these transcript types:
- mRNA, transcript
- ncRNA, lnc_RNA, miRNA, rRNA, tRNA, snRNA, snoRNA
- primary_transcript, pseudogenic_transcript

## Comparison Metrics and Classification

### Transcript-Level Metrics
- **Jaccard_exon**: Overlapping exon bp ÷ union of exon bp
- **Jaccard_CDS_phase**: Overlapping CDS bp (same strand, same codon phase) ÷ union of CDS bp
- **Junction_F1_all**: F1-score over all splice junctions
- **Junction_F1_CDS**: F1-score over CDS junctions only
- **Strand_agree**: Boolean flag (false = antisense)
- **Monoexonic**: True if transcript has no introns

### Gene-Level Classification
Transcript pairs are rolled up to gene-level summaries with these classifications:

**Green (High Confidence):**
- Coding genes: CDS Jaccard ≥ 98% or Junction F1 ≥ 95%
- Non-coding genes: Jaccard ≥ 90% or Junction F1 ≥ 95%
- Monoexonic: Jaccard ≥ 95% and length ratio delta ≤ 10%

**Yellow (Significant but Different):**
- Coding: CDS Jaccard ≥ 40% or Junction F1 ≥ 50%
- Non-coding: Jaccard ≥ 50% or Junction F1 ≥ 50%

**Red (Weak Evidence):**
- Exon Jaccard < 10% or antisense overlaps

**NotMapped:**
- No overlapping partner found

Thresholds are defined in `within_assembly_compare.py:48-66`.

### Split/Merge Detection
- Genes with multiple partners across annotations are flagged
- Affects confidence scoring for downstream use in Salmobase

## Data Organization

```
data/
├── genomes/          # Full genome FASTA and GFF downloads
└── toy-assemblies/   # HoxC A cluster test dataset
    ├── get_gff_subset.sh          # Extract GFF regions
    └── get_all_gff_subsets.sh     # Batch extraction

experiments/
├── liftoff_*/                      # Liftoff cross-assembly results
├── liftofftools_*/                 # Variant/synteny/cluster analysis
├── gffcompare_*/                   # GffCompare comparison tests
├── parseval_*/                     # ParsEval CDS-aware comparison
├── comparison_runs/                # within_assembly_compare.py outputs
└── within_assembly_compare/        # Prototype comparison results

docs/
├── tool-survey.md                  # Survey of candidate tools
├── experiments.md                  # Running log of experiments
├── comparison-script-design.md     # Design spec for comparison script
├── comparison-script-experiments.md # Testing documentation
└── AI-usage.md                     # Documentation of AI usage in project

notebooks/
└── ens_liftoff_vs_native.ipynb    # Gene-level comparison analysis
```

## Common Workflows

### Testing with Toy Dataset
The toy dataset contains HoxC A cluster regions from Atlantic salmon (Ssal_v3.1, ICSASG_v2) and Rainbow trout (Omyk):

```bash
# Scripts are in data/toy-assemblies/
cd data/toy-assemblies
./get_gff_subset.sh       # Extract single region
./get_all_gff_subsets.sh  # Batch extract all regions
```

See `docs/creating-toy-dataset.md` for details.

### Running Full Comparison Pipeline

1. **Sort GFF files:**
```bash
./gff_block_sort.py source_annotation.gff3 -o source.sorted.gff3
./gff_block_sort.py target_annotation.gff3 -o target.sorted.gff3
```

2. **Run Liftoff (cross-assembly mapping):**
```bash
# Liftoff is not in the conda environment for osx-arm64
# Use Docker or Linux environment
```

3. **Run LiftoffTools (variant/synteny analysis):**
```bash
# Also requires Docker/Linux for ARM Mac
```

4. **Compare within target assembly:**
```bash
./within_assembly_compare.py lifted.sorted.gff3 native.sorted.gff3 \
    -o comparison_output.tsv
```

### Analysis in Jupyter
```bash
jupyter lab
# Open notebooks/ens_liftoff_vs_native.ipynb for example analysis
```

## Integration with Salmobase

Outputs feed into the Salmobase database/API:
- Cross-assembly links use confidence metadata (Green/Yellow/Red)
- Related issues: #112 (ID conversion tool), #111 (Cross assembly annotation links), #79 (Cross annotation links improvements)
- See [Salmobase repo](https://gitlab.com/elixir-nmbu/salmobase) for integration details

## Processing Pattern for Large Files

The comparison script uses streaming to handle genome-scale files efficiently:

```python
# Stream processing pattern (from within_assembly_compare.py)
def process_gff(file):
    for gene in stream_genes(file):
        # Process one gene at a time
        # Only keep overlapping partner genes in memory
```

This approach:
- Keeps memory usage constant regardless of genome size
- Requires sorted input to identify when overlap windows close
- Processes genes in order, buffering only active loci

## Design Decisions

### Separate Coding/Non-Coding Flows
- Coding genes use CDS-aware metrics (phase-aligned overlap, CDS junctions)
- Non-coding genes use exon-based metrics only
- Both transcript and gene-level evidence are considered

### Streaming vs. In-Memory
- Initial prototypes loaded full files into memory
- Current implementation streams genes to handle multi-GB genome annotations
- Trade-off: requires sorted input but scales to any genome size

### Why Not Use Existing Tools?
- **GffCompare**: Not CDS-aware, misses important annotation differences
- **ParsEval**: CDS-aware but output format difficult to parse programmatically
- **Custom script**: Combines CDS awareness with structured TSV output tailored for Salmobase integration

See `docs/tool-survey.md` and `docs/experiments.md` for detailed tool evaluations.
