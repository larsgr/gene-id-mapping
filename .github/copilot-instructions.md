# AI Agent Instructions for gene-id-mapping

This project develops a bioinformatics workflow for mapping gene annotations between genome assemblies, with a focus on generating mapping tables for Salmobase. The key goal is to accurately map gene IDs across assemblies while providing quality metrics and confidence scores.

## Project Architecture

The workflow follows a two-step approach:

1. Cross-assembly mapping using Liftoff/LiftoffTools:
   - Lifts gene annotations from source to target assembly
   - Generates variants and synteny information
   - Key files: Outputs in `experiments/liftoff_*/` directories

2. Within-assembly comparison using custom script:
   - Compares lifted annotations with native target annotations 
   - Classifies mappings as Green/Yellow/Red based on overlap metrics
   - Key file: `within_assembly_compare.py`

## Critical Workflows

### Development Environment

```bash
conda env create -f environment.yml
conda activate idmap
```

### Processing Pipeline

1. Sort GFF files using `gff_block_sort.py` (preserves gene blocks)
2. Run Liftoff to map annotations between assemblies
3. Run LiftoffTools for variant/synteny analysis
4. Use `within_assembly_compare.py` for overlap analysis

## Key Conventions

### GFF Processing
- All GFF files must be sorted by `(seqid, gene start)` maintaining gene blocks
- Use `gff_block_sort.py` before comparison operations
- Gene IDs follow format `ID=gene:ENSSASG00000000001`

### Comparison Classifications
- **Green**: High confidence matches (e.g. CDS overlap ≥98% or junction F1 ≥95%)
- **Yellow**: Significant overlap but with differences
- **Red**: Weak evidence or antisense overlaps
- See `docs/comparison-script-design.md` for detailed criteria

### Analysis Conventions
- Group transcript pairs by gene pair before classification
- Consider coding/non-coding status separately
- Calculate both gene-level and transcript-level metrics

## Common Patterns

### Processing Large Files
```python
# Stream processing pattern (used in within_assembly_compare.py)
def process_gff(file):
    for gene in stream_genes(file):
        # Process one gene at a time
        # Only keep overlapping partner genes in memory
```

### Metric Calculation
```python
# Common metrics used across the project
metrics = {
    'Jaccard_exon': overlap_bp / union_bp,
    'Junction_F1': 2 * shared / (tx1_junctions + tx2_junctions),
    'CDS_phase_match': phase_aligned_bp / total_cds_bp
}
```

## Integration Points

- Outputs feed into Salmobase database/API
- Cross-assembly links use metadata for confidence scoring
- See Salmobase repo for integration details

## Project-Specific Notes

- Focus on streaming architecture for memory efficiency
- Maintain separate flows for coding vs non-coding genes
- Consider both transcript and gene-level evidence
- Document AI usage patterns in `docs/AI-usage.md`