# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project develops a bioinformatics workflow for mapping gene annotations across and within genome assemblies. The primary goal is to convert gene IDs from one annotation to another (e.g., from ICSASG_v2 to Ssal_v3.1) with accompanying metadata on mapping confidence and annotation similarity. The workflow will ultimately generate mapping tables for [Salmobase](https://salmobase.org).

## Workflow Structure

### Session Start Protocol
**IMPORTANT:** At the start of each session, automatically check git status:
- Run `git status` to check for uncommitted changes
- If uncommitted changes exist, display them and warn the user:
  ```
  ⚠️ Uncommitted changes detected:
  [show git status output]

  Manual edits should be committed before starting AI work.
  Would you like to commit these changes first? [Yes/No]
  ```
- If Yes: Guide user through creating a commit
- If No: Acknowledge and proceed with the session
- User should commit manual edits from previous sessions before starting new AI work

### Documentation Hierarchy
All work follows this three-tier structure:

1. **README.md → Roadmap Section**
   - High-level tasks using emoji status:
     - ✅ Completed
     - 🔄 In progress (may span multiple sessions/commits)
     - ⬜ Not started
   - Each task links to a detailed log file in docs/

2. **docs/[task-name].md → Detailed Task Log**
   - One log file per Roadmap task
   - Each AI session adds a new dated section with:
     - What was done (commands, files created/modified, analysis performed)
     - Key findings and conclusions
     - Problems encountered and solutions
   - Session section names should match:
     - The entry in docs/AI-usage.md
     - The git commit message for that session

3. **experiments/[experiment-name]/ → Experiment Scripts and Data**
   - Single-use analysis/experiment scripts stored in subdirectories
   - Each subdirectory contains:
     - Analysis scripts (Python, R, bash, etc.)
     - Output data (with .gitignore for files >1MB)
     - Brief README.md if complex
   - Referenced from the detailed task log

### Session End Protocol
**IMPORTANT:** When work appears complete, automatically remind the user to document:
- Detect potential session end (e.g., user says "thanks", "that's all", work is complete)
- Prompt: "This session's work should be documented. Would you like to run /session-end now?"
- The `/session-end` skill will:
  - Gather session summary information via prompts
  - Update docs/AI-usage.md with new entry
  - Update or create detailed task log in docs/
  - Suggest commit message matching session name
  - Verify all documentation is complete
- **Do not consider a session complete until documentation is updated**

## File Organization Rules

### Experiment Data
- Single-use experiments: `experiments/[descriptive-name]/`
- Always create .gitignore for generated files >1MB
- Include script used to generate outputs for reproducibility

### Analysis Reports
- **Rmarkdown reports**: Store in `notebooks/` directory
- Render as both .md (commit) and .html (add to .gitignore)
- Final polished reports can be moved to `docs/` if needed
- Initial/exploratory analyses stay in `notebooks/`

### Tool Scripts
- Reusable tools go in repo root (e.g., `gff_block_sort.py`, `within_assembly_compare.py`)
- Scripts that are part of the final workflow, not one-off analyses
- Should be documented in CLAUDE.md under "Key Scripts" section

### Documentation Files
- `docs/AI-usage.md`: Log of all AI sessions (prompt + summary + reflection)
- `docs/[task-name].md`: Detailed logs per Roadmap task
- `docs/[topic].md`: Design documents, guides, surveys (e.g., tool-survey.md)

## Documentation Practices

**IMPORTANT:** After completing any significant work, you must update the following documentation files:

### 1. README.md
Update when:
- New analysis outputs are generated (e.g., comparison runs, visualizations)
- New scripts, notebooks, or tools are created
- Repository structure changes significantly
- New workflows or pipelines are established

Include:
- Brief description of what was added
- Location of outputs/files
- Purpose and how it fits into the overall project

### 2. docs/AI-usage.md
**Always document your work** by adding an entry with:
- **Prompt**: The user's original request (in a code block)
- **What it did**: Summary bullet list of actions taken, including:
  - Commands executed
  - Files created/modified
  - Analysis performed
  - Any problems encountered and how they were resolved
- **Reflection**: Your assessment of what worked well and what could be improved
  - Note any technical challenges
  - Document error recovery strategies
  - Suggest improvements for future work
  - Comment on the overall workflow effectiveness

### 3. docs/comparison-script-design.md
If you modify `within_assembly_compare.py`, **you MUST update this design document** to reflect:
- Algorithm changes
- New metrics or thresholds
- Classification rule modifications
- Performance improvements

### 4. Session-End Workflow

**MANDATORY:** At the end of each session, run `/session-end` slash command which will:
1. Prompt for session summary information
2. Update docs/AI-usage.md with new entry
3. Update relevant detailed log under docs/
4. Suggest commit message matching session name
5. Verify all documentation is complete

**Do not complete a session without running `/session-end`.**

### Documentation Format Example
```markdown
### (Claude Code) Brief description of task (YYYY-MM-DD)

**Prompt:**
\`\`\`
[User's original request]
\`\`\`

**What it did:**
* Action 1
* Action 2
* Created file X with purpose Y
* Encountered error Z, resolved by doing W

**Reflection:**
[Your assessment of the work, what worked well, areas for improvement]
```

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
- Full liftoff outputs stored in `data/genomes/{species}/liftoff/{source_asm}_to_{target_asm}_{annotation}.gff[.gz]`
  - Example: `data/genomes/AtlanticSalmon/liftoff/ICSASG_v2_to_Ssal_v3.1_Ens.gff3`
  - Annotation suffix: "Ens" = Ensembl, "NCBI" = NCBI
- Experimental/test runs stored in `experiments/liftoff_*/` and `experiments/liftofftools_*/` directories

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
│   └── {species}/
│       └── liftoff/  # Full liftoff outputs: {source_asm}_to_{target_asm}_{annotation}.gff[.gz]
└── toy-assemblies/   # HoxC A cluster test dataset
    ├── get_gff_subset.sh          # Extract GFF regions
    └── get_all_gff_subsets.sh     # Batch extraction

experiments/
├── liftoff_*/                      # Liftoff test/experimental runs
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
# Output should be saved to: data/genomes/{species}/liftoff/{source_asm}_to_{target_asm}_{annotation}.gff[.gz]
```

3. **Run LiftoffTools (variant/synteny analysis):**
```bash
# Also requires Docker/Linux for ARM Mac
```

4. **Compare within target assembly:**
```bash
# Use the liftoff output from data/genomes/{species}/liftoff/
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
