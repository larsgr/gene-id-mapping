# Cross-Assembly Gene ID Mapping Workflow

## Overview

This project aims to develop a bioinformatics workflow for mapping gene annotations across and within genome assemblies. The primary goal is to convert gene IDs from one annotation to another (e.g., from an older to a newer assembly) with accompanying metadata on mapping confidence and annotation similarity.

This workflow will ultimately be used to generate mapping tables for [**Salmobase**](https://salmobase.org), but this repository also serves as a lab notebook documenting tool evaluation, experiments, and methodology.

---

### Goals

- Automatically generate gene ID mappings between genome assemblies.
- Report similarity metrics between gene models (e.g., % CDS overlap, % exon overlap).
- Detect and warn about ambiguous mappings (e.g., non-1:1 relationships).
- Output metadata for confidence scoring and visualization (e.g., green/yellow/red icons).
- Document tool evaluations, findings, and rationale behind design decisions.

---

### Related Issues & Context

This project is related to following issues in the [Salmobase repo](https://gitlab.com/elixir-nmbu/salmobase) 
- [#112 ID conversion tool](https://gitlab.com/elixir-nmbu/salmobase/-/issues/112)
- [#111 Cross assembly annotation links](https://gitlab.com/elixir-nmbu/salmobase/-/issues/111)
- [#79 Cross annotation links improvements](https://gitlab.com/elixir-nmbu/salmobase/-/issues/79)

Initial strategies from these issues include:
- Compare exon/CDS coordinates for overlap-based similarity.
- Provide mapping confidence metadata (tooltips, color indicators).
- Highlight non-1:1 matches for manual review or exclusion.

---


### Roadmap / Progress

1. ✅ Research candidate tools ([tool survey](docs/tool-survey.md)).
2. ✅ Create toy dataset for initial testing ([creating toy dataset](docs/creating-toy-dataset.md)).
3. ✅ Experiment with Liftoff/LiftoffTools and document findings ([experiments: Liftoff](docs/experiments.md#liftoff)).
4. ✅ Test GffCompare on the Liftoff results ([experiments: GffCompare](docs/experiments.md#GffCompare)).
5. ✅ Try ParsEval (CDS‑aware gene model comparison) on the same inputs ([experiments: ParsEval](docs/experiments.md#ParsEval)).
6. 🔄 Develop and test mapping workflow (currently prototyping custom comparison logic in [`within_assembly_compare.py`](within_assembly_compare.py)).
7. ⬜ Integrate with Salmobase outputs.
8. ⬜ Generate final documentation and publish results.

---


### Repository Structure

```
.
├── README.md              # Project overview and progress log
├── environment.yml        # Conda env (gffcompare, gffread, agat, bedtools)
├── docs/                  # Notes, tool reviews, experiments
│   ├── AI-usage.md              # How AI assistants are used here
│   ├── tool-survey.md           # Survey of candidate tools and selection notes
│   ├── creating-toy-dataset.md  # Build steps for the HoxC A toy dataset
│   └── experiments.md           # Running log of experiments and findings
├── experiments/           # Scripts and results for exploratory testing (per-experiment dirs)
│   ├── gffcompare_multi/       # Multi-query gffcompare run building a consensus annotation
│   ├── gffcompare_test/        # Baseline gffcompare reference-vs-Liftoff comparison outputs
│   ├── liftoff_test/           # Liftoff liftover GFF for HoxC A region (ICSASG_v2 → Ssal_v3.1)
│   ├── liftofftools_test/      # LiftoffTools variant-effect summaries plus R inspection script
│   ├── liftofftools_test2/     # Follow-up LiftoffTools run capturing alternate variant effects
│   ├── parseval_html/          # ParsEval-generated HTML dashboards for CDS comparisons
│   ├── parseval_test/          # ParsEval CLI outputs (text/CSV metrics)
│   └── within_assembly_compare/# Prototype Ensembl vs Liftoff on-assembly comparison results
├── data/                  # Toy datasets and annotation files
├── within_assembly_compare.py  # Prototype for on-assembly GFF3 comparison metrics
└── workflow/              # Pipeline (e.g., Snakemake or Nextflow) for production use
```


---

## Overview of results

### The toy-assembly dataset (data/toy-assemblies)

The both copies of the HoxC A cluster from two versions of Atlantic salmon reference genome (Ssal_v3.1 and ICSASG_v2) and Rainbow trout (Omyk) was extracted as a small dataset to use for testing. This includes th corresponding gene annotations (gff) from Ensembl and NCBI for these regions.

The gff extraction scripts are in `data/toy-assemblies`. See [docs/creating-toy-dataset.md](docs/creating-toy-dataset.md) for more details.

### Tool Candidates

Gemini deep research was used to identify relevant existing tools ([full report here](docs/tool-survey.md))

Existing relevant tools identified:

| Tool | Description |
|------|-------------|
| **GffCompare** | Annotation comparison within a single assembly. (note: not CDS aware) |
| **Liftoff / LiftoffTools** | Cross-assembly mapping with a promising "variants" module. |
| **ExOrthist** | Identifies orthologous exons across species/assemblies. (not tested) |
| **ParsEval** | Compares alternative annotations. (note: hard to parse output format) |
| **GeneOverlapAnnotator** | Part of SVAnnotator, returns gene/transcript overlap stats. (not tested) |

We now have full experiment notes for **Liftoff**, **GffCompare**, and **ParsEval** in [docs/experiments.md](docs/experiments.md). ParsEval's reporting confirmed CDS-conserving vs CDS-changing differences but remains cumbersome to parse, motivating development of a custom within-assembly comparison script (`within_assembly_compare.py`).

---

## Getting Started

To begin experimenting:

1. **Clone the repository**:

   ```bash
   git clone https://github.com/larsgr/gene-id-mapping.git
   cd gene-id-mapping
   ```

2. **Set up the development environment**:

   Create the conda environment from the provided file:

   ```bash
   conda env create -f environment.yml
   conda activate idmap
   ```

3. **Explore experiments**:

   See [docs/experiments.md](docs/experiments.md)

---
