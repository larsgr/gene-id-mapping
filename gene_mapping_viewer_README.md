# Gene Mapping Viewer

An interactive R Shiny app for visualizing gene mappings between ICSASG_v2 and Ssal_v3.1 genome assemblies.

## Features

1. **Chromosome Pair Matrix**: A heatmap showing the number of genes mapping between each pair of chromosomes from the two assemblies
   - Only shows main chromosomes (ICSASG_v2: ssa01-ssa29, Ssal_v3.1: 1-29)
   - Scaffolds and unplaced sequences are filtered out
   - Updates dynamically based on active filters

2. **Gene Mapping Categories**:
   - **Liftoff mapped (all)**: All genes from the within-assembly comparison (liftoff results)
   - **Liftoff mapped (ID match only)**: Subset of liftoff results where gene IDs match between assemblies
   - **Stable ID only**: Genes with matching stable IDs in both assemblies that were NOT found by liftoff

3. **Interactive Filters**: Toggle visibility of different gene categories:
   - Check/uncheck to show or hide each category
   - If both "all" and "ID match only" are checked, all liftoff genes are shown
   - If only "ID match only" is checked, only liftoff genes with matching IDs are shown
   - Counts update in real-time as filters change

4. **Interactive Dot Plot**: Click on any chromosome pair to see a detailed view of gene positions with color-coding by mapping quality:
   - **Green**: High-confidence mappings
   - **Yellow**: Significant but different mappings
   - **Red**: Weak evidence mappings
   - **Gray**: Not mapped
   - **Blue**: Stable ID only (not in liftoff comparison)

## Installation

### Required R Packages

```r
install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "readr", "stringr"))
```

Or install from the R console when prompted during first run.

## Usage

### Starting the App

From the command line:

```bash
Rscript gene_mapping_viewer.R
```

Or from within R:

```r
source("gene_mapping_viewer.R")
```

The app will open in your default web browser.

### Using the App

1. **Load Data**: The default file paths are pre-filled. Click "Load Data" to load the gene mapping data. This may take a minute for large genome files.

2. **Apply Filters**: Use the checkboxes to show/hide different gene categories:
   - Start with all filters checked to see the complete picture
   - Uncheck categories to focus on specific subsets
   - The status line shows counts for each category

3. **View Matrix**: The matrix plot shows all chromosome pairs. Darker blue indicates more gene mappings.

4. **Explore Details**: Click on any cell in the matrix to view the detailed dot plot for that chromosome pair.

5. **Interpret the Dot Plot**:
   - X-axis: Position on ICSASG_v2 chromosome (Mb)
   - Y-axis: Position on Ssal_v3.1 chromosome (Mb)
   - Colors indicate mapping quality (Green/Yellow/Red/Blue)
   - Diagonal patterns suggest conserved gene order
   - Scattered points suggest rearrangements
   - Blue dots (Stable ID only) show genes that have the same ID in both assemblies but weren't found by liftoff

### Interpreting Results

**Discrepancies between Liftoff and Stable IDs:**

- **Blue dots on diagonal**: Genes with matching stable IDs at similar positions that liftoff missed - may indicate liftoff sensitivity issues
- **Blue dots off diagonal**: Genes with matching stable IDs at different chromosome positions - suggests potential assembly differences or annotation inconsistencies
- **Green/Yellow/Red dots where gene IDs differ**: Liftoff found an overlap but mapped to a different gene ID - investigate whether the stable IDs changed or if there's a true biological difference

**Use the filters to investigate:**

1. Check only "Stable ID only" to see genes missed by liftoff
2. Check only "Liftoff mapped (ID match only)" to see where liftoff and stable IDs agree
3. Check only "Liftoff mapped (all)" then uncheck "ID match only" to see where they disagree

### Custom Data Files

To use different data files, update the file paths in the sidebar:
- **ICSASG_v2 GFF**: Original source assembly annotation
- **Lifted GFF**: Liftoff output (source genes on target assembly)
- **Ssal_v3.1 Native GFF**: Native target assembly annotation
- **Comparison TSV**: Output from `within_assembly_compare.py`

## Data Requirements

The app expects:
1. GFF3 files with gene features containing `ID=gene:GENE_ID` in the attributes column
2. A comparison TSV file from `within_assembly_compare.py` with columns:
   - `feature`, `annA`, `geneA`, `txA`, `annB`, `geneB`, `txB`, `stats`
   - The `stats` column must contain `class=Green|Yellow|Red|NotMapped`

## Performance

- Initial data loading may take 1-2 minutes for full genome files
- The matrix plot is generated once after loading
- Dot plots are generated on-demand when clicking chromosome pairs
- For very large chromosome pairs (>5000 genes), rendering may take a few seconds

## Troubleshooting

**Error: Cannot find GFF file**
- Check that file paths are correct relative to the working directory
- Use absolute paths if needed

**Error: Package not installed**
- Install missing packages: `install.packages("package_name")`

**App is slow**
- Large GFF files take time to parse
- Consider filtering GFF files to include only gene features before loading

**No genes in dot plot**
- This is expected for chromosome pairs with no mappings
- Check the matrix to see which pairs have gene mappings
