---
name: gff-structure
description: Reference for Atlantic salmon GFF3 file structure (Ensembl and NCBI annotations on Ssal_v3.1). Use this when working with GFF files, parsing gene attributes, or writing scripts that process annotation data.
version: 1.0.0
---

# GFF3 File Structure: Atlantic Salmon Annotations

## File Locations
- Ensembl: `data/genomes/AtlanticSalmon/Ssal_v3.1_Ens.gff3` (plain text)
- NCBI: `data/genomes/AtlanticSalmon/Ssal_v3.1_NCBI.gff3` (plain text)

## Ensembl GFF3

### Gene-level feature types
- `gene` — protein-coding genes (source: `ensembl`)
- `ncRNA_gene` — non-coding RNA genes (source: `ensembl` or `ncrna`)
- `pseudogene` — pseudogenes (source: `ensembl`)

### Attributes format (col 9)
```
ID=gene:ENSSSAG00000110109;Name=ube2d3;biotype=protein_coding;description=ubiquitin-conjugating enzyme E2D 3 [Source:ZFIN%3BAcc:ZDB-GENE-030131-551];gene_id=ENSSSAG00000110109;logic_name=ensembl;version=1
```

### Key attributes
- `ID` — always prefixed with `gene:` (e.g., `gene:ENSSSAG00000110109`)
- `Name` — **optional**, only ~44% of genes have it
- `biotype` — e.g., `protein_coding`, `lncRNA`, `pseudogene`, `snRNA`, `snoRNA`, `miRNA`, `tRNA`, `rRNA`
- `description` — URL-encoded, contains source in bracket suffix: `[Source:ZFIN%3BAcc:ZDB-GENE-030131-551]`
- `gene_id` — same as ID without `gene:` prefix
- `logic_name` — `ensembl` or `ncrna`

### Description sources (extracted from `[Source:XXX;Acc:YYY]`)
- ZFIN — zebrafish nomenclature (most common)
- RFAM — RNA families
- NCBI gene
- HGNC — human gene nomenclature
- miRBase

### Naming patterns
- If a gene has `Name=`, it always also has `description=`
- No genes have `description=` without `Name=`
- lncRNA genes are essentially never named

## NCBI GFF3

### Gene-level feature type
- Only `gene` (no ncRNA_gene/pseudogene distinction at feature level)

### Attributes format (col 9)
```
ID=gene-LOC123729278;Dbxref=GeneID:123729278;Name=LOC123729278;description=lethal(3)malignant brain tumor-like protein 2;gbkey=Gene;gene=LOC123729278;gene_biotype=protein_coding
```

### Key attributes
- `ID` — prefixed with `gene-` (e.g., `gene-LOC123729278` or `gene-wdr32`)
- `Name` — **always present**, but often auto-assigned LOC IDs
- `Dbxref` — contains `GeneID:NNNNN` (NCBI numeric gene ID)
- `gene_biotype` — (note: NOT `biotype` like Ensembl) e.g., `protein_coding`, `lncRNA`, `tRNA`, `rRNA`, `snRNA`, `snoRNA`, `misc_RNA`
- `description` — plain text, URL-encoded (e.g., `%2C` for comma)
- `gene` — same as Name
- `gbkey` — always `Gene`

### Naming patterns
- ALL genes have `Name=` but most are auto-assigned `LOC` IDs (start with "LOC")
- Real gene symbols are lowercase (e.g., `wdr32`, `manba`)
- LOC genes: description may be informative or "uncharacterized LOCnnn"

### Distinguishing real names from LOC IDs
- Check if Name starts with "LOC" followed by digits
- `is_loc = name.startswith("LOC")`

## Key Differences Between Ensembl and NCBI

| Aspect | Ensembl | NCBI |
|--------|---------|------|
| Gene ID prefix | `gene:` | `gene-` |
| Gene feature types | `gene`, `ncRNA_gene`, `pseudogene` | `gene` only |
| Biotype attribute | `biotype` | `gene_biotype` |
| Name attribute | Optional (~44% have it) | Always present (but often LOC) |
| Numeric ID | In `gene_id` attribute | In `Dbxref=GeneID:` |

## Seqid Values
- Both use chromosome numbers: `1`, `2`, ..., `29`
- Plus unplaced scaffolds

## NCBI Gene Orthologs (ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz)
- Pairs stored directionally: "reference" organism in `tax_id`, "other" in `Other_tax_id`
- Atlantic salmon (tax_id=8030) is always in `Other_tax_id` column (never in `tax_id`)
- Human-salmon pairs: filter `tax_id == 9606 & Other_tax_id == 8030`
- `GeneID` = human gene, `Other_GeneID` = salmon gene
- Header starts with `#tax_id` — needs renaming after read
