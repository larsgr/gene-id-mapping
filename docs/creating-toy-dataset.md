Task: create a small "toy" dataset for testing

# roadmap:

1. (done) Extract sequences manually using JBrowse2 at salmobase
2. (partially done without alias handling) Make scripts to extract gff subset 

## Extracting HOXC cluster as a toy-assembly

I choose to extract the both Ss4R copies of one of the HOX clusters (HOXC A) as a subset. Specifically from gene cbx5 to hoxc13a. I know there are some inconsistencies in the salmonid hox cluster gene annotations and since the genes are very conserved across species we have a good idea of what should be the correct gene structure. By keeping both ohnolog regions we also get to test the effect of having duplicated genomes which can complicate matters.

These regions were extract from two versions of Atlantic salmon reference genome (Ssal_v3.1 and ICSASG_v2) so that we can test conversions between those assemblies. The corresponding region in Rainbow trout (Omyk) was also extracted in case we want to test between species as well.

Here are the regions that were extracted (a/b specifies which of the Ss4R duplicates):

Ssal_v3.1_hoxca.fa sequences:

* 13:41179690-41398770 (a)
* 15:92946040-93192783 (b)


ICSASG_v2_hoxca.fa sequences:
* ssa13:35572607-35788439 (a)
* ssa15:88331383-88570007 (b)

Omyk_hoxca.fa sequences:

* 9:64490983-64774759 (b)
* 16:41440669-41645757 (a)

## Extracting the corresponding gff subset

Now that we have the genomic sequence for the toy dataset we need to extract the correspond gene annotations from the gff files. We will use both Ensembl and NCBI gene annotations.

1. Download gff files (use `curl`)
2. Get regions from the fasta headers (note: these are 1-based coordinates)
3. Handle sequence name alias (NCBI uses different names for chromosomes, the fasta files uses Ensembl naming)
4. Extract features from gff falling within these region (use `bedtools intersect`)
5. Shift the coordinates in the gff file to be relative to the region.

How to organize the script(s) for this?

* data/toy-assemblies/get_gff_subset.sh - Main script that does all the above steps. Takes a path/URL of .gff and the subset fasta file and returns a subset gff file. (Alias handling not yet implemented.)
* data/toy-assemblies/get_all_gff_subsets.sh - Just calls the main script for each assembly/annotatoin pair
* convert_seqname_alias.py - Converts the sequence IDs in a tabular file (e.g. gff) using an alias table (TODO)

TODO: Note that there might be features in the gff file that are partially inside the extracted region (e.g. an exon from a gene that stretches outside the region). We might need to "repair" the gff to not contain such partial genes.

