# gff links:
# Ssal_v3.1:
# https://ftp.ensembl.org/pub/release-106/gff3/salmo_salar/Salmo_salar.Ssal_v3.1.106.gff3.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/237/065/GCF_905237065.1_Ssal_v3.1/GCF_905237065.1_Ssal_v3.1_genomic.gff.gz
#
# ICSASG_v2:
# https://ftp.ensembl.org/pub/release-104/gff3/salmo_salar/Salmo_salar.ICSASG_v2.104.gff3.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.gff.gz
#
# Omyk:
# https://ftp.ensembl.org/pub/release-107/gff3/oncorhynchus_mykiss/Oncorhynchus_mykiss.USDA_OmykA_1.1.107.gff3.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/265/735/GCF_013265735.2_USDA_OmykA_1.1/GCF_013265735.2_USDA_OmykA_1.1_genomic.gff.gz
#
# Fasta files:
# Omyk_hoxca.fa
# Ssal_v3.1_hoxca.fa
# ICSASG_v2_hoxca.fa
#
# alias links:
# https://salmobase.org/datafiles/genomes/AtlanticSalmon/ICSASG_v2/alias_filtered.txt
# https://salmobase.org/datafiles/genomes/AtlanticSalmon/Ssal_v3.1/alias.txt


# filter_complete_genes: Filters a GFF file to retain only complete gene features.
#
# Extracting gff subset may result in orphaned gene features if
# a gene spans the region boundaries (e.g. single exon).
# Uses AGAT to filter based on gene IDs found in the GFF. This
# will keep all features belonging to these genes.
# Note1: Removes non-gene fetures like ncRNA and biological region
# Note2: AGAT will add IDs to every feature that don't have one
filter_complete_genes() {
  local INPUT_GFF=$1

  # get all geneIDs
  grep "\tgene\t" $INPUT_GFF | awk -F'\t' '{print $9}' | \
      sed 's/.*ID=//' | cut -d';' -f1 > $INPUT_GFF.genes.txt
  agat_sp_filter_feature_from_keep_list.pl \
    --gff $INPUT_GFF \
    --keep_list $INPUT_GFF.genes.txt \
    -o $INPUT_GFF.fixed

  # remove temporary files
  rm "${INPUT_GFF}.genes.txt"
  rm "${INPUT_GFF}_report.txt"

  # replace the gff with the fixed
  rm $INPUT_GFF
  mv $INPUT_GFF.fixed $INPUT_GFF
}


./get_gff_subset.sh ICSASG_v2_hoxca.fa https://ftp.ensembl.org/pub/release-104/gff3/salmo_salar/Salmo_salar.ICSASG_v2.104.gff3.gz ICSASG_v2_hoxca_Ens.gff
./get_gff_subset.sh Ssal_v3.1_hoxca.fa https://ftp.ensembl.org/pub/release-106/gff3/salmo_salar/Salmo_salar.Ssal_v3.1.106.gff3.gz Ssal_v3.1_hoxca_Ens.gff
./get_gff_subset.sh Omyk_hoxca.fa https://ftp.ensembl.org/pub/release-107/gff3/oncorhynchus_mykiss/Oncorhynchus_mykiss.USDA_OmykA_1.1.107.gff3.gz Omyk_hoxca_Ens.gff


mkdir -p temp

# download NCBI gff files and convert seqnames to Ensembl alias
curl -o temp/Ssal_v3.1.alias https://salmobase.org/datafiles/genomes/AtlanticSalmon/Ssal_v3.1/alias.txt
curl -s -o - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/237/065/GCF_905237065.1_Ssal_v3.1/GCF_905237065.1_Ssal_v3.1_genomic.gff.gz |\
  gunzip -c | python3 convert_seqname_alias.py temp/Ssal_v3.1.alias > temp/Ssal_v3.1_NCBI.gff
./get_gff_subset.sh Ssal_v3.1_hoxca.fa temp/Ssal_v3.1_NCBI.gff Ssal_v3.1_hoxca_NCBI.gff

curl -o temp/ICSASG_v2.alias https://salmobase.org/datafiles/genomes/AtlanticSalmon/ICSASG_v2/alias.txt
curl -s -o - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.gff.gz |\
  gunzip -c | python3 convert_seqname_alias.py temp/ICSASG_v2.alias > temp/Ssal_ICSASG_v2_NCBI.gff
./get_gff_subset.sh ICSASG_v2_hoxca.fa temp/Ssal_ICSASG_v2_NCBI.gff ICSASG_v2_hoxca_NCBI.gff



# Apply filter_complete_genes to all generated GFF subset files
filter_complete_genes ICSASG_v2_hoxca_Ens.gff
filter_complete_genes Ssal_v3.1_hoxca_Ens.gff
filter_complete_genes Omyk_hoxca_Ens.gff
filter_complete_genes Ssal_v3.1_hoxca_NCBI.gff
filter_complete_genes ICSASG_v2_hoxca_NCBI.gff