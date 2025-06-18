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

./get_gff_subset.sh ICSASG_v2_hoxca.fa https://ftp.ensembl.org/pub/release-104/gff3/salmo_salar/Salmo_salar.ICSASG_v2.104.gff3.gz ICSASG_v2_hoxca_Ens.gff
./get_gff_subset.sh Ssal_v3.1_hoxca.fa https://ftp.ensembl.org/pub/release-106/gff3/salmo_salar/Salmo_salar.Ssal_v3.1.106.gff3.gz Ssal_v3.1_hoxca_Ens.gff
./get_gff_subset.sh Omyk_hoxca.fa https://ftp.ensembl.org/pub/release-107/gff3/oncorhynchus_mykiss/Oncorhynchus_mykiss.USDA_OmykA_1.1.107.gff3.gz Omyk_hoxca_Ens.gff
