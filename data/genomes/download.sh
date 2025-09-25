
curl -O https://ftp.ensembl.org/pub/release-104/fasta/salmo_salar/dna/Salmo_salar.ICSASG_v2.dna_sm.toplevel.fa.gz
curl -O https://ftp.ensembl.org/pub/release-104/gff3/salmo_salar/Salmo_salar.ICSASG_v2.104.gff3.gz
curl -O https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz
curl -O https://ftp.ensembl.org/pub/release-106/gff3/salmo_salar/Salmo_salar.Ssal_v3.1.106.gff3.gz

mkdir -p AtlanticSalmon
mv Salmo_salar.ICSASG_v2.dna_sm.toplevel.fa.gz AtlanticSalmon/ICSASG_v2.fa.gz
mv Salmo_salar.ICSASG_v2.104.gff3.gz AtlanticSalmon/ICSASG_v2_Ens.gff3.gz
mv Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz AtlanticSalmon/Ssal_v3.1.fa.gz
mv Salmo_salar.Ssal_v3.1.106.gff3.gz AtlanticSalmon/Ssal_v3.1_Ens.gff3.gz
#gunzip AtlanticSalmon/*.gz
