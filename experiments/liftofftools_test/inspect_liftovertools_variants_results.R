library(tidyverse)
setwd("~/Dropbox/github/larsgr/gene-id-mapping")
gff_v3 <- read_tsv("data/toy-assemblies/Ssal_v3.1_hoxca_Ens.gff",
         comment = "#",col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
gff_v2 <- read_tsv("data/toy-assemblies/ICSASG_v2_hoxca_Ens.gff",
                   comment = "#",col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

comparison <- read_tsv("experiments/liftofftools_test/variant_effects",
                       col_names=c("reference_transcript_id", "target_transcript_id", "dna_identity",
                                   "protein_identity", "variant_effect" ))



genes_v2 <- 
  gff_v2 %>% filter(type=="gene") %>% 
  transmute(gene=sub(".*ID=([^;]+).*","\\1",attributes),GeneName=sub(".*Name=([^;]+).*","\\1",attributes))

tx_v2 <- 
  gff_v2 %>% filter(type=="mRNA") %>% 
  transmute(tx=sub(".*ID=([^;]+).*","\\1",attributes),gene=sub(".*Parent=([^;]+).*","\\1",attributes)) %>% 
  left_join(genes_v2,by="gene")

comparison %>% 
  left_join(tx_v2, by=c("reference_transcript_id"="tx")) %>% 
  arrange(gene) %>% View()

