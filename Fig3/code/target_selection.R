library(tidyverse)
library(here)
library(conflicted)
library(biomaRt)
library(janitor)

conflict_prefer("filter", "dplyr")
conflict_prefer('select', 'dplyr')


# establish the BioMart that will be used
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# get the paralogs to the known human targets
human_targets <- read_csv(here('Fig3/metadata/tocris_targets.csv'),
                          col_names = 'gene_name')

# biomart parameters
filters <- c("gene_name", 'biotype', 'species_id_1010')
value <- list(human_targets$gene_name, 'protein_coding', 'hsapiens')
attributes <- c("wbps_gene_id", 
                "external_gene_id", 
                "wbps_paralog_gene")

human_target_paralogs <- getBM(mart = mart, 
                               filters = filters,
                               value = value,
                               attributes = attributes) %>%
  clean_names() %>% 
  distinct(wbps_paralog_gene, .keep_all = TRUE)

full_human_list <- bind_rows(dplyr::select(human_target_paralogs, gene_id = wbps_gene_id),
                             dplyr::select(human_target_paralogs, gene_id = wbps_paralog_gene)) %>% 
  distinct()

# get the gene names too
filters <- c("wbps_gene_id", 'species_id_1010')
value <- list(full_human_list$gene_id, 'hsapiens')
attributes <- c("wbps_gene_id", 
                "external_gene_id")

full_human_list_name <- getBM(mart = mart, 
                              filters = filters,
                              value = value,
                              attributes = attributes) %>%
  clean_names() %>% 
  distinct(wbps_gene_id, .keep_all = TRUE) %>%
  write_csv(here('Fig3/metadata/human_target_dataset.csv'))

# get the Bm orthologs
filters <- c("gene_name", 'only_brmalaprjna10729_homologue', 'species_id_1010')
value <- list(full_human_list_name$external_gene_id, TRUE, 'hsapiens')
attributes <- c(
  'external_gene_id',
  "wbps_gene_id", 
  "brmalaprjna10729_gene", 
  "brmalaprjna10729_gene_name", 
  "brmalaprjna10729_homolog_ensembl_peptide", 
  "brmalaprjna10729_homolog_perc_id")

# pull the desired data
tocris_orthologs <- getBM(mart = mart, 
                          filters = filters,
                          value = value,
                          attributes = attributes) %>%
  clean_names() %>% 
  arrange(brmalaprjna10729_gene_name, desc(brmalaprjna10729_homolog_perc_id)) %>% 
  distinct(brmalaprjna10729_gene, .keep_all = TRUE) %>%
  write_csv(here('Fig3/metadata/brugia_target_dataset.csv'))


# assign brugia clusters ---------------------------------------------------------

id_map <- read_tsv(here('Fig3', 'metadata', 'id_map.tab')) %>% 
  janitor::clean_names() %>% 
  dplyr::select(-protein_names) %>% 
  mutate(gene_name = str_remove(gene_names, '\\s.*$'))

clusters <- read_csv(here('Fig3', 'metadata', 'human_clusters.csv')) %>% 
  janitor::clean_names() %>% 
  rename(target = uni_prot_kb) %>% 
  left_join(id_map)

ortholog_cluster <- tocris_orthologs %>% 
  left_join(clusters, by = c('external_gene_id' = 'gene_name')) %>% 
  dplyr::select('wbps_gene_id', contains('brmala'), cluster) %>% 
  write_csv(here('Fig3', 'metadata', 'brugia_clusters.csv'))




