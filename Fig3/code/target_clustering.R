library(tidyverse)
library(here)
library(cowplot)
library(factoextra)


blastout <- read_tsv(here('Fig3', 'data', 'all_vs_all.tsv'),
                     col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'),
                     col_types = 'ccnnnnnnnnnn')

dedup <- blastout %>% 
  distinct(qseqid, sseqid, .keep_all = TRUE)

wide <- dedup %>% 
  select(qseqid, sseqid, bitscore) %>% 
  pivot_wider(names_from = sseqid, values_from = bitscore) %>% 
  mutate(
    across(-qseqid, ~replace_na(.x, 0))
  )

dists <- wide %>%
  column_to_rownames(var = "qseqid") %>%
  dist(method = "euclidean")

clust <- hclust(dists, method = "ward.D2")

ord <- clust$order

hm <- as.data.frame(as.matrix(dists)) %>%
  rownames_to_column(var = "qseqid")

hm$qseqid <- factor(hm$qseqid,
                    levels = c(hm$qseqid[ord])
)

hm <- hm %>%
  pivot_longer(-qseqid, names_to = "sseqid", values_to = "dist")

hm$sseqid <- factor(hm$sseqid,
                    levels = c(hm$sseqid[ord])
)

dendro <- ggdendro::ggdendrogram(clust, leaf_labels = TRUE) +
  theme(
    axis.text.x = element_text(size = 2)
  ) +
  NULL

## manually annotate the dendrogram
