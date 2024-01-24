library(tidyverse)
library(readxl)
library(ggbeeswarm)

### Calculating signature scores in bulk data using GSVA
log_abundance <-read.table("bulkRNA_logTPM.tsv", sep="\t")
metadata <- read.table('bulkRNA_metadata.tsv', sep='\t')
json_data <- jsonlite::fromJSON("c2.cp.reactome.v2023.1.Hs.json")
IFNG_signature <- json_data[['REACTOME_INTERFERON_GAMMA_SIGNALING']][['geneSymbols']] 
geneset_list <- list(
    IFNG_signature=IFNG_signature
)

score <- GSVA::gsva(as.matrix(log_abundance), geneset_list, kcdf="Gaussian", verbose=F, parallel.sz=60)
pdat <- 
score %>% 
t() %>% 
as.data.frame() %>% 
rownames_to_column('Ident') %>% 
left_join(metadata)
pdat <-
pdat %>%
mutate_all(~replace_na(., "")) %>%
mutate(int_lab = str_c(SampleTimePoint, '_', Efficacy))

pdat$int_lab <- factor(pdat$int_lab, levels=c('BL_', 'nRCT_non-CR', 'nRCT_CR'))
p <- ggplot(pdat, aes(x=int_lab, y=IFNG_signature)) +
  geom_boxplot(aes(color=int_lab), fill=NA, width=0.8, outlier.shape=NA) +
  geom_quasirandom(aes(color=int_lab), size=1) +
  ggpubr::geom_signif(
    comparisons=list(c('BL_', 'nRCT_non-CR'), c('BL_', 'nRCT_CR'), c('nRCT_non-CR', 'nRCT_CR')),
    step_increase=0.1,
    map_signif_level=F,
    test=t.test,
    textsize=4,
    tip_length = 0.01
  ) +
  labs(title=NULL, x="xlab", y='IFNG_signature') +
  scale_color_manual(values=c('#224b5e', '#3C5488FF', '#00A087FF')) +
  theme_classic() +
  theme(legend.position="none")
ggsave("BL_RCT_IFNG_signature.pdf", p, width=3.5, height=4.5, limitsize=F)

pdat$Efficacy <- factor(pdat$Efficacy, levels=c('CR', 'non-CR'))
roc.list <- pROC::roc(Efficacy ~ IFNG_signature, data = pdat)
p <- pROC::ggroc(roc.list, lwd = 0.65, legacy.axes = TRUE) +
  theme_classic() +
  labs(x = 'False positive rate', y='True positive rate')
ggsave("RCT_CRvsnonCR_IFNG_signature_AUROC.pdf", p, width=3, height=3, limitsize=F)

### Calculating signature scores in single-cell data using scanpy.tl.score_genes
# running these in python
# import scanpy as sc
# for signature in gene_list:
#     sc.tl.score_genes(adata, gene_list[signature], score_name=signature)
# adata.obs.to_csv('signature_score.tsv', sep='\t')
signature_score <- read.table('signature_score.tsv', sep='\t')
gene <- 'IFNG_signature'
input1 <- signature_score
input1$expression <- input1[, gene]

input2 <- 
  input1 %>%
  group_by(group, SampleTimePoint) %>%
  dplyr::summarise(m = mean(expression)) %>%
  ungroup() %>%
  filter(SampleTimePoint %in% c('BL', 'nCT', 'nRT', 'nRCT'))

input3 <- input2 %>% pivot_wider(names_from='SampleTimePoint', values_from='m')

pdat <- as.matrix(input3[,-1])
row.names(pdat) <- input3[[1]]
pdat <- pdat %>% t() %>% scale() %>% t()
pdat <- pdat[c("B", "plasma", "CD4T", "CD8T", "gdT", "ILC", "Macrophage", "DC", "Neutrophil", "Monocyte", 'Platelet'),]

col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#0F7B9F", "white", "#D83215")
)

p <- ComplexHeatmap::Heatmap(pdat,
                            name = gene,
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            row_title = 'Cell type',
                            column_title = "SampleTimePoint",
                            right_annotation = NULL,
                            show_row_names = TRUE,
                            row_names_gp = grid::gpar(fontsize = 10),
                            row_names_side = 'left',
                            border = FALSE, 
                            column_names_side = "bottom",
                            column_names_rot = 45, 
                            column_names_gp = grid::gpar(fontsize = 10),
                            rect_gp = grid::gpar(col="white", lwd=0.5),
                            col = col_fun,
                            row_gap = unit(3, "mm")
                       )

fig.height <- 2 + 0.3*nrow(pdat)
fig.width <- 2 + 0.55*ncol(pdat)
pdf('signature_score_heatmap.pdf', width=fig.width, height=fig.height)
print(p)
dev.off()
