library("plyr")
library("dplyr")
library("ggpubr")
library("ggsci")
library("tidyverse")
library("survival")
library("survminer")
library("data.table")

# run these in python
# CM_info = pd.read_csv('CM_info.csv', index_col=0)
# adata = sc.read_h5ad('tmp.h5ad')
# adata_int = adata[adata.obs['subCluster'].isin(list(CM_info['cellanno']))]
# adata_int = adata_int[adata_int.obs['Tissue']=='Rectum_T']

# sc.tl.rank_genes_groups(adata_int, 'subCluster', use_raw=False, pts=True, method='wilcoxon', key_added=f"rank_genes_groups.subCluster")
# DE_result = sc.get.rank_genes_groups_df(adata_int, key=f"rank_genes_groups.subCluster", pval_cutoff=0.05, log2fc_min=0.25, group=None)
# DE_result = DE_result.sort_values(by=['group','scores'], ascending=[True,False])
# DE_result = DE_result.sort_values(['group', 'logfoldchanges'], ascending=[True, False])
# DE_result = DE_result[(DE_result['pct_nz_group']>0.25) & (DE_result['pct_nz_reference']<0.25)]

# def get_top_n_rows(group, n=10):
#     return group.sort_values(by='scores', ascending=False).head(n)
# top_n_per_group = DE_result.groupby('group').apply(get_top_n_rows)
# top_n_per_group = top_n_per_group.reset_index(drop=True)
# res = pd.merge(CM_info, top_n_per_group, how='left', left_on='cellanno', right_on='group')
# res_dict = res.groupby('CM')['names'].apply(set).to_dict()
# df = res[['CM', 'names']]
# df.columns=['geneset', 'genesymbol']
# # Remove duplicated entries
# df = df[~df.duplicated(['geneset', 'genesymbol'])]
# len(df['geneset'].unique())
# df['geneset'] = 'CM'+df['geneset'].astype('str')
# df.to_csv('CM_signature.csv')

TCGA_data <- readRDS('TCGA_data.rds')
TCGA_data <- TCGA_data[,TCGA_data$sampleType=="Primary.Solid.Tumor"]
assay(TCGA_data,'exprs') <- log1p(assay(TCGA_data,'tpm'))
# geneset_list: loaded from CM_signature.csv
score <- gsva(assay(TCGA_data,'exprs'), geneset_list, method="gsva", kcdf="Gaussian", verbose=F, parallel.sz=40)
score <- score %>% t() %>% data.frame()
score$sample <- TCGA_data$sample
score <- score %>% left_join(as.data.frame(colData(TCGA_data)), by='sample')

sample_info <- read.table(sprintf("TCGA_surdf.csv",dir), header=T, stringsAsFactors=F, check.names=F, sep=",")
sample_info <- sample_info[, c("sample","X_PATIENT","cancer.type.abbreviation","age_at_initial_pathologic_diagnosis", "gender", "race", "ajcc_pathologic_tumor_stage", "OS","OS.time","sampleType","OS.5","OS.3","OS.time.5","OS.time.3" )]
dat <- left_join(score, sample_info, by='sample')

dat = dat[!is.na(dat$OS.time),]
dat = dat[!is.na(dat$age_at_initial_pathologic_diagnosis),]
dat$age = dat$age_at_initial_pathologic_diagnosis
dat$cancerType = dat$cancer.type.abbreviation
#dat$OS.time = round(dat$OS.time/30, 2)
dat = dat[!is.na(dat$ajcc_pathologic_tumor_stage),]
dat$stage = dat$ajcc_pathologic_tumor_stage
unique(dat$stage)
dat = dat[grepl("Stage I",dat$stage),]
unique(dat$stage)
dat$stage = gsub("[ABCD]$","",dat$stage)
unique(dat$stage)

int=c('CM2', 'CM3', 'CM6', 'CM7')
for (i in int){
    pdat <-dat
    res.cut <- surv_cutpoint(pdat, time = "OS.time.use", 
                             event = "OS.use", 
                             variables = i)
    res.cat <- surv_categorize(res.cut)
    pdat$group <- res.cat[,i] 
    pdat$group = factor(pdat$group)

    fit = survfit(Surv(OS.time.use,OS.use) ~ group, data=pdat)

    p = ggsurvplot(fit, pdat, size=0.3, vlegend.labs=unique(pdat$group),
                   surv.median.line="none", pval=T, conf.int=F,
                   palette=c("#00726A", "grey")) + 
                   xlab("Years")

    ggsave(filename = sprintf('%s.KM_plot.pdf', i), p, width = 3,height = 4)
}