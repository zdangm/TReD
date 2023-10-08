options(warn = -1)

work_dir = '~/TReS/'
setwd(work_dir)

library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(data.table)
library(hdf5r)
library(dplyr)
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))


## load compound signatures
gene_list = as.data.frame(fread('row.csv'))
compound_signatures = as.data.frame(fread('compound_signatures/exp_mat.csv'))
row.names(compound_signatures) = gene_list$rid

## load cell name
cell_name = as.character(sapply(colnames(compound_signatures), function(x) strsplit(x, "[_]")[[1]][2]))
immu_cell = as.data.frame(fread('beta_immuncelline_revised.txt'))
compound_signatures = compound_signatures[,which(cell_name %in% immu_cell$cell_iname)] #select immune cell lines

datasets = c('blood', 'lung', 'lymphocytes', 'spleen', 'ALV', 'EXP', 'BALF')
for (dataset_index in 1:7){
  dataset = datasets[dataset_index]
  print(dataset)

  if (dataset %in% c('blood', 'lung', 'lymphocytes', 'spleen')){
    ## load TWAS rs
    load(paste0('/TWAS_rs/twas_', dataset, '_new.Rdata'))
    twas_rs = as.data.frame(twas_rs)
    
    #remove duplicated genes
    twas_rs = twas_rs[order(twas_rs$padj),]
    twas_rs = twas_rs[!duplicated(twas_rs$gene),]
    
    twas_rs$twas_rank = rank(-twas_rs$zscore) / (dim(twas_rs)[1]) #ratio of rank
    twas_sig_gene = subset(twas_rs, padj < 0.05, select = c('GeneID', 'twas_rank'))
  } else{
    ## load dge sig gene list, to reuse twas's code, so named twas_sig_gene
    load(paste0('/DGE_rs/', dataset, '_new.Rdata'))
    dge_rs = as.data.frame(dge_rs)
    
    if(dataset == "BALF"){
      gene_name_col = 'symbol'
      fold_change_col = 'log2FoldChange'
    }else{
      gene_name_col = 'GeneName'
      fold_change_col = 'log2FoldChange'
    }
    
    #remove duplicated genes
    dge_rs = dge_rs[order(abs(dge_rs[,fold_change_col]), decreasing = T),]
    dge_rs = dge_rs[!duplicated(dge_rs[,gene_name_col]),]
    
    twas_sig_gene = subset(dge_rs, select = c('GeneID', 'dge_rank'))
    colnames(twas_sig_gene)[2] = 'twas_rank'
  }
}

#scale the rank to equal step series
twas_sig_gene$twas_rank = rank(twas_sig_gene$twas_rank) / nrow(twas_sig_gene)                                




