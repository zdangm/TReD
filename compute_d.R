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
  ## load TWAS results or DGE results
  dataset = datasets[dataset_index]
  print(dataset)

  if (dataset %in% c('blood', 'lung', 'lymphocytes', 'spleen')){
    twas_rs = as.data.frame(fread(paste0('TWAS_rs/twas_', dataset, '.txt')))
    # convert gene name into entrez ID
    name_index = c()
    for (i in 1:length(twas_rs$genename)){
      try({
        mget(x=as.character(twas_rs$genename[i]), envir=org.Hs.egALIAS2EG)
        name_index = append(name_index, i)
      })
    }
    twas_rs = twas_rs[name_index,]
    twas_rs$GeneID <- mget(x=as.character(twas_rs$genename), envir=org.Hs.egALIAS2EG)
    twas_rs$GeneID <- gsub('[c()"]', '', twas_rs$GeneID) #format text of genes with multiple IDs
    twas_rs <- separate_rows(twas_rs, GeneID)
    twas_rs$padj = p.adjust(twas_rs$pvalue, method = 'fdr')
    save(twas_rs, file = paste0('TWAS_rs/twas_', dataset, '_new.Rdata')) 
  }
   
  
}






