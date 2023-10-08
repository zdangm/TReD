options(warn = -1)
setwd('~/TReS')

library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(data.table)

datasets = c('blood', 'lung', 'lymphocytes', 'spleen', 'ALV', 'EXP', 'BALF')
for (dataset in datasets){
  if (dataset %in% c('blood', 'lung', 'lymphocytes', 'spleen')){
    twas_rs = as.data.frame(fread(paste0('TWAS_rs/twas_', dataset, '.txt')))
    # convert name into entrez ID
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



