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






