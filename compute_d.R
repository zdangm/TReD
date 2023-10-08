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

### get data
## load compound signatures
gene_list = as.data.frame(fread('row.csv'))
compound_signatures = as.data.frame(fread('compound_signatures/exp_mat.csv'))
row.names(compound_signatures) = gene_list$rid

## load TWAS results and DGE results
args = as.numeric(commandArgs(TRUE)) #which dataset
datasets = c('blood', 'lung', 'lymphocytes', 'spleen', 'ALV', 'EXP', 'BALF')
dataset = datasets[args]
twas_rs = as.data.frame(fread(paste0('TWAS_rs/twas_', dataset, '.txt')))

