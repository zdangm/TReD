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



