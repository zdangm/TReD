options(warn = -1)
setwd('~/TReS/compound_d_immune/')

library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)


annotation_df = as.data.frame(fread("annotation_df.txt"))
siginfo_beta = as.data.frame(fread('siginfo_beta.txt'))

datasets = c('blood', 'lung', 'lymphocytes', 'spleen', 'ALV', 'EXP', 'BALF')
for (i in 1:7){
  dataset = datasets[i]
  print(dataset)

  immune_instances_rs = as.data.frame(fread())
}
