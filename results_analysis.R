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

  immune_instances_rs = as.data.frame(fread(paste0(dataset, '/reversal_d.txt')))
  immune_instances_rs = immune_instances_rs[, -which(colnames(immune_instances_rs) == 'V1')]
  colnames(immune_instances_rs)[which(colnames(immune_instances_rs) == 'pert_id')] = 'sig_id'
  immune_instances_rs = merge(immune_instances_rs, siginfo_beta[, c('sig_id', 'pert_id')], by = 'sig_id')
  immune_drug_instances_rs = merge(immune_instances_rs, unique(annotation_df[, c('pert_id', 'common_name')]))
  write.table(immune_drug_instances_rs, file = paste0('immune_drug_instances/', dataset, '.csv'), row.names = F, sep = ',')

  ### subset d>mean+2sd & p<0.05
  sub_immune_drug_instances_rs = immune_drug_instances_rs[which(immune_drug_instances_rs$drug_d > (mean(immune_drug_instances_rs$drug_d) + 2*sd(immune_drug_instances_rs$drug_d))),]
  sub_immune_drug_instances_rs = sub_immune_drug_instances_rs[which(sub_immune_drug_instances_rs$permutation_p < 0.05),]
  
  ## match pert_idose, pert_iname, cell_iname
  sub_immune_drug_instances_rs = merge(sub_immune_drug_instances_rs, siginfo_beta[, c('sig_id', 'pert_idose', 'pert_itime', 'cell_iname')], by = 'sig_id')
  ## whether effecticive under multi-cell
  pert_ids = sub_immune_drug_instances_rs$pert_id
  sub_immune_drug_instances_rs$cell_bool = sapply(1:length(pert_ids), function(index){
    temp_subrs = sub_immune_drug_instances_rs[which(pert_ids == pert_ids[index]), ]
    bool_value = length(unique(temp_subrs$cell_iname)) >= 2
    return(bool_value)
  })

  ### summary 
  if (i == 1){
    all_sub_immune_drug_instances_rs = sub_immune_drug_instances_rs
  } else{
    all_sub_immune_drug_instances_rs = rbind(all_sub_immune_drug_instances_rs, sub_immune_drug_instances_rs)
  }
  sub_immune_drug_instances_cell = sub_immune_drug_instances_rs[sub_immune_drug_instances_rs$cell_bool, ]
  sub_immune_drug_rs_cell = sub_immune_drug_instances_cell %>% group_by(pert_id) %>% dplyr::slice(which.max(drug_d)) %>% as.data.frame

  if (i == 1){
    all_sub_immune_drug_rs_cell = sub_immune_drug_rs_cell
  } else{
    all_sub_immune_drug_rs_cell = rbind(all_sub_immune_drug_rs_cell, sub_immune_drug_rs_cell)
  }
}
