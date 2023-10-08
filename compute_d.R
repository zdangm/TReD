options(warn = -1)

work_dir = '~/TReS/'
setwd(work_dir)
output = 'data/compound_d_immune/'

library(tidyr)
library(data.table)
library(hdf5r)
library(dplyr)
library(foreach)
library(doMC)
registerDoMC(cores=max(detectCores() - 1, 1))


###define function
##twas_sig_rs: a dataframe consisting of twas sig genes' geneid and rank, this rank is a ratio which represents the rank
##drug_sig_rs: a dataframe of drug same as twas_sig_rs
##return: reversal distance, d
compute_d = function(twas_sig_rs, drug_sig_rs){
  twas_rank = twas_sig_rs$twas_rank - 0.5 #vector of twas rank
  drug_rank = drug_sig_rs$drug_rank - 0.5 #vector of drug rank
  twas_rank_norm = sqrt(sum(twas_rank**2)) #norm of twas rank
  drug_rank_norm = sqrt(sum(drug_rank**2)) #norm of drug rank
  cosine_seta = sum(twas_rank*drug_rank) / (twas_rank_norm*drug_rank_norm) #cos(seta)
  d = drug_rank_norm*(-cosine_seta) #the d we need
  return(d)
}

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

  #scale the rank to equal step series
  twas_sig_gene$twas_rank = rank(twas_sig_gene$twas_rank) / nrow(twas_sig_gene)                                
  
  ### start computing rand_drug_d  
  print('start computing rand_drug_d')
  N_PERMUTATIONS <- 1000
  rand_drug_d = foreach(drug_id = 1:ncol(compound_signatures)) %dopar% {
    ## get drug gene rank
    drug_rs = as.data.frame(subset(compound_signatures, select = drug_id))
    drug_rs[1] = rank(-drug_rs[1]) / (dim(drug_rs)[1])  #uniformly distributed
    drug_rs = cbind(gene_list, drug_rs)
    colnames(drug_rs) = c('GeneID', 'drug_rank')
    temp_sig_rs = merge(twas_sig_gene, drug_rs, by = 'GeneID')
    
    ## get null distribution
    #print(paste0('Now we are computing null distribution of ', drug_id))
    null_drug_d_rand1000 = sapply(1:N_PERMUTATIONS, function(i){
      set.seed(i)
      temp_sig_rs$drug_rank = rank(temp_sig_rs$drug_rank) / nrow(temp_sig_rs)
      temp_sig_rs$drug_rank = sample(temp_sig_rs$drug_rank, length(temp_sig_rs$drug_rank))
      twas_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'twas_rank'))
      drug_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'drug_rank'))
      return(compute_d(twas_sig_rs, drug_sig_rs))
    })
    return(null_drug_d_rand1000)
  }
  
  ### start computing drug_d
  print('start computing drug_d')
  drug_d = sapply(1:ncol(compound_signatures), function(drug_id){
    #print(paste0('Now we are computing d of ', drug_id))
    drug_rs = as.data.frame(subset(compound_signatures, select = drug_id))
    drug_rs[1] = rank(-drug_rs[1]) / (dim(drug_rs)[1])
    drug_rs = cbind(gene_list, drug_rs)
    colnames(drug_rs) = c('GeneID', 'drug_rank')
    temp_sig_rs = merge(twas_sig_gene, drug_rs, by = 'GeneID')
    temp_sig_rs$drug_rank = rank(temp_sig_rs$drug_rank)/nrow(temp_sig_rs)
    twas_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'twas_rank'))
    drug_sig_rs = subset(temp_sig_rs, select = c('GeneID', 'drug_rank'))
    return(compute_d(twas_sig_rs, drug_sig_rs))
  })
  
  ### summary                               
  result = data.frame(pert_id = colnames(compound_signatures),
                      drug_d = drug_d,
                      data_set = dataset)
  result$permutation_p = sapply(1:ncol(compound_signatures),function(i) {
    d = result$drug_d[i]
    length(which(rand_drug_d[[i]] >= d)) / length(rand_drug_d[[i]])
  })
  
  if(!dir.exists(paste0(output, dataset, '/'))){
    dir.create(paste0(output, dataset, '/'))
  }
  write.table(result, file = paste0(output, dataset, '/reversal_d.txt'), sep = '\t')                                
}

