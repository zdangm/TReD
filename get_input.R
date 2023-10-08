options(warn = -1)
setwd('~/TReS/data/')

library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(data.table)


### get TWAS rs and DGE rs
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
  } else if (dataset == "ALV"){
    dz_signature <- read.csv("/DGE_rs/ALV_DE.csv")
    dz_signature = na.omit(dz_signature)
    dz_signature$dge_rank = rank(-dz_signature$log2FoldChange) / length(dz_signature$log2FoldChange)
    dz_signature <- dz_signature[which(dz_signature$padj < 0.05),]
    dz_signature <- dz_signature[order(dz_signature$log2FoldChange),]
    dz_signature = dz_signature[-which(dz_signature$GeneName == 'LOC284454'),]
    dz_signature$GeneID <- mget(x=as.character(dz_signature$GeneName), envir=org.Hs.egALIAS2EG) #convert to entrez ID
    dz_signature$GeneID <- gsub('[c()"]', '', dz_signature$GeneID) #format text of genes with multiple IDs
    dz_signature <- separate_rows(dz_signature, GeneID) #genes with multiple IDs, map to all of them as individual rows
    dge_rs = dz_signature
    save(dge_rs, file = paste0('DGE_rs/', dataset, '_new.Rdata'))
  } else if (dataset == "EXP"){
    dz_signature <- read.csv("/DGE_rs/EXP_DE.csv")
    colnames(dz_signature)[c(4,10)] <- c("log2FoldChange", "GeneID") #Relabel columns
    dz_signature = na.omit(dz_signature)
    dz_signature$dge_rank = rank(-dz_signature$log2FoldChange) / length(dz_signature$log2FoldChange)
    dz_signature <- dz_signature[which(dz_signature$padj < 0.05),]
    dz_signature <- dz_signature[which(abs(dz_signature$log2FoldChange) > 2),]
    dz_signature <- dz_signature[order(dz_signature$log2FoldChange),]
    dz_signature <- dz_signature[!is.na(dz_signature$GeneID),] # keep genes with valid entrez ID
    dz_signature <- dz_signature[which(dz_signature$GeneID %in% gene_list$V1),] # keep genes in cmap
    dge_rs = dz_signature
    save(dge_rs, file = paste0('DGE_rs/', dataset, '_new.Rdata'))
  } else if (dataset == "BALF"){
    dz_signature <- read.csv("/DGE_rs/BALF_DE.csv")
    colnames(dz_signature)[c(4,9)] <- c("log2FoldChange", "GeneID") #Relabel columns: estimate is log2FC
    dz_signature = na.omit(dz_signature)
    dz_signature$dge_rank = rank(-dz_signature$log2FoldChange) / length(dz_signature$log2FoldChange)
    dz_signature <- dz_signature[which(dz_signature$p.adjusted < 0.05),]
    dz_signature <- dz_signature[which(abs(dz_signature$log2FoldChange) > 4),]
    dz_signature <- dz_signature[order(dz_signature$log2FoldChange),]
    dz_signature <- dz_signature[!is.na(dz_signature$GeneID),] # keep genes with valid entrez ID
    dz_signature <- dz_signature[which(dz_signature$GeneID %in% gene_list$V1),] # keep genes in cmap
    dge_rs = dz_signature
    save(dge_rs, file = paste0('DGE_rs/', dataset, '_new.Rdata'))
  }  
}

### get drug annotations
compound_info_beta = as.data.frame(fread('/LINCS2_beta/compoundinfo_beta.txt'))

drugbank_info = fread("/drugbank/drugbank vocabulary_5.1.10.csv")
drugbank_info$'Common name'<-toupper(drugbank_info$'Common name')
drugbank_info$Synonyms<-toupper(drugbank_info$Synonyms)

annotation_df = compound_info_beta[, c('pert_id', 'inchi_key', 'cmap_name', 'compound_aliases', 'target', 'moa')]
## set "not" to no nchikey
annotation_df$inchi_key = as.character(annotation_df$inchi_key)
annotation_df$inchi_key = ifelse(nchar(annotation_df$inchi_key) > 0, annotation_df$inchi_key, "not")
## upper
annotation_df$cmap_name = toupper(annotation_df$cmap_name)
annotation_df$compound_aliases = toupper(annotation_df$compound_aliases)
## annotate
annotation_df$inchikey_annot = drugbank_info$'DrugBank ID'[match(annotation_df$inchi_key, drugbank_info$'Standard InChI Key')]
annotation_df$cmap_name_annota2Common_name = drugbank_info$'DrugBank ID'[match(annotation_df$cmap_name, drugbank_info$'Common name')]
annotation_df$compound_aliases_annota2Common_name = drugbank_info$'DrugBank ID'[match(annotation_df$compound_aliases, drugbank_info$'Common name')]
annotation_df$drugbank_id = NA
annotation_df = as.data.frame(unique(annotation_df))
## combine result
for(j in 1:nrow(annotation_df)){
  tmp<-unname(annotation_df[j, c("inchikey_annot","cmap_name_annota2Common_name","compound_aliases_annota2Common_name")])
  tmp<-tmp[!is.na(tmp)]
  tmp<-tmp[tmp!="NULL"]
  tmp<-unique(tmp)
  if(length(tmp)>0){
    annotation_df$drugbank_id[j]<-tmp
  }
}
temp_drugbank_info = drugbank_info[, c('DrugBank ID', 'Common name')]
colnames(temp_drugbank_info) = c('drugbank_id', 'common_name')
annotation_df = merge(annotation_df, temp_drugbank_info, by = 'drugbank_id')
annotation_df = annotation_df[, -match(c('inchikey_annot', 'cmap_name_annota2Common_name', 'compound_aliases_annota2Common_name'), colnames(annotation_df))]
write.table(annotation_df, file = '/LINCS2_beta/annotation_df.txt', sep = '\t', row.names = F)

