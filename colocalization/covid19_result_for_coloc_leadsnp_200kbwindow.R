
library(data.table)
library(coloc) ####### must install latest coloc version 5.2.1 from https://github.com/chr1swallace/coloc
library(stringr)
library(dplyr)

setwd("/home/hecf/")
########## eqtl is built with GRCh37
eqtl_raw_data<-fread("/data/share/eQTL_all/Whole_Blood_cis_eqtl_with_snp_info.txt.gz")
####### nrow(unique(eqtl_raw_data[,'Gene'])
#######  range(eqtl_raw_data[,'MAF'])

############## gwas is built with GRCh38
gwas_raw_data<-fread("/home/hecf/COVID19_HGI_B2_ALL_leave_23andme_20220403_GRCh37.tsv.gz")
colnames(gwas_raw_data)[which(colnames(gwas_raw_data)=='#CHR')]<-'chr'
########## nrow(gwas_raw_data)
i_gene_gwas_case<-44986
i_gene_gwas_allsample<-2356386
i_gene_eqtl_allsample<-573


########### get position annotation for genes
gene_pos_annotation<-fread("/home/hecf/coloc/gencode.v32.GRCh37.txt")

########### get twas significant gene for coloc
twas_blood_sig<-fread("/home/hecf/coloc/sig_fdr_twas_blood.txt")
eqtl_gene<-twas_blood_sig[,'gene']

#i_gene<-2

set.seed(1)
gene_coloc.res_results<-c()
gene_coloc.res_summary<-c()

i_gene<-which(eqtl_gene=="ENSG00000183625")

for(i_gene in 1:nrow(eqtl_gene))
{
  ############### get eqtl data for the number i gene's lead snp 200kb window
  i_eqtl_gene_data<-eqtl_raw_data[Gene==eqtl_gene[i_gene,],c("Gene","chr","rsid","pos","eff_allele","ref_allele","MAF","N_sample","SE","beta",'pval')]
  i_eqtl_gene_data<-i_eqtl_gene_data[pos>=i_eqtl_gene_data[pval==min(i_eqtl_gene_data[,pval]),]$pos-200000 & pos<=i_eqtl_gene_data[pval==min(i_eqtl_gene_data[,pval]),]$pos+200000,]
  colnames(i_eqtl_gene_data)<-paste0(colnames(i_eqtl_gene_data),'.eqtl')
  
  ################### select tmp gwas data for quick merge to get common snp for gwas and eqtl data in a specific gene' lead snp 200kb window
  tmp_select_gwas_data<-gwas_raw_data[chr==as.numeric(i_eqtl_gene_data[1,'chr.eqtl']) & POS>=min(i_eqtl_gene_data[,'pos.eqtl']) & POS<=max(i_eqtl_gene_data[,'pos.eqtl']) ,c("rsid","POS","ALT","REF","all_inv_var_meta_sebeta","all_inv_var_meta_beta","all_inv_var_meta_cases","all_inv_var_meta_controls","all_meta_AF")]
  colnames(tmp_select_gwas_data)<-paste0(colnames(tmp_select_gwas_data),'.gwas')
  
  #################### get common snp for gwas and eqtl data
  merge_gwas_and_eqtl<-merge(tmp_select_gwas_data,i_eqtl_gene_data,by.x=c('rsid.gwas'),by.y=c('rsid.eqtl'))
  i_gene_data<-merge_gwas_and_eqtl
  
  if(nrow(i_gene_data)>0)
  {
    ########################### flip
    ############ make sure the beta direction consistence for eqtl, gwas and LD matrix (r direction)
    i_gene_data<-as.data.frame(i_gene_data)
    i_gene_data[,'beta.eqtl']<-ifelse(i_gene_data[,'eff_allele.eqtl']==i_gene_data[,'ALT.gwas'],i_gene_data[,'beta.eqtl'],-1*i_gene_data[,'beta.eqtl'])
    
    
    ########### calculate sdY for gwas
    i_gene_gwas_p<-i_gene_gwas_case/(i_gene_gwas_allsample)
    i_gene_gwas_sdY<- (i_gene_gwas_p*(1-i_gene_gwas_p))^0.5
    
    
    ################# build gwas data for coloc
    i_gene_gwas<-list(beta=i_gene_data[,"all_inv_var_meta_beta.gwas"],varbeta=i_gene_data[,"all_inv_var_meta_sebeta.gwas"]^2,snp=i_gene_data[,"rsid.gwas"],position=i_gene_data[,"POS.gwas"],type="cc",N=i_gene_gwas_allsample,sdY=i_gene_gwas_sdY)
    ################# build eqtl data for coloc
    i_gene_eqtl<-list(beta=i_gene_data[,"beta.eqtl"],varbeta=i_gene_data[,"SE.eqtl"]^2,snp=i_gene_data[,"rsid.gwas"],position=i_gene_data[,"POS.gwas"],type="quant",N=i_gene_eqtl_allsample,MAF=ifelse(i_gene_data[,"all_meta_AF.gwas"]<0.5,i_gene_data[,"all_meta_AF.gwas"],1-i_gene_data[,"all_meta_AF.gwas"]))
    
    #check_dataset(i_gene_gwas,suffix='.gwas',req=c('snp'),warn.minp=1e-06)
    #check_dataset(i_gene_eqtl,suffix='.eqtl',req=c('snp'),warn.minp=1e-06)
    
    
    ############### get result for coloc.abf
    i_gene_coloc.res<-coloc.abf(i_gene_gwas,i_gene_eqtl)
    
    if(length(i_gene_coloc.res)>0)
    {
      i_gene_coloc.res_results<-i_gene_coloc.res$results%>% filter(SNP.PP.H4 > 0.5)
      i_gene_coloc.res_summary<-as.data.frame(t(i_gene_coloc.res$summary))
      if(nrow(i_gene_coloc.res_results)>0)
      {
        i_gene_coloc.res_results$genename<-as.character(gene_pos_annotation[geneid==eqtl_gene[i_gene,],'genename'])
        gene_coloc.res_results<-rbind(gene_coloc.res_results,i_gene_coloc.res_results)
        
      }
      
      i_gene_coloc.res_summary$genename<-as.character(gene_pos_annotation[geneid==eqtl_gene[i_gene,],'genename'])
      gene_coloc.res_summary<-rbind(gene_coloc.res_summary,i_gene_coloc.res_summary)
      
    }
  }
  print(paste0("done coloc for:",eqtl_gene[i_gene,]," genename:",as.character(gene_pos_annotation[geneid==eqtl_gene[i_gene,],'genename'],"_end\n")))
}
write.table(gene_coloc.res_results,file="/home/hecf/coloc/blood/200kbwindow_result_for_blood.txt",row.names = F,sep="\t")
write.table(gene_coloc.res_summary,file="/home/hecf/coloc/blood/200kbwindow_result_for_blood_summary.txt",row.names = F,sep="\t")