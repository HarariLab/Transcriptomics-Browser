library(VennDiagram)
DEfunction <- function(model_name){
  status = FALSE
  if (model_name == "Sex+Age+Status1"){
    mend_sunshine.clin_final <- mend_sunshine.clin_final[match(colnames(cts_sub1), mend_sunshine.clin_final$Sample_ID),]#making sure the colnames and sample_ids match 
    all(mend_sunshine.clin_final$Sample_ID == colnames(cts_sub1))##this must be TRUE to continue
    model_file <- paste0(model_name, ".csv")#filename that will be written at the end 
    #DESeq args 
    counts = cts_sub1 
    cols = mend_sunshine.clin_final
    designs = ~ Sex + AOD + Status
    status = TRUE 
  }
  if (model_name == "Sex+Age+PMI+Status1"){
    model_file <- paste0(model_name, ".csv")
    mend_sunshine.clin_PMI4 <- mend_sunshine.clin_PMI4[match(colnames(cts_sub_PMI2), mend_sunshine.clin_PMI4$Sample_ID),]
    all(mend_sunshine.clin_PMI4$Sample_ID == colnames(cts_sub_PMI2))
    counts = cts_sub_PMI2
    cols = mend_sunshine.clin_PMI4
    designs = ~ Sex + AOD + PMI + Status
    status = TRUE
  }
  if (model_name == "Sex+Age+PMI+Cohort+Status1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_PMI4$Sample_ID == colnames(cts_sub_PMI2))
    counts = cts_sub_PMI2
    cols = mend_sunshine.clin_PMI4
    designs = ~ Sex + AOD + PMI + Cohort + Status
    status = TRUE
  }
  if (model_name == "Sex+Age+Cohort+Status1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_final$Sample_ID == colnames(cts_sub1))
    counts = cts_sub1
    cols = mend_sunshine.clin_final
    designs = ~ Sex + AOD + Cohort + Status
    status = TRUE
  }
  if (model_name == "Sex+Age+CDR1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_final$Sample_ID == colnames(cts_sub1))
    counts = cts_sub1
    cols = mend_sunshine.clin_final
    designs = ~ Sex + AOD + CDRe
  }
  if (model_name == "Sex+Age+PMI+CDR1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_PMI4$Sample_ID == colnames(cts_sub_PMI2))
    counts = cts_sub_PMI2
    cols = mend_sunshine.clin_PMI4
    designs = ~ Sex + AOD + PMI + CDRe
  }
  if (model_name == "Sex+Age+PMI+Cohort+CDR1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_PMI4$Sample_ID == colnames(cts_sub_PMI2))
    counts = cts_sub_PMI2
    cols = mend_sunshine.clin_PMI4
    designs = ~ Sex + AOD + PMI + Cohort + CDRe
  }
  if (model_name == "Sex+Age+Cohort+CDR1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_final$Sample_ID == colnames(cts_sub1))
    counts = cts_sub1
    cols = mend_sunshine.clin_final
    designs = ~ Sex + AOD + Cohort + CDRe
  }
  if (model_name == "Sex+Age+Braak1"){
    model_file <- paste0(model_name, ".csv")
    mend_sunshine.clin_Braak4 <- mend_sunshine.clin_Braak4[match(colnames(cts_sub_Braak2), mend_sunshine.clin_Braak4$Sample_ID),]
    all(mend_sunshine.clin_Braak4$Sample_ID == colnames(cts_sub_Braak2))
    counts = cts_sub_Braak2
    cols = mend_sunshine.clin_Braak4
    designs = ~ Sex + AOD + BraakTau
  }
  if (model_name == "Sex+Age+PMI+Braak1"){
    model_file <- paste0(model_name, ".csv")
    mend_sunshine.clin_PMI_Braak4 <- mend_sunshine.clin_PMI_Braak4[match(colnames(cts_sub_Braak_PMI3), mend_sunshine.clin_PMI_Braak4$Sample_ID),]
    all(mend_sunshine.clin_PMI_Braak4$Sample_ID == colnames(cts_sub_Braak_PMI3))
    counts = cts_sub_Braak_PMI3
    cols = mend_sunshine.clin_PMI_Braak4
    designs = ~ Sex + AOD + PMI + BraakTau
  }
  if (model_name == "Sex+Age+PMI+Cohort+Braak1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_PMI_Braak4$Sample_ID == colnames(cts_sub_Braak_PMI3))
    counts = cts_sub_Braak_PMI3
    cols = mend_sunshine.clin_PMI_Braak4
    designs = ~ Sex + AOD + PMI + Cohort + BraakTau
  }
  if (model_name == "Sex+Age+Cohort+Braak1"){
    model_file <- paste0(model_name, ".csv")
    all(mend_sunshine.clin_Braak4$Sample_ID == colnames(cts_sub_Braak2))
    counts = cts_sub_Braak2
    cols = mend_sunshine.clin_Braak4
    designs = ~ Sex + AOD + Cohort + BraakTau
  }
  dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = cols,
                               design = designs)
  dds <- DESeq(dds)
  
  ##to assign the results of the DESeq to a variable that can be referenced
  temp1 <- assign(paste0(model_name, ".dds"),0) 
  temp1 <- paste0(model_name, ".dds")
  assign(temp1, dds, envir = .GlobalEnv)
  
  if (status == TRUE){
    res <- results(dds, contrast = c("Status", "Neuro_AD", "Neuro_CO"), alpha = 0.05)
  }
  if (status == FALSE){
    res <- results(dds)
  }
  
  ##writing DESeq results to csv 
  gg <- as.data.frame(str_split_fixed(rownames(res), "_", 1))
  colnames(gg) <- 'GeneID'
  res <- cbind(res, gg)
  res <- as.data.frame(merge(meta.cols, res, by = "GeneID"))
  write.csv(as.data.frame(res),
            file = model_file, quote = FALSE, row.names = FALSE)
  
  ##change these for whichever dataset you're working with at the time
  # df_sex_age_PMI_CDR <- read.csv('/home/sohn/sex_age_PMI_CDR.csv', header =T, row.names =1, sep=",", stringsAsFactors = F, check.names = F)
  #df_sex_age_PMI_CDR <- merge(df_sex_age_PMI_CDR, meta.cols)
  rownames(res) <- res$GeneID
  res_padj <- subset(res, padj<0.05)
  res_pval <- subset(res, pvalue<0.05)
  pval_genesInterest_res <- intersect(res_pval$GeneName, genesInterest)
  padj_genesInterest_res <- intersect(res_padj$GeneName, genesInterest)
  
}

DEfunction("Sex+Age+Braak1")

