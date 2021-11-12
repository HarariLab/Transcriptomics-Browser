library(DESeq2)
## read technical dat
## removed MAP_60180 and MAP_60007 due to irregularity as a control sample 
mend.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2020_09_11_WashU_MendalianVsSporadics_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
sunshine.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_09_16_WashU_SUNSHINE_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
sunshine.tech_KnightADRC <- subset(sunshine.tech, Study == "Knight_ADRC")
sunshine.tech_KnightADRC <- sunshine.tech_KnightADRC[-312, ] # got rid of the unknown sample
sunshine.tech_KnightADRC_sample_names <- sunshine.tech_KnightADRC$Sample_Name
mend.tech_sample_names <- mend.tech$Sample_ID
mend.tech_KnightADRC <- subset(mend.tech, Study == "Knight_ADRC")
mend.tech_KnightADRC_sample_names <- mend.tech_KnightADRC$Sample_ID
mend_sunshine_sample_names <- c(mend.tech_KnightADRC_sample_names, sunshine.tech_KnightADRC_sample_names)
mend_sunshine_sample_names1 <- c(mend_sunshine_sample_names, "GeneID")
asd1 <- mend.tech_KnightADRC[, c("Subj_ID", "Sample_ID")]
asd2 <- sunshine.tech_KnightADRC[, c("Subj_ID", "Sample_Name")]
names(asd2)[2] <- "Sample_ID"
mend_sunshine_sample_name_subj_id <- rbind(asd1,asd2)

## read counts
mend.cts <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_count_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F)
sunshine.cts <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_count_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F)
rownames(mend.cts) <- mend.cts$GeneID
rownames(sunshine.cts) <- sunshine.cts$GeneID
mend.cts_KnightADRC <- mend.cts[ ,mend.tech_KnightADRC_sample_names1]
mend_sunshine.cts <- merge(mend.cts, sunshine.cts, all = TRUE)
mend_sunshine.cts <- mend_sunshine.cts[ ,mend_sunshine_sample_names1]
rownames(mend_sunshine.cts) <- mend_sunshine.cts$GeneID
mend_sunshine.cts <- mend_sunshine.cts[-446] 
mend.tpm <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_norm.tsv',header=T, check.names=F, stringsAsFactors = F)
sunshine.tpm <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_norm.tsv',header=T, check.names=F, stringsAsFactors = F)
mend.fpkm <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix_2nd_read_strand.tsv',header=T, check.names=F, stringsAsFactors = F)
sunshine.fpkm <- read.table('/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F)
mend_sunshine.fpkm <- merge(mend.fpkm, sunshine.fpkm, by = "GeneName")
mend_sunshine.fpkm <- mend_sunshine.fpkm[ ,mend_sunshine_sample_names]
mend_sunshine.tpm <- merge(mend.tpm, sunshine.tpm)
mend_sunshine.tpm <- mend_sunshine.tpm[,mend_sunshine_sample_names1]
rownames(mend_sunshine.tpm) <- mend_sunshine.tpm$GeneID
mend_sunshine.tpm <- mend_sunshine.tpm[-446]
meta.cols <- sunshine.cts[, c('GeneID', 'GeneName', 'GeneBiotype')]  

## read clinical data
mend.clin <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2020_09_11_WashU_MendalianVsSporadics_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F) 
mend.clin['Cohort'] = "MendVSporadic"
sunshine.clin <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_10_13_WashU_SUNSHINE_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F) 
sunshine.clin['Cohort'] <- "Sunshine"
sunshine.clin_KnightADRC <- subset(sunshine.clin, study == "Knight_ADRC")
mend.clin_KnightADRC <- subset(mend.clin, study == "Knight_ADRC" )
mend_sunshine.clin <- bind_rows(mend.clin_KnightADRC, sunshine.clin_KnightADRC)
mend_sunshine.clin$Status <- factor(mend_sunshine.clin$Status)
mend_sunshine.clin$Cohort <- factor(mend_sunshine.clin$Cohort)
mend_sunshine.clin$Sex <- factor(mend_sunshine.clin$Sex)
mend_sunshine.clin$BraakTau <- factor(mend_sunshine.clin$BraakTau)
mend_sunshine.clin <- mend_sunshine.clin[mend_sunshine.clin$Subj_ID != "MAP_60007", ]
mend_sunshine.clin <- mend_sunshine.clin[mend_sunshine.clin$Subj_ID != "MAP_60180", ]
mend_sunshine.clin_PMI <- merge(mend_sunshine_sample_name_subj_id, mend_sunshine.clin, by = "Subj_ID")
mend_sunshine.clin_sample_name <- mend_sunshine.clin_PMI[!duplicated(mend_sunshine.clin_PMI$Sample_ID), ]
mend_sunshine.clin_PMI <- mend_sunshine.clin_sample_name[complete.cases(mend_sunshine.clin_sample_name$PMI),]
mend_sunshine.clin_PMI1 <- mend_sunshine.clin_sample_name[!complete.cases(mend_sunshine.clin_sample_name$PMI),]
PMI_na_values <- mend_sunshine.clin_PMI1$Sample_ID
mend_sunshine.clin_Braak <- merge(mend_sunshine_sample_name_subj_id, mend_sunshine.clin, by = "Subj_ID")
mend_sunshine.clin_sample_name <- mend_sunshine.clin_Braak[!duplicated(mend_sunshine.clin_Braak$Sample_ID), ]
mend_sunshine.clin_Braak <- mend_sunshine.clin_sample_name[complete.cases(mend_sunshine.clin_sample_name$BraakTau),]
mend_sunshine.clin_Braak_PMI <- mend_sunshine.clin_Braak[complete.cases(mend_sunshine.clin_Braak$PMI),]
mend_sunshine.clin_Braak1 <- mend_sunshine.clin_sample_name[!complete.cases(mend_sunshine.clin_sample_name$BraakTau),]
mend_sunshine.clin_Braak_PMI1 <- mend_sunshine.clin_Braak[!complete.cases(mend_sunshine.clin_Braak$PMI),]
Braak_na_values <- mend_sunshine.clin_Braak1$Sample_ID
Braak_PMI_na_values <- mend_sunshine.clin_Braak_PMI1$Sample_ID
mend_sunshine.clin_CDR <- merge(mend_sunshine_sample_name_subj_id, mend_sunshine.clin, by = "Subj_ID")
mend_sunshine.clin_sample_name <- mend_sunshine.clin_CDR[!duplicated(mend_sunshine.clin_CDR$Sample_ID), ]
mend_sunshine.clin_CDR <- mend_sunshine.clin_sample_name[complete.cases(mend_sunshine.clin_sample_name$CDRe),]
mend_sunshine.clin_CDR1 <- mend_sunshine.clin_sample_name[!complete.cases(mend_sunshine.clin_sample_name$CDRe),]
CDR_na_values <- mend_sunshine.clin_CDR1$Sample_ID

mend_cts_norm <- read.table("/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_norm.tsv", header = T, row.names = 1)
mend_cts_norm <- mend_cts_norm[,-c(1,2)]
sunshine_cts_norm <- read.table("/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_norm.tsv", header = T, row.names = 1)
sunshine_cts_norm <- sunshine_cts_norm[,-c(1,2)]
mend.tech_KnightADRC_sample_names1 <- str_replace_all(mend.tech_KnightADRC_sample_names, "-", ".")
mend_cts_norm_KnightADRC <- mend_cts_norm[ ,mend.tech_KnightADRC_sample_names1]
sunshine.tech_KnightADRC_sample_names1 <- str_replace_all(sunshine.tech_KnightADRC_sample_names, "-", ".")
sunshine_cts_norm_KnightADRC <- sunshine_cts_norm[ ,sunshine.tech_KnightADRC_sample_names1]
mend_sunshine.norm <- merge(sunshine_cts_norm_KnightADRC, mend_cts_norm_KnightADRC, by = 0)
row.names(mend_sunshine.norm) <- mend_sunshine.norm$Row.names
mend_sunshine.norm <- mend_sunshine.norm[-1]

all((mend_sunshine.clin2$Sample_ID) %in% colnames(mend_sunshine.cts))
all((mend_sunshine.clin2$Sample_ID) == colnames(mend_sunshine.cts))
cts_sub1 <- mend_sunshine.tpm[,colnames(mend_sunshine.tpm)%in%(mend_sunshine.clin2$Sample_ID)]
all(rownames(mend_sunshine.clin) %in% colnames(cts_sub))
all(rownames(mend_sunshine.clin) == colnames(cts_sub))

all(rownames(mend_sunshine.clin) %in% colnames(cts_norm))
all(rownames(mend_sunshine.clin) == colnames(cts_norm))
cts_norm_sub1 <- mend_sunshine.tpm[,colnames(mend_sunshine.cts)%in%mend_sunshine_sample_names]
all((mend_sunshine.clin2$Sample_ID) %in% colnames(cts_norm_sub1))
all((mend_sunshine.clin2$Sample_ID) == colnames(cts_norm_sub1))
mend_sunshine.clin1 <- mend_sunshine.clin[complete.cases(mend_sunshine.clin$PMI),]

filter.shreshold1 <- floor(length(colnames(cts_norm_sub1))*0.5)  ## data.cols is a vector of all samples in the matrix
keep1 <- rowSums(cts_norm_sub1 >= 0.1) >= filter.shreshold
cts_norm_sub1 <- cts_norm_sub1[keep1,]
cts_sub1 <- cts_sub1[rownames(cts_sub1) %in% rownames(cts_norm_sub1),]
cts_sub1 <- cts_sub1[,!(names(cts_sub1) %in% c("H_VY-60007_S1511546_I.D.21","H_VY-60180_S1511776_VI.N"))]
cts_sub1 <- cts_sub1[rownames(cts_sub1) %in% rownames(cts_norm_sub1),]

cts_sub_PMI1 = cts_sub1[,!(names(cts_sub1) %in% PMI_na_values)] #dropping the samples that are being excluded for having NA values 
cts_sub_Braak1 = cts_sub1[,!(names(cts_sub1) %in% Braak_na_values)]  
cts_sub_Braak_PMI1 = cts_sub_Braak1[,!(names(cts_sub_Braak1) %in% Braak_PMI_na_values)] 
cts_sub_CDR1 = cts_sub1[,!(names(cts_sub1) %in% CDR_na_values)] 
mend_sunshine.clin_PMI1 <- subset(mend_sunshine.clin_PMI1, select = -2)

dds <- DESeqDataSetFromMatrix(countData = cts_sub_CDR, 
                              colData = mend_sunshine.clin_CDR,     
                              design= ~ Sex + AOD + Cohort + Status)  
dds <- DESeq(dds)
sex_age_Cohort_CDR <- results(dds)#log2 fold change for Neuro_AD vs Neuro_CO 
sex_age_Cohort_CDR1 <- subset(sex_age_Cohort_CDR, pvalue<0.05) #filtering for padj<0.1
sex_age_Cohort_Status <- results(dds, contrast = c("Status", "Neuro_AD", "Neuro_CO"))#log2 fold change for Neuro_AD vs Neuro_CO 
sex_age_Cohort_Status1 <- subset(sex_age_Cohort_Status, pvalue<0.05) #filtering for padj<0.1
genesInterest <- c("S100B", "PLD3", "ABCA7", "MMP9", "TMEM106B", "CHI3L1", "TARDBP", "PSEN1", 
              "SNCA", "PSEN2", "APP", "PARK2", "MAPT", "NPTX2", "LRRK2", "GBA", "C9orf72", "SORL1", "NRGN", "GFAP", "GRN")
## add annotation
gg <- as.data.frame(str_split_fixed(rownames(sex_age_Cohort_Status), "_", 1))
colnames(gg) <- 'GeneID'
sex_age_Cohort_Status1 <- cbind(sex_age_Cohort_Status, gg)
sex_age_Cohort_Status1 <- as.data.frame(merge(meta.cols, sex_age_Cohort_Status1, by = "GeneID"))
write.csv(as.data.frame(sex_age_Cohort_Status1), 
                    file = "sex_age_Cohort_Status.csv", quote = FALSE)

##change these for whichever dataset you're working with at the time 
df_sex_age_PMI_Braak <- read.csv('/home/sohn/sex_age_PMI_Braak.csv', header =T, row.names =1, sep=",", stringsAsFactors = F, check.names = F)
df_sex_age_PMI_Braak <- merge(df_sex_age_PMI_Braak, meta.cols1, by=0)
rownames(df_sex_age_PMI_Braak) <- df_sex_age_PMI_Braak$Row.names
df_padj_sex_age_PMI_Braak <- subset(df_sex_age_PMI_Braak, padj<0.05)
df_pval_sex_age_PMI_Braak <- subset(df_sex_age_PMI_Braak, pvalue<0.05)
pval_genesInterest_sex_age_PMI_Braak <- intersect(df_pval_sex_age_PMI_Braak$GeneName, genesInterest)
padj_genesInterest_sex_age_PMI_Braak <- intersect(df_padj_sex_age_PMI_Braak$GeneName, genesInterest)
