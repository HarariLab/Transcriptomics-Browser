library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gggenes)
library(dplyr)
library(shinyalert)
#library(soilprofile2)
library(shinydashboard)

load('/home/sohn/projects/transcriptomics_browser/MEND_AD_vs_CO/de_analysis.RData')
modelName <- c("Sex+Age+Status", "Sex+Age+PMI+Status", "Sex+Age+PMI+Cohort+Status","Sex+Age+Cohort+Status",
               "Sex+Age+Braak", "Sex+Age+PMI+Braak", "Sex+Age+PMI+Cohort+Braak", "Sex+Age+Cohort+Braak",
               "Sex+Age+CDR", "Sex+Age+PMI+CDR", "Sex+Age+PMI+Cohort+CDR", "Sex+Age+Cohort+CDR")
modelStatus <- c("Sex+Age+Status", "Sex+Age+PMI+Status", "Sex+Age+PMI+Cohort+Status","Sex+Age+Cohort+Status")
modelBraak <- c("Sex+Age+Braak", "Sex+Age+PMI+Braak", "Sex+Age+PMI+Cohort+Braak", "Sex+Age+Cohort+Braak")
modelCDR <- c("Sex+Age+CDR", "Sex+Age+PMI+CDR", "Sex+Age+PMI+Cohort+CDR", "Sex+Age+Cohort+CDR")

getTable <- function(d, l, p, t) {
  if (t %in% modelName){
    path <- paste0("/home/sohn/DE_results/all_results/", t, ".csv")
    de_sig_res1 <- read.csv(path, header=T, stringsAsFactors = F, check.names = F, sep=",")
    de_sig_res1 <- subset(de_sig_res1, padj<0.7)
    rownames(de_sig_res1) <- NULL
  }
  de_sig_res1$baseMean <- signif(de_sig_res1$baseMean, digits = 5)
  de_sig_res1$log2FoldChange <- signif(de_sig_res1$log2FoldChange, digits = 4)
  de_sig_res1$lfcSE <- signif(de_sig_res1$lfcSE, digits = 4)
  de_sig_res1$stat <- signif(de_sig_res1$stat, digits = 4)
  de_sig_res1$pvalue <- signif(de_sig_res1$pvalue, digits = 4)
  de_sig_res1$padj <- signif(de_sig_res1$padj, digits = 4)
  names(de_sig_res1)[names(de_sig_res1) == 'log2FoldChange'] <- "FC"
  names(de_sig_res1)[names(de_sig_res1) == 'padj'] <- "padjvalue"
  de_sig_res1$pvalue <- scientific(de_sig_res1$pvalue, digits = 3)
  de_sig_res1$padjvalue <- scientific(de_sig_res1$padjvalue, digits = 3)
  de_sig_res1 <- de_sig_res1[-c(1,4,6,7)]
  if (d == "All"){
    s1 <- de_sig_res1
  }
  if (d != "All"){
    s1 <- subset(de_sig_res1, de_sig_res1$GeneBiotype == d)
  }
  if(l != 0){
    s1 <- subset(s1, abs(s1$FC) > l)
  }
  if(p != 0){
    s1 <- subset(s1, as.numeric(s1$padjvalue) < p)
  }
  return(s1)
  #return (de_sig_res)
}

### read the DE results and prepare data for the ShinyApp
plotGeneFunction <- function(g, t) {
  if (!is.na(g)){
    if (t == "Sex+Age+CDR"){
      dds <- Sex_Age_CDR.dds
      g.groups <- as.data.frame(colData(dds)[10])##look at colData to make sure they are getting the correct items.
      dds_res <- Sex_Age_CDR.res
    }
    if (t == "Sex+Age+PMI+CDR"){
      dds <- Sex_Age_PMI_CDR.dds
      g.groups <- as.data.frame(colData(dds)[11])
      dds_res <- Sex_Age_PMI_CDR.res
    }
    if (t == "Sex+Age+PMI+Cohort+CDR"){
      dds <- Sex_Age_PMI_Cohort_CDR.dds
      g.groups <- as.data.frame(colData(dds)[11])
      dds_res <- Sex_Age_PMI_Cohort_CDR.res
    }
    if (t == "Sex+Age+Cohort+CDR"){
      dds <- Sex_Age_Cohort_CDR.dds
      g.groups <- as.data.frame(colData(dds)[11])
      dds_res <- Sex_Age_Cohort_CDR.res
    }
    if (t == "Sex+Age+Status"){
      dds <- Sex_Age_Status.dds
      g.groups <- as.data.frame(colData(dds)[3])
      dds_res <- Sex_Age_Status.res
    }
    if (t == "Sex+Age+PMI+Status"){
      dds <- Sex_Age_PMI_Status.dds
      g.groups <- as.data.frame(colData(dds)[4])
      dds_res <- Sex_Age_PMI_Status.res
    }
    if (t == "Sex+Age+PMI+Cohort+Status"){
      dds <- Sex_Age_PMI_Cohort_Status.dds
      g.groups <- as.data.frame(colData(dds)[4])
      dds_res <- Sex_Age_PMI_Cohort_Status.res
    }
    if (t == "Sex+Age+Cohort+Status"){
      dds <- Sex_Age_Cohort_Status.dds
      g.groups <- as.data.frame(colData(dds)[4])
      dds_res <- Sex_Age_Cohort_Status.res
    }
    if (t == "Sex+Age+Braak"){
      dds <- Sex_Age_Braak.dds
      g.groups <- as.data.frame(colData(dds)[12])
      dds_res <- Sex_Age_Braak.res
    }
    if (t == "Sex+Age+PMI+Braak"){
      dds <- Sex_Age_PMI_Braak.dds
      g.groups <- as.data.frame(colData(dds)[12])
      dds_res <- Sex_Age_PMI_Braak.res
    }
    if (t == "Sex+Age+PMI+Cohort+Braak"){
      dds <- Sex_Age_PMI_Cohort_Braak.dds
      g.groups <- as.data.frame(colData(dds)[12])
      dds_res <- Sex_Age_PMI_Cohort_Braak.res
    }
    if (t == "Sex+Age+Cohort+Braak"){
      dds <- Sex_Age_Cohort_Braak.dds
      g.groups <- as.data.frame(colData(dds)[12])
      dds_res <- Sex_Age_Cohort_Braak.res
    }
    #dds <- paste0("`", t, ".dds`")
    gene.counts <- counts(dds, normalized=FALSE)
    gene.norm.exp <- log2(counts(dds, normalized=FALSE)+1)
    #plotCounts(dds, gene=which(de_sig_res$GeneName==g), intgroup="Status")
    gene.id <- dds_res$GeneID[dds_res$GeneName == g]
    gene.cts <- melt(gene.counts[rownames(gene.counts) == gene.id,])
    gene.exp <- melt(gene.norm.exp[rownames(gene.norm.exp) == gene.id,])
    gene.cts$sample <- rownames(gene.cts)
    gene.exp$sample <- rownames(gene.exp)
    rownames(gene.cts) <- NULL
    rownames(gene.exp) <- NULL
    colnames(gene.cts) <- c('raw_counts', 'sample')
    colnames(gene.exp) <- c('norm_counts', 'sample')
    gene.data <- merge(gene.cts, gene.exp)
    ## extract groups 
    if (t %in% modelStatus){
      g.groups$sample <- rownames(g.groups)
      rownames(g.groups) <- NULL
      ## merge all 
      gene.data <- merge(gene.data, g.groups, sort =F)
      gene.data <- gene.data[gene.data$Status %in% c('Neuro_CO', 'Neuro_AD'), ]
      gene.data$Status <- factor(gene.data$Status, levels = c('Neuro_CO', 'Neuro_AD', as.character(unique(gene.data$Status)[!unique(gene.data$Status) %in% c('Neuro_CO', 'Neuro_AD')])))
      ### plot the expression 
      p <- ggplot(gene.data, aes(x=Status, y=norm_counts)) + geom_boxplot(aes(fill=Status)) + theme_bw() + geom_jitter(size=2) +
        theme(axis.text.x = element_text(size=12, face="bold"), axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12),
              plot.title = element_text(size=14, hjust =0.5, face ="bold"),
              legend.position = "none") +
        labs (x="", y="Log2(Expression+1)", title = paste0(g)) +
        #scale_fill_viridis_d() +
        scale_x_discrete(labels=c(paste0('Controls\n(n=', nrow(gene.data[gene.data$Status=="Neuro_CO",]),')'), 
                                  paste0('AD\n(n=', nrow(gene.data[gene.data$Status=="Neuro_AD",]),')')))
      print(p)
    }
    if (t %in% modelBraak){
      g.groups$sample <- rownames(g.groups)
      rownames(g.groups) <- NULL
      ## merge all 
      gene.data <- merge(gene.data, g.groups, sort =F)
      gene.data <- gene.data[gene.data$BraakTau %in% c(0,1,1.5,2,3,4,5,6 ), ]
      gene.data$BraakTau <- factor(gene.data$BraakTau, levels = c(0,1,1.5,2,3,4,5,6, as.character(unique(gene.data$BraakTau)[!unique(gene.data$BraakTau) %in% c(0,1,1.5,2,3,4,5,6)])))
      ### plot the expression 
      p <- ggplot(gene.data, aes(x=BraakTau, y=norm_counts)) + geom_boxplot(aes(fill=BraakTau)) + theme_bw() + geom_jitter(size=2) +
        theme(axis.text.x = element_text(size=12, face="bold"), axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12),
              plot.title = element_text(size=14, hjust =0.5, face ="bold"),
              legend.position = "none") +
        labs (x="", y="Log2(Expression+1)", title = paste0(g)) +
        #scale_fill_viridis_d() +
         scale_x_discrete(labels=c(paste0('Braak 0\n(n=', nrow(gene.data[gene.data$BraakTau==0,]),')'), 
                                   paste0('Braak 1\n(n=', nrow(gene.data[gene.data$BraakTau==1,]),')'),
                                   paste0('Braak 1.5\n(n=', nrow(gene.data[gene.data$BraakTau==1.5,]),')'),
                                   paste0('Braak 2\n(n=', nrow(gene.data[gene.data$BraakTau==2,]),')'),
                                   paste0('Braak 3\n(n=', nrow(gene.data[gene.data$BraakTau==3,]),')'),
                                   paste0('Braak 4\n(n=', nrow(gene.data[gene.data$BraakTau==4,]),')'),
                                   paste0('Braak 5\n(n=', nrow(gene.data[gene.data$BraakTau==5,]),')'),
                                   paste0('Braak 6\n(n=', nrow(gene.data[gene.data$BraakTau==6,]),')')))
      print(p)
    }
    if (t %in% modelCDR){
      g.groups$sample <- rownames(g.groups)
      rownames(g.groups) <- NULL
      ## merge all 
      gene.data <- merge(gene.data, g.groups, sort =F)
      gene.data <- gene.data[gene.data$CDRe %in% c(0,0.5,1,2,3), ]
      gene.data$CDRe <- factor(gene.data$CDRe, levels = c(0,0.5,1,2,3, as.character(unique(gene.data$CDRe)[!unique(gene.data$CDRe) %in% c(0,0.5,1,2,3)])))
      ### plot the expression 
      p <- ggplot(gene.data, aes(x=CDRe, y=norm_counts)) + geom_boxplot(aes(fill=CDRe)) + theme_bw() + geom_jitter(size=2) +
        theme(axis.text.x = element_text(size=12, face="bold"), axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12),
              plot.title = element_text(size=14, hjust =0.5, face ="bold"),
              legend.position = "none") +
        labs (x="", y="Log2(Expression+1)", title = paste0(g)) +
        #scale_fill_viridis_d() +
        scale_x_discrete(labels=c(paste0('CDR 0\n(n=', nrow(gene.data[gene.data$CDRe==0,]),')'), 
                                  paste0('CDR 0.5\n(n=', nrow(gene.data[gene.data$CDRe==0.5,]),')'),
                                  paste0('CDR 1\n(n=', nrow(gene.data[gene.data$CDRe==1,]),')'),
                                  paste0('CDR 2\n(n=', nrow(gene.data[gene.data$CDRe==2,]),')'),
                                  paste0('CDR 3\n(n=', nrow(gene.data[gene.data$CDRe==3,]),')')))
      print(p)
    }
  }
}

plotVol <- function(d, l, p, g, t) {
  x=0
  if (t %in% modelName){
    path <- paste0("/home/sohn/DE_results/all_results/", t, ".csv")
    de_sig_res <- read.csv(path, header=T, stringsAsFactors = F, check.names = F, sep=",")
    rownames(de_sig_res) <- NULL
  }
  de_sig_res$direction <- 'NC'
  de_sig_res[de_sig_res$log2FoldChange > 0 & !is.na(de_sig_res$padj) & de_sig_res$padj < 0.05, 'Direction'] <- 'Up'
  de_sig_res[de_sig_res$log2FoldChange < 0 & !is.na(de_sig_res$padj) & de_sig_res$padj < 0.05, 'Direction'] <- 'Down'
  s1 <- de_sig_res
  #print(s1)
  if (l != 0){
    #s1 <- de_sig_res$GeneName[abs(de_sig_res$log2FoldChange) > as.numeric(l)]
    s1 <- subset(s1, abs(s1$log2FoldChange) > l)
    x=1
  }
  if (p != 0){
    #s1 <- de_sig_res$GeneName[abs(de_sig_res$padj) < as.numeric(p)]
    s1 <- subset(s1, as.numeric(s1$padj) < p)
    x=1
  }
  ### plot the data
  v <- ggplot(de_sig_res, aes(x=log2FoldChange, y=-log10(pvalue))) + theme_bw()
  if (length(g) != 0){
    v <- v + geom_text_repel(fontface="bold", size=4, aes(label = ifelse(GeneName==g, GeneName, '')))
  }
  if (x==1) {
    v <- v + geom_point(aes(color = ifelse(GeneName %in% s1$GeneName, 'blue', 'gray90')))
    v <- v + scale_color_manual(name="", values=c("blue"="blue","gray90"="gray90"), labels =c("blue"=paste0("LogFC > abs(",l,")"),"gray90"="Others"))
  } else {
    v <- v + geom_point(aes(color=Direction))
    v <- v + scale_color_manual(name="", values=c("Up"="orange2","Down"="skyblue3"), labels =c("up"="Up","dn"="Down"))
  }
  #v <- v + geom_text(aes(label=ifelse(de_sig_res$pvalue < 0.055, de_sig_res$Gene, '')), color="black", size=2, vjust=-1, fontface="plain")
  v = v + theme(axis.text.x=element_text(size=12, color="black"),
                axis.text.y=element_text(size=12, color="black"),
                axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
                #plot.title = element_text(size = 12, hjust=0.5, color="black", face="bold"),
                legend.position="bottom", legend.key.size = unit(5,'mm'),
                legend.text=element_text(size=12),
                panel.border = element_rect(linetype='solid', color='black'))
  if (l != 0){
    v <- v + geom_vline(xintercept = l, linetype="dashed")
    v <- v + geom_vline(xintercept = -l, linetype="dashed")
  }
  # if (p != 0){
  #   v <- v + geom_hline(yintercept = -log10(p), linetype="dashed")
  # }
  v <- v + guides (color = guide_legend(override.aes = list(size = 8)))
  v
}

shinyServer(function(input, output, session){
  modelData <- reactive({
    return(input$modelDE)
  })
  FCText <- reactive({
    return(input$FCThresholdText)
  })
  FCData <- reactive({
    return(input$FCThresholdSlider)
  })
  pvalueData <- reactive({
    return(input$pvalSlider)
  })
  
  # observeEvent(input$FCThresholdText,{
  #   if(is.na(as.numeric(as.character(input$FCThresholdText)))) {
  #     shinyalert("Please enter a valid number", type = "error")
  #     print("Getting here or wat")
  #     stop()
  #   }
  # })

  observeEvent(input$pvalText, {
    # if(!is.numeric(input$pvalText)) {
    #   shinyalert("Please enter a valid padjusted value between 0 and 1", type = "error")
    #   #stop()
    # }
    if ((as.numeric(input$pvalText) != input$pvalSlider) &
        input$pvalText != "" &  input$pvalSlider != "")
    {
      updateSliderInput(
        session = session,
        inputId = 'pvalSlider',
        value = input$pvalText
      )
    } else {
      if (input$pvalText == "") {
        updateSliderInput(session = session,
                          inputId = 'pvalSlider',
                          value = 0)
      }
    }
  })
  observeEvent(input$pvalSlider, {
    if ((as.numeric(input$pvalText) != input$pvalSlider) &
        input$pvalSlider != "" & input$pvalText != "")
    {
      updateTextInput(
        session = session,
        inputId = 'pvalText',
        value = input$pvalSlider
      )
    }
  })
  
  observeEvent(input$FCThresholdText, {
    if ((as.numeric(input$FCThresholdText) != input$FCThresholdSlider) &
        input$FCThresholdText != "" &  input$FCThresholdSlider != "")
    {
      updateSliderInput(
        session = session,
        inputId = 'FCThresholdSlider',
        value = input$FCThresholdText
      )
    } else {
      if (input$FCThresholdText == "") {
        updateSliderInput(session = session,
                          inputId = 'FCThresholdSlider',
                          value = 0)
      }
    }
  })
  observeEvent(input$FCThresholdSlider, {
    if ((as.numeric(input$FCThresholdText) != input$FCThresholdSlider) &
        input$FCThresholdSlider != "" & input$FCThresholdText != "")
    {
      updateTextInput(
        session = session,
        inputId = 'FCThresholdText',
        value = input$FCThresholdSlider
      )
    }
  })

  output$mytable = DT::renderDataTable(getTable(input$selectedGeneBiotype, FCData(), pvalueData(), modelData()), server = FALSE, selection=list(mode = "single"))
  # highlight selected rows in the scatterplot
  output$gene_plot = renderPlot({
    g = input$mytable_rows_selected
    g = getTable(input$selectedGeneBiotype, FCData(), pvalueData(), modelData())[g, "GeneName"]
    if (length(g) && nzchar(g)) plotGeneFunction(g, modelData()) 
  })
  
  output$vol_plot <- renderPlot({ 
    g = input$mytable_rows_selected
    g = getTable(input$selectedGeneBiotype, FCData(), pvalueData(), modelData())[g, "GeneName"]
    plotVol(input$selectedGeneBiotype, FCData(), pvalueData(), g, modelData()) 
  })
  
})  
