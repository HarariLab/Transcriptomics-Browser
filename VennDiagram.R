vennDiagramFunction <- function(model_name, direction){
    dataFileMend <- paste0("/home/sohn/", model_name, "_MEND1.csv") #csv file with deseq results, must have GeneName column 
    dataFileSunshine <- paste0("/home/sohn/", model_name, "_SUNSHINE1.csv")
    dataFileBoth <- paste0("/home/sohn/", model_name, "1.csv")
    mendVal <- read.csv(dataFileMend, header =T, sep=",", stringsAsFactors = F, check.names = F)
    mendVal <- subset(mendVal, padj < 0.05)
    mendValUp <- subset(mendVal, log2FoldChange > 0)
    mendValDown <- subset(mendVal, log2FoldChange < 0)
    sunshineVal <- read.csv(dataFileSunshine, header =T, sep=",", stringsAsFactors = F, check.names = F)
    sunshineVal <- subset(sunshineVal, padj < 0.05)
    sunshineValUp <- subset(sunshineVal, log2FoldChange > 0)
    sunshineValDown <- subset(sunshineVal, log2FoldChange < 0)
    bothVal <- read.csv(dataFileBoth, header =T, sep=",", stringsAsFactors = F, check.names = F)
    bothVal <- subset(bothVal, padj < 0.05)
    bothValUp <- subset(bothVal, log2FoldChange > 0)
    bothValDown <- subset(bothVal, log2FoldChange < 0)
    if (direction == "Down"){
      vennDiagramFile <- paste0(model_name, "_Down.png") 
      vennDiagramTitle <- paste0(model_name, " Downregulated") 
      vennDiagramVal <- venn.diagram(
        x = list(mendValDown$GeneName, sunshineValDown$GeneName, bothValDown$GeneName),
        category.names = c(paste("MEND", nrow(mendValDown)) , paste("SUNSHINE", nrow(sunshineValDown)) , paste("MEND+SUNSHINE", nrow(bothValDown))),
        filename = vennDiagramFile,
        main = vennDiagramTitle,
        cat.just = list(c(0.6,1), c(0.6,1), c(0,0)),
        #output=TRUE
      )
    }
    if (direction == "Up"){
      vennDiagramFile <- paste0(model_name, "_Up.png") 
      vennDiagramTitle <- paste0(model_name, " Upregulated") 
      vennDiagramVal <- venn.diagram(
        x = list(mendValUp$GeneName, sunshineValUp$GeneName, bothValUp$GeneName),
        category.names = c(paste("MEND", nrow(mendValUp)) , paste("SUNSHINE", nrow(sunshineValUp)) , paste("MEND+SUNSHINE", nrow(bothValUp))),
        filename = vennDiagramFile,
        main = vennDiagramTitle,
        cat.just = list(c(0.6,1), c(0.6,1), c(0,0)),
        #output=TRUE
      )
    }
}

vennDiagramFunction("Sex+Age+Status", "Down")
vennDiagramFunction("Sex+Age+Status", "Up")
vennDiagramFunction("Sex+Age+PMI+Status", "Down")
vennDiagramFunction("Sex+Age+PMI+Status", "Up")
vennDiagramFunction("Sex+Age+Neuro+Astro+Status", "Down")
vennDiagramFunction("Sex+Age+Neuro+Astro+Status", "Up")

vennDiagramFunction("Sex+Age+CDR", "Down")
vennDiagramFunction("Sex+Age+CDR", "Up")
vennDiagramFunction("Sex+Age+PMI+CDR", "Down")
vennDiagramFunction("Sex+Age+PMI+CDR", "Up")
vennDiagramFunction("Sex+Age+Neuro+Astro+CDR", "Down")
vennDiagramFunction("Sex+Age+Neuro+Astro+CDR", "Up")

vennDiagramFunction("Sex+Age+Braak", "Down")
vennDiagramFunction("Sex+Age+Braak", "Up")
vennDiagramFunction("Sex+Age+PMI+Braak", "Down")
vennDiagramFunction("Sex+Age+PMI+Braak", "Up")
vennDiagramFunction("Sex+Age+Neuro+Astro+Braak", "Down")
vennDiagramFunction("Sex+Age+Neuro+Astro+Braak", "Up")
