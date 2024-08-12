## Gene-level differential expression analysis using DESeq2 ##
## Based on the DE-expression analysis from HBC: https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/
#Written by Emil Karpinski (2023)

## Setup
#Installing BiocManager and Bioconducter + CRAN libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # BiocManager::install("DESeq2")
# # install.packages("tidyverse")
# # install.packages("RColorBrewer")
# # install.packages("pheatmap")
# # BiocManager::install("DEGreport")
# # BiocManager::install("tximport")
# # install.packages("ggplot2")
# # install.packages("ggrepel")

### Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)

#Getting stuff passed from the command line
CL_Args <- commandArgs(trailingOnly = TRUE)

#Printing to stdout to monitor progress/quality control
CL_Args

## List all directories containing data
#full.names gives the full path to each file, and pattern will only return those ending in the first Argument
samples <- list.files(path="../Salmon/", full.names = T, pattern = CL_Args[2])

#Getting the second file path and adding it to the vector
samples <- c(samples,list.files(path = "../Salmon/", full.names = T, pattern = CL_Args[3]))

#Printing to stdout to monitor progress/quality control
samples

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "../Salmon//","") %>% str_replace("_Salmon", "")

# Loading the annotation table I generated for Elemax1 from the rna.fna file
tx2gene <- read.delim("../../../DEGAnnotationTable/EleMax1_Transcript2GeneID.txt")

#Reading in the Salmon quant.sf files here, and assigning each Transcript ID to the Genesymbol directly.
#Tested this for ~3 loci and it worked identically to the version in the old script.
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("TranscriptID", "GeneSymbol")], countsFromAbundance="lengthScaledTPM")


######## Reading in a metadata file instead of creating one. Seems unnecessary and looks like it ouputs the same.
METADATA <- read.csv("Metadata.txt",header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE)
meta <- droplevels(METADATA[names(files),])

#DESeq2 DE analysis
## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Set)
## Run analysis
dds <- DESeq(dds)

#Looping through all contrasts relative to the Wild Type (which will be the H2O Control)
#There's probably a better way to do this, but this is easy.
#Will loop through all the levels in meta
for(i in unique(meta[,1])){
  if(i != "WT"){
    #Retaining the sample name for downstream titles
    SampleName <- i
    
    #Setting the comparison
    Comp_Contrast <- c("Set", i, "WT")
    
    ## Extract results for comparison above. 
    #This alpha is for gene-filtering if DESeq2. I'm just gonna leave this at 0.05, but there's an explanation here in the "Genes with a low mean normalized counts" section (https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html)
    Res_Table <- results(dds, contrast=Comp_Contrast, alpha = 0.05)
    
    #
    #Extracting the significant DE Genes
    #
    ### Set thresholds for the adjusted P-value after FDR correction
    padj.cutoff <- 0.05
    
    # Create a tibble of results
    Res_Table_tb <- Res_Table %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
    
    # Subset the tibble to keep only significant genes and sorting by fold-change
    SigDEGs <- Res_Table_tb %>% dplyr::filter(padj < padj.cutoff)
    SigDEGs <- SigDEGs %>% dplyr::arrange(log2FoldChange)

    #Generating the name of the output file
    OutputTableName <- paste(i,"_vs_WT_SigDEGs.txt",sep="")
    write.table(SigDEGs, file=OutputTableName, sep="\t", quote=F, row.names = FALSE)
    
    #
    ####VOLCANO PLOTS###
    #Adds a column indicating whether the gene passes the adjusted P-value threshold and is significant
    Res_Table_tb <- Res_Table_tb %>% dplyr::mutate(Threshold = padj < 0.05)
    
    #Setting the title and filename for the plots
    Plotname <- paste(i,"_vs_WT",sep="")
    PlotFilename <- paste(Plotname,".pdf",sep="")
    
    #Starting the plot printing.
    pdf(PlotFilename)
    
    ## Create an empty column to indicate which genes to label
    Res_Table_tb <- Res_Table_tb %>% dplyr::mutate(genelabels = "")
    
    #Making a vector with the TP53 downstream gene labels
    GeneList <- c("CDKN1A","SFN","RPRM","GADD45A","GTSE1","ZNF385A","FAS","PIDD1","BAX","SIVA1","TP53I3","SHISA5","AIFM2","IGFBP3","SERPINE1","DDB2","PTEN","STEAP3","MDM2","TP53")

    #Looping through and labelling the genes above.
    for(g in GeneList){
      Res_Table_tb$genelabels[match(g,Res_Table_tb$gene)] <- g
    }

    #Making the plot. 
    #Have to save to a variable then print since it's in a loop.
    VolcanoPlot <- ggplot(Res_Table_tb, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(colour = Threshold)) +
      geom_text_repel(aes(label = genelabels),max.overlaps = 1000000, force=10) +
      ggtitle(Plotname) +
      xlab("log2 fold change") + 
      ylab("-log10 adjusted p-value") +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25))) 
    print(VolcanoPlot)
    
    #Stopping the printing
    dev.off()
    
  }
}

##Plotting some general graphs to check quality and sense of data

#Making the PCA plot with the two main PCA vectors
#Starting the plot printing.
PCA_Title <- paste(SampleName,"_PCA.pdf",sep="")
pdf(PCA_Title)

### Transform counts for data visualization
##The blind=TRUE argument is to make sure that the rlog() function does not take our sample groups into account - i.e. does the transformation in an unbiased manner. 
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
#By default the ntop takes the top 500 most variable genes. But this can be increased by changing it as shown.
PCAPlot <- plotPCA(rld, intgroup="Set",ntop = 500)
print(PCAPlot)

#Stopping the printing
dev.off()


##Plotting Hierarchical Clustering

#Starting the plot printing.
SHC_Title <- paste(SampleName,"_HierCluster.pdf",sep="")
pdf(SHC_Title)

### Extract the rlog matrix from the object
rld_mat <- assay(rld)    

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

### Plot heatmap using the correlation matrix and the metadata object
HClusterPlot <- pheatmap(rld_cor, annotation = meta)

#Stopping the printing
dev.off()


##Plotting the data dispersion
#Starting the plot printing.
Disp_Title <- paste(SampleName,"_Dispersion.pdf",sep="")
pdf(Disp_Title)

##Look at the dispersion of the data
#Should be a nice decreasing curved line - with dispersion decreasing as expression increases. 
DisPlot <- plotDispEsts(dds)

#Stopping the printing
dev.off()


