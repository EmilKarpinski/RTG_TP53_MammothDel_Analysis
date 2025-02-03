#Importing stuff
library(Biostrings)
library(TFBSTools)
suppressMessages(library(JASPAR2014))

#Defining parameters for JASPAR retrieval
#Not sure about this all versions = True thing. 
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["all_versions"]] <- TRUE

#Retrieving the set from JASPAR2014. 
#The number here is 263 with the above opts parameters (205 with all versions set to False). This is way less than JASPAR2022 on the website (n=1205), but not sure if that's accessible.
#Appears it's missing many of the variations for these TFs (e.g. has MA0003.1, but not MA0003.2, MA0003.3, MA0003.4 from the web tool)
PFMatrixList <- getMatrixSet(JASPAR2014, opts)
PFMatrixList

#Converting the PFMatrixList to a PWMatrixList
PWMatrixList <- toPWM(PFMatrixList, pseudocounts=0.8)

#Setting the working directory to the folder on dropbox.
#setwd("C:/Users/Emil/Dropbox\ (HMS)/Emil_Workbooks/TranscriptionFactorPredicition")

#Catching the command line arguments
myargs = commandArgs(trailingOnly=TRUE)
DataFile = myargs[1]

#Reading in the data csv
Data <- read.csv(DataFile, header = FALSE)
#Deleting an empty row
#Data <- Data[-3,]

#Preforming the TFBS predicition on the current file. 
sitesetList <- searchSeq(PWMatrixList, DNAString(Data[1,4]), seqname=Data[1,1], min.score="60%", strand=Data[1,5])

#Getting P-values for the predicted locations. 
#Using the sampling type here since "TFMPValue" takes literally forever and crashes.
PVal <- pvalues(sitesetList, type="sampling")
#Converting the list into a single dataframe
PVal_List <- stack(PVal)

#Storing the results in a variable for easier filtering. 
GFF3 <- writeGFF3(sitesetList)

#Filtering for a p-value of less than or equal to 0.0001. 
FilteredGFF3 <- GFF3[which(PVal_List[,1] <= 0.0001),]

#Creating the output file and outputting the data
Output <- paste(Data[1,6],"_p0.0001.csv",sep="")
write.csv(FilteredGFF3,file = Output)










