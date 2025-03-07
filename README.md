# Summary
Scripts used in the analysis of TP53 and TP53 retrogene (RTG) CRISPR-Cas9 knockout data, as well identifying and analyzing woolly mammoth noncoding deletions. Includes scripts to filter sgRNAs for off-targets, clean and analyze amplicon nanopore data, and preform Asian elephant bulk RNA-seq analysis.

# File Summary
## Cas9Filtering
### Cas9_MMFiltering_V1_1.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Filters Guidescan2 output (https://github.com/pritykinlab/guidescan-cli) following mismatch estimation. Searches the ouput files for guides with the same name and ranks their mismatches. If 1 or more off-targets are predicted to be targetable, then that guide and all off-targets are removed from the output file. Targettability is deemed viable when the sum of the mismatch scores along the sgRNA for two or more sites is less than 1. Targetting scope is editable by changing the MM_Scores array to change the mismatch score weights in the 5' to 3' direction. 

## Nanopore
### CigarDel.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Examines the presence of absence of particular deletions within end-to-end Nanopore amplicon data, aligned to an amplicon reference in SAM format. Requires users to predefine a Region File containg the length of the amplicon, the proportion of missing data required for a deletion to be marked as present, and the deletion start and end coordinates (see ExampleRegionFile.txt and RegionFile_Template.txt). Multiple deletion regions can be assessed in the same analysis. Output is a plain text document containing each read name, and the deletion intervals that were detected at the request proportion parameter. 

### FastqFilter.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Filters raw sequencing read data to minimum and maximum sizes, as a first step in nanopore analysis. Reads along an input fastq file and outputs all reads that are between the specified minimum and maximum lengths provided as command line arguments. 

## RNAseq
### DEG_Analysis_O2_Pairwise.TP53.R
<ins>Language</ins>: R

<ins>Summary</ins>: Analyzes bulk RNA-seq data in a pairwise fastion (i.e. one sample condition against a control). Requires Salmon has been previously run on raw sequencing data and quant-sf files are produced. Also requires a tab-delimited metadata file (Metadata.txt) to be present in the same directory. All control samples need to be denoted with a "WT" in the second column (case-sensitive) (see ExampleMetadata.txt and Metadata_Template.txt). Designed for use with the Asian elephant reference genome (EleMax1; GCF_024166365.1) annotations combined with the TP53 retrogene loci added. A file assinging each transcript ID to genes is provided for EleMax1 with TP53 RTGs added. Paths to the Salmon output directories and the transcript-gene data file are hardcoded and may need to be changed. 

### DEG_Analysis_O2_Pairwise.LogFoldChange.R
<ins>Language</ins>: R

<ins>Summary</ins>: Identical to the TP53 RNAseq analysis script above, except labels the top 10 genes with the highest absolute log2 fold change instead of those related to the TP53 pathway. Used during mammoth deletion analysis. 

## TFBSDeletions
### GenerateWindows.awk
<ins>Language</ins>: Awk

<ins>Summary</ins>: Takes a fasta and generates windows along it of a given size, with a given slide (e.g. 25 bp windows starting every 10 bp). It also filters out windows if they contain more than the specified amount of missing data - via AmpProp=X where X is a decimal of how many N's are allowed in the given window. This runs better if you add "LANG=C" before the command (e.g. LANG=C ~/bin/GenerateWindows.awk -v WinSize=25 WinSlide=10 AmbProp=0.2 [input fasta]). 

### SplitDepth.awk
<ins>Language</ins>: Awk

<ins>Summary</ins>: Splits a genome-wide depth file into individual contig depth files. 

### CoverageCalc.awk
<ins>Language</ins>: Awk

<ins>Summary</ins>: Takes a tab delimited file containing windows (akin to a bed file - [Contig Name] [Start] [End]), and a samtools depth file for that contig, and outputs ranges which are some standard deviation threseholds below the mean coverage for the contig. If the value of StDevMultiplier is set really high (e.g. 99), it will restrict windows to only those which have no coverage. 

### ConcatenateWindows.awk
<ins>Language</ins>: Awk

<ins>Summary</ins>: Takes a tab delimited file containing windows and concatonates overlap or contigous ones into a single interval for easier downstream processing.

### TFBSTools_JASPAR2014.R
<ins>Language</ins>: R

<ins>Summary</ins>: Takes in a CSV containing some information on the contig and the sequence, and predicts transcription factor binding sites within the region. 

# Source Information
Scripts were generated as part of Karpinski et al. (2025). If used as part of your work, please cite the script used and the source publication: [DOI]
