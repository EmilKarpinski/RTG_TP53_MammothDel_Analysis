# Summary
Scripts used in the analysis of TP53 and TP53 retrogene (RTG) CRISPR-Cas9 knockout data. Includes scripts to filter sgRNAs for off-targets, clean and analyze amplicon nanopore data, and preform Asian elephant bulk RNA-seq analysis.

# File Summary
## Cas9Filtering
### Cas9_MMFiltering_V1_1.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Filters Guidescan2 output (https://github.com/pritykinlab/guidescan-cli) following mismatch estimation. Searches the ouput files for guides with the same name and ranks their mismatches. If 1 or more off-targets are predicted to be targetable, then that guide and all off-targets are removed from the output file. Targettability is deemed viable when the sum of the mismatch scores along the sgRNA for two or more sites is less than 1. Targetting scope is editable by changing the MM_Scores array to change the mismatch score weights in the 5' to 3' direction. 




# Source Information
