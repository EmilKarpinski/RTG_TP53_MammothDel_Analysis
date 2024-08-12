#!/usr/bin/env python3

#Written by: Emil Karpinski (2024-03-01)
#Curent Version: V1.1
#Changelog - 1.1: Changed the output summary document to be based on the name of the input file instead of the general "Read_Deletion_Summary" as I think this was tripping up parralelized slurm jobs


#Program that will identify which deletions are present in a read using the cigar string and a tab delimited region file
#Usage: CigarDel.py <Input SAM File> <Region File> <Output File>

#Importing some packages
#sys to catch command line arguments
#re to do the cigar string splitting
#math to do one line of rounding during the thresehold stuff

import sys
import re
import math

#Cathing command line input
#infile = the sam to parse
#regionfile - the tab delimited region file. Also contains the info on how long the ref is and the % of bases that need to be deleted to call the deletion on the first line.
#outfile - the new final sam to ouptut

infile = sys.argv[1]
regfile = sys.argv[2]
outfile = sys.argv[3]


#Declaring and initalizing some variables
count = 0
RefLength = 0
Regions = []
MissProp = 0

#Making the name of the output summary file.
SummaryFile = outfile.rsplit(".")[0] + ".Read_Deletion_Summary.txt"


#Quick hack to delete the output so it doesn't keep appending to the end.
with open(outfile, "w") as Output:
    Output.close()
with open(SummaryFile,"w") as RDSummary:
    RDSummary.close()


#This little code block reads the RegionFile and parses it to get the relevant data.
with open(regfile, "r+") as RegionFile:
    for line in RegionFile:
        if (count == 0):
            info = line.split()
            RefLength = int(info[0])
            MissProp = float(info[1])
            count = 1
        else:
            l = line.split()
            Regions.append(l)

RegionFile.close()

#This code block does the processing of the sam file
with open(infile, "r+") as SamFile:
    #Everytime through we wipe the old data (and initialize the variables the first time)
    for line in SamFile:
        ReadName =""
        Pattern = ""
        Start = 0
        CIGAR = ""
        SummaryString = ""
        #Since Sam files have the header, need to check if the line starts with an @ indicating the sam header so we just write it to the output sam.
        if (line.startswith("@")):
            with open(outfile, "a+") as Output:
                Output.write(line)
        #For all non header lines we grab some info, then reconstruct the reference as a full length CIGAR string, then parse each region in order along that string to see if there's a deletion that meets the specified criteria.
        else:
            CurrRead = line.split('\t')
            ReadName = CurrRead[0]
            SummaryString += ReadName
            Start = int(CurrRead[3])
            #Note: this outputs whitespace in CIGAR[0] since the cigar string always starts with a digit. Functionall this means CIGAR is 1-indexed, not 0-indexed.
            CIGAR = re.split('(\d+)', CurrRead[5])

            #If the read does not start on position 1 of the reference (i.e. it starts mapping in the middle of the reference) we add "-"'s until we get to the start site.
            #Basically this just makes indexing later on easier.
            if (Start > 1):
                for i in range(0,Start):
                    Pattern += "-"
            
            #Reconstructs the CIGAR string for that molecule.
            #Note this ignores all bases that are not in the reference. So no insertions and nothing that would be masked
            for i in range(2,len(CIGAR), 2):
                if (CIGAR[i] != "S" and CIGAR[i] != "I" and CIGAR[i] != "P" and CIGAR[i] != "H"):
                    for c in range(0,int(CIGAR[i-1])):
                        Pattern += CIGAR[i]

            #Like at the start if the read doesn't go all the way to the end of the molecule this fills it in with dashes "-"'s
            if (len(Pattern)<RefLength):
                Add=RefLength-len(Pattern)
                for i in range(0,Add):
                    Pattern += "-"
            
            #Loops through all the regions and counts the number of deletions ("D" in the CIGAR string) within the deletion. Then checks to see if it passes the thresehold and if yes then that region is deleted in that read.
            #The thresehold is set by the length of the region and the float multiplier that's the second column in line 1 of the region file.
            for reg in Regions:
                RegStart = int(reg[0])+1
                RegEnd = int(reg[1])+1
                DelCounter = 0
                Thresehold = math.ceil((RegEnd - RegStart + 1)*MissProp)

                for i in range(RegStart, RegEnd+1):
                    if (Pattern[i]=="D"):
                        DelCounter+=1

                #Checks to see if the region passes the deletion thresehold. If it does then it appends the deletion to the name of that read in the same file, and creates a line for the summary table.
                if (DelCounter >= Thresehold):
                    ReadName = ReadName + "_" + reg[0] + "-" + reg[1]
                    SummaryString = SummaryString +'\t' + reg[0] + "-" + reg[1]
            
            #Modifies the original name of the read and prints the output same and summary table.
            CurrRead[0] = ReadName
            with open(outfile, "a+") as Output:
                Output.write('\t'.join(CurrRead[0:]))
            with open(SummaryFile,"a+") as RDSummary:
                RDSummary.write(SummaryString+'\n')




