#!/usr/bin/env python3

#Written by: Emil Karpinski (2023-03-22)
#Curent Version: V1.1

#Program for use with GuideScan2 output to filter out guides using a scoring system as opposed to MM based alone. 
#Program runs through the output file searching for guides with the same name and ranks them and their mismatches. If 1 or more off-targets are found that are predicted to be biologically viable based on the scoring system below it removes those guides from the file.
#Scoring is based on an adhoc consensus score based on the results from Anderson, E.M., et al (2015). doi: 10.1016/j.jbiotec.2015.06.427; Lee, H.J., et al. (2021). doi: 10.3390/ijms22126457; Bravo, J.P.K., et al. (2022). doi: 10.1038/s41586-022-04470-1; Ricci, C.G., et al. (2019). doi: 10.1021/acscentsci.9b00020 ; Fu, B.X.H., et al. (2016). doi: 10.1093/nar/gkw417; Fu, Y., et al. (2013). doi: 10.1038/nbt.2623 
#V1 - one scoring scheme and argument parsing just based on position. 
#V1.1 - Changed the line where the guidename is stored in memory to inside the Guidescore <1 loop as oppossed to before it. Important when dealing with non-target genomes as you can end up with a name, but no lines stored and this changes the off-target filtering.
#Usage Cas9_MMFiltering.py <Inputfile - from Guidescan2> <Outputfile - new reduced guidescan2 file>

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

#Defining the scoring list. Note this is adhoc based on some papers (above) and goes in the 5' to 3' direction (i.e. MM_Scores[1] is the furtherest position from the PAM; MM_Scores[20] is the PAM adjacent base). 
#Since zero-indexed starting with an empty zero here for easier looping.
MM_Scores=[0,0.33,0.33,0.33,0.4,0.4,0.5,0.5,0.5,0.5,0.5,1,1,1,1,1,1,1,1,1,1]

#For printing the first line
FirstLineFlag=0

#UniqueFlag for tracking if that guide is unique or not. And storing the variable for unique guide line storage
UniqueFlag=0
StoredLines=""


#Function that does the mismatch math and returns the cumulative mismatch score
def GuideScore_Calc(guide):
    #Resetting some variables
    count=0
    score=0

    #Loops through each character (tracking the position with count) and checks if its lower. If it is, the adds the current positions MM penalty to the scores variable.
    for char in guide:
        count+=1
        if char.islower() and count <= 20:
            #print(char, count)
            score+=MM_Scores[count]
    return(score)


#Opens the input file and loops through each line
with open(infile, "r+") as InputFile:
    for line in InputFile:      

        #Splits the current line of input on commas since the default input is a csv.
        CurrLine=line.split(',')

        #These if-statments go thrlsough the four possible cases that may arise.
        #Case 1: the first line in the file. Necessary to print the headers.
        if FirstLineFlag == 0:
            FirstLineFlag = 1
            GuideName = CurrLine[0]
            with open(outfile, "w+") as OutputFile:
                OutputFile.write(str(line))


        #Case 2: The name of the current guide is the same as the one stored in memory and the unique flag has not yet been tripped (i.e. there have been either no other match sequences for this guide or none that have passed the scoring filter)
        #If true then we calculate the score for the new mismatch sequence and if its an off-target (score<1) trip the UniqueFlag and remove any guides stored for printing.
        elif CurrLine[0] == GuideName and UniqueFlag == 0:
            Curr_Score = GuideScore_Calc(CurrLine[6])

            if Curr_Score < 1:
                UniqueFlag = 1
                StoredLines=""

        #Case 3: A shortcut for Case2 where if one other locus has already failed and tripped the UniqueFlag there's no reason to redo everything. 
        elif CurrLine[0] == GuideName and UniqueFlag == 1:
            pass

        #Case 4: The case where the current guide being examined isn't the same as the previous one. So need to print the previous stuff and store the new stuff.
        elif CurrLine[0] != GuideName:

            #Print previous lines to file. Checking just that this isn't blank here so that we can open Outputfile as above with w+ (write mode as opposed to append more here), to overwrite the output.
            if StoredLines != "":
                with open(outfile, "a+") as OutputFile:
                    OutputFile.write(str(StoredLines))

            #Reseting the line storage 
            StoredLines=""
            UniqueFlag = 0

            #Getting the score of the current match. I think this should always be zero (i.e. Guidescan always outputs the perfect match first), but leaving this in in case there is no perfect match in the genome and guidescan finds an offtarget (though maybe that should be used as input since mismatches will now be off)    
            Curr_Score = GuideScore_Calc(CurrLine[6])
            if Curr_Score < 1:
                #Storing the first instance of the guide and it's name
                StoredLines=line
                GuideName = CurrLine[0]


#Printing anything stored in StoredLines after done reading from the file
with open(outfile, "a+") as OutputFile:
    OutputFile.write(str(StoredLines))

#Closing Output and Inputfiles
OutputFile.close()
InputFile.close()