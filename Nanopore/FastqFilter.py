#!/usr/bin/env python3

#Written by: Emil Karpinski (2024-02-29)
#Curent Version: V1.0

#Program that will filter fastq by read length
#Usage python3 FastqFilter.py <Input Fastq> <Output Fastq> <Min Length> <Max Length>


#Cathing command line input
#Infile - the original fastq 
#Outfile - the fastq to output
#MinLength - the minimum fragment length
#Also importing gzip to allow reading and writing of compressed fastqs
import sys
import gzip

infile = sys.argv[1]
outfile = sys.argv[2]
MinLength = int(sys.argv[3])
MaxLength = int(sys.argv[4]) + 1
Flag = 0

#Small function that will just loop through for writing. 
#Need to add the encode line here to write it in bytes since it's required for gzip files.
def WriteLines():
    for l in lines:
        OutputFile.write(l.encode())

#Open the file and creates a line variable
#Need the while True here to keep reading after the first four lines. However this also necessitates a break 
with gzip.open(infile , 'rt') as InputFile:
    while True:
        lines = []

        #Loops through the first four lines and puts them into the lines array
        for i in range(4):
            try:
                lines.append(InputFile.readline())
            except StopIteration:
                break    


        #Checks if the length of the thrid line is greater than or equal to the minimum length if so prints it.
        #Proably should simplify this and just write an empty file at the start, but atm using this flag to just check if I need to open and write to the file or append to it. 
        #Note: these are greater than since the length function counts the \n as a character so it adds +1 to the length
        if len(lines[1]) > MinLength and len(lines[1]) <= MaxLength and Flag==0:
            Flag=1
            with gzip.open(outfile, "wb") as OutputFile:
                WriteLines()
		
        elif len(lines[1]) > MinLength and len(lines[1]) <= MaxLength and Flag==1:
            with gzip.open(outfile, "ab") as OutputFile:
                WriteLines()

        #Need this line here to eventually break out of the while loop
        elif lines[1]=="":
            break

#Closing the output and input
InputFile.close()
OutputFile.close()
