#!/usr/bin/awk -f

#This program calculates the mean coverage and standard deviation for the contig and flags windows below it. 
#This takes two files as input:
    #A tab delimited file of windows formatted as [contig] [window start] [window end]
    #A samtools depth file for the contig outputting all positions. 
#Usage: CoverageCalc.awk -v StdevMultiplier=XXX [Input Windows] [Input Samtools depth]
#Written by: Emil Karpinski (2022-06-09)

BEGIN {
    #Sets the default print seperator in awk to be a tab
    OFS = "\t"; 

    #Setting the preliminary Flag
    Flag = 0;

    #Intitalizing some counters
    Count = 0;
    WindowCount = 0;

}
#FNR > 1 skips the header and NR == FNR loops I think until the file is finished. 
#Essentially this reads in the first file and creates an associative array such that every position in one of the windows is set to 0. The exact values here don't matter, so it doesn't matter if it gets overwritten.
FNR>1 && NR==FNR{
    #Storing the window start position in the Window Start array (note this one starts at 2 since FNR tracks the number of records we processed from this file and FNR == 1 is the headers).
    WindowStart[FNR] = $2;

    #Set values of all of the positions in our window to 0. 
    for(i = $2; i<=$3; i++){
        Depth[i]=0;
    }

    if(Flag == 0){
        Name = $1;
        Output = Name "" "_LowCovWindows.txt"
        Winsize = $3 - $2; 
        Flag = 1;

        #Print the header information
        print("Name","WinStart","WinEnd") > Output;
    }
    #Counting to keep track of the number of unique windows.
    WindowCount+=1;
    next;
}

#Sets the value of each position in the depth array.
#Also including the name check here in case I don't want to split the original bam, so I don't overwrite the depth once we switch scaffolds. 
FNR>1{
    if(Depth[$2] == 0 && $1 == Name){
        Depth[$2] = $3;
    }
    next;
}

END { 
    #Calculating some values for subsequent filtering. 
    #Need to also count the number of unique values here. 
    for(i in Depth){
        Sum+=Depth[i];
        SumSq+=Depth[i]*Depth[i];
        Count+=1;

    }
    #Calculating the mean and standard deviation. 
    OverallMean = Sum/Count; 
    StDev = sqrt(SumSq/Count - (Sum/Count)**2);
    TempThres = (OverallMean-(StDev*StdevMultiplier));

    print("Overall Mean: ", OverallMean);
    print("StDev: ", StDev);
    print("Coverage Thresehold: ", TempThres);

    if(TempThres< 0){
        CovThres = 0;
        print("Coverage Thresehold less than 0. Outputing all zero cov positions instead.");
    }
    else{
        CovThres = TempThres;
    }

    #Calculating the window mean coverage.
    for(i = 2; i<= WindowCount; i++){
        #Resetting the sum. 
        Sum=0;

        #Looping though all depths in the window and summing the coverage, then calculating the mean. 
        for(s = WindowStart[i]; s<=(WindowStart[i]+Winsize); s++){
            Sum+=Depth[s];
        }
        Mean = Sum/(Winsize+1);

        #If the mean is less than or equal to [StdevMultiplier] x StDev below the Overall Mean for all windows; flag it and print that window. 
        if(Mean <= CovThres){
            print(Name, WindowStart[i], WindowStart[i]+Winsize, Mean) >> Output;
        }
    }
}
