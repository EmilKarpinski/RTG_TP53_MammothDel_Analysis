#!/usr/bin/awk -f

#This program splits a large samtools depth file into smaller depth files by scaffold/chromosome  
#Takes as input:
    #A samtools depth file
#Usage: SplitDepth.awk [Input Samtools Depth File]
#Written by: Emil Karpinski (2022-06-17)

BEGIN {
    #Sets the default print seperator in awk to be a tab
    OFS = "\t"; 

    #Intitalizing some variables
    Name = "";
    Output = "";
    Line = "";

}

#Catching first line and making the output file name variable
FNR == 1{
    Name = $1;
    Line = $0;
    Output = Name "" "_Depth.txt";
    FirstLineFlag = 0;
}

#Main body of program.  
FNR>1{
    #Checking if next line is from the same contig as the previous. 
    if($1 == Name){
        if(FirstLineFlag == 0){
            print(Line)>Output;
            FirstLineFlag = 1;
        }
        else{
            print(Line)>>Output;
        }
        Line = $0; 
    }
    else if($1 != Name){
        #Printing the final line of the output.
        print(Line)>>Output;

        #Reseting the variables for the new contig. 
        Name = $1;
        Line = $0;
        Output = Name "" "_Depth.txt";
        FirstLineFlag = 0;
    }
}
#Printing the final window
END { 
    print(Line)>>Output;
}

