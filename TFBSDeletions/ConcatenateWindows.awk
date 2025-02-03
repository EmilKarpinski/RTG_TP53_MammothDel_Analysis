#!/usr/bin/awk -f

#This program concates contigous windows into larger regions that can more easily be converted to a bed file for further filtering.  
#Takes as input:
    #A tab delimited text file of windows formatted as [contig] [window start] [window end] [anything else]
#Usage: ConcatonateWindows.awk [Input Windows]
#Written by: Emil Karpinski (2022-06-15)

BEGIN {
    #Sets the default print seperator in awk to be a tab
    OFS = "\t"; 

    #Intitalizing some variables
    Name = "";
    Lower = 0;
    Upper = 0;

}

#Catching headers from the first record
FNR == 1{
    print($1, $2, $3);
}

#Main body of program. Printing and doing math. 
FNR>1{
    #Checking if the lower bound of the window is within the previous window (for sliding windows) or equal to it (for discrete windows) 
    if($2 >= Lower && $2 <= Upper){
        Name = $1; 
        Upper = $3;
    }
    #If not then print the new previous larger window and store the new window variables. 
    else{
        #Not printing the first line. 
        if(Lower != 0){
            print(Name, Lower, Upper);
        }
        Name = $1;
        Lower = $2;
        Upper = $3;
    }
}
#Printing the final window
END { 
    print(Name, Lower, Upper);
}

