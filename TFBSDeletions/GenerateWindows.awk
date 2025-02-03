#!/usr/bin/awk -f

#This program takes a fasta and generates windows of a specified size. 
#Usage: GenerateWindows.awk -v WinSize=XXX WinSlide=XXX AmbProp=XXX TestData_GenerateWindows.fasta
#Written by: Emil Karpinski (2022-06-09)

BEGIN {
    #Sets the default print seperator in awk to be a tab
    OFS = "\t"; 	
	#Initializes an empty string to store the Name name.
	Name = "";
	#Intializes an empty string to store the Sequence. 
	Seq = "";

    ##Uses the split function which splits the string "A,C,T,G,N" into the nucleotides array based on the "," seperator. 
    ###I'm not sure what the purpose of N_nuc here is but it contains the value 5? Maybe this corresponds to the number of nucleotides in the analysis such that this could be extended. 
    N_nuc = split("A,C,T,G,N",nucleotides,",");
    #Setting a flag to print the header the first time around. 
    Flag = 0;
}

#This function counts each nucleotide and returns a tab delimited
#string of counts (including 0s in the order A,C,T,G,N)
##Initializes the CountWindow function and catches the passed variable in a new one called window (typically contains the contents of the Seq variable I think)
function CountWindow(window){

    ##This for loop loops from 1 to 5, acting as an index for the nucleotides array and setting the count of each nucleotide to 0. Necessary in case one of the characers isn't present 
    ##This count it stored in the Count array which has a BP index (i.e. Count["A"]; Count["C"])
    for(nuc=1;nuc<=N_nuc;nuc++){
        Count[nucleotides[nuc]] = 0;
    }
    ##As above splits the string stored in the window variable into the array chars, based on nothing (so should be every characters).
    ##Also generates a variable n which stores the length of the chars array (i.e. the number of characters). 
    n=split(window,chars,"");
    ##This loop loops through every character in the new chars array, and for each one does two things:
    ##First it converts it to upper case (in case the fasta has lower case bases in it. 
    ##Second it updates the corresponding entry in the Count array if the base is a A, C, T, or G, else it updates the count of the N entry (a general catch all for N's and all other ambigious bases).
    for(c=1;c<=n;c++){
        char = toupper(chars[c]);
        if(char ~ /[ACTG]/){
            Count[char]++;
        } else {Count["N"]++}
    }
    ##Creates a new variable called str (which is used below in creating the output), and sets it equal to the number of A's in the window. 
    ##This appears to be just a shortcut. You could probably do this in the loop as well, but it means more code. 
    str = Count["A"];
    ##Loops from 2 to 5 (the value of N_nuc), and appends the value of the Count[] element consistent with the current base ot the value of str.
    ###I'm pretty sure you could make str an empty string above and loop through 5 times as well. 
    for(nuc=2;nuc<=N_nuc;nuc++){
        str=str"\t"Count[nucleotides[nuc]];
    }

    ##This returns the value of string to the line that called it. 
    return str;
}
##This is a really truncated if-statment which basically says if you encounter a ">" at the start of line "^" do the code in the block.
###Functionally I think this means that 
/^>/{ 
    if(Flag == 0){
        ##Generates a substring from the line with everything from the first character
        Name=substr($0,2);
        Output = Name "" ".txt";
        FilteredOutput = Name "" "_AmbFiltered.txt";
        
        #Print the headers
        print("Name","WinStart","WinEnd","A","C","T","G","N") > Output;
        print("Name","WinStart","WinEnd","A","C","T","G","N") > FilteredOutput;
        Flag = 1;
    }

    ##Generates a substring from the line with everything from the first character
    Name=substr($0,2);
    Output = Name "" ".txt";
    FilteredOutput = Name "" "_AmbFiltered.txt";

    ##Sets the window start position to 1q
    winpos=1;
    ##Immediately stops processing the current line and goes to the next one. 
    ###I suspect this is just a time saver, and is unnecessary. Just stops the program from evaulating the next block of code
    next;
}
#Each line of Sequence is appended to the current known Sequence
{Seq = Seq $0 }
{
    #Once at least a windows worth of Sequence is loaded, count the window
    #Then strip the window off the front of the loaded Sequence
    while(length(Seq) >= WinSize){
        ##Generates the substring of the window from the Seq starting at the first position and equal to the window size
        window=substr(Seq,1,WinSize);
        ##Calculates the final base in the window
        winend = winpos + WinSize - 1;

        ##Prints the name of the Name, the start and end position of the window and then runs the CountWindow function (above) passing it the window substring. This returns the string containing the number of A, T, C, and Gs seperated by tabs already. 
        print(Name, winpos, winend, CountWindow(window)) >> Output;
                #Also printing the window if it is less than or equal to the ambiguity filter. 
        if(Count["N"]/WinSize <= AmbProp){
            print(Name, winpos, winend, str) >> FilteredOutput;
        }
        ##Updates the winpos integer by adding the WinSize produces a new integer corresponding to the rightmost position outside the window. 
        winpos+=WinSlide;
        ##Produces a new substrating from the Sequence starting at the rightmost base outside the window. 
        ##This uses WinSize+1 as opposed to winpos since it effectively truncates the Sequence by removing the old window. So winpos wouldn't work. 
        Seq = substr(Seq,WinSlide+1);
    }

} 
END { #Finish processing any remaining Sequence
    ##As above checks if there's any Sequence length (i.e. length != 0) and 
    if(length(Seq)>WinSize){
        ##Calculates the final position in the Sequence, except here will create a window from the last tested position to the end of the Sequence (i.e. a smaller window comprising only the remainder of the Seq). 
        winend = winpos + length(Seq) - 1;

        ##Prints the name of the Name, the start and end position of the window and then runs the CountWindow function (above) passing it the final, partial, substring. This returns the string containing the number of A, T, C, and Gs seperated by tabs already. 
        print(Name, winpos, winend, CountWindow(Seq)) >> Output;
        
        #Also printing the window if it is less than or equal to the ambiguity filter. 
        if(Count["N"]/WinSize <= AmbProp){
            print(Name, winpos, winend, str) >> FilteredOutput;
        }
    }
}

