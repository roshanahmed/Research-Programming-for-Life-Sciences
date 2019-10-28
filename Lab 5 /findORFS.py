#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Name: Roshan Ahmed (roahmed)
# Group Members: Irvin Ramos (iraramos) Kourosh Kouhmareh (kkouhmar)

from sequenceAnalysis import OrfFinder, FastAreader

"""
Pseudo Code Outline 

__init__: 
    sequence = DNA from fastA file 
    readingFrame = list of open reading frames found
    
readingFrame:
    
    for frames 1,2,3 
        for :
            codon = splice sequence to get codon 
            if codon is startCodon:
                save codon to list 
                start = ORF list [0]
            if codon is stop codon and we found startCodon:
                start = ""
                stop = start + 3
            if codon is stop 
            
reverseframe:
    
    same as reading frame, just do math in reverse.
            
storeFrame:
    append reading frame list based on start codons found 
    
reverseFrame:
    take reverse complement 
    use dictionary and reversed function
    return
    
    
"""



########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(inCL=None):
    '''
    Find some genes.  
    '''
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    
###### replace the code between comments.
    print (myCommandLine.args)
        # myCommandLine.args.inFile has the input file name
        # myCommandLine.args.outFile has the output file name
        # myCommandLine.args.longestGene is True if only the longest Gene is desired
        # myCommandLine.args.start is a list of start codons
        # myCommandLine.args.minGene is the minimum Gene length to include
        #
#######
        
    myReader = FastAreader(myCommandLine.args.inFile)
    
    with open(myCommandLine.args.outFile, "w") as textFile:
        for head, sequence in myReader.readFasta():                            #uses fastA reader to get arguments for text file
            textFile.write(str(myCommandLine.args) + "\n")                     #convert command line arguments to string to print 
            textFile.write(head + "\n")                                        #write header and new line to text file 
            myOrf = OrfFinder(sequence)                                        #sequence from OrfFinder
            myOrf.findOrfs()                                                   #output of saved leading strand from OrfFinder
            myOrf.reverseOrfs()                                                #output of saved lagging strand from OrfFinder
            
            sortedGenes = sorted(myOrf.orfList, key = lambda orf: orf[3], reverse = True)  #sorted list of orfs 
            
            for frame, start, stop, length in sortedGenes:                     #format the frame, start, stop, and length in sorted 
                textFile.write('{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(frame, start, stop, length))  
        
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN
    
    
    