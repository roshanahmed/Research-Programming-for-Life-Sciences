#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 18:51:57 2019

@author: roshanahmed
"""

#!/usr/bin/env python3
# Name: Roshan Ahmed (roahmed)
# Group Members: Kourosh Kouhmareh (kkouhmar) Ana Paula Kitos Vasconcelos (akitosva)

from sequenceAnalysis import FastAreader

class findUnique:
    """
    Program takes a trna Sequence and first cleans it removing unwanted characters. 
    Then using a powerset function creates a powerset of each sequence. 
    Compares subsests to find unique sequences. 
    Uses a find essentials function to take unique sequences and remove larger substring. 
    Prints out ordered essentials by their starting position.
    """
    
    powerList = []  
    def __init__ (self, header, sequence):
        """ 
        Constructor function initializes all lists and sets. 
        Gets cleaned sequence from the cleanSeq fucntion. 
        """
        
        sequence = self.cleanSeq(sequence)
        findUnique.powerList.append(self)
        self.uniqueList = []                            
        self.trnaSet = set()
        self.uniqueSet = set()
        self.essentialSet = set()
        self.header = header 
        self.sequence = sequence
        self.trnaSet = self.powerSet(sequence)
        
    def cleanSeq(self, sequence):
        """
        Cleans sequence by replacing all the dashes and underscores.
        """
        
        newSeq = sequence.replace("-", "") 
        newSeq = sequence.replace ("_", "")
        return(newSeq)
        
    def powerSet(self, sequence):
        """
        Creates power set and returns a set of substrings for each tRNA sequence. 
        """
        
        length = len(self.sequence)
        for index in range(0, length):
            for char in range(index + 1, length + 1):
                self.trnaSet.add(self.sequence[index: char])
        return(self.trnaSet)
        
    def uniqueSeq(self):
        """
        Finds unique sequences by comparing the power sets and finding the difference 
        between values in the set. 
        """
        
        self.uniqueSet = self.trnaSet.copy()
        for trnaSeq in findUnique.powerList:
            if trnaSeq.trnaSet is not self.trnaSet:
                self.uniqueSet = self.uniqueSet.difference(trnaSeq.trnaSet)
        return self.uniqueSet
            
    def essentialSeq(self):
        """
        Takes uniques and removes superstrings. 
        Returns the minimal form of each substring which are the essentials. 
        """
        
        self.essentialSet = self.uniqueSet.copy()
        for sub1 in self.uniqueSet:
            for sub2 in self.uniqueSet:
                if sub1 in sub2 and sub1 is not sub2:
                    self.essentialSet.discard(sub2)
        return self.essentialSet
            
    
    def outputFormat(self):
        """
        Aligns the sequences so output is formatted correctly.
        """
        
        formatList = []
        for val in self.essentialSet:
            index = self.sequence.find(val)
            while index is not -1:
                formatList.append(''.join(['.'*index,val]))
                index = self.sequence.find(val, index + 1)
                
        formatList.sort(key = len)
        return formatList
        
########################################################################
# Main
# Here is the main program
# 
########################################################################

def main():
    """
    Prints out headers and formats list. 
    Calls other functions in class findUnique to run. 
    """
    
    trnaReader = FastAreader()
    trnaSeq= []
    
    for header, sequence in trnaReader.readFasta():trnaSeq = findUnique(header, sequence)
    findUnique.powerList.sort(key = lambda x:x.header)              #sorts powerList 
    for index in range (0, len(findUnique.powerList)):              
        headLine = findUnique.powerList[index].header.replace(" ","")
        seqLine = findUnique.powerList[index].sequence
        
        print(headLine)
        print(seqLine)
        
        uniques = findUnique.powerList[index].uniqueSeq()           #calls powerList function
        essentials = findUnique.powerList[index].essentialSeq()     #calls essential function
        aligned = findUnique.powerList[index].outputFormat()        #calls outputFormat function
        for sequenceAlign in aligned:print(sequenceAlign)           #prints formatted list 
        
        
if __name__ == "__main__":
    main()  
