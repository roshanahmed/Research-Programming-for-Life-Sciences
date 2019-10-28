#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:33:12 2019

@author: roshanahmed

# Name: Roshan Ahmed (roahmed)
# Group Members: Irvin Ramos (iraramos) Kourosh Kouhmareh (kkouhmar)
"""

class OrfFinder:
    
    """
    Find open reading frames from a sequence. 
    Use for loops to parse through the reading frames.
    Use start flags to check if a start codon is found. 
    Returns list of Orfs. 
    Below is code for dangled ends, was not producing right output.
    """

#                if codon in self.stopCodons:                                    #check to see if we have a stop codon
#                    stopFlag = 1
                    
#                    
#                if startFlag == 0 and stopFlag == 1:                           #if we have a stop codon with no start codon
#                    start = 1
#                    stop = position + 3 
#                    length = (stop - start) + 1
#                    self.saveOrfs(((frame % 3) + 1), start, stop, length)
#                    
                    
#            if startFlag == 1 and stopFlag == 0:    #if no stop codon was found 
#                
#                start = self.startList[0] + 1
#                stop = len(self.sequence)
#                length = (stop - start) + 1 
#                self.saveOrfs(((frame % 3) + 1), start, stop, length)

#           elif codon in self.stopCodons:                                    #check to see if we have a stop codon
#                stopFlag = 1
        
#                elif startFlag == 0 and stopFlag == 1:                         #if we have a stop codon with no start codon
#                    start = (len(reverseSequence)) - (pos + 2)
#                    stop = (len(reverseSequence))
#                    length = (stop - start) 
#                    self.saveOrfs(-1 * ((frame % 3) + 1), start, stop, length)
#                    
#                    
#            if startFlag == 1 and stopFlag == 0:                               #if no stop codon was found 
#                
#                start = self.startList[0] + 1
#                stop = len(self.sequence)
#                length = (stop - start) + 1 
#                self.saveOrfs([frame, start, stop, length])
#    

    
    def __init__(self, sequence):
        
        """
        Initializes list of orfs and list of acceptables codons 
        """
        
        self.sequence = sequence
        self.orfList = []
        self.startCodons = ["ATG"]                                             #list of start codons
        self.stopCodons = ["TAA", "TAG", "TGA"]                                #list of stop codons
        
    def findOrfs(self):
        """
        Function takes in sequence, creates list of start codons, and returns values of:
        frame, start, stop, and length to the save orfs function. 
        """
        startList = []                                                         #list of start codons
        for frame in range(0,3):                                               #for all three reading frames
            startFlag = 0                                                      #set start Flag to 0 
            for position in range(frame, len(self.sequence), 3):               #for each reading frame parse sequence
                codon = self.sequence[position: position + 3]                  #slice codon to check 
                
                if codon in self.startCodons:                                  #check to see if we have a start codon
                    startList.append(position)                                 #store position of start
                    startFlag = 1                                              #increment flag
                
                elif startFlag == 1 and codon in self.stopCodons:              #if we have both a start and stop codon 
                    if frame == 1:                                             #if we are in the first frame
                        start = startList[0] + 1                               #set start to first codon in list 
                        stop = position + 3                                    #stop is position in at end of codon read

                    else: 
                        start = startList[0] + 1                               #same as above
                        stop = position + 3
                
                    length = (stop - start) + 1                                #calculate length given stop and start
                    if length > 100:                                           #checks length of output before returning to Orfs
                        self.saveOrfs(((frame % 3) + 1), start, stop, length)
                        
                    startList.clear()                                          #clear list for next reading frame
                    startFlag = 0                                              #re-set flag to 0
                         
        return(self.orfList)
        
    
    def reverseOrfs(self):
        """
        Same as the previous function but for the lagging strand.
        For the calculations have to use the reverse end, because the location
        on the original strand as at the opposite end. 
        """
        
        reverseSequence = self.reverseComplement()                             #sequence from the reverse complement 
        startList = []                                                         #create empty start list

        for frame in range(0,3):
            startFlag = 0 
            for position in range(frame, len(reverseSequence), 3):
                codon = reverseSequence[position: position + 3]
                
                if codon in self.startCodons:                                   #check to see if we have a start codon
                    startList.append(position)                                  #store position of start
                    startFlag = 1            
            
                elif startFlag == 1 and codon in self.stopCodons:               #if we have both a start and stop codon 
                    start = (len(reverseSequence)) - (position + 2)
                    stop = (len(reverseSequence)) - startList[0]
                    length = (stop - start) + 1
                    if length > 100:
                        self.saveOrfs(-1 * ((frame % 3) + 1), start, stop, length)
                        
                    startList.clear()
                    startFlag = 0
                    
        return(self.orfList)
        
        
    def reverseComplement(self):
        """
        Returns the reverse complement of the leading strand. 
        Using a dictionary of bases and their complements. 
        """
        complementDict = {"A": "T", "T": "A", "G": "C", "C": "G"}                      #dictionary of bases and thair complements
        reverseSeq = (list(reversed(self.sequence)))                                   #reverses sequence and casts as a list
        reverseComp = "".join(complementDict.get(base, base) for base in reverseSeq)   #eliminates whitespace and joins string 
        return reverseComp                                            
    
    
    def saveOrfs(self, frame, start, stop, length):
        """
        function that saves Orfs 
        """
        self.orfList.append([frame, start, stop, length])                              #append list of orfs in that order

                    
class NucParams:
    """
    Create Dictionaries to store information from Fasta Files.
    Dictionary 1: Amino Acid Composition 
    Dictionary 2: Nucleotide Composition 
    Dictionary 3: Codon Composition

    Includes list of valid nucleotides to check input 
    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    validNucleotides = ["A", "G", "C", "T", "U", "N"]     #list of valid Nucleotides


    def __init__ (self, inString=''):
        """
        Intializes dictionaries with proper keys and sets values to zero. 
        """
        
        self.aaComp = {}                                        
        self.nucComp = {}
        self.codonComp = {}
        
        for aa in self.rnaCodonTable.values():            #Intialize dictionary 
            self.aaComp[aa] = 0 
            
        for nuc in self.validNucleotides:                 #only take ACTGUN
            self.nucComp[nuc] = 0 
                 
        for codon in self.rnaCodonTable.keys():           #Intialize dictionary 
            self.codonComp[codon] = 0 
                
    
    def addSequence (self, inSeq):
        """
        Puts sequences in dictionaries made in __init__ method 
        """
        
        rnaSeq = inSeq.replace("T", "U")                  #replace T with U for RNA seq
        
        for nuc in rnaSeq:
            if nuc in NucParams.validNucleotides:         #if valid nucleotide 
                self.nucComp[nuc] += 1                    #count how many nucleotides
                
        for group in range(0, len(rnaSeq), 3):            #use range and 3 step to look at codons
            codon = rnaSeq[group: group + 3]              
            if codon in self.rnaCodonTable.keys():
                self.codonComp[codon] += 1                #adds sequence to dictionary and count
                aminoAcid = self.rnaCodonTable[codon]
                self.aaComp[aminoAcid] += 1               #adds amino acid to dicionary and count 
                
            
    def aaComposition(self):
        """
        Returns amino acid comp dictionary and counts. 
        """
        
        return self.aaComp
    
    
    def nucComposition(self):
        """
        Returns nucleotide dictionary from addSequence.
        """
        
        return self.nucComp
    
    
    def codonComposition(self):
        """
        Returns codon composition dictionary and counts. 
        """
        
        return self.codonComp
    
    
    def nucCount(self):
        """
        Returns sum of all the valid nucleotides. 
        """
        
        return sum(self.nucComp.values())
    
      
import sys
class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
        

class ProteinParam :
    """
    Program calculates the properties of peptide and outputs composition, pI and other properties. 
    """
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    
    def __init__ (self, protein):
        """
        Takes list of string and creates string object 
        """
        
        bigAA = protein.upper()
        self.aminoDict = {}
        
        for aa in ProteinParam.aa2mw.keys(): 
            self.aminoDict[aa] = 0 
            for value in bigAA:
                if value is aa:
                    self.aminoDict[aa] += 1 
                
    
    def aaCount (self):
        return sum(self.aminoDict.values())
    
    
    def pI (self):
        """
        isoelectric point using mypH values, mypH is incremented till charge is as close to 0 
        """
        
        bestCharge = 10000000                                                 #large number used to keep if running
        bestpH = 0 
        mypH = 0 
        
        while mypH < 14.01 :                                
            aaCharge = abs((self._charge_(mypH)))                             #checks charge at given pH
            if aaCharge < bestCharge:                                         #continuing loop 
                bestCharge = aaCharge      
                bestpH = mypH
            mypH += 0.01                                                      #adds .01 to pH value and continues 
        return bestpH                                                         #return pH at which charge is 0 
        
            
    def aaComposition (self) :
        """
        returns dictionary created in _init_ method 
        
        """
        return (self.aminoDict)

    
    def _charge_ (self, pH):
        """
        Calculates net charge of peptide chain at any pH. Uses pKa values of charged amino acids.
        input: pH from pI method 
        output: net charge of peptide at given pH
        """
        
        positive = 0 
        for aa in self.aa2chargePos:     #for all the positive amino acids 
            positive += ((self.aminoDict[aa]) * (10 ** self.aa2chargePos[aa])/((10 ** self.aa2chargePos[aa]) + (10 ** pH))) 
            nTerm = (10 ** self.aaNterm)/((10 ** self.aaNterm) + (10 ** pH))
        positive += nTerm                #adds charge of n-terminus 
            
        negative = 0 
        for aa in self.aa2chargeNeg:     #for all the negative amino acids
            negative += ((self.aminoDict[aa]) * (10 ** pH)/((10 ** self.aa2chargeNeg[aa]) + (10 ** pH)))
            cTerm = (10 ** pH)/((10 ** self.aaCterm) + (10 ** pH))
        negative += cTerm                #adds charge of c-terminus 
            
        return (positive - negative)     #returns charge difference 
    

    def molarExtinction (self):
        """
        Molar Extinction calculation based on number of tyrosines, tryptophans, and cysteines in peptide.
        Uses extinction coeffcients at 280nm, returns peptide extinction coefficent. 
        """
        
        tyrosine = self.aminoDict.get("Y") * self.aa2abs280.get("Y")          #multiplies number of tyrosines by coefficents 
        tryptophan = self.aminoDict.get("W") * self.aa2abs280.get("W")        #multiplies number of tryptophans by coefficents 
        cysteine = self.aminoDict.get("C") * self.aa2abs280.get("C")          #multiplies number of cysteines by coefficents 
        return (tyrosine + tryptophan + cysteine)                             #sum of all extinction coefficents

    
    def massExtinction (self):
        """
        Calculates mass extinction, divides molar extinction coefficent by molecular weight .
        """
        
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

   
    def molecularWeight (self):
        """
        Calculates total molecular weight of peptide given sequence.
        Sums weight of all amino acids in peptide. 
        """
        
        totalWeight = 0 
        pBond = (self.mwH2O * (self.aaCount() - 1))                           #caclulates weight of h2o released during pbond formation
        
        for aa in self.aminoDict:                                             #for every amino acid in peptide sequence
            if self.aminoDict.get(aa) is not 0:
                weight = self.aa2mw.get(aa) * self.aminoDict.get(aa)          #multiply weight of aa by number of aa
                totalWeight += weight                                         #add to previous weights 
                
        return (totalWeight - pBond)                                          #subtract weight of h2o released 