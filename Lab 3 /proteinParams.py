#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 13:40:33 2019

@author: roshanahmed
"""

#!/usr/bin/env python3
# Name: Roshan Ahmed (roahmed)
# Group Members: None

"""
Program takes a protein sequence and using a dictionary of molecular weights, 
as well as a function for pI, molar/molecular extinction, caluclates properties 
of given sequence. 

Input: A protein sequence 
Output: number of amino acids, total molecular weight, 
        Molar/Mass extinction coefficient,
        theoretical isoelectric point (pI), amino acid composition. 

"""
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
            
# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
            
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()