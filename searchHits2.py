#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 16:47:45 2021

@author: klaus
"""


import os
from Bio.Seq import Seq
import pandas as pd
import re
import regex
import os
import numpy as np
import translatePattern as tp
import gregexprindex as gr
from Bio.Seq import Seq
import py2bit
import writeHits2 as wrtHts2


pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)



def _preprocess_me2(gRNAs, max_mismatch):
    if len(gRNAs) <= 10:
        return(False)
    gRNA_min_width = min(gRNAs["width"])
    if gRNA_min_width < 6 or gRNA_min_width // (max_mismatch + 1) < 3:
        return(False)
    other = []
    for i in range(len(gRNAs)):
        if len(gRNAs["seq"][i]) == gRNAs["seq"][i].count("A") + gRNAs["seq"][i].count("C") + gRNAs["seq"][i].count("G") + gRNAs["seq"][i].count("T"):
            other.append(0)
        else:
            other.append(len(gRNAs["seq"][i]) - gRNAs["seq"][i].count("A") + gRNAs["seq"][i].count("C") + gRNAs["seq"][i].count("G") + gRNAs["seq"][i].count("T"))      
    if sum(other) != 0:
        return(False)
    return True


def _searchHitsInOneSeq2(gRNAs, seq, seqname, BSgenomeName, PAM, PAM_pattern, PAM_size, max_mismatch,
                         allowed_mismatch_PAM, outputDir, PAM_location = "3prime", 
                         baseEditing = False, targetBase = "C", editingWindow = [3, 4, 5, 6, 7, 8]):
    
    if _preprocess_me2(gRNAs, max_mismatch) == True:
        patterns = PDict(gRNAs, max_mismatch = max_mismatch)
    else:
        patterns = gRNAs

    all_plus_matches = []
  
    for i in range(len(patterns)):
        #seq12 = regex.findall("(" + patterns.loc[i, ("seq")] + "){s<=" + str(max_mismatch) + "}", seq, overlapped = True)
        seqmatch = []
        seq12 = regex.finditer("(" + patterns.loc[i, ("seq")] + "){s<=" + str(max_mismatch) + "}", seq, overlapped = True)

        for j in seq12:         
            
            seqmatch.append(j)
         
        all_plus_matches.append(pd.Series(seqmatch))
        
      
     
    revseq = str(Seq(seq).reverse_complement())
    all_minus_matches = []
    
    for i in range(len(patterns)):
       
        revseqmatch = []
       
        revseq12 = regex.finditer("(" + patterns.loc[i, ("seq")] + "){s<=" + str(max_mismatch) + "}", revseq, overlapped = True)

        for j in revseq12:
           
            revseqmatch.append(j)
     
            
        all_minus_matches.append(pd.Series(revseqmatch))

   
    for i in range(len(gRNAs)):
     
        patternID = re.sub("'", "", gRNAs.loc[i, "names"])
        if len(patternID) < 1:
            patternID = "pattern" + i
        pattern = Seq(gRNAs.loc[i, ("seq")])
        # by default PAM is NGG or NAG
        plus_matches = all_plus_matches[i]
      
       
        if len(plus_matches) > 0:
            
             names_plus_matches = [patternID] * len(plus_matches)
             wrtHts2.writeHits2(gRNA = pattern, seqname = seqname, matches = plus_matches, strand = "+",
                     file = outputDir + "outfile.txt", gRNA_size = len(pattern), PAM = PAM, PAM_pattern = PAM_pattern, 
                     max_mismatch = max_mismatch, chrom_len = len(seq), appendd = True, 
                     PAM_location = PAM_location, PAM_size = PAM_size,  mismatch_names = names_plus_matches,
                     allowed_mismatch_PAM = allowed_mismatch_PAM, BSgenomeName = BSgenomeName,
                     baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        if pattern.reverse_complement() != pattern:
            
            minus_matches = all_minus_matches[i]
            if len(minus_matches) > 0:
               
                names_minus_matches = [patternID] * len(minus_matches)
                wrtHts2.writeHits2(gRNA = pattern, seqname = seqname, matches = minus_matches, strand = "-",
                       file = outputDir + "outfile.txt", gRNA_size = len(pattern), PAM = PAM, PAM_pattern = PAM_pattern, 
                       max_mismatch = max_mismatch, chrom_len = len(seq), appendd = True, 
                       PAM_location = PAM_location, PAM_size = PAM_size, mismatch_names = names_minus_matches,
                       allowed_mismatch_PAM = allowed_mismatch_PAM, BSgenomeName = BSgenomeName,
                       baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow) 
    
   
    

def searchHits2(gRNAs, BSgenomeName, outputDir, chromToSearch = ["all"], chromToExclude = "",
                max_mismatch = 3, PAM_size = 3, gRNA_size = 20, PAM = "NGG", PAM_pattern = "NRG$",
                allowed_mismatch_PAM = 1, PAM_location = "3prime", 
                baseEditing = False, targetBase = "C", editingWindow = [4, 8]):
 
    if type(gRNAs) != pd.core.frame.DataFrame:
        raise Exception("gRNAs is required as a DNAStringSet object!")
    if BSgenomeName == None:
        raise Exception("BSgenomeName is required as BSgenome object!")
    
    max_mismatch = max_mismatch
    seqnames = list(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms().keys())

    if chromToSearch != ["all"]:
        seqnames = set(seqnames).intersection(chromToSearch)
    if len(chromToExclude) > 0:
        seqnames = set(seqnames).difference(chromToExclude)
    appendd = False
    if len(gRNAs) < 1:
        return gRNAs

    if gRNAs["width"][0] == gRNA_size + PAM_size:
        if PAM_location == "3prime":
            for i in range(len(gRNAs)):
                gRNAs.loc[i, ("seq")] = gRNAs.loc[i, ("seq")][0:gRNA_size]
                gRNAs.loc[i, ("width")] = len(gRNAs.loc[i, ("seq")])
        else:
            for i in range(len(gRNAs)):
                gRNAs.loc[i, ("seq")] = gRNAs.loc[i, ("seq")][PAM_size:gRNA_size + PAM_size]
                gRNAs.loc[i, ("width")] = len(gRNAs.loc[i, ("seq")])
    elif gRNAs["width"][0] != gRNA_size:
        raise Exception("the gRNA length needs to be equal to the specified gRNA.size (or gRNA.size plus PAM.size\n)")
    for seqname in seqnames:
        print(">>> Finding all hits in sequence", seqname, "...\n")
        subject = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(seqname)
        _searchHitsInOneSeq2(gRNAs = gRNAs, seq = subject, seqname = seqname, PAM = PAM,
                             PAM_pattern = PAM_pattern, PAM_size = PAM_size, max_mismatch = max_mismatch,
                             allowed_mismatch_PAM = allowed_mismatch_PAM, outputDir = outputDir,
                             PAM_location = PAM_location, BSgenomeName = BSgenomeName,
                             baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        print(">>> DONE searching\n")   
    
    if os.path.isfile(outputDir + "outfile.txt") == True:
        hits = pd.read_table(outputDir + "outfile.txt")
        os.remove(outputDir + "outfile.txt")
        return hits
    else:
        print("No matching found, please check your input sequence, and make sure you are using the right genome. You can also alter your search criteria such as increasing max.mismatch!")
        return pd.DataFrame()
    
    
