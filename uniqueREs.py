# -*- coding: utf-8 -*-
"""
Created on Thu May 14 00:47:50 2020

@author: Klaus
"""

import os
import numpy as np
import pandas as pd
import py2bit
from Bio.Seq import Seq
import biostrings as bio
import re
import pyreadr
import rdata
import isPatternUnique as isp


def character(length):
    a =[]
    for i in range(0,length):
        a.append("")
    return a




def uniqueREs(REcutDetails, summary, offTargets, BSgenomeName, scanUpstream = 100, scanDownstream= 100):
    
    a = 0
    if "REpattern" not in REcutDetails.columns:
        REcutDetails["REpattern"] = "NA"
        REcutDetails["REname"] = "NA"
        a += 1
    REwithName = pd.concat([pd.DataFrame(REcutDetails["REpattern"]), pd.DataFrame(REcutDetails["REname"])],axis=1).drop_duplicates()
    
    REs = character(np.shape(summary)[0])
    summary["id"] = summary["names"]
    summary["id"] = summary["names"].str.cat(summary["forViewInUCSC"], sep="-")
    offTargets["id"] = offTargets["name"]
    offTargets["id"] = offTargets["name"].str.cat(offTargets["forViewInUCSC"], sep="-")
    summary = pd.merge((summary),(offTargets[["id", "chrom", "chromStart", "chromEnd", "strand"]]))
    REs = character(np.shape(summary)[0])
    if np.shape(summary)[0]>0:
        Start = summary["chromStart"].astype(int)
        End = summary["chromEnd"].astype(int)
        strand = summary["strand"]
        chr = summary["chrom"]
        
        Start = Start - scanUpstream
        End = End + scanDownstream
        
        
       
        width = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms().values()
       
        chromosome = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms().keys()
        chromosome = pd.DataFrame(list(chromosome))
        width = pd.DataFrame(list(width))
        
       
        widthchromosome = pd.concat([width, chromosome], axis=1)
        widththisChr = []
        for i in chr:
            for j in range(len(chromosome)):
                if i == widthchromosome.iloc[j, 1]:
                    widththisChr.append(widthchromosome.iloc[j, 0])
        widththisChr = pd.Series(widththisChr)
      
        for i in range(len(Start)):
            a = 0
            thisChr = chr[i]
            if not pd.isnull(thisChr) and thisChr != "":
                thisEnd = min(End[i], widththisChr[i])
                thisStart = max(1, Start[i]) - 1
                thisStrand = str(strand[i])
                if thisStrand == "+":
                    scanSequence = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(thisChr, int(thisStart), int(thisEnd))
                else:
                    scanSequence = str(Seq(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(thisChr, int(thisStart), int(thisEnd)))[::-1].reverse_complement())[::-1]

        
                REnames = (str(summary["REname"][i])).split()
               
                REpatterns = []
                for j in range(len(REcutDetails)):
                    if str(REcutDetails["REname"][j]) in REnames:
                        if REcutDetails["REpattern"][j] not in REpatterns:
                            REpatterns.append(REcutDetails["REpattern"][j])
              
                REnames = []
                for k in range(len(REwithName)):
                    if REwithName.iloc[k, 0] in REpatterns:
                        REnames.append(REwithName.iloc[k, 1])
                
                
                for t in range(len(REpatterns)):
                    if isp.isPatternUnique(scanSequence, bio.DNAStringSet(REpatterns))[t] == "Yes":
                        if a == 0:
                            REs[i] = REnames[t]
                            a += 1
                        else:
                            REs[i] = REs[i] + " " + REnames[t]
            
                        

            else:
                REs[i] = ""
        if a != 1:
            REcutDetails.drop(["REname", "REpattern"], axis = 1, inplace = True)
  
    return REs
   
                