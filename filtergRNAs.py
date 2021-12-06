#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 08:37:02 2021

@author: klaus
"""


import os
import re
import pandas as pd
import numpy as np
from findgRNAs import findgRNAs
from Bio.Seq import Seq
import translatePattern as tP
import gregexprindex as gr
import biostrings as bio

pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)


def filtergRNAs(all_gRNAs, pairOutputFile = None, findgRNAsWithREcutOnly = False, format = "fasta",
                REpatternFile = os.getcwd() + os.sep + os.sep.join(["extdata", "NEBenzymes.fa"]),
                minREpatternSize = 4, overlap_gRNA_positions = [17, 18], overlap_allpos = True):
    
    
    
    if len(all_gRNAs) == 0:
        raise Exception("all_gRNAs contains no gRNAs!")
    if REpatternFile == None:
        raise Exception("REpatternFile containing the restriction enzyme cut pattern is required!")
    if os.path.isfile(REpatternFile) == False:
        raise Exception("REpatternFile specified as ", REpatternFile, " does not exists!")
    if format != "fasta" and format != "fastq":
        raise Exception("format needs to be either fasta or fastq!")
    
    patterns = pd.read_csv(REpatternFile, header = None)
    seq = []
    names = []
    width = []
    for i in range(len(patterns)):
        if i % 2 == 0:
            names.append(patterns.iloc[i][0].replace(">", ""))
        else:
            seq.append(patterns.iloc[i][0])
            width.append(len(patterns.iloc[i][0]))
    
    patterns = pd.concat([pd.DataFrame(width, columns=(["width"])), pd.DataFrame(seq, columns=(["seq"])), pd.DataFrame(names, columns=(["names"]))], 1).drop_duplicates(subset=["seq"]).reset_index(drop=True)
    patterns = patterns[patterns["width"] >= minREpatternSize ].reset_index(drop=True)
    countRE = pd.concat([pd.DataFrame(patterns["names"]), pd.DataFrame(patterns["names"]).count(axis="columns")], 1)
    
    names1 = []
    width1 = []
    seq1 = []
    names2 = []
    width2 = []
    seq2 = []
    for i in range(len(patterns)):
        if countRE[0][i] == 1:
            if patterns["names"][i] in list(countRE["names"]):
                names1.append(patterns["names"][i])
                width1.append(patterns["width"][i])
                seq1.append(patterns["seq"][i])
        else:
            if patterns["names"][i] in list(countRE["names"]):
                names2.append(patterns["names"][i])
                width2.append(patterns["width"][i])
                seq2.append(patterns["seq"][i])
            
    singleRE = pd.concat([pd.DataFrame(width1, columns=(["width"])), pd.DataFrame(seq1, columns=(["seq"])), pd.DataFrame(names1, columns=(["names"]))], 1)
    duplicateRE = pd.concat([pd.DataFrame(width2, columns=(["width"])), pd.DataFrame(seq2, columns=(["seq"])), pd.DataFrame(names2, columns=(["names"]))], 1)
    duplicateRE = duplicateRE.sort_values(by=["names"])
    if len(duplicateRE) > 0:
         print("More than one pattern found for the following REs and these RE will be skipped! If you want to include them in the search, please correct the REpattern file and run this again!")
         print(duplicateRE.drop_duplicates(subset=["names"]).reset_index(drop=True))
    patterns = singleRE
    seqs = list(all_gRNAs["seq"])
    seq_names = list(all_gRNAs["names"])
    min_pStart_plus = min(overlap_gRNA_positions)-1
    max_pStart_plus = max(overlap_gRNA_positions)-1
    gRNAs_RE = pd.concat([pd.DataFrame([], columns=(["gRNAPlusPAM"])), pd.DataFrame([], columns=(["REcutgRNAName"])),
                          pd.DataFrame([], columns=(["REname"])), pd.DataFrame([], columns=(["REpattern"])),
                          pd.DataFrame([], columns=(["REcutStart"])), pd.DataFrame([], columns=(["REcutEnd"]))], 1)
    
    for j in range(len(patterns)):
        pattern_name = re.sub("'", "", patterns["names"][j])
        pattern = patterns["seq"][j]
        this_pattern_size = len(pattern)
        revpattern = str(Seq(pattern).reverse_complement())
       
        if revpattern != pattern:
            revpattern = tP.translatePattern(revpattern)
            def minus_gRNAs(all_gRNAs):
                minus_gRNAs = []
                for i in range(len(all_gRNAs)):
                    res1 = gr.gregexpr_index(revpattern, seqs[i]) 
                    if res1 == []:
                        res1 = [-1]
                    for k in range(len(res1)):
                        if res1[k] != -1:
                            res1[k] += 1
                        if res1[k] >= 0 and overlap_allpos != True and ((min_pStart_plus >= res1[k] and min_pStart_plus <= (res1[k] + this_pattern_size -1)) or (max_pStart_plus >= res1[k] and max_pStart_plus <= (res1[k] + this_pattern_size -1))):
                            res1[k] -= 1
                            minus_gRNAs.append([seqs[i], seq_names[i], pattern_name, str(Seq(patterns["seq"][j]).reverse_complement()), this_pattern_size - 1 + res1[k], res1[k]])
                        elif res1[k] >= 0 and overlap_allpos == True and min_pStart_plus >= res1[k] and max_pStart_plus <= (res1[k] + this_pattern_size -1):
                            res1[k] -= 1
                            minus_gRNAs.append([seqs[i], seq_names[i], pattern_name, str(Seq(patterns["seq"][j]).reverse_complement()), this_pattern_size - 1 + res1[k], res1[k]])
                           
                return pd.DataFrame(minus_gRNAs)
            
            minus_gRNAs = minus_gRNAs(all_gRNAs)
           
            if len(minus_gRNAs) > 1:
                #print(minus_gRNAs)
                #print(j)
                #print(np.shape(minus_gRNAs))
                minus_gRNAs = minus_gRNAs.rename(columns={0: "gRNAPlusPAM", 1: "REcutgRNAName", 2: "REname", 3: "REpattern", 4: "REcutStart", 5:"REcutEnd"})
                gRNAs_RE = pd.concat([gRNAs_RE,  minus_gRNAs], 0).reset_index(drop=True)
                  
        pattern = tP.translatePattern(pattern)
        
        def pos_gRNAs(all_gRNAs):
            pos_gRNAs = []
            for i in range(len(all_gRNAs)):
                res1 = gr.gregexpr_index(pattern, seqs[i])
            
                if res1 == []:
                    res1 = [-1]   
                for k in range(len(res1)):
                    if res1[len(res1) -1] <= min_pStart_plus:
                        if res1[k] != -1:
                            res1[k] += 1
                        if res1[k] >= 0 and overlap_allpos != True and ((min_pStart_plus >= res1[k] and min_pStart_plus <= (res1[k] + this_pattern_size -1)) or (max_pStart_plus >= res1[k] and max_pStart_plus <= (res1[k] + this_pattern_size -1))):
                            res1[k] -= 1
                            pos_gRNAs.append([seqs[i], seq_names[i], pattern_name, str(patterns["seq"][j]), res1[k]+1, this_pattern_size  + res1[k]])
                        elif res1[k] >= 0 and overlap_allpos == True and min_pStart_plus >= res1[k] and max_pStart_plus <= (res1[k] + this_pattern_size -1):
                            res1[k] -= 1
                            pos_gRNAs.append([seqs[i], seq_names[i], pattern_name, str(patterns["seq"][j]), res1[k]+1, this_pattern_size  + res1[k]])
          
            return pd.DataFrame(pos_gRNAs)
      
        pos_gRNAs = pos_gRNAs(all_gRNAs)
       
        if len(pos_gRNAs) > 0:
            pos_gRNAs = pos_gRNAs.rename(columns={0: "gRNAPlusPAM", 1: "REcutgRNAName", 2: "REname", 3: "REpattern", 4: "REcutStart", 5:"REcutEnd" })
            gRNAs_RE = pd.concat([gRNAs_RE,  pos_gRNAs], 0).reset_index(drop=True)
   
    if pairOutputFile != None and os.path.isfile(pairOutputFile) != None:
         gRNAs_RE_plus = gRNAs_RE.rename(columns={"gRNAPlusPAM": "ForwardgRNAPlusPAM", "REcutgRNAName": "ForwardREcutgRNAName", "REname": "ForwardREname", "REpattern": "ForwardREpattern", "REcutStart": "ForwardREcutStart", "REcutEnd":"ForwardREcutEnd"})
         gRNAs_RE_minus = gRNAs_RE.rename(columns={"gRNAPlusPAM": "ReversegRNAPlusPAM", "REcutgRNAName": "ReverseREcutgRNAName", "REname": "ReverseREname", "REpattern": "ReverseREpattern", "REcutStart": "ReverseREcutStart", "REcutEnd":"ReverseREcutEnd"})
         pairgRNAs = pd.read_csv(pairOutputFile, sep = "\t",)
         ann_gRNAs = pairgRNAs.merge(gRNAs_RE_plus, how = "left", on = "ForwardgRNAPlusPAM")
         ann_gRNAs = ann_gRNAs.merge(gRNAs_RE_minus, how = "left", on = "ReversegRNAPlusPAM")
         ann_gRNAs = ann_gRNAs.fillna("NA")
         #withRE = ann_gRNAs[["ReversegRNAPlusPAM", "ReversegRNAName", "ForwardgRNAPlusPAM", "ForwardgRNAName"]].drop_duplicates().reset_index(drop=True)
         withRE = ann_gRNAs.iloc[:,0:4].drop_duplicates().reset_index(drop=True)
       
         if np.shape(withRE)[0] == 0 and findgRNAsWithREcutOnly == True:
             raise Exception("No pairs with RE sites!")
         gRNAs = bio.DNAStringSet((list(withRE["ForwardgRNAPlusPAM"])+list(withRE["ReversegRNAPlusPAM"]))).reset_index(drop=True)
         gRNAs = pd.concat([gRNAs, pd.DataFrame(list(withRE["ForwardgRNAName"])+list(withRE["ReversegRNAName"]), columns=(["names"])).reset_index(drop=True)],1)
    else:
        if np.shape(gRNAs_RE)[0] == 0 and findgRNAsWithREcutOnly == True:
            raise Exception("No gRNAs with RE sites!")
        temp = pd.concat([gRNAs_RE["gRNAPlusPAM"], gRNAs_RE["REcutgRNAName"]], 1)
        gRNAs = pd.concat([bio.DNAStringSet(list(temp["gRNAPlusPAM"])).reset_index(drop=True), temp["REcutgRNAName"]],1).rename(columns={"REcutgRNAName": "names"})
        ann_gRNAs = gRNAs_RE
   
    gRNAs_withRE = gRNAs.drop_duplicates().reset_index(drop=True)
    gRNAREcutDetails = ann_gRNAs
    
 
    return gRNAs_withRE, gRNAREcutDetails


