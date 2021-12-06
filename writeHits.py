#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:53:59 2021

@author: klaus
"""

from Bio.Seq import Seq
import pandas as pd
import re
import os
import numpy as np
import translatePattern as tp
import gregexprindex as gr

pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",50)
pd.set_option("display.max_colwidth",50)
pd.set_option("display.width",None)



def writeHits(gRNA, seqname, matches, strand, file, chrom_len, seqs, gRNA_size = 20, 
            PAM = "NGG", PAM_pattern = "NRG", max_mismatch = 4, 
            appendd = False,  PAM_location = "3prime", PAM_size = 3,
            allowed_mismatch_PAM = 1, targetBase = "C"):
   
    if not gRNA or type(gRNA) != Seq:
        raise Exception("gRNA is required as a Bio.Seq.Seq object!")
    if not seqname:
        raise Exception("seqname is required as character!")
    if matches == None or type(matches) != re.Match: 
        raise Exception("matches is required as re.Match object!")
    if strand != "+" and strand != "-": 
        raise Exception("strand is required as + or - !")
    if os.path.isfile(file) == appendd and appendd != False:
        raise Warning("existing file ", file, " will be overwritten with 'appendd = False'")
    if os.path.isfile(file) != appendd and appendd == True:
        appendd = False
    
    Lmismatch = []
    TF = []
    for i in range(len(gRNA)):
        if m[0][i] == gRNA[i]:
            TF.append(False)
        else:
            TF.append(True)
    Lmismatch.append(TF)
    
    for i in range(len(Lmismatch[0])):
        if Lmismatch[0][i] == False and pd.isna(Lmismatch[0])[i] == False:
            Lmismatch[0][i] = 0
        else:
            Lmismatch[0][i] = 1
    
    if len(np.shape(Lmismatch)) == 2:
        Lmismatch = pd.DataFrame(Lmismatch)
    
    if PAM_location == "3prime":
        Lmismatch = Lmismatch.iloc[:,:gRNA_size]

    elif np.shape(Lmismatch)[1] == gRNA_size + PAM_size:
        start_pos = PAM_size
        end_pos = PAM_size + gRNA_size
        Lmismatch = Lmismatch.iloc[:,start_pos:end_pos]

    n_mismatch = sum(Lmismatch.iloc[0][:])
    
    colnames = []
    for i in range(len(Lmismatch.columns)):
        colnames.append("IsMismatch.pos" + str(i+1))
    
    Lmismatch.columns = colnames
    
    if PAM_location == "3prime":
        gRNAPlusPAM = str(gRNA)+str(PAM)
    else:
        gRNAPlusPAM = str(PAM)+str(gRNA)
   
    old_start = matches.span(0)[0]
    old_end = matches.span(0)[1]
    
    if strand == "-":
       new_start1 = chrom_len - old_end
       new_end1 = chrom_len - old_start
       if PAM_location == "3prime":
           new_start1 = new_start1 - PAM_size
       else:
           new_end1 = new_end1 + PAM_size
           
    if PAM_location == "3prime":
        new_start = old_start
        new_end = old_end + PAM_size
    else:
        new_start = old_start - PAM_size
        new_end = old_end
    
    starts = max(new_start ,0)
    ends = min(new_end, chrom_len)  
    # sequence fetch is the same for plus and minus strand
    # because, revcomplement of the sequences (seqs) are used for minus strand
    
    OffTargetSequence = seqs[starts:ends]
 
    # coordinate needs to be changed for minus strand
    
    if strand == "-":
        starts = max(new_start1, 0)
        ends = min(new_end1, chrom_len)
    
    
    hits = pd.DataFrame([[strand, seqname, starts+1, ends, m_names, str(gRNAPlusPAM), str(OffTargetSequence), n_mismatch,
        chrom_len]], columns=(["strand", "chrom", "starts", "ends", "name", "gRNAPlusPAM", "OffTargetSequence", "n.mismatch", "chrom.len"]))
    
    hits = pd.concat([Lmismatch, hits], axis=1)
  
    
    if int(hits["n.mismatch"]) <= max_mismatch and int(hits["ends"])+1 - int(hits["starts"]) == gRNA_size + PAM_size:
        hits = hits[hits["n.mismatch"] <= max_mismatch]
    else:
        hits = hits[:0]
    
    PAM_pattern = tp.translatePattern(PAM_pattern)
    
    if np.shape(hits)[0] > 0:
        for i in range(np.shape(hits)[0]):
            pos_plus = gr.gregexpr_index(PAM_pattern, hits.iloc[i]["OffTargetSequence"])
            if len(pos_plus) > 0:
                containPAM = 1
            else:
                containPAM = 0
           
        if containPAM == 1:
            hits = hits
        else:
            hits = hits[:0]
        
        if np.shape(hits)[0] > 0:
            if baseEditing == True:
                n_targetBase = []
                for i in range(np.shape(hits)[0]):
                    n_targetBase.append(hits.loc[i, ("OffTargetSequence")][min(editingWindow):max(editingWindow)])
                
                n_targetBase1 = []
                for i in range(len(n_targetBase)):
                    n_targetBase1.append(n_targetBase[i].count(targetBase))
                
                hits = pd.concat([hits, pd.DataFrame(n_targetBase1, columns = (["n.targetBase1"]))], 1)
                hits = hits[hits["n.targetBase1"] > 0]
                hits.drop("n.targetBase1", axis = 1, inplace = True) 
        
        if np.shape(hits)[0] > 0:
            if PAM_location == "3prime":
                PAM_sequence = hits["OffTargetSequence"][0][gRNA_size:gRNA_size + PAM_size]
                
            else:
                PAM_sequence = hits["OffTargetSequence"][0][0:PAM_size]
        
        
            if PAM_sequence == PAM[:len(PAM_sequence)]:
                n_PAM_mismatch = 0
            else:
                n_PAM_mismatch = 1
            
            if n_PAM_mismatch <= allowed_mismatch_PAM:
                hits = hits
            else:
                hits = hits[:0]
           
            forViewInUCSC = hits["chrom"]
            score = np.repeat(100, np.shape(hits)[0])
            del hits["chrom.len"]
            hits["forViewInUCSC"] = forViewInUCSC
            hits["score"] = pd.Series(score)
            hits["forViewInUCSC"] = hits["forViewInUCSC"].map(str) + ":" + hits["starts"].map(str) + "-" + hits["ends"].map(str)
  
            
            hits.to_csv(file, sep = "\t", index = False)
                        
         
    return hits 
    
    
    
