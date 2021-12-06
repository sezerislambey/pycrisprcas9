#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 16:31:47 2021

@author: klaus
"""

from Bio.Seq import Seq
import pandas as pd
import re
import os
import numpy as np
import translatePattern as tp
import gregexprindex as gr
from Bio.Seq import Seq
import py2bit

pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",50)
pd.set_option("display.max_colwidth",50)
pd.set_option("display.width",None)
pd.set_option('expand_frame_repr', False)



def writeHits2(gRNA, seqname, matches, strand, file, chrom_len, mismatch_names, BSgenomeName,
               gRNA_size = 20, PAM = "NGG", PAM_pattern = "NRG$", max_mismatch = 4, appendd = False,
               PAM_location = "3prime", PAM_size = 3, allowed_mismatch_PAM = 1, targetBase = "C",
               baseEditing = False, editingWindow = [4, 8]):
    
    if not gRNA or type(gRNA) != Seq:
        raise Exception("gRNA is required as a Bio.Seq.Seq object!")
    if not seqname:
        raise Exception("seqname is required as character!")

    if strand != "+" and strand != "-": 
        raise Exception("strand is required as + or - !")
    if os.path.isfile(file) == appendd and appendd != True:
        raise Warning("existing file ", file, " will be overwritten with 'appendd = False'")
    if os.path.isfile(file) != appendd and appendd == True:
        appendd = False
    
    Lmismatch = []
    for j in range(len(matches)):
        TF = []
        for i in range(len(gRNA)):
            if matches[j][0][i] == gRNA[i]:
                TF.append(False)
            else:
                TF.append(True)
        Lmismatch.append(TF)
   
    for i in range(len(Lmismatch)):
        for j in range(len(Lmismatch[i])):
            if Lmismatch[i][j] == False and pd.isna(Lmismatch[i])[j] == False:
                Lmismatch[i][j] = 0
            else:
                Lmismatch[i][j] = 1

    if len(np.shape(Lmismatch)) == 2:
        Lmismatch = pd.DataFrame(Lmismatch)
    
    n_mismatch = []
    for i in range(len(Lmismatch)):
        n_mismatch.append(sum(Lmismatch.iloc[i][:]))
    
    colnames = []
    for i in range(len(Lmismatch.columns)):
        colnames.append("IsMismatch.pos" + str(i+1))
    
    Lmismatch.columns = colnames
   
    if PAM_location == "3prime":
        gRNAplusPAM = str(gRNA)+str(PAM)
    else:
        gRNAplusPAM = str(PAM)+str(gRNA)
    
    old_start = []
    old_end = []
    for i in range(len(matches)):
        old_start.append(matches[i].span()[0] + 1)
        old_end.append(matches[i].span()[1])
    new_start = []
    new_end = []
    if strand == "-":
        for i in range(len(old_start)): 
            new_start.append(chrom_len - old_end[i]+1)
            new_end.append(chrom_len - old_start[i]+1)
            if PAM_location == "3prime":
                new_start[i] = new_start[i] - PAM_size
            else:
                new_end[i] = new_end[i] + PAM_size
    
    elif PAM_location == "3prime":
        for i in range(len(old_start)):
            new_start.append(old_start[i])
            new_end.append(old_end[i] + PAM_size)
    else:
        for i in range(len(old_start)):
            new_start.append(old_start[i] - PAM_size)
            new_end.append(old_end[i])
    starts = []
    ends = []
    for i in range(len(new_start)):
        starts.append(max(new_start[i], 0))
        ends.append(min(new_end[i], chrom_len))
    
    OffTargetSequence = []
    if strand == "+":
        for i in range(len(starts)):
            OffTargetSequence.append(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(seqname, starts[i]-1, ends[i]))
    if strand == "-":
        for i in range(len(starts)):
            OffTargetSequence.append(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(seqname)[::-1])
            OffTargetSequence[i] = str(Seq(OffTargetSequence[i]).reverse_complement())[starts[i]-1:ends[i]][::-1]
            
   
    hits = pd.concat([pd.DataFrame([strand] * len(starts), columns=(["strand"])),
                      pd.DataFrame([seqname] * len(starts), columns=(["chrom"])),
                      pd.DataFrame(starts, columns=(["chromStart"])),
                      pd.DataFrame(ends, columns=(["chromEnd"])),
                      pd.DataFrame(mismatch_names, columns=(["name"])),
                      pd.DataFrame([gRNAplusPAM] * len(starts), columns=(["gRNAPlusPAM"])),
                      pd.DataFrame(OffTargetSequence, columns=(["OffTargetSequence"])),
                      pd.DataFrame(n_mismatch, columns=(["n.mismatch"])),
                      pd.DataFrame([chrom_len] * len(starts), columns=(["chrom.len"]))], 1)
    hits = pd.concat([Lmismatch, hits], axis=1)

    hits =  hits[(hits["n.mismatch"] <= max_mismatch) & (hits["chromEnd"] + 1 - hits["chromStart"] == gRNA_size + PAM_size)].reset_index(drop=True)
  
    

    PAM_pattern = tp.translatePattern(PAM_pattern)
    containPAM = []
    if np.shape(hits)[0] > 0:
        for i in range(np.shape(hits)[0]):
           
            for j in range(len(PAM_pattern)):
                pos_plus = re.findall(PAM_pattern[j], hits.loc[i, ("OffTargetSequence")])
                if len(pos_plus) > 0:
                    break
            if len(pos_plus) > 0:
                containPAM.append(1)
            else:
                containPAM.append(0)
        
        
        hits["containPAM"] = pd.DataFrame(containPAM)
        hits = hits[hits["containPAM"] == 1].reset_index(drop=True)
        hits.drop("containPAM", axis = 1, inplace = True)
      
      
       
        if np.shape(hits)[0] > 0:
            if baseEditing == True:
                n_targetBase = []
                for i in range(np.shape(hits)[0]):
                    n_targetBase.append(hits.loc[i, ("OffTargetSequence")][min(editingWindow)-1:max(editingWindow)])
                
                n_targetBase1 = []
                for i in range(len(n_targetBase)):
                    n_targetBase1.append(n_targetBase[i].count(targetBase))
                
                hits = pd.concat([hits, pd.DataFrame(n_targetBase1, columns = (["n.targetBase1"]))], 1)
                hits = hits[hits["n.targetBase1"] > 0].reset_index(drop=True)
                hits.drop("n.targetBase1", axis = 1, inplace = True) 
        
        if np.shape(hits)[0] > 0:
            PAM_sequence = []
            if PAM_location == "3prime":
                for i in range(np.shape(hits)[0]):
                    PAM_sequence.append(hits.loc[i, ("OffTargetSequence")][gRNA_size:gRNA_size + PAM_size])
            else:
                for i in range(np.shape(hits)[0]):
                    PAM_sequence.append(hits.loc[i, ("OffTargetSequence")][0:PAM_size])

            n_PAM_mismatch = []
            for i in range(len(PAM_sequence)):
                n_PAM_mismatch.append(PAM.count(PAM_sequence[i]))
    
          
            hits["n_PAM_mismatch"] = pd.DataFrame(n_PAM_mismatch)
            hits = hits[hits["n_PAM_mismatch"] <= allowed_mismatch_PAM]
            hits.drop("n_PAM_mismatch", axis = 1, inplace = True)
    
            if np.shape(hits)[0] > 0:
                forViewInUCSC = hits["chrom"]
                score = np.repeat(100, np.shape(hits)[0])
                del hits["chrom.len"]
                hits["forViewInUCSC"] = forViewInUCSC
                hits["score"] = pd.Series(score)
                hits["forViewInUCSC"] = hits["forViewInUCSC"].map(str) + ":" + hits["chromStart"].map(str) + "-" + hits["chromEnd"].map(str)
                if os.path.isfile(file) == False:
                    hits.to_csv(file, sep = "\t", index = False)
                else:
                    hits.to_csv(file, sep = "\t", index = False, mode = "a+", header=None)
   
    #return hits 
    
    
