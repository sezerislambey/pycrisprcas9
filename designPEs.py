#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 09:58:26 2021

@author: klaus
"""


import numpy as np
import pandas as pd
import biostrings as bio
from Bio.Seq import Seq
import re


def _getSeq(inputSeq, starts, length, strand):
     seq_len = len(inputSeq)
     n_rec = len(starts)
     ends = []
     seq = []
     for i in range(n_rec):
         if strand =="minus":
             ends.append(max(1, min(seq_len, starts[i]) +1))   
             starts[i] = max(1, starts[i] - length[i] + 1)
         else:
             ends.append(max(1, min(seq_len, (starts[i] + length[i]))))
             starts[i] = max(1, starts[i])
         seq.append(str(inputSeq[starts[i]:ends[i]]))
     
     seq = bio.DNAStringSet(seq)
     return seq


test = 0
if test == 1:
    potential_PEs = designPEs(inputSeq, PAM_size = PAM_size, PAM_location = PAM_location, 
                gRNA_size = gRNA_size, overlap_gRNA_positions = overlap_gRNA_positions,
                PBS_length = PBS_length, paired_gRNAs = paired_gRNAs, 
                RT_template_length = RT_template_length, RT_template_pattern = RT_template_pattern,
                corrected_seq = corrected_seq, targeted_seq_length_change = targeted_seq_length_change,
                bp_after_target_end = bp_after_target_end, target_start = target_start,
                target_end = target_end, primeEditingPaired_output =  primeEditingPaired_output)


def designPEs(inputSeq, target_start, target_end, paired_gRNAs, RT_template_length, RT_template_pattern, 
              corrected_seq, targeted_seq_length_change, PAM_size = 3, PAM_location = "3prime", 
              gRNA_size = 20, overlap_gRNA_positions = [17, 18], PBS_length = 13, bp_after_target_end = 15,
              primeEditingPaired_output = "pairedgRNAsPE.xls", appendd = False, col_names = True):
    
    ###Hsap_GATA1_ex2_gR7r or Hsap_GATA1_ex2_gR17f
    ###Obtain the cut_start from the names of gRNAs
    ### insert just before target_start
    ### mutate/delete all bases between target_start and target_end inclusive
    
    def cutSiteIn5primeOfTargetSite(cutDistanceFromTarget):
        for i in range(len(cutDistanceFromTarget)):
            if cutDistanceFromTarget[i] < 0:
                cutDistanceFromTarget[i] = True
            else:
                cutDistanceFromTarget[i] = False
        
        return cutDistanceFromTarget
    

    cut_start_plus = []
    cut_start_minus = []
    gRNAs_plus_cut_d_target_start = []
    gRNAs_minus_cut_d_target_end = []
    actual_RT_template_length_plus = []
    actual_RT_template_length_minus = []
    ForwardgRNA_PAM_start = []
    ForwardgRNA_PAM_end = []
    ReversegRNA_PAM_start = []
    ReversegRNA_PAM_end = []
    for i in range(len(paired_gRNAs)):
        cut_start_plus.append(paired_gRNAs["ForwardgRNAName"][i][paired_gRNAs["ForwardgRNAName"][i].find("_gR")+3:])
        #cut_start_plus.append(re.sub("NA_gR", "", paired_gRNAs["ForwardgRNAName"][i]))
        cut_start_plus[i] = int(re.sub("f", "", cut_start_plus[i]))
        #cut_start_minus.append(re.sub("NA_gR", "", paired_gRNAs["ReversegRNAName"][i]))
        cut_start_minus.append(paired_gRNAs["ReversegRNAName"][i][paired_gRNAs["ReversegRNAName"][i].find("_gR")+3:])
        cut_start_minus[i] = int(re.sub("r", "", cut_start_minus[i]))
    
        gRNAs_plus_cut_d_target_start.append(cut_start_plus[i] - target_start)
        gRNAs_minus_cut_d_target_end.append(target_end - cut_start_minus[i])
       
        actual_RT_template_length_plus.append(target_start - cut_start_plus[i]  - 1 +  bp_after_target_end) 
        actual_RT_template_length_minus.append(cut_start_minus[i] - target_end -1  + bp_after_target_end)
        ForwardgRNA_PAM_start.append(cut_start_plus[i] + gRNA_size - min(overlap_gRNA_positions) + 1)
        ForwardgRNA_PAM_end.append(ForwardgRNA_PAM_start[i] + PAM_size - 1)
        ReversegRNA_PAM_start.append(cut_start_minus[i]- gRNA_size + min(overlap_gRNA_positions) - 1)
        ReversegRNA_PAM_end.append(ReversegRNA_PAM_start[i] - PAM_size + 1)
    

    paired_gRNAs = pd.concat([paired_gRNAs, pd.DataFrame(cut_start_minus, columns=(["ReversegRNA.cut.start"])), 
                        pd.DataFrame(cut_start_plus, columns=(["ForwardgRNA.cut.start"])),	
                        pd.DataFrame(ReversegRNA_PAM_start, columns=(["ReversegRNA.PAM.start"])),
                        pd.DataFrame(ReversegRNA_PAM_end,  columns=(["ReversegRNA.PAM.end"])),
                        pd.DataFrame(ForwardgRNA_PAM_start, columns=(["ForwardgRNA.PAM.start"])),
                        pd.DataFrame(ForwardgRNA_PAM_end, columns=(["ForwardgRNA.PAM.end"])),
                        pd.DataFrame(gRNAs_minus_cut_d_target_end, columns=(["ReversegRNA.cut.start.d.targetEnd"]))*-1,
                        pd.DataFrame(gRNAs_plus_cut_d_target_start, columns=(["ForwardgRNA.cut.start.d.targetStart "])),
                        pd.DataFrame(cutSiteIn5primeOfTargetSite(gRNAs_minus_cut_d_target_end), columns=(["ReversegRNA.cut.5prime.targetEnd"])),
                        pd.DataFrame(cutSiteIn5primeOfTargetSite(gRNAs_plus_cut_d_target_start), columns=(["ForwardgRNA.cut.5prime.targetStart"])),
                        pd.DataFrame(actual_RT_template_length_minus, columns=(["ReversegRNA.RT.template.length"])),
                        pd.DataFrame(actual_RT_template_length_plus, columns=(["ForwardgRNA.RT.template.length"]))], 1)
    paired_gRNAs0 = []
    for i in range(len(paired_gRNAs)):
        if (paired_gRNAs["ReversegRNA.cut.5prime.targetEnd"][i] == True and 
            paired_gRNAs["ReversegRNA.RT.template.length"][i]  <= max(RT_template_length) and 
            paired_gRNAs["ReversegRNA.RT.template.length"][i] >= min(RT_template_length)-1 and 
            (paired_gRNAs["ForwardgRNA.PAM.end"][i] < target_start or 
             paired_gRNAs["ForwardgRNA.PAM.start"][i] > target_end)) or (paired_gRNAs["ForwardgRNA.cut.5prime.targetStart"][i] == True and 
         paired_gRNAs["ForwardgRNA.RT.template.length"][i] <= max(RT_template_length) and
         paired_gRNAs["ForwardgRNA.RT.template.length"][i] >= min(RT_template_length)-1 and 
         (paired_gRNAs["ReversegRNA.PAM.start"][i] < target_start or 
          paired_gRNAs["ReversegRNA.PAM.end"][i] > target_end)):
            paired_gRNAs0.append(paired_gRNAs.loc[i])
    paired_gRNAs = pd.DataFrame(paired_gRNAs0).reset_index(drop=True)
    
    PBS_plus = _getSeq(inputSeq, starts = list(paired_gRNAs["ForwardgRNA.cut.start"] - PBS_length),
                       length = np.repeat(PBS_length, len(paired_gRNAs["ForwardgRNA.cut.start"])),
                       strand = "plus").reset_index(drop=True)
    for i in range(len(PBS_plus)):
        PBS_plus.loc[i, ("seq")] =  str(Seq(PBS_plus.loc[i, ("seq")]).reverse_complement())
    
    
    PBS_minus = _getSeq(inputSeq, starts = list(paired_gRNAs["ReversegRNA.cut.start"] - 1),
                       length = np.repeat(PBS_length, len(paired_gRNAs["ReversegRNA.cut.start"])),
                       strand = "plus").reset_index(drop=True)
    for i in range(len(PBS_minus)):
        PBS_minus.loc[i, ("seq")] =  PBS_minus.loc[i, ("seq")][::-1]
    
  
    actual_RT_template_plus_seq1 = _getSeq(inputSeq, starts = list(paired_gRNAs["ForwardgRNA.cut.start"]), 
                                          length = target_start - paired_gRNAs["ForwardgRNA.cut.start"] - 1,
                                          strand = "plus")
    
    actual_RT_template_minus_seq1 = _getSeq(inputSeq, starts = list(paired_gRNAs["ReversegRNA.cut.start"] - 2), 
                                          length = paired_gRNAs["ReversegRNA.cut.start"] - target_end - 1,
                                          strand = "minus")
    
    if targeted_seq_length_change < 0:
        actual_RT_template_plus_seq2 = list(_getSeq(inputSeq, starts = [target_end], 
                                               length = [bp_after_target_end], strand = "plus")["seq"])[0]
        actual_RT_template_minus_seq2 = list(_getSeq(inputSeq, starts = [target_start - 2],
                                                length = [bp_after_target_end], strand = "minus")["seq"])[0]
    elif targeted_seq_length_change == 0:
        paired_gRNAs["ReversegRNA.RT.template.length"] = paired_gRNAs["ReversegRNA.RT.template.length"] + target_end - target_start + 1
        paired_gRNAs["ForwardgRNA.RT.template.length"] = paired_gRNAs["ForwardgRNA.RT.template.length"] + target_end - target_start + 1
        actual_RT_template_plus_seq2 = corrected_seq + list(_getSeq(inputSeq, starts = [target_end], length = [bp_after_target_end], strand = "plus")["seq"])[0]
        actual_RT_template_minus_seq2 = corrected_seq + list(_getSeq(inputSeq, starts = [target_start - 2], length = [bp_after_target_end], strand = "minus")["seq"])[0]
        
    else:
        paired_gRNAs["ReversegRNA.RT.template.length"] = paired_gRNAs["ReversegRNA.RT.template.length"] + targeted_seq_length_change
        paired_gRNAs["ForwardgRNA.RT.template.length"] = paired_gRNAs["ForwardgRNA.RT.template.length"] + targeted_seq_length_change
        actual_RT_template_plus_seq2 = corrected_seq + list(_getSeq(inputSeq, starts = [target_start - 1], length = [bp_after_target_end], strand = "plus")["seq"])[0]
        actual_RT_template_minus_seq2 = corrected_seq + list(_getSeq(inputSeq, starts = [target_end - 1], length = [bp_after_target_end], strand = "minus")["seq"])[0]
    
    actual_RT_template_seq_plus = []
    actual_RT_template_seq_minus = []
    for i in range(len(actual_RT_template_plus_seq1)):
        actual_RT_template_seq_plus.append(str(Seq(list(actual_RT_template_plus_seq1["seq"])[i] + str(actual_RT_template_plus_seq2)).reverse_complement()))
        actual_RT_template_seq_minus.append(str(Seq(str(actual_RT_template_minus_seq2) + list(actual_RT_template_minus_seq1["seq"])[i]))[::-1])
    actual_RT_template_seq_plus = bio.DNAStringSet(actual_RT_template_seq_plus)
    actual_RT_template_seq_minus = bio.DNAStringSet(actual_RT_template_seq_minus)
    
    paired_gRNAs = pd.concat([paired_gRNAs, pd.DataFrame(list(actual_RT_template_seq_minus["seq"]), columns=(["ReversegRNA.RT.template.seq"])),
                        pd.DataFrame(list(actual_RT_template_seq_plus["seq"]), columns=(["ForwardgRNA.RT.template.seq"])), 
                        pd.DataFrame(list(PBS_minus["seq"]), columns=(["ReversegRNA.PBS"])),
                        pd.DataFrame(list(PBS_plus["seq"]), columns=(["ForwardgRNA.PBS"]))], 1)
    
    #list(paired_gRNAs, actual_RT_template_plus_seq1, actual_RT_template_plus_seq2, actual_RT_template_minus_seq1, actual_RT_template_minus_seq2) 
 
    paired_gRNAs.to_csv(primeEditingPaired_output, sep = "\t", index = False)
    
    return paired_gRNAs