#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 13:40:10 2021

@author: klaus
"""

import os
import numpy as np
import pandas as pd
import biostrings as bio
from Bio.Seq import Seq
from translatePattern import translatePattern
import re
import designPEs
import gregexprindex as gr
import character as ch
from calculategRNAEfficiency import calculategRNAEfficiency


pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)



def getgRNA_cut_sites(subject, subjectname, PAM ="NGG", gRNA_pattern = "", gRNA_size = 20, cut_site = 17,
                               PAM_size = 3, calculategRNAEfficacy = True, baseBeforegRNA = 4,
                               baseAfterPAM = 3, reverse_subject = False, PAM_location = "3prime",
                               baseEditing = False, targetBase = "C", editingWindow = [4, 8]):
    seq_len = len(subject)
    pos_PAMs0 = gr.gregexpr_index(PAM, subject)
    pos_PAMs0.sort()
    pos_PAMs = []

    if PAM_location == "3prime":
    
        for i in range(len(pos_PAMs0)):
            if pos_PAMs0[i] >= gRNA_size:
                pos_PAMs.append(pos_PAMs0[i])
               
    else:
        for i in range(len(pos_PAMs0)):
            if pos_PAMs0[i] != -1:
                pos_PAMs.append(pos_PAMs0[i])
                
    if len(pos_PAMs) > 0:
        if PAM_location == "3prime":
            starts_gRNA0 = []
            for i in range(len(pos_PAMs)):
                starts_gRNA0.append(pos_PAMs[i] - gRNA_size)
        else:
            starts_gRNA0 = []
            for i in range(len(pos_PAMs)):
                starts_gRNA0.append(pos_PAMs[i] + PAM_size)
       
        ends_gRNA = []
        starts_gRNA = []
       
        for i in range(len(pos_PAMs)):
           
            if starts_gRNA0[i] <= seq_len - gRNA_size + 1 and starts_gRNA0[i] >= 0:
                starts_gRNA.append(starts_gRNA0[i])
                ends_gRNA.append(starts_gRNA[i] + gRNA_size - 1)
                
    
        if gRNA_pattern != "":
            gRNA_seqs = []
            for i in range(len(starts_gRNA)):
                gRNA_seqs.append(str(subject[starts_gRNA[i]:ends_gRNA[i]+1]))
           
            pos_set2 = []
            for i in range(len(gRNA_seqs)):
               for j in range(len(gRNA_seqs[i])):
                   pos_set2.append(j)
          
            for i in range(len(pos_PAMs)):
                if pos_set2[i] != 0:
                    pos_PAMs[i] = None
            
            pos_PAMs = [i for i in pos_PAMs if i is not None]
        
        if baseEditing and len(pos_PAMs) > 0:
            if PAM_location == "3prime":
                starts_gRNA= []
                for i in range(len(pos_PAMs)):
                    starts_gRNA.append(pos_PAMs[i] - gRNA_size)
               
            else:
                starts_gRNA = []
                for i in range(len(pos_PAMs)):
                    starts_gRNA.append(pos_PAMs[i] + PAM_size)
            
            ends_gRNA = []
            for i in range(len(starts_gRNA)):
                ends_gRNA.append(starts_gRNA[i] + gRNA_size)
            
            for i in range(len(pos_PAMs)):
                if ends_gRNA[i] <= len(subject):
                    pos_PAMs[i] = pos_PAMs[i]
                    starts_gRNA[i] = starts_gRNA[i]
                    ends_gRNA[i] = ends_gRNA[i]
                else:
                    pos_PAMs[i] = None
                    starts_gRNA[i] = None
                    ends_gRNA[i] = None
            pos_PAMs = [i for i in pos_PAMs if i is not None]
            starts_gRNA = [i for i in starts_gRNA if i is not None]
            ends_gRNA = [i for i in ends_gRNA if i is not None]
            gRNA_seqs = []
            for i in range(len(starts_gRNA)):
                gRNA_seqs.append(str(subject[starts_gRNA[i]:ends_gRNA[i]]))
         
            n_targetBase = []
            for i in range(len(gRNA_seqs)):
                n_targetBase.append(gRNA_seqs[i][min(editingWindow)-1:max(editingWindow)])
           
            n_targetBase1 = []
            for i in range(len(n_targetBase)):
                n_targetBase1.append(n_targetBase[i].count(targetBase))
        
            
            if 1 not in n_targetBase1:
                pos_PAMs = []
                starts_gRNA = []
                ends_gRNA = []
            else:
                for i in range(len(n_targetBase1)):
                    if n_targetBase1[i] == 1:
                        pos_PAMs = [pos_PAMs[i]]
                        starts_gRNA = [starts_gRNA[i]]
                        ends_gRNA = [ends_gRNA[i]]

        seq = []
        if len(pos_PAMs) > 0:
            if PAM_location == "3prime":
                for i in range(len(starts_gRNA)):
                    seq.append(str(subject[starts_gRNA[i]:ends_gRNA[i] + PAM_size + 1]))
            else:
                 for i in range(len(starts_gRNA)):
                    seq.append(str(subject[starts_gRNA[i] - PAM_size:ends_gRNA[i] + 1]))
            
            extendedSequence = seq
                 
            if calculategRNAEfficacy == True:
                extended_starts = []
                for i in range(len(starts_gRNA)):
                    extended_starts.append(starts_gRNA[i] - baseBeforegRNA)
                    if extended_starts[i] < 0:
                        extended_starts[i] = 0
                        
                extended_ends = []
                if PAM_location == "3prime":
                    for i in range(len(starts_gRNA)):
                        extended_ends.append(ends_gRNA[i] + PAM_size + baseAfterPAM)
                else:
                    for i in range(len(starts_gRNA)):
                        extended_ends.append(starts_gRNA[i] + baseAfterPAM - 1)
                
                for i in range(len(extended_ends)):
                    if extended_ends[i] > len(subject):
                        extended_ends[i] = len(subject)
                
                extendedSequence = []
                for i in range(len(extended_ends)):
                    extendedSequence.append(str(subject[extended_starts[i]:extended_ends[i] + 1]))
            
             
            if reverse_subject == True:
                for i in range(len(seq)):
                    if len(seq) != len(subjectname):
                        subjectname = subjectname*len(seq)
                
                gRNAs_cut = pd.concat([pd.DataFrame(seq, columns=(["seq"])), pd.DataFrame(subjectname, columns=(["name"]))["name"].str.cat("_gR"+(seq_len-pd.DataFrame(starts_gRNA)-cut_site+1).astype(str)+"r"), seq_len-pd.DataFrame(starts_gRNA, columns=(["starts_gRNA"])), pd.DataFrame(["-"]*len(seq), columns=(["strand"])), pd.DataFrame(extendedSequence, columns=(["extendedSequence"]))], 1).dropna(thresh=2)
               
            else:
                for i in range(len(seq)):
                    if len(seq) != len(subjectname):
                        subjectname = subjectname*len(seq)
                    
                gRNAs_cut = pd.concat([pd.DataFrame(seq, columns=(["seq"])), pd.DataFrame(subjectname, columns=(["name"]))["name"].str.cat("_gR"+(pd.DataFrame(starts_gRNA)+cut_site).astype(str)+"f"), pd.DataFrame(starts_gRNA, columns=(["starts_gRNA"]))+1, pd.DataFrame(["+"]*len(seq), columns=(["strand"])), pd.DataFrame(extendedSequence, columns=(["extendedSequence"]))], 1).dropna(thresh=2)
                
        else:
            gRNAs_cut = ""
    else:
        gRNAs_cut = ""
    
    return gRNAs_cut

## Very inefficient (perform in quadratic time).





def compute_pair_index3(plus_start, minus_start, min_gap, max_gap):
    subject = pd.concat([pd.DataFrame(minus_start, columns = (["start"])), pd.DataFrame(minus_start, columns = (["end"])), pd.DataFrame([1]*len(minus_start), columns = (["width"]))], 1)
    query = pd.concat([pd.DataFrame(plus_start, columns = (["start"])), pd.DataFrame(plus_start, columns = (["end"])), pd.DataFrame([1]*len(plus_start), columns = (["width"]))], 1)
    queryHits = []
    subjectHits = []
    d = []
    gap = []
    
    for j in range(len(query)):
        for i in range(len(subject)):
            if abs(subject.iloc[i][0] - query.iloc[j][0]) <= max_gap:    
                d.append(query.iloc[j][0] - subject.iloc[i][0])
                
                if min_gap <= query.iloc[j][0] - subject.iloc[i][0] <= max_gap:
                    queryHits.append(j)
                    subjectHits.append(i)
                    gap.append(query.iloc[j][0] - subject.iloc[i][0])
    hits = pd.concat([pd.DataFrame(queryHits, columns = (["queryHits"])), pd.DataFrame(subjectHits, columns=(["subjectHits"]))], 1)
   
    return [list(hits["queryHits"]), list(hits["subjectHits"]), list(gap)]
    

    
def findgRNAs(inputFilePath, pairOutputFile, efficacyFile, chrom_acc = None, corrected_seq = None,
              featureWeightMatrixFile = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "DoenchNBT2014.csv"])),
              targeted_seq_length_change = None, target_start = None, target_end = None,
              rule_set = "Root_RuleSet1_2014", baseEditing = False, targetBase = "C", format = "fasta",
              PAM = "NGG", PAM_size = 3, findPairedgRNAOnly = False, annotatePaired = True,
              paired_orientation = "PAMout", enable_multicore = False, n_cores_max = 6, gRNA_pattern = "",
              gRNA_size = 20, overlap_gRNA_positions = [17, 18], primeEditing = False, PBS_length = 13,
              RT_template_length = [8, 28], 
              RT_template_pattern = "D$", bp_after_target_end = 15, editingWindow = [4, 8],
              primeEditingPaired_output = "pairedgRNAsForPE.xls", min_gap = 0, max_gap = 20, name_prefix = "", 
              baseBeforegRNA = 4, baseAfterPAM = 3, calculategRNAEfficacy = False, PAM_location = "3prime"):
    
    
    def DNAStringSet(seq):
        
        width = []
        for i in range(len(seq)):
            width.append(len(seq[i]))
                
        name = ["NA"]
        DNAStringSet = pd.concat([pd.DataFrame(width, columns=["width"], index=range(1, len(width)+1)), pd.DataFrame(seq, columns=["seq"],index=range(1, len(seq)+1)), pd.DataFrame(name, columns=["name"],index=range(1, len(seq)+1))],axis=1)
    
        return DNAStringSet

    if type(inputFilePath) == list:
        width = []
        for i in range(len(inputFilePath)):
            width.append(len(inputFilePath[i]))
                
        name = ["NA"]
        subjects= pd.concat([pd.DataFrame(width, columns=["width"], index=range(1, len(width)+1)), pd.DataFrame(inputFilePath, columns=["seq"],index=range(1, len(inputFilePath)+1)), pd.DataFrame(name, columns=["name"],index=range(1, len(inputFilePath)+1))],axis=1)
        
    else:
        name = []
        width = []
        chrom = []
        for i in range(np.shape(inputFilePath)[0]):
            width.append(len(inputFilePath.iloc[0][i]))
            chrom.append(inputFilePath.iloc[0][i].upper())
            name.append(inputFilePath["names"][i])
        subjects = pd.concat([pd.DataFrame(width, columns=["width"], index=range(1, len(width)+1)), pd.DataFrame(chrom, columns=["seq"],index=range(1, len(inputFilePath)+1)), pd.DataFrame(name, columns=["name"],index=range(1, len(inputFilePath)+1))],axis=1)
   
    cut_site = min(overlap_gRNA_positions)
    def match_arg(x, lst):
        return [el for el in lst if x in el]
    
    
    PAM = translatePattern(PAM)
    
    gRNA_pattern = translatePattern(gRNA_pattern)
    if gRNA_pattern == [""]:
        gRNA_pattern = ""

    min_subject = gRNA_size + PAM_size 
    subjects = subjects[subjects["width"] >= min_subject]
    if len(subjects) == 0:
        raise Exception("The input file contains no sequence! This could be caused by wrong format of the file. If file is created in mac, you could reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the command line")
    
    toAppend = False
    colNames = True
    subjects["name"] = re.sub("'", "", subjects["name"].values[0])
    subjects["name"] = re.sub(" ", "", subjects["name"].values[0])
    subjects["name"] = re.sub("\t", ":", subjects["name"].values[0])
    
    
    
    subject = Seq(str(subjects.iloc[0][1]))
    subjectname = list(subjects["name"])
   
    revsubject = subject.reverse_complement()
    plus_gRNAs = getgRNA_cut_sites(subject, subjectname, PAM = PAM, 
                         gRNA_pattern = gRNA_pattern, 
                         gRNA_size = gRNA_size,
                         cut_site = cut_site,
                         PAM_size = PAM_size, 
                         calculategRNAEfficacy = calculategRNAEfficacy,
                         baseBeforegRNA = baseBeforegRNA,
                         baseAfterPAM = baseAfterPAM,
                         PAM_location = PAM_location,
                         baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
    
    minus_gRNAs = getgRNA_cut_sites(revsubject, subjectname, PAM = PAM,
                         gRNA_pattern = gRNA_pattern,
                         gRNA_size = gRNA_size,
                         cut_site = cut_site,
                         PAM_size = PAM_size,
                         calculategRNAEfficacy = calculategRNAEfficacy,
                         baseBeforegRNA = baseBeforegRNA,
                         baseAfterPAM = baseAfterPAM,
                         reverse_subject = True,
                         PAM_location = PAM_location,
                         baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
   
    if len(plus_gRNAs) > 0:
        n_plus_gRNAs = np.shape(plus_gRNAs)[0]
    else:
        n_plus_gRNAs = 0
    if len(minus_gRNAs) > 0:
        n_minus_gRNAs = np.shape(minus_gRNAs)[0]
    else:
        n_minus_gRNAs = 0
    
    if n_minus_gRNAs == 0 and n_plus_gRNAs == 0:
        print("No gRNAs found in the input sequence", subjectname)
        all_gRNAs = []
    elif findPairedgRNAOnly == True and n_minus_gRNAs * n_plus_gRNAs == 0:
        print("No gRNAs found in the input sequence", subjectname)
        all_gRNAs = []
    elif annotatePaired == True or findPairedgRNAOnly == True:
        temp = pd.DataFrame(columns = (["ReversegRNAPlusPAM", "ReversegRNAName", "ForwardgRNAPlusPAM", "ForwardgRNAName", "gap"]))
        if len(subjects) > 1:
            toAppend = True
            colNames = False
        
        if n_minus_gRNAs > 0 and n_plus_gRNAs > 0:
            plus_start = list(plus_gRNAs["starts_gRNA"])
            minus_end = list(minus_gRNAs["starts_gRNA"])
            plus_end = []
            minus_start = []
            for i in range(len(plus_start)):
                plus_end.append(plus_start[i] + gRNA_size + PAM_size - 1)
            for i in range(len(minus_end)):
                minus_start.append(minus_end[i] - gRNA_size - PAM_size + 1)
            if paired_orientation == "PAMout":
                pair_index = compute_pair_index3(plus_start, minus_end, min_gap, max_gap)
                plus_index = pair_index[0]
                minus_index = pair_index[1]
                gap = pair_index[2]
            else:
                pair_index = compute_pair_index3(minus_start, plus_end, min_gap, max_gap)
                plus_index = pair_index[1]
                minus_index = pair_index[0]
                gap = pair_index[2]
            #print("plus_index:", plus_index, "minus_index:", minus_index, "gap:", gap);
            #print("minus_start:", minus_start, "plus_start:", plus_start, "minus_end:", minus_end, "plus_end:", plus_end);
          
            if findPairedgRNAOnly != True:
                all_gRNAs = bio.DNAStringSet(list(plus_gRNAs["seq"])+list(minus_gRNAs["seq"])).reset_index(drop=True)   
                all_gRNAs = pd.concat([all_gRNAs, pd.DataFrame(list(plus_gRNAs["name"])+list(minus_gRNAs["name"]), columns=(["names"]))], 1)
                forEffi = pd.concat([plus_gRNAs, minus_gRNAs], 0).reset_index(drop=True)
            
            if len(minus_index) > 0 and len(plus_index) > 0:
                paired = pd.concat([pd.DataFrame(minus_gRNAs["seq"][minus_index]).rename(columns={"seq":"ReversegRNAPlusPAM"}).reset_index(drop=True), pd.DataFrame(minus_gRNAs["name"][minus_index]).rename(columns={"name":"ReversegRNAName"}).reset_index(drop=True), pd.DataFrame(plus_gRNAs["seq"][plus_index]).rename(columns={"seq":"ForwardgRNAPlusPAM"}).reset_index(drop=True), pd.DataFrame(plus_gRNAs["name"][plus_index]).rename(columns={"name":"ForwardgRNAName"}).reset_index(drop=True), pd.DataFrame(gap, columns = (["gap"]))], 1)
                
                if findPairedgRNAOnly == True:
                    plus_index = np.unique(plus_index)
                    minus_index = np.unique(minus_index)
                    all_gRNAs = bio.DNAStringSet(list(plus_gRNAs["seq"][plus_index])+list(minus_gRNAs["seq"][minus_index])).reset_index(drop=True)   
                    all_gRNAs = pd.concat([all_gRNAs, pd.DataFrame(list(plus_gRNAs["name"][plus_index])+list(minus_gRNAs["name"][minus_index]), columns=(["names"]))], 1)
                    forEffi = pd.concat([plus_gRNAs.loc[plus_index], minus_gRNAs.loc[minus_index]], 0).reset_index(drop=True)
                
                if primeEditing == True:
                    if target_start != None or target_end != None or targeted_seq_length_change != None or corrected_seq != None:
                        if type(target_start) == int and type(target_end) == int and type(targeted_seq_length_change) == int:
                            paired = designPEs.designPEs(subject,
                                                PAM_size = PAM_size,
                                                gRNA_size = gRNA_size,
                                                overlap_gRNA_positions = overlap_gRNA_positions,
                                                PBS_length = PBS_length,
                                                paired_gRNAs = pd.DataFrame(paired),
                                                RT_template_length = RT_template_length,
                                                RT_template_pattern = RT_template_pattern,
                                                corrected_seq = corrected_seq,
                                                targeted_seq_length_change = targeted_seq_length_change,
                                                bp_after_target_end = bp_after_target_end,
                                                target_start = target_start,
                                                target_end = target_end,
                                                primeEditingPaired_output =  primeEditingPaired_output, 
                                                col_names = colNames, appendd = toAppend)
                        else:
                            if type(target_start) != int:
                                raise Exception("target_start != int")
                            elif target_end != inr:
                                raise Exception("target_end != int")
                            else:
                                raise Exception("targeted_seq_length_change != int")
                            
                    else:
                        if target_start == None:
                            raise Exception("target_start = None")
                        elif target_end == None:
                            raise Exception("target_end = None")
                        elif targeted_seq_length_change == None:
                            raise Exception("targeted_seq_length_change = None")
                        else:
                            raise Exception("corrected_seq = None")
                    
                    paired = paired.iloc[:,0:5]
                    all_gRNAs = pd.concat([bio.DNAStringSet(np.unique(list(paired.iloc[:,2]) + list(paired.iloc[:,0]))).reset_index(drop=True), pd.DataFrame(np.unique(list(paired.iloc[:,3]) + list(paired.iloc[:,1])), columns=(["name"])).reset_index(drop=True)], 1)
                    forEffi0 = []
                    
                    for i in range(len(forEffi)):
                        if forEffi.iloc[:,1][i] in list(all_gRNAs["name"]):
                            forEffi0.append(i)
                    forEffi = forEffi.iloc[forEffi0,:].reset_index(drop=True)
                    
                if np.shape(paired)[0] == 1:
                    paired.to_csv(pairOutputFile, index = False) 
                else:
                   
                    paired.sort_values(by=["ForwardgRNAName"]).to_csv(pairOutputFile, sep = "\t", index = False)
                        
                ## if paired found
                
            elif findPairedgRNAOnly == True:
                print("No paired gRNAs found for sequence", subjectname)
                all_gRNAs = []
                forEffi = ""
                
                ### if plus_gRNAs and minus_gRNAs not empty
                
        elif findPairedgRNAOnly != True and n_plus_gRNAs > 0:
            all_gRNAs = pd.concat([bio.DNAStringSet(list(plus_gRNAs.iloc[:,0])).reset_index(drop=True), pd.DataFrame(list(plus_gRNAs.iloc[:,1]), columns=(["name"]))], 1)
            forEffi = plus_gRNAs
        elif findPairedgRNAOnly != True and n_minus_gRNAs > 0:
            all_gRNAs = pd.concat([bio.DNAStringSet(list(minus_gRNAs.iloc[:,0])).reset_index(drop=True), pd.DataFrame(list(minus_gRNAs.iloc[:,1]), columns=(["name"]))], 1)
            forEffi = minus.gRNAs
            
            ### annotatePaired and (no paired found or findPairedOnly)
        
    elif n_minus_gRNAs > 0 and n_plus_gRNAs > 0:

        all_gRNAs = pd.concat([bio.DNAStringSet(list(plus_gRNAs.iloc[:,0]) + list(minus_gRNAs.iloc[:,0])).reset_index(drop=True), pd.DataFrame(list(plus_gRNAs.iloc[:,1]) + list(minus_gRNAs.iloc[:,1]), columns=(["name"]))], 1)
        forEffi = pd.concat([plus_gRNAs, minus_gRNAs], 0).reset_index(drop=True)
       
    elif n_minus_gRNAs > 0:
        all_gRNAs = pd.concat([bio.DNAStringSet(list(minus_gRNAs.iloc[:,0])).reset_index(drop=True), pd.DataFrame(list(minus_gRNAs.iloc[:,1]), columns=(["name"]))] ,1)
        forEffi = minus_gRNAs
    elif n_plus_gRNAs > 0:
        all_gRNAs = pd.concat([bio.DNAStringSet(list(plus_gRNAs.iloc[:,0])).reset_index(drop=True), pd.DataFrame(list(plus_gRNAs.iloc[:,1]), columns=(["name"]))] ,1)
        forEffi = plus_gRNAs
    else:
        all_gRNAs = []
        
        #if len(all_gRNAs) == 0 and findPairedgRNAOnly != True:
            #print("No gRNAs found in the input sequence", subjectname))
   
    
    if len(all_gRNAs) > 0:
        all_gRNAs_df = forEffi
    else:
        return  pd.DataFrame(all_gRNAs, columns=(["names"]),)
    if calculategRNAEfficacy == True and np.shape(all_gRNAs_df)[0]*np.shape(all_gRNAs_df)[1] > 4:
        
        featureWeightMatrix = featureWeightMatrixFile
        if rule_set == "Root_RuleSet1_2014":
            
            effi = calculategRNAEfficiency(all_gRNAs_df.iloc[:,4], 
                                      baseBeforegRNA = baseBeforegRNA, 
                                      featureWeightMatrix = featureWeightMatrix, 
                                      enable_multicore = enable_multicore,
                                      n_cores_max = n_cores_max,
                                      gRNA_size = gRNA_size)
        elif rule_set == "Root_RuleSet2_2016":
            effi = calculategRNAEfficiency2(all_gRNAs_df.iloc[:,4])
        elif rule_set == "CRISPRscan":
            effi = calculategRNAEfficiencyCRISPRscan(all_gRNAs_df.iloc[:,4], featureWeightMatrix = featureWeightMatrix)
        
        extendedSequences = pd.concat([all_gRNAs_df, pd.DataFrame(effi, columns=(["gRNAefficacy"]))], 1)
        extendedSequences = extendedSequences.rename(columns={"seq": "gRNAplusPAM", "starts_gRNA": "start"})
       
        if PAM_location == "3prime":
            for i in range(len(extendedSequences)):
                if len(extendedSequences.iloc[i,4]) < baseBeforegRNA + gRNA_size + PAM_size + baseAfterPAM:
                    
                    extendedSequences.iloc[i,5] = "extended sequence too short"
        else:
            for i in range(len(extendedSequences)):
                if len(extendedSequences.iloc[i,4]) < baseBeforegRNA + baseAfterPAM:
                    extendedSequences.iloc[i,5] = "extended sequence too short"
        extendedSequences.to_csv(efficacyFile, sep = "\t", index = False)
    
    all_gRNAs = bio.DNAStringSet(list(all_gRNAs_df.iloc[:,0])).reset_index(drop=True)
    all_gRNAs = pd.concat([all_gRNAs, pd.DataFrame(list(all_gRNAs_df.iloc[:,1]), columns=(["names"]))] ,1)
  
    return all_gRNAs
             
