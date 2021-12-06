#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 10:55:10 2021

@author: klaus
"""



import os
import pandas as pd
import numpy as np
import numeric as nm
import annotateOffTargets
import py2bit
from Bio.Seq import Seq
import calculategRNAEfficiency
from calculategRNAEfficiencyCRISPRscan import *


def filterOffTarget(scores, BSgenomeName, outputDir, txdb, orgAnn, chrom_acc = None, 
                    min_score = 0.01, topN = 200, topN_OfftargetTotalScore = 20,
                    annotateExon = True,  ignore_strand = True, oneFilePergRNA = False,
                    fetchSequence = True, upstream = 200, downstream = 200,
                    baseBeforegRNA = 4, baseAfterPAM = 3, gRNA_size = 20, PAM_location = "3prime", PAM_size = 3,
                    featureWeightMatrixFile = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "DoenchNBT2014.csv"])),
                    rule_set = "Root_RuleSet1_2014", calculategRNAefficacyForOfftargets = True):
    
    
    if fetchSequence == True and BSgenomeName == None:
        raise Exception("To fetch sequences, BSgenomeName is required as BSgenome object!")
    
    scores = scores[scores["score"] >= min_score]
   
    scores1 = pd.DataFrame()
    savethis_score = pd.DataFrame()
    for i in range(len(scores.columns)):
        e = []
        if "IsMismatch.pos" not in scores.columns[i]:
            e.append(scores.columns[i])
            scores1 = pd.concat([scores1, scores[e]],axis=1)
    scores = scores1
    OfftargetFile = outputDir + "OfftargetAnalysis.xls"
    OfftargetSummary = outputDir + "Summary.xls"
    gRNAsPlusPAM = scores["name"].drop_duplicates().reset_index(drop=True)
    names = list(gRNAsPlusPAM)
    top5OfftargetTotalScore = nm.numeric(len(names))
    topNOfftargetTotalScore = top5OfftargetTotalScore
    temp = pd.concat([pd.DataFrame(names, columns=(["names"])), pd.DataFrame(gRNAsPlusPAM), pd.DataFrame(top5OfftargetTotalScore, columns=(["top5OfftargetTotalScore"])), pd.DataFrame(topNOfftargetTotalScore, columns=(["topNOfftargetTotalScore"]))],1).rename(columns={"name":"gRNAsPlusPAM"})
    mismatch_distance2PAM = np.array([ [None] * 11] * len(names))
    appendd = False
  
    for i in range(len(gRNAsPlusPAM)):
        this_score = scores[scores["name"] == list(gRNAsPlusPAM)[i]]
       
        this_score = this_score.sort_values(by=["score", "n.mismatch"], ascending=(False)).reset_index(drop=True)
        maxN = min(topN+1, np.shape(this_score)[0])
        this_score = this_score[0:maxN]
        maxN_totalScore = min(maxN, (topN_OfftargetTotalScore + 1))
        if this_score["n.mismatch"].values[0] == 0 and this_score["NGG"].values[0] == 1:
            start_ind = 2
            end_ind = min(maxN, 6)
            end_forSummary = 11
        else:
            start_ind = 1
            maxN = maxN - 1
            maxN_totalScore = maxN_totalScore - 1
            end_forSummary = 10
            end_ind = min(maxN, 5)

        temp.iloc[i, 2] = sum(this_score["score"][start_ind-1:end_ind])
        
        #if maxN < 6:
        #    temp.iloc[i,3] = sum(this_score["score"][2-1:maxN])
        #else:
        #    temp.iloc[i,3] = sum(this_score["score"][2-1:6])
        if maxN < maxN_totalScore:
            temp.iloc[i, 3] = sum(this_score["score"][start_ind-1:maxN])
        else:
            temp.iloc[i, 3] = sum(this_score["score"][start_ind-1:maxN_totalScore])
        temp.iloc[i, 1] = this_score["gRNAPlusPAM"].drop_duplicates()
        if this_score["mismatch.distance2PAM"].iloc[0] == "" or this_score["mismatch.distance2PAM"].iloc[0] == ["-", "-", "-","-"]:
            for y in range(np.shape(mismatch_distance2PAM)[1]):
                mismatch_distance2PAM[i][y] = "NMM"
        else:
            for y in range(np.shape(mismatch_distance2PAM)[1]):
                mismatch_distance2PAM[i][y] = "perfect match not found"
     
        forSummary = this_score.iloc[start_ind-1:end_forSummary,:]
        forSummary = forSummary.sort_values(by=["score"], ascending=False).reset_index(drop=True)
       
        for y in range(1, 11):
            if forSummary.empty == True:
                mismatch_distance2PAM[i][y] = "NA"
            elif np.shape(forSummary)[0] == 10:
                for y in range(1, 11):
                    mismatch_distance2PAM[i][y] = forSummary["mismatch.distance2PAM"][y-1]
            else:
                for y in range(1, np.shape(forSummary)[0]+1):
                    mismatch_distance2PAM[i][y] = forSummary["mismatch.distance2PAM"][y-1]
                    r = y + 1
        for y in range(1, 11):
            if  mismatch_distance2PAM[i][y] == "NMM":
                mismatch_distance2PAM[i][y] = ""
     
        this_score = pd.concat([this_score["name"], this_score["gRNAPlusPAM"],this_score["OffTargetSequence"],
                          this_score["score"], this_score["n.mismatch"], this_score["mismatch.distance2PAM"], 
                          this_score["alignment"],this_score["NGG"], this_score["forViewInUCSC"], 
                          this_score["strand"], this_score["chrom"], this_score["chromStart"],
                          this_score["chromEnd"]], 1)
        
        if oneFilePergRNA == True and np.shape(this_score)[0] > 0:
            this_score.to_csv(outputDir + "OfftargetAnalysis-" + str(temp.iloc[i][0]) + ".xls", sep = "\t", index = False)
        if i == 0 and np.shape(this_score)[0] > 0:
            savethis_score = pd.concat([savethis_score, this_score], 0).reset_index(drop=True)
            appendd = True
        elif np.shape(this_score)[0] > 0:
            savethis_score = pd.concat([savethis_score, this_score], 0).reset_index(drop=True)
            appendd = True
    savethis_score.to_csv(OfftargetFile, sep ="\t", index = False)
    
    for p in range(len(temp)):
        if temp.loc[p, ("top5OfftargetTotalScore")] == 0 :
            temp.loc[p, ("top5OfftargetTotalScore")] = "NA"
        if temp.loc[p, ("topNOfftargetTotalScore")] == 0:
            temp.loc[p, ("topNOfftargetTotalScore")] = "NA"
 
    
    temp = pd.concat([temp, pd.DataFrame(mismatch_distance2PAM)], 1)
    temp = temp.rename(columns ={temp.columns[4]: "top1Hit(onTarget)MMdistance2PAM", temp.columns[3]: "top" + str(topN_OfftargetTotalScore) + "OfftargetTotalScore"})
    for i in range(5,len(temp.columns)):
        temp = temp.rename(columns ={temp.columns[i]: "topOfftarget" + str(i-4) + "MMdistance2PAM"}) 
    
    Offtargets = pd.read_csv(OfftargetFile, sep = "\t", header = 0, keep_default_na=False)
  
    if annotateExon == True:
        txdb = os.getcwd() + os.sep + os.sep.join(["extdata", txdb, txdb + ".sqlite"])
        orgAnn = os.getcwd() + os.sep + os.sep.join(["extdata", orgAnn, orgAnn +".sqlite"])
        Offtargets = annotateOffTargets.annotateOffTargets(Offtargets, txdb, orgAnn, ignore_strand)
    ontargets = Offtargets[Offtargets["n.mismatch"] == 0]
    
    if calculategRNAefficacyForOfftargets != True and np.shape(ontargets)[0] > 0:
        chr = list(ontargets["chrom"].astype(str))
        strand = list(ontargets["strand"])
        Start = []
        End = []
        for i in range(len(strand)):
            if strand[i] == "-":
                Start.append(ontargets.loc[i, ("chromStart")] - baseAfterPAM)
                End.append(ontargets.loc[i, ("chromEnd")] + baseBeforegRNA)
            else:
                Start.append(ontargets.loc[i, ("chromStart")] - baseBeforegRNA)
                End.append(ontargets.loc[i, ("chromEnd")] + baseAfterPAM)
    
    elif calculategRNAefficacyForOfftargets == True and np.shape(Offtargets)[0] > 0:
        chr = list(Offtargets["chrom"].astype(str))
        strand = list(Offtargets["strand"])
        Start = []
        End = []
        if PAM_location == "3prime":
            for i in range(len(strand)):
                if strand[i] == "-":
                    Start.append(Offtargets.loc[i, ("chromStart")] - baseAfterPAM)
                    End.append(Offtargets.loc[i, ("chromEnd")] + baseBeforegRNA)
                else:
                    Start.append(Offtargets.loc[i, ("chromStart")] - baseBeforegRNA)
                    End.append(Offtargets.loc[i, ("chromEnd")] + baseAfterPAM)
        else:
            for i in range(len(strand)):
                if strand[i] == "-":
                    Start.append(Offtargets.loc[i, ("chromStart")] - baseAfterPAM + gRNA_size)
                    End.append(Offtargets.loc[i, ("chromEnd")] + baseBeforegRNA - PAM_size)
                else:
                    Start.append(Offtargets.loc[i, ("chromStart")] - baseBeforegRNA + PAM_size)
                    End.append(Offtargets.loc[i, ("chromEnd")] + baseAfterPAM - gRNA_size)
    if (calculategRNAefficacyForOfftargets == True and np.shape(Offtargets)[0] > 0) or (calculategRNAefficacyForOfftargets != True and np.shape(ontargets)[0] > 0):
        starts = []
        ends = []
        extendedSequence = []
      
        for i in range(len(Start)):
            starts.append(max(Start[i], 0))
            ends.append(min(End[i], py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms(chr[i])))
            if strand[i] == "-":
                extendedSequence.append(str(Seq(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i])[::-1]).reverse_complement())[starts[i]-1:ends[i]][::-1])
            else:
                extendedSequence.append(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i])[starts[i]-1:ends[i]])
        
    
        if rule_set == "Root_RuleSet1_2014":
            gRNAefficiency = calculategRNAEfficiency.calculategRNAEfficiency(extendedSequence, baseBeforegRNA = baseBeforegRNA, featureWeightMatrix = featureWeightMatrixFile)
        elif rule_set == "Root_RuleSet2_2016":
            gRNAefficiency = calculategRNAEfficiency2(extendedSequence)
        elif rule_set == "CRISPRscan":
            gRNAefficiency = calculategRNAEfficiencyCRISPRscan(extendedSequence, featureWeightMatrix = featureWeightMatrixFile)
        if calculategRNAefficacyForOfftargets != True and np.shape(ontargets)[0] > 0:
            ontargets = pd.concat([ontargets,  pd.DataFrame(extendedSequence, columns=(["extendedSequence"])), pd.DataFrame(gRNAefficiency, columns=(["gRNAefficacy"]))], 1)
            Offtargets = Offtargets.merge(ontargets, how = "inner", sort=True)
        else:
            Offtargets  = pd.concat([Offtargets,  pd.DataFrame(extendedSequence, columns=(["extendedSequence"])), pd.DataFrame(gRNAefficiency, columns=(["gRNAefficacy"]))], 1)
    
    if fetchSequence == True:
        strand = list(Offtargets["strand"])
        chr = list(Offtargets["chrom"])
        Start = []
        End = []
        for i in range(len(strand)):
            if strand[i] == "-":
                Start.append(Offtargets.loc[i, ("chromStart")] - downstream)
                End.append(Offtargets.loc[i, ("chromEnd")] + upstream)
            else:
                Start.append(Offtargets.loc[i, ("chromStart")] - upstream)
                End.append(Offtargets.loc[i, ("chromEnd")] + downstream)
        starts = []
        ends = []
        seq = []
        for i in range(len(Start)):
            starts.append(max(Start[i], 0))
            ends.append(min(End[i], py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms(chr[i])))
            if strand[i] == "-":
                seq.append(str(Seq(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i])[::-1]).reverse_complement())[starts[i]-1:ends[i]][::-1])
            else:
                seq.append(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i])[starts[i]-1:ends[i]])
        Offtargets = pd.concat([Offtargets, pd.DataFrame(seq, columns=(["flankSequence"]))], 1)
    if "NGG" in Offtargets.columns:
        Offtargets = Offtargets.rename(columns={"NGG": "isCanonicalPAM"})
    
   
    temp.to_csv(OfftargetSummary, index = False)
    Offtargets.sort_values(by = ["name", "score", "OffTargetSequence"]).to_csv(OfftargetFile, index = False)
    offtargets = Offtargets.drop_duplicates()
    summary = temp
   
    return offtargets, summary
      
