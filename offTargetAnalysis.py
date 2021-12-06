#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 11:39:35 2021

@author: klaus
"""

from multiprocessing import Pool
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import re
import os
import numpy as np
import math
import translatePattern as tp
import gregexprindex as gr
import biostrings as bio
import getSeqFromBed as gSFB
import findgRNAs as f
import filtergRNAs as filtergRNAs
import searchHits2
import writeHits2 as wrtHts2
import buildFeatureVectorForScoring as bFVFS
import getOfftargetScore
import getOfftargetScore2
import filterOffTarget
import uniqueREs



pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", None)
pd.set_option('expand_frame_repr', True)


 

def offTargetAnalysis(inputFilePath, outputDir, BSgenomeName=None, txdb=None, orgAnn=None, chrom_acc = None, 
                      corrected_seq = None, target_start = None, target_end = None, 
                      targeted_seq_length_change = None, gRNAoutputName = None, format = "fasta", 
                      header = False, findgRNAs = True, exportAllgRNAs = "all", 
                      findgRNAsWithREcutOnly = False, annotatePaired = True, 
                      REpatternFile = os.getcwd() + os.sep + os.sep.join(["extdata", "NEBenzymes.fa"]),
                      minREpatternSize = 4, overlap_gRNA_positions = [17, 18], findPairedgRNAOnly = False, 
                      paired_orientation = "PAMout", enable_multicore = False, n_cores_max = 6,
                      min_gap = 0, max_gap = 20, gRNA_name_prefix = "", PAM_size = 3, gRNA_size = 20,
                      PAM = "NGG", chromToSearch = ["all"],
                      chromToExclude = ["chr17_ctg5_hap1", "chr4_ctg9_hap1", "chr6_apd_hap1",
                                        "chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4",
                                        "chr6_mcf_hap5", "chr6_qbl_hap6","chr6_ssto_hap7"],
                      max_mismatch = 3, PAM_pattern = "NNG$", allowed_mismatch_PAM = 1, gRNA_pattern = "",
                      baseEditing = False, targetBase = "C", editingWindow = [4, 8],
                      editingWindow_offtargets = [4, 8], primeEditing = False, PBS_length = 13,
                      RT_template_length = [8, 28],
                      RT_template_pattern = "D$", bp_after_target_end = 15,
                      primeEditingPaired_output = "pairedgRNAsForPE.xls", min_score = 0, topN = 1000,
                      topN_OfftargetTotalScore = 10, annotateExon = True, ignore_strand = True,
                      fetchSequence = True, upstream = 200, downstream = 200, upstream_search = 0,
                      downstream_search = 0, weights = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
                                                        0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 
                                                        0.804, 0.685, 0.583], 
                      baseBeforegRNA = 4, baseAfterPAM = 3, featureWeightMatrixFile = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "DoenchNBT2014.csv"])), 
                      useScore = True, useEfficacyFromInputSeq = False, outputUniqueREs = True, 
                      gRNA_backbone = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
                      temperature = 37, overwrite = False, scoring_method = "Hsu-Zhang",
                      subPAM_activity = {"AA" : 0,
                                         "AC" : 0,
                                         "AG" : 0.259259259,
                                         "AT" : 0,
                                         "CA" : 0,
                                         "CC" : 0,
                                         "CG" : 0.107142857,
                                         "CT" : 0,
                                         "GA" : 0.069444444,
                                         "GC" : 0.022222222,
                                         "GG" : 1,
                                         "GT" : 0.016129032,
                                         "TA" : 0,
                                         "TC" : 0,
                                         "TG" : 0.038961039,
                                         "TT" : 0},
                      subPAM_position = [22, 23], PAM_location = "3prime", rule_set = "Root_RuleSet1_2014",
                      calculategRNAefficacyForOfftargets = True,
                      mismatch_activity_file = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "NatureBiot2016SuppTable19DoenchRoot.csv"])),
                      method_indelFreq = "Lindel", baseBeforegRNA_indelFreq = 13, baseAfterPAM_indelFreq = 24):
    
    if type(inputFilePath) == str:
        if  os.path.isfile(inputFilePath) == True:
            inputFilePath = pd.read_csv(inputFilePath, sep = ">", header = None)
            inputFilePathseq = inputFilePath[0].dropna().reset_index(drop=True)
            inputFilePathname = inputFilePath[1].dropna().reset_index(drop=True)
            inputFilePath = pd.concat([inputFilePathseq, inputFilePathname],1).rename(columns={0:"seq", 1:"names"})
        else:
            inputFilePath = inputFilePath.upper()
            if inputFilePath.count("A") + inputFilePath.count("T") + inputFilePath.count("G") + inputFilePath.count("C") == len(inputFilePath):
               inputFilePath = pd.concat([pd.DataFrame([inputFilePath], columns=(["seq"])), pd.DataFrame(["NA"], columns = (["names"]))], 1)
            else:
                raise Exception("inputFilePath has to be a DNA string of characters!")
    else:
        raise Exception("inputFilePath has to be a DNA string or string of characters!")
    

    print("Validating input ...\n")
  
    
    if rule_set == "DeepCpf1":
        baseBeforegRNA = 8
        baseAfterPAM = 26
        if scoring_method == "CFDscore" and subPAM_activity["TT"] < 1:
            subPAM_activity = { "AA" : 0,
                                "AC" : 0,
                                "AG" : 0,
                                "AT" : 0.1,
                                "CA" : 0,
                                "CC" : 0,
                                "CG" : 0,
                                "CT" : 0.05,
                                "GA" : 0,
                                "GC" : 0,
                                "GG" : 0,
                                "GT" : 0.05,
                                "TA" : 0.2,
                                "TC" : 0.1,
                                "TG" : 0.1,
                                "TT" : 1}
    elif rule_set in ["Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan"]:
        if PAM_location == "3prime":
            baseBeforegRNA = 4
            baseAfterPAM = 3
        else:
            baseBeforegRNA = 4 + PAM_size
            baseAfterPAM = 3 + gRNA_size
    mismatch_activity = mismatch_activity_file
    required_col = ["Mismatch.Type", "Position", "Percent.Active"]
    
    if scoring_method == "CFDscore":
        mismatch_activity = mismatch_activity_file
        required_col = ["Mismatch.Type", "Position", "Percent.Active"]
        if len(mismatch_activity.columns.intersection(required_col)) != len(required_col):
            raise Exception("Please rename the mismatch activity file column to contain at least these 3 column names: Mismatch.Type, Position, Percent.Active\n")
    elif scoring_method == "Hsu-Zhang":
        if len(weights) !=  gRNA_size:
            raise Exception("Please make sure the size of weights vector equals to the gRNA_size!\n")
    
    if findgRNAsWithREcutOnly == True and findgRNAs == True and os.path.isfile(os.getcwd() + os.sep + os.sep.join(["extdata", "NEBenzymes.fa"])) != True:
        raise Exception("Please specify an REpattern file as fasta file with restriction enzyme recognition sequences!")
    
    pairOutputFile = outputDir + os.sep+ "pairedgRNAs.xls"
    REcutDetailFile = outputDir + os.sep + "REcutDetails.xls"
    bedFile = outputDir + os.sep + "gRNAspycrisprcas9.bed"  
  
    if gRNAoutputName == None:
        gRNAoutputName = "inputseq"
    if format =="bed":
        if os.path.isfile(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName + ".2bit"])) != True:
            raise Exception("BSgenomeName is required as BSgenome object when input file is in bed format!")
        inputFilePath = gSFB.getSeqFromBed(inputFilePath, header = header, BSgenomeName = BSgenomeName)
        #### format for filtergRNAs
        format = "fasta"
    if findgRNAs == True:
        print("Searching for gRNAs ...\n")
        efficacyFile = outputDir + os.sep + "gRNAefficacy.xls"
        if chromToSearch == "" or useEfficacyFromInputSeq == True:
            potential_gRNAs = f.findgRNAs(inputFilePath,
               overlap_gRNA_positions = overlap_gRNA_positions,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
               primeEditing = primeEditing,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired_orientation = paired_orientation,
               pairOutputFile = pairOutputFile, PAM = PAM,
               PAM_location = PAM_location,
               gRNA_pattern = gRNA_pattern, PAM_size = PAM_size,
               gRNA_size = gRNA_size, min_gap = min_gap,
               max_gap = max_gap, name_prefix = gRNA_name_prefix,
               format = format, featureWeightMatrixFile = featureWeightMatrixFile, 
               baseBeforegRNA = baseBeforegRNA, baseAfterPAM = baseAfterPAM,
    	       calculategRNAEfficacy = True, efficacyFile = efficacyFile, corrected_seq = "T",
               targeted_seq_length_change = 0, target_start = 46,  target_end = 46,
               rule_set = rule_set, chrom_acc = chrom_acc)
        else:
            potential_gRNAs = f.findgRNAs(inputFilePath,
               overlap_gRNA_positions = overlap_gRNA_positions,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
               primeEditing = primeEditing,
               PBS_length = PBS_length,
               RT_template_length = RT_template_length,
               RT_template_pattern = RT_template_pattern,
               corrected_seq = corrected_seq,
               targeted_seq_length_change = targeted_seq_length_change,
               bp_after_target_end = bp_after_target_end,
               target_start = target_start,
               target_end = target_end,
               primeEditingPaired_output =  outputDir + primeEditingPaired_output,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired_orientation = paired_orientation,
               enable_multicore = enable_multicore,
               n_cores_max = n_cores_max,
               pairOutputFile = pairOutputFile, PAM = PAM,
	       gRNA_pattern = gRNA_pattern, PAM_size = PAM_size,
               PAM_location = PAM_location,
               gRNA_size = gRNA_size, min_gap = min_gap,
               max_gap = max_gap, name_prefix = gRNA_name_prefix,
               efficacyFile = efficacyFile, featureWeightMatrixFile = featureWeightMatrixFile, 
               format = format,  rule_set = rule_set, chrom_acc = chrom_acc)
        potential_gRNAs["names"] = potential_gRNAs["names"].str.replace(">", "")
        if len(potential_gRNAs) == 0:
            return(print("no gRNAs found!"))
        if len(potential_gRNAs) > 0 and exportAllgRNAs == "fasta" or exportAllgRNAs == "all":
            saveFasta = open(outputDir + "allgRNAs.fa", "w+")
            for i in range(len(potential_gRNAs)):
                if ">" in potential_gRNAs["names"][i]:
                    saveFasta.write(potential_gRNAs["names"][i] +"\n")
                else:
                    saveFasta.write(">" + potential_gRNAs["names"][i] +"\n")
                saveFasta.write(potential_gRNAs["seq"][i]+"\n")
            saveFasta.close()
        
        if len(potential_gRNAs) > 0 and exportAllgRNAs == "genbank" or exportAllgRNAs == "all":
            if type(inputFilePath) == pd.core.frame.DataFrame:
                subjects = []
                subjects.append(inputFilePath.iloc[0][0].upper())
                subjects = bio.DNAStringSet(subjects)
            else:
                subjects = bio.DNAStringSet(inputFilePath.upper())
    
            subjects["names"] = inputFilePath["names"]
            locuses = list(subjects["names"])
            names_gRNA = list(potential_gRNAs["names"])
           
            for i in range(len(locuses)):
                thisLocus = re.sub("'", "", locuses[i])
                thisLocus = re.sub(" ", "", thisLocus)
                thisLocus = re.sub(">", "", thisLocus)
                thisSeq = subjects.iloc[i][1].lower()
                n_bp = len(thisSeq)
                temp = []
                
                for j in range(len(names_gRNA)):
                    temp.append(re.sub(">", "", names_gRNA[j]).split(thisLocus + "_gR"))
                
                locus = "LOCUS       " + thisLocus + "                     " + str(n_bp) + " bp    dna     linear   UNK"
                definition = "DEFINITION  pycrisprcas9 output for " + gRNAoutputName +" sequence"
                accession = "ACCESSION   unknown"
                features = "FEATURES             Location/Qualifiers"
                header = pd.DataFrame([locus, definition, accession, features], index=(["locus", "definition", "accession", "features"]))
                found_gRNA = 0
                
                for j in range(len(temp)):
                    if len(temp[j]) > 1:
                        found_gRNA = found_gRNA + 1
                        if found_gRNA == 1:
                            thisFile = thisLocus + ".gbk"
                            savegbk = open(outputDir + thisFile, "w+")
                            for k in range(len(header)):
                                savegbk.write(header[0][k]+"\n")
                            savegbk.close()
                        temp1 = []
                        for l in range(len(temp[j])):
                            if  len(gr.gregexpr_index("f", temp[j][l])) > 0:   
                                temp1.append(re.sub("f", "", temp[j][l]))
                                isForward = True
                            else:
                                temp1.append(re.sub("r", "", temp[j][l]))
                                isForward = False
                        feature = temp1[0]
                        location = temp1[1]
                        if isForward ==True:
                            Start = location
                            End = int(Start) + max(overlap_gRNA_positions) -  min(overlap_gRNA_positions)
                            savegbk = open(outputDir + thisFile, "a+")
                            savegbk.write("     misc_bind       " + str(Start) + ".." + str(End) + "\n")
                            savegbk.write("                     /note=\"gRNAf" + str(feature) + "\"" + "\n")
                            savegbk.close()
                        else:
                            End = location
                            Start = int(End) - max(overlap_gRNA_positions) + min(overlap_gRNA_positions)
                            savegbk = open(outputDir + thisFile, "a+")
                            savegbk.write("     misc_bind       complement(" + str(Start) + ".." + End + ")" +"\n")
                            savegbk.write("                     /note=\"gRNAr" + feature + "\"" + "\n")
                            savegbk.close()
                if found_gRNA > 0:
                    savegbk = open(outputDir + thisFile, "a+")
                    savegbk.write("ORIGIN" + "\n")
                    savegbk.close()
                    seq_lines = math.floor(len(thisSeq) / 60) + 1
                    for k in range(1, seq_lines+1):
                        line_start = (k - 1) * 60 + 1
                        line_end = min(line_start + 59, len(thisSeq))
                        n_leading_spaces = 9 - len(str(line_start))
                        leading_spaces = " " * n_leading_spaces
                        seq_thisLine = thisSeq[line_start-1:line_end]
                        len_thisLine = len(seq_thisLine)
                        n_seg = math.floor(len_thisLine /10) + 1
                        for m in range(1, n_seg+1):
                             seg_start = (m -1) * 10 + 1
                             seg_end = min(seg_start + 9, len_thisLine)
                             if m == 1:
                                 seq_thisLine_formatted = seq_thisLine[seg_start-1:seg_end]
                             else:
                                 seq_thisLine_formatted = seq_thisLine_formatted + " " + seq_thisLine[seg_start-1:seg_end]
                        savegbk = open(outputDir + thisFile, "a+")
                        savegbk.write(leading_spaces + str(line_start) + " " + seq_thisLine_formatted +"\n")
                        savegbk.close()
                    
                    savegbk = open(outputDir + thisFile, "a+")
                    savegbk.write("//" + "\n")
                    savegbk.close()
       
        if findPairedgRNAOnly == True and len(potential_gRNAs) >0:
            gRNAs_RE = filtergRNAs.filtergRNAs(potential_gRNAs, pairOutputFile = pairOutputFile, 
                findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, REpatternFile = REpatternFile, 
                format = format,  minREpatternSize = minREpatternSize,
                overlap_gRNA_positions = overlap_gRNA_positions)
            REcutDetails = gRNAs_RE[1]
            REcutDetails.sort_values('ForwardgRNAName').to_csv(REcutDetailFile, index = False)
            
        elif len(potential_gRNAs) >0:
            gRNAs_RE = filtergRNAs.filtergRNAs(potential_gRNAs, findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
                                REpatternFile = REpatternFile, format = format, 
                                minREpatternSize = minREpatternSize,
                                overlap_gRNA_positions = overlap_gRNA_positions)
            REcutDetails = gRNAs_RE[1]
            REcutDetails.sort_values('REcutgRNAName').to_csv(REcutDetailFile, index = False)
      
        if findgRNAsWithREcutOnly == True:
            gRNAs = gRNAs_RE[0]
            
        else:
            gRNAs = potential_gRNAs
     
        if annotatePaired == True or findPairedgRNAOnly == True:
            pairedInformation = pd.read_csv(pairOutputFile, sep="\t")
        
    else:
        if type(inputFilePath) != pd.core.frame.DataFrame:
            if os.path.isfile(outputDir + inputFilePath) != True:
                raise Exception("inputfile specified as ", inputFilePath, " does not exists!")
            if format == "fasta" or format == "fastq":
                potential_gRNAs = pd.read_csv(inputFilePath)
            else:
                raise Exception("format needs to be either fasta,fastq or bed!")
        else:
            potential_gRNAs = inputFilePath
            if len(potential_gRNAs.columns) == 0:
                for i in range(len(potential_gRNAs)):
                    potential_gRNAs = potential_gRNAs.rename(columns={potential_gRNAs.columns[i]:"gRNAs" + str(i+1)})
        gRNAs_RE = filtergRNAs.filtergRNAs(potential_gRNAs,  REpatternFile = REpatternFile, format = format, 
                                  minREpatternSize = minREpatternSize, 
                                  overlap_gRNA_positions = overlap_gRNA_positions)
        REcutDetails = gRNAs_RE[1]
        REcutDetails.sort_values("REcutgRNAName").to_csv(REcutDetailFile, index = False)
       
        if findgRNAsWithREcutOnly == True:
            gRNAs = gRNAs_RE[0]
        else:
            gRNAs = potential_gRNAs
        pairedInformation = ""
    
    if len(chromToSearch) == 1 and chromToSearch == [""]:
        print("Done. Please check output files in directory ", outputDir, "\n")
        return gRNAs
    if len(gRNAs) != 0:
        for i in range(len(gRNAs)):
            gRNAs.loc[i, ("names")] = re.sub("\t", "", gRNAs.loc[i, ("names")])
            gRNAs.loc[i, ("names")] = re.sub("\n", "", gRNAs.loc[i, ("names")])
            gRNAs.loc[i, ("names")] = re.sub(" ", "", gRNAs.loc[i, ("names")])
  
    hits = searchHits2.searchHits2(gRNAs = gRNAs, PAM = PAM, PAM_pattern = PAM_pattern, outputDir = outputDir,
                                   BSgenomeName = BSgenomeName, chromToSearch = chromToSearch,
                                   chromToExclude = chromToExclude, max_mismatch = max_mismatch,
                                   PAM_size = PAM_size, gRNA_size = gRNA_size, allowed_mismatch_PAM = allowed_mismatch_PAM,
                                   PAM_location = PAM_location, baseEditing = baseEditing, targetBase = targetBase,
                                   editingWindow = editingWindow_offtargets)

    if np.shape(hits)[0] > 0:
        
        print("Building feature vectors for scoring ...\n")
        featureVectors = bFVFS.buildFeatureVectorForScoring(hits = hits, canonical_PAM = PAM, gRNA_size = gRNA_size, 
            subPAM_position = subPAM_position, PAM_location = PAM_location, PAM_size = PAM_size)
        print("Calculating scores ...\n")
        
        if scoring_method == "CFDscore":
            scores = getOfftargetScore2.getOfftargetScore2(featureVectors, subPAM_activity = subPAM_activity, mismatch_activity_file = mismatch_activity_file)
        else:
            scores = getOfftargetScore.getOfftargetScore(featureVectors, weights = weights)
       
        print("Annotating, filtering and generating reports ...\n")
        offTargets = filterOffTarget.filterOffTarget(scores = scores, outputDir = outputDir,
            BSgenomeName = BSgenomeName, fetchSequence = fetchSequence, txdb = txdb,
            orgAnn = orgAnn, ignore_strand = ignore_strand, min_score = min_score, topN = topN, 
            topN_OfftargetTotalScore = topN_OfftargetTotalScore, upstream = upstream, 
            downstream = downstream, annotateExon = annotateExon, baseBeforegRNA = baseBeforegRNA, 
            baseAfterPAM = baseAfterPAM, gRNA_size = gRNA_size, PAM_location = PAM_location, 
            PAM_size = PAM_size, featureWeightMatrixFile = featureWeightMatrixFile, rule_set = rule_set, 
            chrom_acc = chrom_acc, calculategRNAefficacyForOfftargets = calculategRNAefficacyForOfftargets)
        print("Done annotating\n")
        
        summary = pd.read_csv(outputDir + "Summary.xls", keep_default_na=False)
        
        if np.shape(summary)[1] == 1:
            summary1 = []
            for i in range(np.shape(offTargets[1])[1]):
                summary0= []
                for j in range(np.shape(offTargets[1])[0]):
                    summary0.append(offTargets[1].iloc[j][i])
                summary1.append(summary0)
            summary = pd.DataFrame(summary1, index=offTargets[1].columns)

        for i in range(len(summary.columns)):
            if "topOfftarget" in summary.columns[i]:
                y = summary.loc[:, summary.columns[i]]
                for j in range(len(y)):
                    if y[j] == "NA":
                        y[j] = ""
                        summary.loc[j, summary.columns[i]] = y[j]
       
        if findgRNAs == True and (annotatePaired == True or findPairedgRNAOnly == True):
            print("Add paired information...\n")
            PairedgRNAName = []
            for i in range(np.shape(summary)[0]):
                a = ""
                for j in range(len(pairedInformation[pairedInformation["ForwardgRNAName"] == summary.loc[i, ("names")]]["ReversegRNAName"].drop_duplicates().values)):
                    a = a + pairedInformation[pairedInformation["ForwardgRNAName"] == summary.loc[i, ("names")]]["ReversegRNAName"].drop_duplicates().values[j]
                    if len(pairedInformation[pairedInformation["ForwardgRNAName"] == summary.loc[i, ("names")]]["ReversegRNAName"].drop_duplicates().values) > 1:
                        a = a + " "
                for k in range(len(pairedInformation[pairedInformation["ReversegRNAName"] == summary.loc[i, ("names")]]["ForwardgRNAName"].drop_duplicates().values)):
                    a = a + pairedInformation[pairedInformation["ReversegRNAName"] == summary.loc[i, ("names")]]["ForwardgRNAName"].drop_duplicates().values[k]
                    if len(pairedInformation[pairedInformation["ReversegRNAName"] == summary.loc[i, ("names")]]["ForwardgRNAName"].drop_duplicates().values) > 1 :
                        a = a + " "
                PairedgRNAName.append(a)
        print("Add RE information...\n")
        if findPairedgRNAOnly == True and findgRNAs == True:
            REname = []
            for i in range(np.shape(summary)[0]):
                a = ""
                for j in range(len(REcutDetails[REcutDetails["ForwardREcutgRNAName"] == summary.loc[i, ("names")]]["ForwardREname"].drop_duplicates().values)):
                    a = a + REcutDetails[REcutDetails["ForwardREcutgRNAName"] == summary.loc[i, ("names")]]["ForwardREname"].drop_duplicates().values[j]
                    if len(REcutDetails[REcutDetails["ForwardREcutgRNAName"] == summary.loc[i, ("names")]]["ForwardREname"].drop_duplicates().values) > 1 and len(REcutDetails[REcutDetails["ForwardREcutgRNAName"] == summary.loc[i, ("names")]]["ForwardREname"].drop_duplicates().values) != j+1:
                        a = a + " "
                for k in range(len(REcutDetails[REcutDetails["ReverseREcutgRNAName"] == summary.loc[i, ("names")]]["ReverseREname"].drop_duplicates().values)):
                    a = a + REcutDetails[REcutDetails["ReverseREcutgRNAName"] == summary.loc[i, ("names")]]["ReverseREname"].drop_duplicates().values[k]
                    if len(REcutDetails[REcutDetails["ReverseREcutgRNAName"] == summary.loc[i, ("names")]]["ReverseREname"].drop_duplicates().values) > 1 and len(REcutDetails[REcutDetails["ReverseREcutgRNAName"] == summary.loc[i, ("names")]]["ReverseREname"].drop_duplicates().values) != k+1:
                        a = a + " "
                REname.append(a)
            summary = pd.concat([summary, pd.DataFrame(PairedgRNAName, columns=(["PairedgRNAName"])), pd.DataFrame(REname, columns=(["REname"]))], 1)
        else:
            REname = []
            for i in range(np.shape(summary)[0]):
                a = ""
                for j in range(len(REcutDetails[REcutDetails["REcutgRNAName"] == summary.loc[i, ("names")]]["REname"].values)):
                    a = a + re.sub("NA", "", REcutDetails[REcutDetails["REcutgRNAName"] == summary.loc[i, ("names")]]["REname"].values[j])
                    if len(REcutDetails[REcutDetails["REcutgRNAName"] == summary.loc[i, ("names")]]["REname"].values) > 1 and len(REcutDetails[REcutDetails["REcutgRNAName"] == summary.loc[i, ("names")]]["REname"].values) != j+1:
                        a = a + " "   
                REname.append(a)   
            summary = pd.concat([summary, pd.DataFrame(REname, columns=(["REname"]))], 1)
       
        seq = list(summary["gRNAsPlusPAM"])
        print("write gRNAs to bed file...\n")
        on_target = offTargets[0]
     
        #for i in range(len(on_target)):
        on_target = on_target[on_target.loc[:, ("n.mismatch")] == 0].drop_duplicates()
       
        on_target = on_target[on_target.loc[:, ("isCanonicalPAM")] == 1].drop_duplicates().sort_values("name", ascending=False).reset_index(drop=True)
     
        if np.shape(on_target)[0] > 0:
            gRNA_bed = pd.concat([on_target["chrom"], on_target["chromStart"],
                                 on_target["chromEnd"], on_target["name"], 
                                 on_target["gRNAefficacy"] * 1000, 
                                 on_target["strand"], 
                                 on_target["chromStart"], 
                                 on_target["chromEnd"]], 1).drop_duplicates().sort_values("name", ascending=False).reset_index(drop=True)
            if useScore != True:
                gRNA_bed = pd.concat([gRNA_bed, pd.DataFrame(["255,0,0"]*np.shape(gRNA_bed)[0])], 1)
                for i in range(np.shape(gRNA_bed)[0]):
                    if gRNA_bed.loc[i, ("strand")] == "-":
                        gRNA_bed.iloc[i, 8] = "0,255,0"
            for i in range(np.shape(gRNA_bed)[0]):
                gRNA_bed.iloc[i, 1] = gRNA_bed.iloc[i, 1].astype(int) - 1
                gRNA_bed.iloc[i, 2] = gRNA_bed.iloc[i, 2].astype(int)
                if gRNA_bed.iloc[i, 5] == "+":
                    gRNA_bed.iloc[i, 6] = gRNA_bed.iloc[i, 1] + min(overlap_gRNA_positions) - 1
                    gRNA_bed.iloc[i, 7] = gRNA_bed.iloc[i, 1] + max(overlap_gRNA_positions)
                else:
                    gRNA_bed.iloc[i, 6] = gRNA_bed.iloc[i, 2] - max(overlap_gRNA_positions)
                    gRNA_bed.iloc[i, 7] = gRNA_bed.iloc[i, 2] - min(overlap_gRNA_positions) + 1
            savebed = open(bedFile, "w+")
            savebed.write("track name=\"gRNA sites\" description=\"pycrisprcas9\" visibility=2 useScore=1 itemRgb=\"On\""+"\n")
            savebed.close()
            gRNA_bed.to_csv(bedFile, mode = "a+", sep = "\t", header = None, index = False)
            on_target = pd.concat([on_target["name"], on_target["forViewInUCSC"], on_target["extendedSequence"], on_target["gRNAefficacy"]], 1).drop_duplicates().rename(columns={"name": "names"})
            
            if useEfficacyFromInputSeq == True:
                on_target = on_target.iloc[:,:2]
                inputEfficacy = pd.read_csv(efficacyFile, sep="\t")
                inputEfficacy = inputEfficacy[["name","extendedSequence","gRNAefficacy"]]
                on_target = on_target.merge(inputEfficacy.rename(columns={"name":"names"}), how="inner", on="names").sort_values("names", ascending=False).reset_index(drop=True)
          
            if np.shape(on_target)[0] > 0:
                summary = on_target.merge(summary, how="inner", on="names").drop_duplicates().sort_values("names").reset_index(drop=True)
            summary.to_csv(outputDir + "Summary.xls", index = False)
            print("Scan for REsites in flanking region...\n")
            if outputUniqueREs == True:
                REs_isUnique100 = uniqueREs.uniqueREs(REcutDetails = REcutDetails, summary = summary, offTargets=offTargets[0],  BSgenomeName = BSgenomeName, scanUpstream = 100, scanDownstream = 100)
                REs_isUnique50 = uniqueREs.uniqueREs(REcutDetails = REcutDetails, summary = summary, offTargets=offTargets[0],  BSgenomeName = BSgenomeName, scanUpstream = 50, scanDownstream = 50)
                summary = pd.concat([summary, pd.DataFrame(REs_isUnique100, columns=(["uniqREin200"])), pd.DataFrame(REs_isUnique50, columns=(["uniqREin100"]))], 1)
            else:
                REs_isUnique100 = ""
                REs_isUnique50 = ""
        else:
            print("No on-target found for the input gRNAs with your search criteria!")
            gRNA_bed = ""
            REs_isUnique100 = ""
            REs_isUnique50 = ""
        if "id" in summary.columns:
            summary.drop("id", axis = 1, inplace = True)
        if "id" in offTargets[0].columns:
            offTargets[0].drop("id", axis = 1, inplace = True)
                
        gRNAs_notInGenome = set(gRNAs["names"]).difference(set(summary["names"]))
   
        if len(gRNAs_notInGenome) > 0:
            dat2 = pd.DataFrame(np.array([ ["NA"] * np.shape(summary)[1]] * len(gRNAs_notInGenome)))
            for i in range(np.shape(summary)[1]):
                dat2 = dat2.rename(columns={i:summary.columns.values[i]})
            for i in range(len(dat2["names"])):
                dat2.loc[i, ("names")] = list(gRNAs_notInGenome)[i]
                for j in range(len(gRNAs["names"])):
                    if gRNAs.loc[j, ("names")] == dat2.loc[i, ("names")]:
                        dat2.loc[i, ("gRNAsPlusPAM")] = gRNAs.loc[j, ("seq")][0:gRNA_size] + PAM  
                    
            summary = pd.concat([summary, dat2], 0).reset_index(drop = True)
            
        if np.shape(on_target)[0] == 0:
            summary.sort_values("names").reset_index(drop=True).to_csv(outputDir + "Summary.xls", index = False)
        else:
            summary.sort_values("forViewInUCSC").reset_index(drop=True).to_csv(outputDir + "Summary.xls", index = False)
       
        print("Done. Please check output files in directory \n", outputDir, "\n")
        return {"on_target": on_target, "summary": summary, "offtarget": offTargets[0], "gRNAs.bedFormat": gRNA_bed, "REcutDetails": REcutDetails, "REs_isUnique100": REs_isUnique100, "REs_isUnique50": REs_isUnique50}
    
    else:
        x = []
        if PAM_location == "3prime":
            for i in range(np.shape(gRNAs)[0]):
                x.append(gRNAs.loc[i, ("seq")][0:gRNA_size] + PAM)
        else:
            for i in range(np.shape(gRNAs)[0]):
                x.append(PAM + gRNAs.loc[i, ("seq")][0:gRNA_size])
        
        summary = pd.concat([gRNAs["names"], pd.DataFrame(x, columns=(["gRNAsPlusPAM"])), pd.DataFrame(["NA"]*len(gRNAs), columns=(["top5OfftargetTotalScore"])),
                            pd.DataFrame(["NA"] * len(gRNAs), columns=(["top10OfftargetTotalScore"])),
                            pd.DataFrame(["NA"] * len(gRNAs), columns=(["top1Hit.onTarget.MMdistance2PAM "]))], 1)
        
        summary.to_csv(outputDir + "Summary.xls", index = False)

   
        return summary

