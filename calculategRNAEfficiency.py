#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 10:10:05 2021

@author: klaus
"""

import pandas as pd
import numpy as np
import re
from joblib import Parallel, delayed
import multiprocessing


pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",50)
pd.set_option("display.max_colwidth",50)
pd.set_option("display.width",None)
pd.set_option('expand_frame_repr', True)


"""
extendedSequence = ["TGGATTGTATAATCAGCATGGATTTGGAAC","TCAACGAGGATATTCTCAGGCTTCAGGTCC",
"GTTACCTGAATTTGACCTGCTCGGAGGTAA","CTTGGTGTGGCTTCCTTTAAGACATGGAGC",
"CATACAGGCATTGAAGAAGAATTTAGGCCT","AGTACTATACATTTGGCTTAGATTTGGCGG",
"TTTTCCAGATAGCCGATCTTGGTGTGGCTT","AAGAAGGGAACTATTCGCTGGTGATGGAGT"]
"""

def calculategRNAEfficiency(extendedSequence, baseBeforegRNA, featureWeightMatrix, gRNA_size = 20, enable_multicore = False, n_cores_max = 6):
    featureWeightMatrix.iloc[:, 0] = featureWeightMatrix.iloc[:, 0].str.upper()
    features = []
    fWeights = []
    featureNames = []
    featureStart = []
    featureEnd = []
    for i in range(0, len(featureWeightMatrix.iloc[:, 0])):
        if featureWeightMatrix.iloc[i, 0] == "INTERCEPT":
            efficiency = featureWeightMatrix.iloc[i,1]
   
        elif featureWeightMatrix.iloc[i, 0] == "GC_LOW":
            GClow= featureWeightMatrix.iloc[i,1]
        
        elif featureWeightMatrix.iloc[i, 0] == "GC_HIGH":
            GChigh= featureWeightMatrix.iloc[i,1]
         
        else:
            features.append(featureWeightMatrix.iloc[i,0])
            fWeights.append(featureWeightMatrix.iloc[i, 1])
            featureNames.append(re.sub(r'[0123456789]',"", featureWeightMatrix.iloc[i,0]))
            featureStart.append(int(re.sub(r"[ACGT]", "", featureWeightMatrix.iloc[i,0])))
            featureEnd.append(int(re.sub(r"[ACGT]", "", featureWeightMatrix.iloc[i,0])) + len(re.sub(r'[0123456789]',"", featureWeightMatrix.iloc[i,0])) - 1)
    featureWeights = pd.concat([pd.DataFrame(pd.Series(featureStart),columns=["featureStart"]),pd.DataFrame(pd.Series(featureEnd), columns=["featureEnd"]),pd.DataFrame(pd.Series(featureNames),columns=["featureNames"]),pd.DataFrame(pd.Series(fWeights),columns=["fWeights"])],axis=1)
    featureWeights.index = featureWeights.index + 1
    
    for i in range(len(featureWeights)):
        featureWeights.iloc[i,0] = str(featureWeights.iloc[i,0])
        featureWeights.iloc[i,1] = str(featureWeights.iloc[i,1])
    
    featureWeights = featureWeights.sort_values(by = ["featureStart", "featureEnd"])
    
    n_cores = multiprocessing.cpu_count() - 1
    n_cores = min(n_cores, len(extendedSequence))
    n_cores = min(n_cores, n_cores_max)
    
    
    n_C = []
    n_G = []
    GC_content = []
    GClow_coef = []
    GChigh_coef = []
    for i in range(len(extendedSequence)):
        
        n_C.append(list(extendedSequence[i][baseBeforegRNA:baseBeforegRNA+gRNA_size]).count("C"))
        n_G.append(list(extendedSequence[i][baseBeforegRNA:baseBeforegRNA+gRNA_size]).count("G"))
        GC_content.append(n_C[i]+n_G[i])
        if GC_content[i] < 10:
            GClow_coef.append(1)
            GChigh_coef.append(0)
        elif GC_content[i] > 10:
            GClow_coef.append(0)
            GChigh_coef.append(1)
        else:
            GClow_coef.append(0)
            GChigh_coef.append(0)
    
    for i in range(0, len(featureWeightMatrix.iloc[:, 0])):
        if featureWeightMatrix.iloc[i, 0] == "INTERCEPT":
            efficiency = featureWeightMatrix.iloc[i,1]
    
    for i in range(len(featureWeightMatrix)):
        if featureWeightMatrix.iloc[i, 0] == "gc_low": 
            GClow = featureWeightMatrix.iloc[i, 1]
        elif featureWeightMatrix.iloc[i, 0] == "gc_high":
            GChigh = featureWeightMatrix.iloc[i, 1]

   
    efficiency = efficiency + GClow_coef * (np.float64(10) - GC_content) * GClow + GChigh_coef * ( GC_content - np.float64(10)) * GChigh
    
    def numeric(length):
        a = []
        for i in range(0, length):
            a.append(0)
        return a

   
    all_features = []
    for j in range(len(extendedSequence)):
        this_features = []
        for i in range(len(featureWeights)):
            if extendedSequence[j][int(featureWeights.iloc[i, 0])-1:int(featureWeights.iloc[i, 1])] == featureWeights.iloc[i, 2]:
                this_features.append(1) 
            else:
                this_features.append(0)
        all_features.append(this_features)

    efficiency = efficiency + np.dot(all_features, featureWeights.iloc[:,3])
    
    efficiency = 1 / (1 + np.exp(-efficiency))
    
    
    return efficiency


