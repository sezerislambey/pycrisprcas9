#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:45:02 2021

@author: klaus
"""



import os
import pandas as pd
import numpy as np
import re


def getOfftargetScore2(featureVectors,
                       mismatch_activity_file = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "NatureBiot2016SuppTable19DoenchRoot.csv"])),
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
                                          "TT" : 0}):
    
    
    mismatch_activity = mismatch_activity_file
    required_col = ["Mismatch.Type", "Position", "Percent.Active"]
    if len(mismatch_activity.columns.intersection(required_col)) != len(required_col):
        raise Exception("Please rename the mismatch activity file column to contain at least these 3 column names: Mismatch.Type, Position, Percent.Active\n")
    
    position = mismatch_activity["Position"]
    r_nu_d_nu = mismatch_activity["Mismatch.Type"]
    weights = mismatch_activity["Percent.Active"]
    
    #mismatch.activity[mismatch.activity$Mismatch.Type == "rA:dG" & 
    #    mismatch.activity$Position == 10,]$Percent.Active
    ##### by default weights is a column vector
    ##### the mismatch activity  is given as pos 20, 19, 18,....1 distance from PAM,    ##### Position named as 1, 2, 3, .....20 though 
    ##### and the featureVectors is in the same order now
    ##### so no need to reverse any more. weights = rev(weights)
    featureVectors.loc[:, ("score")] = 1
    for i in range(len(featureVectors["score"])):
        featureVectors.loc[i, ("score")] = subPAM_activity[featureVectors.loc[i, ("SubPAM")]]
    fv_geThan1Mismatch = featureVectors[featureVectors["n.mismatch"] >= 1].reset_index(drop=True)
    fv_lessThan1Mismatch = featureVectors[featureVectors["n.mismatch"] < 1].reset_index(drop=True)
    
    
    
    if np.shape(fv_geThan1Mismatch)[0] > 0:
        mismatch = pd.DataFrame()
        mismatch_pos = pd.DataFrame()
        for i in range(len(fv_geThan1Mismatch.columns.values)):
            e =[]
            if "IsMismatch.pos" in fv_geThan1Mismatch.columns[i]:
                e.append(fv_geThan1Mismatch.columns[i])
                mismatch = pd.concat([mismatch,fv_geThan1Mismatch[e]],axis=1)
        
        for j in range(np.shape(mismatch)[0]):
            a = []
            for k in range(np.shape(mismatch)[1]):
                a.append(mismatch.iloc[j][k])
            mismatch_pos=pd.concat([mismatch_pos,pd.DataFrame(a, columns=([mismatch.index.values[j]]))], axis=1)
        
        pos = []
        for y in range(len(fv_geThan1Mismatch.columns)):
            if "IsMismatch.pos" in fv_geThan1Mismatch.columns[y]:
                pos.append(y)
        min_pos = min(pos)
        
        mismatch_index = []
        for x in range(len(fv_geThan1Mismatch.index.values)):
      
            for c in pos:
                if fv_geThan1Mismatch.iloc[x][c]==1:
                    mismatch_index.append(c+1)
        score_new = []
       
        
        for i in range(len(fv_geThan1Mismatch)):
            thisMismatch = fv_geThan1Mismatch.loc[i, ("mismatch.type")]
            score_new1 = fv_geThan1Mismatch.loc[i, ("score")]
            for j in range(len(mismatch_index)):
                a = mismatch_activity[mismatch_activity["Position"] == mismatch_index[j]]
                a = a[a["Mismatch.Type"] == thisMismatch[j]].reset_index(drop=True)
                score_new1 = score_new1 * a.loc[0, ("Percent.Active")]
            score_new.append(score_new1)
        fv_geThan1Mismatch["score"] = score_new
  
    if np.shape(fv_geThan1Mismatch)[0] > 0:
        score = fv_geThan1Mismatch
        if np.shape(fv_lessThan1Mismatch)[0] > 0:
            score = pd.concat([fv_lessThan1Mismatch, fv_geThan1Mismatch], axis=0)
    else:
        score = fv_lessThan1Mismatch
    
    score["score"] = round(score["score"],6)
    score = score.sort_values(by="name", ascending=False)
    
    return score
    
