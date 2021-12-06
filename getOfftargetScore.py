# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:48:06 2020

@author: Klaus
"""


import pandas as pd
import numpy as np
import re
import buildFeatureVectorForScoring as bFVFS


pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)



def getOfftargetScore(featureVectors, weights=[0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
                                               0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]):
    fv_lessThan2Mismatch = featureVectors[featureVectors["n.mismatch"] < 2].reset_index(drop=True)
    
   
    if np.shape(fv_lessThan2Mismatch)[0] > 0:
        mismatch = pd.DataFrame()
        mismatch_pos = pd.DataFrame()
        for i in range(len(fv_lessThan2Mismatch.columns.values)):
            e =[]
            if "IsMismatch.pos" in fv_lessThan2Mismatch.columns[i]:
                e.append(fv_lessThan2Mismatch.columns[i])
                mismatch = pd.concat([mismatch,fv_lessThan2Mismatch[e]],axis=1)
        
        for j in range(np.shape(mismatch)[0]):
            a = []
            for k in range(np.shape(mismatch)[1]):
                a.append(mismatch.iloc[j][k])
            mismatch_pos=pd.concat([mismatch_pos,pd.DataFrame(a, columns=([mismatch.index.values[j]]))],axis=1)
      
        fv_lessThan2Mismatch["score"] = 100*(1-np.dot(weights,mismatch_pos))
    
    
    fv_geThan2Mismatch = featureVectors[featureVectors["n.mismatch"] >= 2].reset_index(drop=True)
  
    if np.shape(fv_geThan2Mismatch)[0] > 0:
        mismatch = pd.DataFrame()
        mismatch_pos = pd.DataFrame()
        for i in range(len(fv_geThan2Mismatch.columns.values)):
            e =[]
            if "IsMismatch.pos" in fv_geThan2Mismatch.columns[i]:
                e.append(fv_geThan2Mismatch.columns[i])
                mismatch = pd.concat([mismatch,fv_geThan2Mismatch[e]],axis=1)
        
        for j in range(np.shape(mismatch)[0]):
            a = []
            for k in range(np.shape(mismatch)[1]):
                a.append(mismatch.iloc[j][k])
            mismatch_pos=pd.concat([mismatch_pos,pd.DataFrame(a, columns=([mismatch.index.values[j]]))],axis=1)
            
        fv_geThan2Mismatch["score"] = 100
        pos = []
        for y in range(len(fv_lessThan2Mismatch.columns)):
            if "IsMismatch.pos" in fv_lessThan2Mismatch.columns[y]:
                pos.append(y)
        min_pos = min(pos)
        
        mismatch_index = []
       
        for x in range(len(fv_geThan2Mismatch.index.values)):
            a = []
            for c in pos:
                if fv_geThan2Mismatch.iloc[x][c]==1:
                    a.append(c+1)
            mismatch_index.append(a)
        
        score_new = fv_geThan2Mismatch["score"].astype(float)  
   
        '''for v in range(np.shape(mismatch_index)[0]):
            for b in range(np.shape(mismatch_index)[1]):
                score_new[score_new.index.values[v]] = score_new[score_new.index.values[v]] * (1 - weights[mismatch_index[v][b]])
        '''
       
        for v in range(np.shape(mismatch_index)[0]):
            for b in mismatch_index[v]:
                score_new[score_new.index.values[v]] = score_new[score_new.index.values[v]] * (1 - weights[b-1])
        
        fv_geThan2Mismatch["score"] = score_new
        fv_geThan2Mismatch["score"] = fv_geThan2Mismatch["score"] / (((19-fv_geThan2Mismatch["mean.neighbor.distance.mismatch"])/19)*4+1)
        fv_geThan2Mismatch["score"] = fv_geThan2Mismatch["score"] / fv_geThan2Mismatch["n.mismatch"]**2
        fv_geThan2Mismatch.loc[:, ("score")][fv_geThan2Mismatch.loc[:, ("score")] < 0] = 1
       
    if np.shape(fv_geThan2Mismatch)[0] > 0:
        score = fv_geThan2Mismatch
        if np.shape(fv_lessThan2Mismatch)[0] > 0:
            score = pd.concat([fv_lessThan2Mismatch,fv_geThan2Mismatch],axis=0)
    else:
        score = fv_lessThan2Mismatch
    
    score["score"] = round(score["score"], 1)
    score = score.sort_values(by=["name", "score"], ascending=False).reset_index(drop=True)
    
    return score

