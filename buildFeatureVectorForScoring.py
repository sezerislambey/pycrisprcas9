# -*- coding: utf-8 -*-
"""
Created on Sun May 31 10:47:34 2020

@author: Klaus
"""

import pandas as pd
import numpy as np
from Bio.Seq import Seq
import translatePattern as tp
from Bio.SeqFeature import SeqFeature, FeatureLocation


pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",50)
pd.set_option("display.max_colwidth",50)
pd.set_option("display.width",None)
pd.set_option('expand_frame_repr', True)




def mismatches_as_IntegerList(mismatches):
       
    mismatch_pos = []
    
    for t in range(np.shape(mismatches)[0]):
        a = []
        for y in range(np.shape(mismatches)[1]):
            
            if mismatches.iloc[t, y] == 1:
                a.append(y+1)
        mismatch_pos.append(a)
    
    mismatch_pos = pd.DataFrame(mismatch_pos)

    if np.shape(mismatch_pos)[0] > 0 and np.shape(mismatch_pos)[1] > 0 :
        for p in range(np.shape(mismatch_pos)[0]):
            for l in range(np.shape(mismatch_pos)[1]):
                try:
                    mismatch_pos.iloc[p, l] = int(mismatch_pos.iloc[p, l])
                except:
                    mismatch_pos.iloc[p, l] = "-"
    else:
        for i in range(np.shape(mismatches)[0]):
            for j in range(2):
                mismatch_pos.loc[i, j] = "-"
  
    return mismatch_pos



def buildFeatureVectorForScoring(hits, gRNA_size=20, canonical_PAM = "NGG", subPAM_position = [22, 23], PAM_size = 3, PAM_location = "3prime"):

    if np.shape(hits)[0] == 0:
        raise("Empty hits!")
    def DNAStringSet(seq):
    
        width = []
        for i in range(len(seq)):
            width.append(len(seq[i]))

        seq = pd.concat([pd.DataFrame(width, columns = ["width"], index=range(len(width))), pd.DataFrame(pd.Series(seq), columns = ["seq"], index=range(len(seq)))], axis=1)

        return seq
  
    subject = DNAStringSet(list(hits["OffTargetSequence"]))
   
    canonical_PAM_list = tp.translatePattern(canonical_PAM)
    isCanonical_PAM = []
   
    if PAM_location == "3prime":
        for i in range(np.shape(subject)[0]):
            if subject.iloc[i,1][gRNA_size:] in canonical_PAM_list:
                isCanonical_PAM.append(1)
            else:
                isCanonical_PAM.append(0)
        
    else:
        for i in range(np.shape(subject)[0]):
            if subject.iloc[i,1][0:len(canonical_PAM_list[0])] in canonical_PAM_list:
                isCanonical_PAM.append(1)
            else:
                isCanonical_PAM.append(0)

    PAM = []
    for j in range(len(subject)):
        PAM.append(subject.iloc[j,1][subPAM_position[0]-1:subPAM_position[1]])
    mismatches = []
    mismatches = pd.DataFrame(mismatches)
    
    for k in range(len(hits.columns.values)):
        if "IsMismatch.pos" in hits.columns.values[k]:
            mismatches = pd.concat([mismatches, hits.iloc[:,k]],axis=1)
   
    mismatch_pos = mismatches_as_IntegerList(mismatches)
   
    if PAM_location == "3prime":
        at = []
        for p in range(np.shape(mismatch_pos)[0]):
            atstart = []
            atend = []
            atwidth = []
            for r in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[p,r] == "-":
                    atstart.append(mismatch_pos.iloc[p,r])
                    atend.append(mismatch_pos.iloc[p,r])
                    atwidth.append(mismatch_pos.iloc[p,r])
            
                else:
                    atstart.append(int(mismatch_pos.iloc[p,r]))
                    atend.append(int(mismatch_pos.iloc[p,r]))
                    atwidth.append(len(range(int(mismatch_pos.iloc[p,r]),int(mismatch_pos.iloc[p,r]+1))))
               
            at.append(pd.concat([pd.DataFrame(atstart,columns=["Start"]),pd.DataFrame(atend,columns=["End"]),pd.DataFrame(atwidth,columns=["Width"])],axis=1))
        at00 = pd.DataFrame()
        for h in range(len(pd.Series(at))):    
            at00 = pd.concat([at00, pd.Series(at)[h]], axis=0)
        at = at00

    else:
        at_old = []
        at = []
        for p in range(np.shape(mismatch_pos)[0]):
            atstart_old = []
            atend_old = []
            atwidth_old = []
            atstart = []
            atend = []
            atwidth = []
            for r in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[p,r] == "-":
                    atstart_old.append(mismatch_pos.iloc[p,r])
                    atend_old.append(mismatch_pos.iloc[p,r])
                    atwidth_old.append(mismatch_pos.iloc[p,r])
                    atstart.append(mismatch_pos.iloc[p,r])
                    atend.append(mismatch_pos.iloc[p,r])
                    atwidth.append(mismatch_pos.iloc[p,r])
            
                else:
                    atstart_old.append(int(mismatch_pos.iloc[p,r]))
                    atend_old.append(int(mismatch_pos.iloc[p,r]))
                    atwidth_old.append(len(range(int(mismatch_pos.iloc[p,r]),int(mismatch_pos.iloc[p,r]+1))))
                    atstart.append(int(mismatch_pos.iloc[p,r]+PAM_size))
                    atend.append(int(mismatch_pos.iloc[p,r]+PAM_size))
                    atwidth.append(len(range(int(mismatch_pos.iloc[p,r]+PAM_size),int(mismatch_pos.iloc[p,r]+PAM_size+1))))
            at_old.append(pd.concat([pd.DataFrame(atstart_old,columns=["Start"]),pd.DataFrame(atend_old,columns=["End"]),pd.DataFrame(atwidth_old,columns=["Width"])],axis=1))
            at.append(pd.concat([pd.DataFrame(atstart,columns=["Start"]),pd.DataFrame(atend,columns=["End"]),pd.DataFrame(atwidth,columns=["Width"])],axis=1))
        at00 = pd.DataFrame()
        at01 = pd.DataFrame()
        for h in range(len(pd.Series(at))):    
            at00 = pd.concat([at00, pd.Series(at)[h]], axis=0).reset_index(drop=True)
            at01 = pd.concat([at01, pd.Series(at_old)[h]], axis=0).reset_index(drop=True)
        at_old = at01
        at = at00
  
    if sum(hits["n.mismatch"]) > 0:
        d_nucleotide = []
        r_nucleotide = []
        a = hits["gRNAPlusPAM"]
        a = DNAStringSet(list(a))
        for r in range(len(subject)):
            d_nucleotide1 = []
            r_nucleotide1 = []
            for j in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[r][j] == "-":
                    d_nucleotide1.append("-")
                else:
                    if PAM_location == "3prime":
                        feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r,j]-1), int(mismatch_pos.iloc[r,j])), type="gene", strand=-1)
                        d_nucleotide1.append(str(Seq(subject.iloc[r, 1]).complement())[feature.location.start:feature.location.end])
                    else:
                        feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r][j]-1+PAM_size), int(mismatch_pos.iloc[r][j]+PAM_size)), type="gene", strand=-1)
                        d_nucleotide1.append(str(Seq(subject.iloc[r, 1]).complement())[feature.location.start:feature.location.end])
                    
            d_nucleotide.append(d_nucleotide1)
            
            for h in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[r][h] == "-":
                    r_nucleotide1.append("-")
                else:
                    if PAM_location == "3prime":
                        feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r,h]-1), int(mismatch_pos.iloc[r,h])), type="gene", strand=-1)
                        r_nucleotide1.append(a.iloc[r, 1][feature.location.start:feature.location.end])
                    else:
                        feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r,h]-1+PAM_size), int(mismatch_pos.iloc[r,h]+PAM_size)), type="gene", strand=-1)
                        r_nucleotide1.append(a.iloc[r,1][feature.location.start:feature.location.end])
            r_nucleotide.append(r_nucleotide1)
        
        for i in range(np.shape(r_nucleotide)[0]):
            for j in range(np.shape(r_nucleotide)[1]):
                if r_nucleotide[i][j] == "T":
                    r_nucleotide[i][j] = "U"
        
        d_nu_r_nu = []
        for b in range(len(r_nucleotide)):
            for v in range(len(r_nucleotide[0])):
                if r_nucleotide[b][v] != "-":
                    d_nu_r_nu.append("r"+r_nucleotide[b][v]+":d"+d_nucleotide[b][v])
        arr_ind = np.array(np.where(mismatches != 0))+1        
        arr_ind = pd.concat([pd.DataFrame(arr_ind[0],index=range(len(arr_ind[0])),columns=["row"]),pd.DataFrame(arr_ind[1],index=range(len(arr_ind[0])),columns=["col"])], axis=1).sort_values(by=["col", "row"]).reset_index(drop=True)
        ind_row = arr_ind["row"].sort_values().reset_index(drop=True)
        
        s = 1
        e = 1
        partitioningwidth = []
        partitioningstart = []
        partitioningend = []
      
        for i in range(1, len(mismatches)+1):
            if i in list(ind_row): 
                partitioningwidth.append(list(ind_row).count(i))
                partitioningstart.append(s)
                e = list(ind_row).count(i) + s - 1
                s = e + 1
                partitioningend.append(e)
            else:
                partitioningwidth.append(0)
                partitioningstart.append(s)
                partitioningend.append(e)
                e = s - 1
        partitioning = pd.concat([pd.DataFrame(partitioningstart,columns=["Start"]),pd.DataFrame(partitioningend,columns=["End"]),pd.DataFrame(partitioningwidth,columns=["Width"])],axis=1)
        
        d_nu_r_nu2 = []
        for i in range(len(partitioning)):
            d = []
            if partitioning.iloc[i,0] > partitioning.iloc[i,1]:
                d_nu_r_nu2.append([""])
            elif partitioning.iloc[i,0] == partitioning.iloc[i,1]:
                d_nu_r_nu2.append(d_nu_r_nu[partitioning.iloc[i,0]])
            else:
                for k in range(partitioning.iloc[i,0]-1, partitioning.iloc[i,1]):
                   
                    d.append(d_nu_r_nu[k])
                d_nu_r_nu2.append(d)
        d_nu_r_nu2 = (pd.Series(d_nu_r_nu2))
        d_nu_r_nu = list(d_nu_r_nu2)
       
    else:
        d_nu_r_nu = ""

    if PAM_location == "3prime":
        mismatch_distance2PAM = []
        for i in range(np.shape(mismatch_pos)[0]):
            m = []
            for j in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[i,j] == "-":
                    m.append("-")
                else:
                    m.append(str(int(gRNA_size + 1 - mismatch_pos.iloc[i,j])))
            mismatch_distance2PAM.append(m)  
    else:
        mismatch_distance2PAM = []
        for i in range(np.shape(mismatch_pos)[0]):
            m = []
            for j in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[i,j] == "-":
                    m.append("-")
                else:
                    m.append(str(int(mismatch_pos.iloc[i,j])))
            mismatch_distance2PAM.append(m)                 

    alignment1 = DNAStringSet(["."]*np.shape(hits)[0])*gRNA_size
    
    if PAM_location == "3prime":
        ex = []
        for r in range(len(subject)):
            ex1 = []
            for j in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[r][j] == "-":
                    ex1.append("-") 
                else:
                    feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r][j]-1), int(mismatch_pos.iloc[r][j])), type="gene", strand=-1)
                    ex1.append(str(Seq(subject.iloc[r, 1]))[feature.location.start:feature.location.end])
            ex.append(ex1)
       
        z = 0
        alignment = []
        for i in range(np.shape(alignment1)[0]):
            bb = []
            cc = list(alignment1.iloc[i,1])
            for j in range(np.shape(ex)[1]):
                #z += j
                if at.iloc[z,0] == "-":
                    pass
                else:
                    cc[int(at.iloc[z,0])-1] = ex[i][j]
                z += 1
           
            rr = ""
            for i in range(len(cc)):
                rr += cc[i]
            alignment.append(rr)
            
    else:
        ex = []
        for r in range(len(subject)):
            ex1 = []
            for j in range(np.shape(mismatch_pos)[1]):
                if mismatch_pos.iloc[r][j] == "-":
                    ex1.append("-") 
                else:
                    feature = SeqFeature(FeatureLocation(int(mismatch_pos.iloc[r][j]-1+PAM_size), int(mismatch_pos.iloc[r][j]+PAM_size)), type="gene", strand=-1)
                    ex1.append(str(Seq(subject.iloc[r, 1]))[feature.location.start:feature.location.end])
            ex.append(ex1)
        
        z = 0
        alignment = []
        for i in range(np.shape(alignment1)[0]):
            bb = []
            cc = list(alignment1.iloc[i,1])
            for j in range(np.shape(ex)[1]):
                #z += j
                if at_old.iloc[z,0] == "-":
                    pass
                else:
                    cc[int(at_old.iloc[z,0])-1] = ex[i][j]
                    
                z += 1
                
            rr = ""
            for i in range(len(cc)):
                rr += cc[i]
            alignment.append(rr)
  
    mean_neighbor_distance_mismatch = []
    '''for i in range(np.shape(mismatch_pos)[0]):
        if mismatch_pos.iloc[i,0] == "-":
            mean_neighbor_distance_mismatch.append("NaN")
        else:
            mean_neighbor_distance_mismatch.append(float(mismatch_pos.iloc[i,1]) - float(mismatch_pos.iloc[i,0]))
    '''
    for i in range(np.shape(mismatch_pos)[0]):
        if  "-" in  list(mismatch_pos.iloc[i, :]):
            if  np.shape(mismatch_pos)[1] == list(mismatch_pos.iloc[i, :]).count("-"):
                mean_neighbor_distance_mismatch.append("NaN")
            else:
                mean_neighbor_distance_mismatch.append(np.mean(np.diff(mismatch_pos.iloc[i, :len(mismatch_pos.iloc[i,:])-list(mismatch_pos.iloc[i,:]).count("-")])))
        else:
            mean_neighbor_distance_mismatch.append(np.mean(np.diff(mismatch_pos.iloc[i,:])))
   
    no_neighbor_idx = []    
    for i in range(len(mean_neighbor_distance_mismatch)):
        if mean_neighbor_distance_mismatch[i] == "NaN":
            no_neighbor_idx.append(True)
        else: 
            no_neighbor_idx.append(mean_neighbor_distance_mismatch[i] <= 1)
    
    for i in range(len(no_neighbor_idx)):
        if no_neighbor_idx[i] == True:
            mean_neighbor_distance_mismatch[i] = gRNA_size

    features = pd.concat([pd.DataFrame(pd.Series(mismatch_distance2PAM),columns=["mismatch.distance2PAM"]),pd.DataFrame(pd.Series(alignment), columns=["alignment"]),pd.DataFrame(pd.Series(isCanonical_PAM),columns=["NGG"]),pd.DataFrame(pd.Series(mean_neighbor_distance_mismatch),columns=["mean.neighbor.distance.mismatch"]),pd.DataFrame(pd.Series(d_nu_r_nu),columns=["mismatch.type"]),pd.DataFrame(pd.Series(PAM),columns=["SubPAM"])],axis=1)
    features["mismatch.type"] = features["mismatch.type"].fillna("")
    for i in range(np.shape(features)[0]):
        if features["mismatch.distance2PAM"].iloc[i] == ["-", "-"]:
            features.loc[i, ("mismatch.distance2PAM")] = ""
    return pd.concat([hits, features], axis = 1)
 
