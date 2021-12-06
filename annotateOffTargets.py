#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 13:28:30 2021

@author: klaus
"""

import buildFeatureVectorForScoring as bFVFS
import pandas as pd
import getOfftargetScore
import pyranges as pr
import sqlite3
import re
import numpy as np
import character as ch



def annotateOffTargets(scores, txdb, orgAnn, ignore_strand = True):
    
    conn = sqlite3.connect(txdb)
    cursor = conn.execute("""SELECT * FROM exon""")
    exons = cursor.fetchall()
    
    conn2 = sqlite3.connect(orgAnn)

    
    score_RD = pd.concat([scores["chrom"], scores["chromStart"], scores["chromEnd"], scores["strand"], scores["forViewInUCSC"]], 1)
 
    seqnames = []
    strand = []
    ranges1 = []
    ranges2 = []
    for i in range(len(exons)):
       seqnames.append(exons[i][2])
       strand.append(exons[i][3])
       ranges1.append(exons[i][4])
       ranges2.append(exons[i][5])
    
    ranges = []
    for j in range(len(ranges1)):
        ranges.append(str(ranges1[j]) + "-" + str(ranges2[j]))
    allExons = pd.concat([pd.DataFrame(seqnames, columns=(["chrom"])), pd.DataFrame(ranges1,columns=(["chromStart"])), pd.DataFrame(ranges2,columns=(["chromEnd"])), pd.DataFrame(strand, columns=(["strand"]))], 1)
  
   
    a = 0
    b = 0
    for i in range(len(allExons["chrom"])):
        if "Chr" in allExons["chrom"][i]:
            a += 1
    for i in range(len(scores["chrom"])):
        if "Chr" in scores["chrom"][i]:
            b += 1

    if a == 0 & b < 0:
        for i in range(len(allExons["chrom"])):
            allExons["chrom"][i] = "Chr" + allExons["chrom"][i]
    
    allExons = allExons[allExons["chrom"].isin(list(scores["chrom"]))]
   
    if ignore_strand:
        score_plus_RD = pd.concat([scores["chrom"], scores["chromStart"], scores["chromEnd"], scores["strand"], scores["forViewInUCSC"]], 1)
        score_plus_RD["strand"] = "+"
        score_minus_RD = pd.concat([scores["chrom"], scores["chromStart"], scores["chromEnd"], scores["strand"], scores["forViewInUCSC"]], 1)
        score_minus_RD["strand"] = "-"
        score_RD = pd.concat([score_plus_RD, score_minus_RD], 0)
        scores = pd.concat([scores, scores], 0)

    
    ann_scores = []
    for i in range(len(score_RD)):
        ann = False
        for j in range(len(allExons)):
            if list(score_RD.iloc[i])[0] == list(allExons.iloc[j])[0] and list(score_RD.iloc[i])[1] >= list(allExons.iloc[j])[1] and list(score_RD.iloc[i])[2] <= list(allExons.iloc[j])[2] and list(score_RD.iloc[i])[3] == list(allExons.iloc[j])[3]:
                ann_scores.append(True)
                ann = True
                break
        if ann == False:
            ann_scores.append(False)
    ann_scores = pd.DataFrame(ann_scores, columns=(["inExon"]))
    
    inExon = pd.concat([pd.DataFrame(list(score_RD["forViewInUCSC"]), columns=(["forViewInUCSC"])), ann_scores], 1)            
    
    
    cursor1 = conn.execute("""SELECT * FROM gene""")
    gene = pd.DataFrame(cursor1.fetchall(), columns=(["gene_id", "tx_id"])).sort_values(["gene_id", "tx_id"])
    
    cursor2 = conn.execute("""SELECT * FROM transcript""")
    trans = pd.DataFrame(cursor2.fetchall(), columns=(["tx_id", "tx_name", "tx_type", "tx_chrom", "tx_strand", "tx_start", "tx_end"]))
    
    genes0 = gene.merge(trans, how ="inner", on = "tx_id")
    
    genes1 = genes0.drop_duplicates(subset=["gene_id", "tx_start"], keep="first")
    genes2 = genes0.drop_duplicates(subset=["gene_id", "tx_start"], keep="last")
    allGenes = genes1.merge(genes2, how = "inner", on = ["gene_id", "tx_chrom", "tx_strand"]).sort_values("gene_id")
    allGenes = allGenes[["gene_id", "tx_chrom",  "tx_start_x", "tx_end_y","tx_strand"]].rename(columns={"tx_start_x": "tx_start", "tx_end_y": "tx_end"})
   
    score_RD = score_RD.reset_index(drop = True)
    queryHits = []
    subjectHits = []

    for i in range(len(score_RD)):
        for j in range(len(allGenes)):
            if score_RD["chrom"][i] == allGenes["tx_chrom"][j] and score_RD["strand"][i] == allGenes["tx_strand"][j] and score_RD["chromStart"][i] >= allGenes["tx_start"][j] and score_RD["chromEnd"][i] <= allGenes["tx_end"][j]:
                queryHits.append(i)
                subjectHits.append(j)
   
    queryHits = pd.DataFrame(queryHits, columns=(["queryHits"]))
    subjectHits = pd.DataFrame(subjectHits, columns=(["subjectHits"]))
    overlapGenes = pd.concat([queryHits, subjectHits], 1)

    entrez_id = ch.character(np.shape(scores)[0])
    symbol = entrez_id
  
    for i in range(len(overlapGenes)):
        entrez_id[overlapGenes["queryHits"][i]] = allGenes["gene_id"][overlapGenes["subjectHits"][i]]
    
    symbol = pd.DataFrame(symbol, columns = (["symbol"]))
    entrez_id = pd.DataFrame(entrez_id, columns = (["entrez_id"]))
 
    scores = pd.concat([scores.reset_index(drop=True), entrez_id, symbol], 1)
   
    if len(overlapGenes["queryHits"]) > 0 and  orgAnn != None:
        cursor = conn2.execute("""SELECT * FROM genes""")
        txdb_geneinfo = pd.DataFrame(conn2.execute("""SELECT * FROM gene_info""").fetchall(), columns=(["_id", "gene_name", "symbol"]))
        txdb_genes = pd.DataFrame(cursor.fetchall(),  columns=(["_id", "gene_id"]))
        egSYMBOL = txdb_genes.merge(txdb_geneinfo, how = "inner", on = "_id")[["gene_id", "symbol"]]
        if "flybase_id" == egSYMBOL.columns[1]:
           for i in range(len(egSYMBOL["gene_id"])):
                for j in range(len(scores)):
                    if  egSYMBOL["gene_id"][i] == scores["entrez_id"][j]:
                        scores["symbol"][j] = egSYMBOL["flybase_id"][i]
        else:
            for i in range(len(egSYMBOL["gene_id"])):
                for j in range(len(scores)):
                    if  egSYMBOL.loc[i, ("gene_id")] == scores.loc[j, ("entrez_id")]:
                        scores.loc[j, ("symbol")] = egSYMBOL.loc[i, ("symbol")]
                        
    scores = pd.concat([scores, inExon["inExon"]], 1)
    
    inIntron = []
    for i in range(len(scores)):
        if scores.loc[i, ("entrez_id")] != "" and scores.loc[i, ("inExon")] == False:
            inIntron.append(True)
        else:
            inIntron.append(False)
            
    scores = pd.concat([scores, pd.DataFrame(inIntron, columns = (["inIntron"]))], 1).sort_values(by=["entrez_id", "name"], ascending=False, na_position = "last")

    scores = pd.concat([scores[scores["entrez_id"] != ""], scores.drop_duplicates("forViewInUCSC")], 0).drop_duplicates(["name", "gRNAPlusPAM", "strand", "forViewInUCSC"]).reset_index(drop=True)

    if "FBgn" in list(entrez_id)[0] and orgAnn != None:
        temp_id = scores["entrez_id"]
        scores["entrez_id"] = scores["symbol"]
        scores["symbol"] = temp_id
        del temp_id
    
    if orgAnn != None:
        scores = pd.concat([scores["name"], scores["gRNAPlusPAM"], scores["OffTargetSequence"],
                            scores["inExon"], scores["inIntron"], scores["entrez_id"], scores["symbol"],
                            scores["score"], scores["n.mismatch"], scores["mismatch.distance2PAM"],
                            scores["alignment"], scores["NGG"], scores["forViewInUCSC"],
                            scores["strand"], scores["chrom"], scores["chromStart"], 
                            scores["chromEnd"]], 1).rename(columns = {"symbol":"gene"}).reset_index(drop = True)
    else:
        scores = pd.concat([scores["name"], scores["gRNAPlusPAM"], scores["OffTargetSequence"],
                            scores["inExon"], scores["inIntron"], scores["entrez_id"], scores["symbol"],
                            scores["score"], scores["n.mismatch"], scores["mismatch.distance2PAM"],
                            scores["alignment"], scores["NGG"], scores["forViewInUCSC"],
                            scores["strand"], scores["chrom"], scores["chromStart"], 
                            scores["chromEnd"]], 1).rename(columns = {"symbol":"gene"}).reset_index(drop = True)
    for i in range(np.shape(scores)[0]):
        if scores.loc[i, ("inExon")] == False:
            scores.loc[i , ("inExon")] = ""
        if scores.loc[i, ("inIntron")] == False:
            scores.loc[i, ("inIntron")] = ""
    
    
    return scores
    
