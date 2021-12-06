# -*- coding: utf-8 -*-
"""
Created on Tue May 19 18:17:55 2020

@author: Klaus
"""

import os
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import py2bit


pd.set_option("display.max_rows", 100)
pd.set_option("display.max_columns", 50)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.width", None)



def getSeqFromBed(inputFilePath, BSgenomeName, header = False, upstream = 0, downstream = 0):
    
    bed = pd.read_csv(inputFilePath, sep="\t", header = header)
   
    
    if np.shape(bed)[1] < 3:
        raise TypeError("inputfile specified as", inputFilePath, " is not a valid bed file!")

    
    if np.shape(bed)[1] >= 6:
        strand = list(bed.iloc[:, 5])
        for i in range(len(strand)):
            if strand[i] != "+" and "-" and "*":
                strand[i] = "*"
    
    else:
        strand = np.repeat("+", np.shape(bed)[0], axis=0)
    
    if np.shape(bed)[1] >= 4:
        names = bed.iloc[:, 3]
    else:
        for i in range(np.shape(bed)[0]):
            for j in range(np.shape(bed)[1]):
                bed.iloc[i, j] = str(bed.iloc[i, j])
    
        names = bed.iloc[:, 0].str.cat(bed.iloc[:, 1].str.cat(bed.iloc[:, 2], sep="-"), sep=":")
    
    chr = bed.iloc[:, 0]
   
    for i in range(np.shape(bed)[0]):
        for j in range(1, np.shape(bed)[1]):
            bed.iloc[i, 1] = int(bed.iloc[i, 1])
            bed.iloc[i, 2] = int(bed.iloc[i, 2])
            
    Start = bed.iloc[:, 1]
    End = bed.iloc[:, 2]
    
    for i in range(len(strand)):
        if strand[i] == "-":
            Start[i] = Start[i] - downstream
            End[i] = End[i] + upstream
        else:
            Start[i] = Start[i] - upstream
            End[i] = End[i] + downstream
   
    k = []
    for i in range(len(Start)):
        k.append(1)
        if i == len(Start)-1:
            starts = list(pd.Series(pd.concat([Start, pd.DataFrame(k)], axis=1).apply(max, 1)))
     
    width = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms().values()
    chromosome = py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).chroms().keys()
    chromosome = pd.DataFrame(list(chromosome))
    width = pd.DataFrame(list(width))
    widthchromosome = pd.concat([width, chromosome], axis=1)
    
    endd = []
    for i in chr:
        for j in range(len(chromosome)):
            if i == widthchromosome.iloc[j, 1]:
                endd.append(widthchromosome.iloc[j, 0])
    endd = pd.DataFrame(endd)
    ends = list(pd.concat([End, endd], axis=1).apply(min, 1)) 
   
    seq = []
    for i in range(len(Start)):
        if strand[i] == "+":
            seq.append(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i], starts[i]-1, ends[i]))
        else:
            seq.append(str(Seq(py2bit.open(os.getcwd() + os.sep + os.sep.join(["extdata", BSgenomeName, BSgenomeName + ".2bit"])).sequence(chr[i], starts[i]-1, ends[i]))[::-1].reverse_complement())[::-1])
    
    seqnames = []
    for i in range(len(names)):
        seqnames.append(seq[i])
        
    def DNAStringSet(seq):
    
        width = []
        for i in range(len(seq)):
            width.append(len(seq[i]))

        seq = pd.concat([pd.DataFrame(width, columns = ["width"], index = range(len(width))), pd.DataFrame(seq, columns = ["seq"], index=range(len(seq)))], axis=1)

        return seq

    names = list(names)
    names = pd.DataFrame(names, columns = ["names"])
    seq = DNAStringSet(seqnames)
    seq = pd.concat([seq, names], axis=1)

    return seq
    

