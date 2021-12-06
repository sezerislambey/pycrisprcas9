#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:32:25 2021

@author: klaus
"""


from annotateOffTargets import *
from buildFeatureVectorForScoring import *
from buildFeatureVectorForScoring2 import *
from calculategRNAEfficiency import *
from calculategRNAEfficiencyCRISPRscan import *
from designPEs import *
from filtergRNAs import *
from filterOffTarget import *
from findgRNAs import *
from getOfftargetScore import *
from getOfftargetScore2 import *
from getSeqFromBed import *
from isPatternUnique import *
from offTargetAnalysis import *
from searchHits2 import *
from translatePattern import *
from uniqueREs import *
from writeHits import *
from writeHits2 import *
import time




inputFilePath = "/home/klaus/spyder/crispr/extdata/inputseq.fa"
featureWeightMatrixFile = pd.read_csv(os.getcwd() + os.sep + os.sep.join(["extdata", "Morenos-Mateo.csv"]))


results = offTargetAnalysis(inputFilePath, annotatePaired = False, scoring_method = "CFDscore",
                            PAM_pattern = "NNG$|NGN$", chromToSearch = ["chrX"], annotateExon = False,
                            max_mismatch = 0, BSgenomeName = "BSgenome.Hsapiens.UCSC.hg19", 
                            outputDir = "/home/klaus/", overwrite = True)

print(results)