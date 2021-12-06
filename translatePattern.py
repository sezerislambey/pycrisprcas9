# -*- coding: utf-8 -*-
"""
Created on Sun May 31 15:32:20 2020

@author: Klaus
"""

import re
ambi_base_map = {'Y':['C', 'T'],
                 'R':['A', 'G'],
                 'S':['C', 'G'],
                 'W':['A', 'T'],
                 'K':['G', 'T'],
                 'M':['A', 'C'],
                 'B':['C', 'G', 'T'],
                 'D':['A', 'G', 'T'],
                 'H':['A', 'C', 'T'],
                 'V':['A', 'C', 'G'],
                 'N':['A', 'C', 'G', 'T']}
 
 

def translatePattern(ambiguous_seq): 
    ambiguous_base = re.findall(r'[RYMKSWHBVDN]', ambiguous_seq) 
    out = [ambiguous_seq] 
    for ambi_base in ambiguous_base: #
        result = []
        for i in ambi_base_map[ambi_base]:
            for j in out:
                result.append(j.replace(ambi_base, i, 1))
        out = result
    return out
 
