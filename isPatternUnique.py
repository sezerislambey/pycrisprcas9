import character as ch
import gregexprindex as gr
import translatePattern as tP
from Bio.Seq import Seq
import biostrings as bio



def isPatternUnique(seq, patterns):

    if seq == None and patterns == None and len(seq) == 0 and len(patterns) == 0:
        IsUnique = "NA"
    else:
        IsUnique = ch.character(len(patterns))
        seq = str(seq)

        for i in range(len(patterns)):
            pattern = patterns.iloc[i, 1]
            revpattern = Seq(pattern).reverse_complement()
            pattern = str(pattern)
            
            if len(pattern) == 0:
                IsUnique[i] = "NA"
                continue

            pattern = tP.translatePattern(pattern)
            revpattern = str(revpattern)
      
            revpattern = tP.translatePattern(revpattern)
            
            resplus = gr.gregexpr_index(pattern, seq)
            nplus = len(resplus)
            
            if revpattern != pattern:
                resminus = gr.gregexpr_index(revpattern, seq)
                nminus = len(resminus)
          
            else:
                nminus = 0
            if nplus + nminus > 1:
                IsUnique[i] = "No"
            elif nplus + nminus == 1:
                IsUnique[i] = "Yes"
            else:
                IsUnique[i] = "Not Found"
    return IsUnique