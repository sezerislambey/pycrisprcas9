import pandas as pd
from Bio.Seq import Seq


def DNAStringSet(seq):

    width = []
    for i in range(len(seq)):
        width.append(len(seq[i]))

    DNAStringSet = pd.concat([pd.DataFrame(width, columns=["width"]), pd.DataFrame(seq, columns=["seq"])],axis=1)

    #for i in range(len(DNAStringSet)):
    #    DNAStringSet.iloc[i, 1] = Seq(DNAStringSet.iloc[i, 1])

    return DNAStringSet
