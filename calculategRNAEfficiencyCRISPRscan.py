import re
import pandas as pd
import numpy as np
import numeric as nm

pd.set_option("display.max_rows",100)
pd.set_option("display.max_columns",50)
pd.set_option("display.max_colwidth",50)
pd.set_option("display.width",None)
pd.set_option('expand_frame_repr', True)



def calculategRNAEfficiencyCRISPRscan(extendedSequence, featureWeightMatrix):

    featureWeightMatrix.iloc[:, 0] = featureWeightMatrix.iloc[:, 0].str.upper()
    features = []
    fWeights = []
  
    for i in range(0, len(featureWeightMatrix.iloc[:, 0])):
        if featureWeightMatrix.iloc[i, 0] == "INTERCEPT":
            efficiency = featureWeightMatrix.iloc[i,1]
        else:
            features.append(featureWeightMatrix.iloc[i, 0])
            fWeights.append(featureWeightMatrix.iloc[i, 1])

    featureNames = []
    featureStart = []
    featureEnd = []

    for i in range(0, len(features)):
        featureNames.append(re.sub("[0-9]+", "", features[i]))

    for i in range(0, len(features)):
        featureStart.append(int(re.sub("[ACGT]+", "", features[i])))

    for i in range(0, len(features)):
        featureEnd.append(featureStart[i] + len(featureNames[i]) - 1)


    featureStart = pd.DataFrame(featureStart,  columns=["featureStart"])
    featureEnd = pd.DataFrame(featureEnd, columns=["featureEnd"])
    featureNames = pd.DataFrame(featureNames, columns=["featureNames"])
    fWeights = pd.DataFrame(fWeights, columns=["fWights"])

    featureWeights = pd.concat([featureStart, featureEnd, featureNames, fWeights], axis=1, names=True)

    for i in range(0,len(featureWeights.iloc[0])):
        for j in range(len(featureEnd)):
            featureWeights.iloc[j, i] = str(featureWeights.iloc[j, i])

    featureWeights = featureWeights.iloc[np.argsort(featureWeights.iloc[:, 0])]


    thisfeatures = nm.numeric(np.shape(featureWeights)[0])

    for i in range(len(thisfeatures)):
        for j in range(len(extendedSequence)):
            thisfeatures[i] = int(extendedSequence[j][int(featureWeights.iloc[i,0])-1:int(featureWeights.iloc[i,1])] == featureWeights.iloc[i, 2])

    allfeatures = []

    for r in range(0, len(extendedSequence)):
        k = []
        allfeatures.append(k)
        for n in range(0, len(thisfeatures)):
            k.append(int(extendedSequence[r][int(featureWeights.iloc[n,0])-1:int(featureWeights.iloc[n,1])] == featureWeights.iloc[n, 2]))


    allfeatures = pd.DataFrame(allfeatures)

    for i in range(len(extendedSequence)):
        for j in range(len(thisfeatures)):
            allfeatures.iloc[i,j] = float(allfeatures.iloc[i, j])

    for i in range(len(featureWeights)):
        featureWeights.iloc[i, 3] = float(featureWeights.iloc[i, 3])

    efficiency = 100 * (efficiency + np.dot(allfeatures.iloc[:], featureWeights.iloc[:,3]))

    return efficiency
