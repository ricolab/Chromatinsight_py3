###########################
### Chromatinsight v3.1 ###
###########################
#
# a set of methods
# used to use a random forest algorithm
# to detect differential patterns between two sets of samples
# analysed with ChIP-seq histone modifications
# and pre-binarised by ChromHMM
#
# Author: Marco Trevisan-Herraz, PhD
# Computational Epigenomics Laboratory
# Newcastle University, UK
# 2019-2021
#
# requirements:
# Python 3.x
# pandas installed, see https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html
# sklearn installed, see https://scikit-learn.org/stable/install.html
# Quick way to install both:
# pip install pandas
# pip install scikit-learn
#

import glob
import os
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.ensemble import RandomForestClassifier
import random


#######################################################

# load any _tab separated values_ file
def load2stringList(fileName, removeCommas=False, splitChar="\t"):
    with open(fileName, "r") as reader:
        fullList = []
        for myRow in reader:
            myRowStrip = myRow.strip()
            if len(myRowStrip) > 0:
                thisRow = myRowStrip.split(splitChar)
                thisRow = [item.strip() for item in thisRow]

                if removeCommas:
                    thisRow = [item[1:len(item) - 1] if item.endswith('"') and item.startswith('"') else item for item
                               in thisRow]

                fullList.append(thisRow)
    return fullList


# ------------------------------------------------------

def stringList2inputDataFile(input, format=['s', 'f', 'f'], fillEmptyPositions=False, emptyFiller=""):
    result = []
    for myRow in input:
        if fillEmptyPositions or len(myRow) >= len(format):
            resultRow = []
            for i in range(len(format)):
                if i > len(myRow) - 1:
                    resultRow.append(emptyFiller)
                else:
                    stringy = myRow[i].strip()
                    if format[i] == 's':
                        resultRow.append(stringy)
                    elif format[i] in ('f', 'i'):
                        if stringy and (stringy[0].isdigit() or (stringy[0] == '-' and stringy[1].isdigit())):
                            resultRow.append(float(stringy) if format[i] == 'f' else int(float(stringy)))
                        else:
                            resultRow = []
                            break
            if resultRow:
                result.append(resultRow)
    return result


# ------------------------------------------------------

def removeHeader(myList):
    if myList:
        myList.pop(0)
    return myList


# ------------------------------------------------------

def saveFile(fileName, list, header=""):
    with open(fileName, "w") as writer:
        if header:
            writer.write(header + "\n")
        for row in list:
            saveRow(writer, row)


# ------------------------------------------------------

def saveRow(writer, rowList):
    line = "\t".join(map(str, rowList)) + "\n"
    writer.write(line)


# ------------------------------------------------------

def mergeRegionFiles(regionFileFolder="", minDistance=1000, regionFileId="", outputFile=""):
    if not outputFile.strip():
        outputFile = os.path.join(regionFileFolder, "mergeRegionFiles_output.bed.txt")

    chromList = [str(x) for x in range(1, 23)] + ["X", "Y"]
    regionFiles = glob.glob(os.path.join(regionFileFolder, "*.bed"))
    allFileData = []

    for regionFile in regionFiles:
        regionList = stringList2inputDataFile(removeHeader(load2stringList(regionFile)), format=["s", "i", "i", "f"])
        allFileData.append(regionList)

    allChromDictionary = {chrom: [] for chrom in chromList}
    for file in allFileData:
        for fileRow in file:
            allChromDictionary[fileRow[0]].extend(fileRow[1:3])

    for chrom in chromList:
        allChromDictionary[chrom].sort()
        n = 0
        while n < len(allChromDictionary[chrom]) - 1:
            if abs(allChromDictionary[chrom][n] - allChromDictionary[chrom][n + 1]) < minDistance:
                allChromDictionary[chrom][n] = int(
                    (allChromDictionary[chrom][n] + allChromDictionary[chrom][n + 1]) / 2)
                del allChromDictionary[chrom][n + 1]
            else:
                n += 1

    finalList = []
    for chrom in chromList:
        n = 0
        while n < len(allChromDictionary[chrom]) - 1:
            finalList.append([chrom, allChromDictionary[chrom][n], allChromDictionary[chrom][n + 1]])
            n += 1

    saveFile(outputFile, finalList, "chr\tstart\tend")

    message = f"\n\nProcessed {len(regionFiles)} files.\nGenerated {len(finalList)} subregions.\nNew file saved at {outputFile}"
    print(message)


# ------------------------------------------------------

def joinData(fileList=[], histmod="ac", verbose=False):
    histmodMap = {"ac": ("H3K27ac", 0), "H3K27ac": ("H3K27ac", 0), "me1": ("H3K4me1", 1), "H3K4me1": ("H3K4me1", 1)}
    histmod, histmodPos = histmodMap.get(histmod, (None, -1))

    myFiles = fileList
    myList = []
    if verbose: print("Loading files...")
    badFileIndices = []
    for i, myFile in enumerate(myFiles):
        badFile = False
        if verbose: print(myFile)
        with open(myFile, "r") as reader:
            myFileContents = [myFile]
            skipThis = True
            counter = 0
            while True:
                counter += 1
                myLine = reader.readline().strip()
                if not myLine: break
                if not skipThis:
                    myLine = myLine.split("\t")
                    if counter == 2:
                        if len(myLine) >= histmodPos + 1:
                            if myLine[histmodPos] != histmod:
                                badFile = True
                                break
                        else:
                            badFile = True
                            break
                    if counter > 2:
                        myValue = int(myLine[histmodPos])
                        if myValue == 2:
                            badFile = True
                            break
                        myFileContents.append(myValue)
                else:
                    skipThis = False
            if not badFile:
                myList.append(myFileContents)
            else:
                if verbose: print(f"Warning: {myFile} is a bad file, skipping.")
                badFileIndices.append(i)

    if verbose: print("Generating main data frame...")
    myListP = pd.DataFrame(myList)
    if verbose: print("Main data frame generated...")
    myListPi = myListP.set_index(0)

    return myListPi, badFileIndices


# ------------------------------------------------------

def testPrediction(groupingFile="",
                   regionFile="",
                   testSize=0.3,
                   totRandomStates=11,
                   chrom="",
                   histmod="ac",
                   verbose=False,
                   interRegionTested=True,
                   binSize=200,
                   outputFolder="",
                   output="output.txt",
                   randomize=False,
                   randomizeMethod="scramble",
                   label_seed=None,
                   RF_seed=None):
    outputFile = os.path.join(outputFolder, output)

    chromList = [f"chr{chrom}" for chrom in range(1, 23)] + ["chrX"] if not chrom else [chrom]
    medianPos = totRandomStates // 2

    if len(chromList) == 1:
        outputFile = outputFile.replace("*", chromList[0])

    if regionFile:
        regionList = stringList2inputDataFile(removeHeader(load2stringList(regionFile)), format=["s", "i", "i"])
    else:
        regionList = [["0", 0, 0]]

    if not groupingFile:
        print("A grouping file indicating the path to the files and a group identifier is needed.")
        print("If an asterisk (*) is included in the filename, it will be replaced by the chromosome.")
        print("Example (there are two tab-separated columns, and no header):")
        print("file_1.txt\tgroupA")
        print("file_2.txt\tgroupA")
        print("...")
        print("file_3.txt\tgroupB")
        print("file_4.txt\tgroupB")
        print("...")
        return

    groupingList = stringList2inputDataFile(load2stringList(groupingFile), format=["s", "s"])
    fileList = [element[0] for element in groupingList]
    sampleLabelList = [element[1] for element in groupingList]
    sampleLabelSet = list(set(sampleLabelList))

    if len(sampleLabelSet) != 2:
        print("Chromatinsight works with two groups of samples,")
        print(f"so the number of unique labels must be exactly 2.")
        print(f"In the grouping file there are {len(sampleLabelSet)}")
        print("Namely:")
        print(sampleLabelSet)
        print()
        return

    myScoreList = []
    for singleChrom in chromList:
        thisChromSampleLabelList = sampleLabelList

        fileList_chromReplaced = [singleFile.replace("*", singleChrom) for singleFile in fileList]
        myData, badFileIndices = joinData(fileList_chromReplaced, histmod=histmod, verbose=verbose)
        for i in reversed(badFileIndices):
            del thisChromSampleLabelList[i]
        if verbose: print("Data joined.")

        if randomize:
            if verbose: print("Randomising labels, as requested...")
            random.seed(label_seed)
            if randomizeMethod == "coin":
                thisChromSampleLabelList = [sampleLabelSet[random.randint(0, 1)] for _ in range(len(myData))]
            elif randomizeMethod == "scramble":
                random.shuffle(thisChromSampleLabelList)

        myScoreChrom = []
        previousRegionEnd = 0
        interTADLabel = "Starting"
        for chromRegion in regionList:
            thisChrom = f"chr{chromRegion[0]}"
            if thisChrom == singleChrom:
                regionID = f"{singleChrom}_{previousRegionEnd}-{chromRegion[1]}_{interTADLabel}"
                regionStart = previousRegionEnd // binSize
                regionEnd = chromRegion[1] // binSize

                if regionStart < regionEnd - 1 and interRegionTested:
                    if verbose: print(f"Getting patterns in inter-region {regionID}")
                    if chromRegion[2] == 0: regionEnd = len(myData.iloc[0, :])
                    regionCoordinates = f"{singleChrom}:{regionStart * binSize}-{regionEnd * binSize}"

                    thisData = myData.iloc[:, regionStart:regionEnd - 1].copy()  # Make an explicit copy
                    thisData.loc[:, "group"] = thisChromSampleLabelList
                    myScores = [getScore(thisData, testSize, randomState, RF_seed) for randomState in
                                range(totRandomStates)]
                    myScoreChrom.append([regionID] + myScores)
                    print(f"{interTADLabel} {regionCoordinates}: {myScores}, median = {sorted(myScores)[medianPos]}")
                    interTADLabel = "interTAD"

                regionID = f"{singleChrom}_{chromRegion[1]}-{chromRegion[2]}_TAD"
                regionStart = chromRegion[1] // binSize
                regionEnd = chromRegion[2] // binSize

                if regionStart < regionEnd - 1:
                    if verbose: print(f"Getting patterns in region {regionID}")
                    if chromRegion[2] == 0: regionEnd = len(myData.iloc[0, :])
                    regionCoordinates = f"{singleChrom}:{regionStart * binSize}-{regionEnd * binSize}"

                    thisData = myData.iloc[:, regionStart:regionEnd - 1].copy()  # Make an explicit copy
                    thisData.loc[:, "group"] = thisChromSampleLabelList

                    myScores = [getScore(thisData, testSize, randomState, RF_seed) for randomState in
                                range(totRandomStates)]
                    myScoreChrom.append([regionID] + myScores)
                    print(f"TAD {regionCoordinates}: {myScores}, median = {sorted(myScores)[medianPos]}")
                previousRegionEnd = chromRegion[2]

        regionStart = previousRegionEnd // binSize
        regionEnd = len(myData.iloc[0, :])
        if regionStart < regionEnd - 1 and interRegionTested:
            regionID = f"{singleChrom}_{previousRegionEnd}-{regionEnd * binSize}_Ending"
            regionCoordinates = f"{singleChrom}:{regionStart * binSize}-{regionEnd * binSize}"

            thisData = myData.iloc[:, regionStart:regionEnd - 1].copy()  # Make an explicit copy
            thisData.loc[:, "group"] = thisChromSampleLabelList

            myScores = [getScore(thisData, testSize, randomState, RF_seed) for randomState in range(totRandomStates)]
            myScoreChrom.append([regionID] + myScores)
            print(f"Ending {regionCoordinates}: {myScores}, median = {sorted(myScores)[medianPos]}")

        myScoreList.append(myScoreChrom)

    saveFile(outputFile, myScoreList[0], header="chrom_init-end_region")
    print(f"Results saved in file {outputFile}")

    return myScoreList


# ------------------------------------------------------

def getScore(myData, testSize, randomState=None, RF_seed=None):
    train_index, test_index = next(
        StratifiedShuffleSplit(test_size=testSize, random_state=randomState).split(myData, myData["group"]))
    myData_train = myData.iloc[train_index, :]
    myData_test = myData.iloc[test_index, :]

    rf = RandomForestClassifier(n_estimators=200, random_state=RF_seed)
    myDataLength = myData.shape[1] - 1
    myFit = rf.fit(myData_train.iloc[:, 0:myDataLength], myData_train["group"])
    myScore = rf.score(myData_test.iloc[:, 0:myDataLength], myData_test["group"])

    return myScore

#######################################################
