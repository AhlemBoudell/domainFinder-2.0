import os
import matplotlib.pyplot as plt
from sklearn import preprocessing
import shutil
import numpy as np

from apyori import apriori
import math

class Domain:
    startPos = 0
    size = 0
    E_value = 0
    all_Feat_counts = []
    all_feats = []
    uniFeatsCount = 0

    def __init__(self, SP=0, S=0, F_list=[], allFs=[]):
        self.startPos = SP
        self.size = S
        self.all_Feat_counts = F_list
        self.all_feats = allFs
        self.uniFeatsCount = 0

    def getStartPos(self):
        return self.startPos

    def getUniFeatsCount(self):
        tempStr = self.getfeatsSum()
        self.uniFeatsCount = tempStr.count("1")
        return self.uniFeatsCount

    def getAnno(self):
        return self.allAnno

    def getAnnoStr(self):
        tempStr = ""
        for a in self.allAnno:
            tempStr += a + ","
        if tempStr == "":
            tempStr = "-"
        return tempStr

    def getfeats(self):
        self.getUniFeatsCount()
        return self.all_feats

    def getfeatsSum(self):
        ftsSum = self.all_feats[0]
        for i in range(1, len(self.all_feats)):
            for j in range(0, self.all_feats[i].__len__()):
                if self.all_feats[i][j] == "1" and ftsSum[j] == "0":
                    ftsSum[j] == "1"
        return ftsSum

    def getSize(self):
        return self.size

    def getFlist(self):
        tempList = []
        for tfc in self.all_Feat_counts:
            tempList.append(int(tfc))
        return tempList

    def setDomain(self, SP, S, EL, Fs, As):
        self.startPos = SP
        self.size = S
        # self.E_value = E
        self.all_E_values = EL
        self.all_feats = Fs
        self.calculateESum()
        self.allAnno = As

    def printDomain(self):
        print("Start Pos: " + self.startPos.__str__() + "   E value: " + self.E_value.__str__())

class DomainClass:
    classSize = 0
    classGroups = []
    classGroupsCounts = []

    def __init__(self, s=0):
        self.classSize = s
        self.classGroups = []
        self.classGroupsCounts = []

    def getSize(self):
        return self.classSize

    def getGroups(self):
        return self.classGroups

    def getOneGroupGroups(self, index):
        str = ""
        self.classGroups[index].sort()
        for g in self.classGroups[index]:
            str += g + ","
        return str[0:str.__len__() - 1]

    def getGroupCounts(self):
        return self.classGroupsCounts

    def getGroupCountsStr(self):
        str = ""
        for g in self.classGroupsCounts:
            str = str + g.__str__() + ","
        return str[0:str.__len__() - 1]

    def getGroupsStr(self):
        str = ""
        for g in self.classGroups:
            subStr = ""
            for sg in g:
                subStr = subStr + sg + ","
            subStr = subStr[0:subStr.__len__() - 1]
            str = str + subStr.__str__() + "_"
        return str[0:str.__len__() - 1]

    def addGroup(self, group):
        if group.__len__() == self.classSize:
            found = False
            for i in range(0, self.classGroups.__len__()):
                if self.AreSameList(self.classGroups[i], group):
                    self.classGroupsCounts[i] += 1
                    found = True
                    break

            if not found:
                self.classGroups.append(group)
                self.classGroupsCounts.append(1)
        else:
            print("Wrong Group size!")

    def AreSameList(self, L1, L2):
        if L1.__len__() != L2.__len__():
            return False
        else:
            for item in L1:
                if not L2.__contains__(item):
                    return False
            return True

    def getDiffGroups(self):
        return self.classGroups.__len__()


# checked!
def wigFilesModifier():
    stages = ["Emb", "L3"]
    chrms = ["I"]
    WS = 200
    Directory = "/home/boudela/mountboudela/HisMod/"
    for stage in stages:
        fileDir = Directory + stage + "/"
        folders = os.listdir(fileDir)
        for chrm in chrms:
            for folder in folders:
                if folder.startswith("modEncode"):
                    print(folder)
                    newDir = fileDir + folder + "/signal_data_files/"
                    files = os.listdir(newDir)
                    for file in files:
                        if file.endswith(".wig"):
                            inputFile = newDir + file
                            outputFile = inputFile.split(".wig")[0] + "_" + chrm + ".txt"
                            # print(inputFile)
                            # print(outputFile)
                            textFileGenerator(inputFile, outputFile, chrm)

    for stage in stages:
        fileDir = Directory + stage + "/"
        folders = os.listdir(fileDir)
        for chrm in chrms:
            for folder in folders:
                if folder.startswith("modEncode"):
                    print(folder)
                    newDir = fileDir + folder + "/signal_data_files/"
                    files = os.listdir(newDir)
                    for file in files:
                        if file.endswith(".txt") and not file.__contains__("AvgPeaks"):
                            fname = newDir + file
                            avgFile = fname.split(".txt")[0] + "AvgPeaks.txt"
                            avgPeakFileGenerator(fname, avgFile, WS)
# checked!
def textFileGenerator(infile, outfile, chrm):
    inData = open(infile, "r")
    outData = open(outfile, "w")
    chrmstr = "chr" + chrm
    # print(infile)
    hline = inData.readline()
    line = inData.readline().strip()
    flag = False
    span = 0

    outData.write("sp \t ep \t peakValue \n")

    while line != "":
        if line.__contains__("chrom"):
            lparts = line.split("chrom=")[1]
            currentChrm = lparts.split()[0]
            if chrm == currentChrm or chrmstr == currentChrm:
                flag = True
                span = int(lparts.split("=")[1])
                line = inData.readline().strip()
            else:
                flag = False
                span = 0
                line = inData.readline().strip()

        elif line.__contains__("span"):
            lparts = line.split("span=")[1]
            span = int(lparts.split()[0])
            line = inData.readline().strip()

        while line.split().__len__() == 2 and flag:
            lparts = line.split()
            sp = int(lparts[0])
            ep = sp + span - 1
            outData.write(sp.__str__() + "\t" + ep.__str__() + "\t" + lparts[1] + "\n")
            line = inData.readline().strip()
        while (line.split().__len__() == 2 and not flag) or (line.startswith("#")) or (line.startswith("track")):
            line = inData.readline().strip()

    inData.close()
    outData.close()
# checked!
def avgPeakFileGenerator(infile, outfile, WS):
    inData = open(infile, "r")
    outData = open(outfile, "w")
    outData.write("WindowSize=" + WS.__str__() + "\n")
    outData.write("WindowSP \t AvgPeakValue \n")
    # print(infile)
    hline = inData.readline()
    line = inData.readline().strip()
    sp = 0
    ep = sp + WS
    count = 0
    totalPeaks = 0
    exitFlag = False
    while line != "":
        lparts = line.split()
        tempSP = int(lparts[0])
        tempEP = int(lparts[1])
        while (tempSP >= sp and tempSP < ep) or (tempEP >= sp and tempEP < ep) or (tempSP < sp and tempEP >= ep):
            count += 1
            totalPeaks += float(lparts[2])
            if tempEP >= ep:
                break
            else:
                line = inData.readline().strip()
                if line == "":
                    exitFlag = True
                    break
                else:
                    lparts = line.split()
                    tempSP = int(lparts[0])
                    tempEP = int(lparts[1])
        if exitFlag:
            break
        else:
            if count > 0:
                outData.write(sp.__str__() + " \t " + float(totalPeaks / count).__str__() + "\n")
            else:
                outData.write(sp.__str__() + " \t " + "-" + "\n")

            sp = ep
            ep = sp + WS
            count = 0
            totalPeaks = 0

    if count > 0:
        outData.write(sp.__str__() + " \t " + float(totalPeaks / count).__str__() + "\n")
    else:
        outData.write(sp.__str__() + " \t " + "-" + "\n")

    inData.close()
    outData.close()
# checked!
def appendHMdata(oldFile, newFile, hmFile, hmName):
    oldReader = open(oldFile, "r")
    hmReader = open(hmFile, "r")
    newWriter = open(newFile, "w")

    for i in range(0, 2):
        oldLine = oldReader.readline()
        hmLine = hmReader.readline()
        newWriter.write(oldLine)

    oldLine = oldReader.readline()
    hmLine = hmReader.readline()

    while oldLine != "":
        if hmLine.__contains__("-"):
            newWriter.write(oldLine)
        else:
            hmSig = hmLine.strip().split()[1]
            hmSigf = float(hmSig)
            if hmSigf < 0:
                newWriter.write(oldLine)
            else:
                tempLine = ""
                oldparts = oldLine.strip().split()
                if oldparts.__len__() == 1:
                    tempLine = oldLine.strip() + "\t" + hmName + ":" + hmSig + ","
                else:
                    tempLine = oldLine.strip() + hmName + ":" + hmSig + ","
                newWriter.write(tempLine + "\n")

        oldLine = oldReader.readline()
        hmLine = hmReader.readline()
        if hmLine == "":
            break

    while oldLine != "":
        newLine = newWriter.write(oldLine)
        oldLine = oldReader.readline()

    oldReader.close()
    hmReader.close()
    newWriter.close()
# checked!
def peakCalling():
    stages = ["Emb", "L3"]
    chrms = ["I"]
    WS = 200

    for stage in stages:
        fileDir = "/home/boudela/mountboudela/HisMod/" + stage + "/testData/"
        folders = os.listdir(fileDir)
        for folder in folders:
            if chrms.__contains__(folder.split("=")[1]):
                newDir = fileDir + folder + "/"
                fname = "Chr_" + folder.split("=")[1] + "_Peaks.txt"
                outData = open(newDir + fname, "w")
                chrmLenFile = open("/home/boudela/mountboudela/chrLen.txt", "r")
                found = False
                chrmLen = 0
                line = chrmLenFile.readline()
                while line != "":
                    lparts = line.split()
                    if lparts[0] == folder.split("=")[1]:
                        found = True
                        chrmLen = int(lparts[1])
                        break
                    line = chrmLenFile.readline()
                chrmLenFile.close()
                if found:
                    outData.write("WindowSize=" + WS.__str__() + "\n")
                    outData.write("StartPosition \t HM:AvgPeak," + "\n")
                    windowSP = 0
                    while windowSP < chrmLen:
                        outData.write(windowSP.__str__() + " \n")
                        windowSP += WS
                    outData.close()
                    files = os.listdir(newDir)
                    for file in files:
                        if file != fname:
                            ''' 
                                1. get HM name
                                2. open a copy file
                                3. read from original, print into file header lines + append new HM info
                                4. delete original file, rename the copy
                            '''
                            hmName = getHMName(file)
                            oldFile = newDir + fname
                            newFile = newDir + "tempFile.txt"
                            HMFile = newDir + file
                            print(hmName)
                            appendHMdata(oldFile, newFile, HMFile, hmName)
                            os.remove(oldFile)
                            os.rename(newFile, oldFile)
# checked!
def standarizeDataFiles():
    stages = ["L3", "Emb"]
    chrms = ["I"]
    WS = 200
    Directory = "/home/boudela/mountboudela/HisMod/"
    for stage in stages:
        fileDir = Directory + stage + "/"
        folders = os.listdir(fileDir)
        for chrm in chrms:
            for folder in folders:
                if folder.startswith("modEncode"):
                    print(folder)
                    newDir = fileDir + folder + "/signal_data_files/"
                    files = os.listdir(newDir)
                    for file in files:
                        if file.endswith("AvgPeaks.txt"):
                            inputFile = newDir + file
                            outputFile = inputFile.split(".txt")[0] + "_Stand.txt"
                            xValues = []
                            yValues = []
                            dashIndecies = []

                            inData = open(inputFile, "r")
                            outData = open(outputFile, "w")
                            header = inData.readline()
                            outData.write(header)
                            header = inData.readline()
                            outData.write(header)
                            line = inData.readline().strip()
                            while line != "":
                                lparts = line.split()
                                xValues.append(lparts[0])
                                if lparts[1] == "-":
                                    dashIndecies.append(xValues.__len__() - 1)
                                else:
                                    yValues.append(float(lparts[1]))
                                line = inData.readline().strip()

                            standardized_Y = preprocessing.scale(yValues)

                            j = 0
                            for i in range(0, xValues.__len__()):
                                if dashIndecies.__contains__(i):
                                    outData.write(xValues[i] + "\t - \n")
                                else:
                                    outData.write(xValues[i] + "\t" + standardized_Y[j].__str__() + "\n")
                                    j += 1

                            inData.close()
                            outData.close()
# checked!
def getHMName(fname):
    hmName = ""
    index = 0
    if fname.upper().__contains__("H2"):
        index = fname.upper().index("H2")
    elif fname.upper().__contains__("H3"):
        index = fname.upper().index("H3")
    elif fname.upper().__contains__("H4"):
        index = fname.upper().index("H4")
    else:
        print("ERROR!!")
        print(fname)
        return "Extra"

    if fname.upper().__contains__("AC"):
        indexLast = fname.upper().index("AC")
        hmName = fname[index:indexLast + 2]
        hmName = hmName.split("_")[0]
    elif fname.upper().__contains__("ME"):
        indexLast = fname.upper().index("ME")
        hmName = fname[index:indexLast + 3]
        hmName = hmName.split("_")[0]
    elif fname.upper().__contains__("TETRA"):
        indexLast = fname.upper().index("TETRA")
        hmName = fname[index:indexLast + 5]
        hmName = hmName.split("_")[0]
        # print(fname)
    else:
        hmName = fname[index:index + 2]
        hmName = hmName.split("_")[0]
    # print(hmName)
    hmName = hmName.split(":")[0]
    return hmName
# checked!
def avgReplicas():
    stages = ["L3", "Emb"]
    chrms = ["I"]
    WS = 200
    Directory = "/home/boudela/mountboudela/HisMod/"
    for stage in stages:
        fileDir = Directory + stage + "/"
        folders = os.listdir(fileDir)
        for chrm in chrms:
            for folder in folders:
                if folder.startswith("modEncode"):
                    # print(folder)
                    newDir = fileDir + folder + "/signal_data_files/"
                    files = os.listdir(newDir)
                    Windows = []
                    wFlag = False
                    SignalLists = []
                    h1 = ""
                    h2 = ""
                    fcount = 0
                    hm = ""
                    modCode = ""
                    for file in files:
                        if file.endswith("Stand.txt"):
                            fcount += 1
                        elif file.endswith("AvgFilesData.txt"):
                            os.remove(newDir + file)
                    files = os.listdir(newDir)
                    if fcount > 1:
                        for file in files:
                            if file.endswith("Stand.txt"):
                                if hm == "" or hm == "Extra":
                                    hm = getHMName(file)
                                    modCode = file.split("_")[0]
                                    print(modCode)
                                    print(hm)
                                inData = open(newDir + file, "r")

                                h1 = inData.readline()
                                h2 = inData.readline()
                                Ys = []
                                line = inData.readline().strip()
                                while line != "":
                                    lparts = line.split()
                                    if not wFlag:
                                        Windows.append(lparts[0])

                                    Ys.append(lparts[1])
                                    line = inData.readline().strip()
                                SignalLists.append(Ys)
                                inData.close()
                                wFlag = True

                        outData = open(newDir + hm + "_" + modCode + "_AvgFilesData.txt", "w")
                        outData.write(h1)
                        outData.write(h2)
                        # making all lists of same size
                        maxL = 0
                        for i in range(0, SignalLists.__len__()):
                            if SignalLists[i].__len__() > maxL:
                                maxL = SignalLists[i].__len__()

                        for i in range(0, SignalLists.__len__()):
                            while SignalLists[i].__len__() < maxL:
                                SignalLists[i].append("-")

                        while Windows.__len__() < maxL:
                            prevValue = float(Windows[Windows.__len__() - 1])
                            Windows.append((prevValue + WS).__str__())

                        for i in range(0, maxL):
                            tempSig = 0
                            count = 0
                            for j in range(0, SignalLists.__len__()):
                                if SignalLists[j][i] != "-":
                                    tempSig += float(SignalLists[j][i])
                                    count += 1
                            if count > 0:
                                outData.write(Windows[i] + "\t" + (tempSig / count).__str__() + "\n")
                            else:
                                outData.write(Windows[i] + "\t - \n")

                        outData.close()
# checked!
def combineHMFiles():
    WS = 200
    stages = ["L3", "Emb"]
    chr = "I"
    for stage in stages:
        maxFts = 0
        print(stage)
        ODir = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/"
        ofilename = ODir + "listOfHMpeaks.txt"
        fileW = open(ofilename, "w")
        fileW.write("WindowSize= " + WS.__str__() + "\n")
        fileW.write("Window \t HMCount \t HM:Peak, \n")
        tempW = 0
        chrmLenFile = open("/home/boudela/mountboudela/chrLen.txt", "r")
        found = False
        chrmLen = 0
        line = chrmLenFile.readline()
        while line != "":
            lparts = line.split()
            if lparts[0] == chr:
                found = True
                chrmLen = int(lparts[1])
                break
            line = chrmLenFile.readline()
        chrmLenFile.close()
        if found:
            while tempW <= (chrmLen + (WS)):
                fileW.write(tempW.__str__() + " \t 0 \n")
                tempW += WS

        fileW.close()

        folders = os.listdir(ODir)
        for folder in folders:
            # print(folder)
            if folder.startswith("H"):
                newDir = ODir + folder + "/"
                files = os.listdir(newDir)
                if files.__len__() > 0:
                    inputR = open(newDir + files[0], "r")
                    sigR = open(ofilename, "r")
                    outputF = open(ODir + "tempFile.txt", "w")
                    hline = inputR.readline()
                    hline = inputR.readline()
                    sigLine = sigR.readline()
                    outputF.write(sigLine)
                    sigLine = sigR.readline()
                    outputF.write(sigLine)

                    sigLine = sigR.readline().strip()
                    hline = inputR.readline().strip()

                    while hline != "":
                        pvalue = hline.split()[1]
                        if pvalue != "-":
                            pvalue = pvalue[0:7]
                            pvalue = float(pvalue)
                            if pvalue > 3:
                                sigLineParts = sigLine.split()
                                ftsCount = int(sigLineParts[1])
                                ftsCount += 1
                                if ftsCount > maxFts:
                                    maxFts = ftsCount

                                if sigLine.split().__len__() > 2:
                                    outputF.write(sigLineParts[0] + " \t " + ftsCount.__str__() + " \t " + sigLineParts[
                                        2] + folder + ":" + pvalue.__str__() + ",\n")
                                else:
                                    outputF.write(sigLineParts[
                                                      0] + " \t " + ftsCount.__str__() + " \t " + folder + ":" + pvalue.__str__() + ",\n")
                            else:
                                outputF.write(sigLine + "\n")
                        else:
                            outputF.write(sigLine + "\n")

                        sigLine = sigR.readline().strip()
                        hline = inputR.readline().strip()
                    while sigLine != "":
                        outputF.write(sigLine + "\n")
                        sigLine = sigR.readline().strip()

                    sigR.close()
                    inputR.close()
                    outputF.close()
                    os.remove(ofilename)
                    os.rename(ODir + "tempFile.txt", ofilename)
        print(maxFts)
# checked!
def groupHMfiles():
    stages = ["L3", "Emb"]
    chrms = ["I"]
    WS = 200

    Directory = "/home/boudela/mountboudela/HisMod/"

    for stage in stages:
        fileDir = Directory + stage + "/"
        folders = os.listdir(fileDir)
        for chrm in chrms:
            for folder in folders:
                if folder.startswith("modEncode"):

                    newDir = fileDir + folder + "/signal_data_files/"
                    files = os.listdir(newDir)
                    for file in files:
                        if file.endswith("AvgFilesData.txt"):
                            print(folder)
                            hmName = getHMName(file.upper())
                            hmDir = fileDir + "HMfiles/" + hmName + "/"
                            if not os.path.exists(hmDir):
                                os.mkdir(hmDir)

                            # os.rename(newDir+file, hmDir+file)
                            shutil.copyfile(newDir + file, hmDir + file)
# checked!
def hmList():
    WS = 200
    stages = ["L3", "Emb"]
    allLists = []
    allHMs = []
    for stage in stages:
        print(stage)
        names = []
        ODir = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/"
        fileW = open(ODir + "listOfHmNames.txt", "w")
        files = os.listdir(ODir)
        for file in files:
            if file.startswith("H"):
                names.append(file.upper())
                if not allHMs.__contains__(file.upper()):
                    allHMs.append(file.upper())
        names.sort()
        fileW.write("hmCount=" + names.__len__().__str__() + "\n")

        for name in names:
            fileW.write(name + "\n")

        allLists.append(names)
        fileW.close()

    print(allLists.__len__())
    print(allLists[0].__len__())
    commonHMs = []

    for hm in allHMs:
        flag = True
        for tempList in allLists:
            if not tempList.__contains__(hm):
                flag = False
                break
        if flag:
            commonHMs.append(hm)
    commonHMs.sort()
    print(commonHMs.__len__())
    print(commonHMs)

    stages = ["L3", "Emb"]
    for stage in stages:
        print(stage)
        ODir = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/"
        fileW = open(ODir + "listOfCommonHmNames.txt", "w")
        for name in commonHMs:
            fileW.write(name + "\n")
        fileW.close()
# checked!
def maxDomainWithNoGaps():
    WS = 200
    stages = ["L3", "Emb"]
    chr = "I"
    for stage in stages:
        print(stage)
        ifile = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMpeaks.txt"
        ofile = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/maxDomains.txt"

        file = open(ifile, 'r')
        result = open(ofile, 'w')

        header = file.readline()
        temp = file.readline()

        result.write(header)
        result.write(temp)

        maxNoGapList = []
        maxLenNoGapDomain = Domain()
        maxHeightNoGapDomain = Domain()
        line = file.readline()
        while (line != ''):
            lparts = line.strip().split()
            if len(lparts) > 2:
                StartPos = lparts[0]
                domainList = []
                domainListFeats = []
                domainSigList = []
                Alist = []

                while (int(lparts[1]) > 0):
                    domainList.append(int(lparts[1]))
                    fts = []
                    sigs = []
                    temp = lparts[2].split(",")
                    temp = temp[0:temp.__len__() - 1]
                    for t in temp:
                        tparts = t.split(":")
                        fts.append(tparts[0])
                        sigs.append(float(tparts[1]))
                    domainListFeats.append(fts)
                    domainSigList.append(sigs)

                    line = file.readline()
                    if (line == ''):
                        break
                    else:
                        lparts = line.strip().split()

                tempDomain = Domain(StartPos, domainList.__len__(), domainSigList, domainList, domainListFeats,
                                    Alist)
                maxNoGapList.append(tempDomain)

            else:
                line = file.readline()

        maxNoGapList.sort(key=lambda x: x.size, reverse=True)

        result.write("Total Number of Domains = " + maxNoGapList.__len__().__str__() + "\n")
        result.write("StartPosition \t DomainSize(inWindows) \t [classSize|ClassCombinations] \t FeatureAnnotations \n")
        for i in range(0, maxNoGapList.__len__()):
            # print("StartPosition: " + maxNoGapList[i].getStartPos().__str__() + "\t Size: " + maxNoGapList[i].getSize().__str__() + "\n" )
            # strCombs = FeatureGroupsFinderForDomains(maxNoGapList[i].getfeats())
            result.write(maxNoGapList[i].getStartPos().__str__() + "\t" + maxNoGapList[i].getSize().__str__() + "\t" +
                         maxNoGapList[i].getfeats().__str__() + "\t" + maxNoGapList[i].getAnnoStr() + "\n")
            # result.write(maxNoGapList[i].getStartPos().__str__() + "\t" + maxNoGapList[i].getSize().__str__() +"\t"+ strCombs+"\t"+ maxNoGapList[i].getAnnoStr()+"\n")
            # result.write(maxNoGapList[i].getStartPos().__str__() + "\t" + maxNoGapList[i].getSize().__str__() +"\t" +maxNoGapList[i].getFlist().__str__()+"\t"+ strCombs+"\n")

        file.close()
        result.close()
# checked!
def Apriori():
    WS = 200
    stages = ["L3", "Emb"]
    for stage in stages:
        print(stage)
        ODir = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/"
        fileW = open(ODir + "AprioriResults.txt", "w")
        records = []
        files = os.listdir(ODir)
        for file in files:
            if file.endswith("listOfHMpeaks.txt"):
                fileR = open(ODir + "listOfHMpeaks.txt", "r")
                tempLine = fileR.readline()
                tempLine = fileR.readline()
                tempLine = fileR.readline().strip()
                while tempLine != "":
                    tempRec = []
                    lparts = tempLine.split()
                    if lparts.__len__() > 2:
                        hms = lparts[2].split(",")
                        line = ""
                        for hm in hms:
                            if hm != "":
                                hm = hm.split(":")[0]
                                tempRec.append(hm)
                    # else:  # no hms in window
                    if tempRec.__len__() > 0:
                        records.append(tempRec)

                    tempLine = fileR.readline().strip()
        print(records.__len__())
        # for r in records:
        #    print(r)
        print("============")

        association_rules = apriori(records, min_support=0.0045, min_confidence=0.2, min_lift=3, min_length=2)
        association_results = list(association_rules)
        print(len(association_results))
        print(association_results[0])
        print("============")

        for item in association_results:
            # first index of the inner list
            # Contains base item and add item
            pair = item[0]
            items = [x for x in pair]
            # print("Rule: " + items[0] + " -> " + items[1])
            fileW.write("Rule: " + str(list(item.ordered_statistics[0].items_base)) + " -> " + str(
                list(item.ordered_statistics[0].items_add)) + "\n")

            # second index of the inner list
            fileW.write("Support: " + str(item[1]) + "\n")

            # third index of the list located at 0th
            # of the third index of the inner list

            fileW.write("Confidence: " + str(item[2][0][2]) + "\n")
            fileW.write("Lift: " + str(item[2][0][3]) + "\n")
            fileW.write("=====================================" + "\n")

        fileW.close()
# checked!
def categoriesDef():
    WS = 200
    stages = ["L3", "Emb"]
    chr = "I"
    for stage in stages:
        print(stage)
        hmfile = open("/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHmNames.txt", "r")
        inputfile = open("/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMpeaks.txt", "r")
        outputfile = open("/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMpeaks_Cat.txt", "w")

        HMsAry = []
        HMsStr = ""
        line = hmfile.readline().strip()
        line = hmfile.readline().strip()
        while line != "":
            HMsAry.append(line)
            HMsStr = HMsStr + line + " "
            line = hmfile.readline().strip()

        hmfile.close()
        outputfile.write("HMs=" + HMsStr + "\n")
        line = inputfile.readline().strip()
        outputfile.write(line + "\n")
        line = inputfile.readline().strip()
        outputfile.write("Window \t HMCount \t HMCatHex  \tHMCat" + "\n")
        line = inputfile.readline().strip()
        while line != "":
            lparts = line.split()
            ftsCount = int(lparts[1])
            if ftsCount == 0:
                outputfile.write(line + "\n")
            else:
                binStr = [0] * HMsAry.__len__()
                count = 0
                fts = lparts[2].split(",")
                for f in fts:
                    if f != "":
                        hm = f.split(":")[0]
                        ind = HMsAry.index(hm)
                        binStr[ind] = 1
                        count += 1
                tempstr = ""
                for b in binStr:
                    tempstr += b.__str__()

                intValue = int(tempstr, 2)
                hexV = hex(intValue)

                # catV = count.__str__()+ ":" +intValue.__str__()+ ":" + tempstr
                catV = hexV[2:] + "\t\t" + tempstr
                # catV = tempstr
                outputfile.write(lparts[0] + " \t " + lparts[1] + "\t\t " + catV + "\n")

            line = inputfile.readline().strip()
# checked!
def isPartialCat(subStr, origStr):
    for i in range(0, subStr.__len__()):
        if subStr[i] == "1" and origStr[i] == "0":
            return False
    return True
# checked!
def listCategories():
    stages = ["L3", "Emb"]
    for stage in stages:
        print(stage)
        inputfile = open("/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMpeaks_Cat.txt", "r")
        outputfile = open("/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMCategories.txt", "w")
        allCats = []
        catsCount = []
        catsPartialCount = []

        line = inputfile.readline().strip()
        outputfile.write(line + "\n")
        line = inputfile.readline().strip()
        line = inputfile.readline().strip()
        line = inputfile.readline().strip()

        while line != "":
            if line.split().__len__() > 2:
                cat = line.split()[3]
                if not allCats.__contains__(cat):
                    allCats.append(cat)
                    catsCount.append(1)
                    catsPartialCount.append(0)
                else:
                    ind = allCats.index(cat)
                    catsCount[ind] += 1

            line = inputfile.readline().strip()

        inputfile.close()

        for i in range(0, allCats.__len__()):
            for j in range(0, allCats.__len__()):
                if i != j and isPartialCat(allCats[i], allCats[j]):
                    catsPartialCount[i] += catsCount[j]

        catList = []
        for c in range(0, allCats.__len__()):
            ftsCount = allCats[c].count("1")
            catList.append([allCats[c], ftsCount, catsCount[c], catsPartialCount[c]])

        catList.sort(key=lambda x: (x[1], x[2] + x[3]), reverse=True)

        outputfile.write("CategoriesCount=" + catList.__len__().__str__() + "\n")
        outputfile.write("Category \t\t\tftsCount \t recurrence \t subRecrrences \t CatInHex " + "\n")

        for c in range(0, catList.__len__()):
            intValue = int(catList[c][0], 2)
            hexV = hex(intValue)
            outputfile.write(
                catList[c][0] + "\t" + catList[c][1].__str__() + "\t\t " + catList[c][2].__str__() + "\t\t " +
                catList[c][3].__str__() + "\t\t " + hexV[2:] + "\n")

        outputfile.close()
# checked!
def maxDomainWithNoGapsCat():
    minDomSize = 2
    minHMCount = 2
    stages = ["L3", "Emb"]
    for stage in stages:
        print(stage)
        ifile = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/listOfHMpeaks_Cat.txt"
        ofile = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/maxDomains_Cat.txt"
        cfile = "/home/boudela/mountboudela/HisMod/" + stage + "/HMfiles/maxDomains_CatCombos.txt"

        file = open(ifile, 'r')
        result = open(ofile, 'w')
        combo = open(cfile, 'w')

        header = file.readline()

        key = header.split("=")[1].split()
        temp = file.readline()

        result.write(header)
        result.write(temp)
        combo.write(header)
        combo.write(temp)
        header = file.readline()

        maxNoGapList = []
        maxLenNoGapDomain = Domain()
        maxHeightNoGapDomain = Domain()
        line = file.readline()
        while (line != ''):
            lparts = line.strip().split()
            if len(lparts) > 2:
                StartPos = lparts[0]
                domainList = []
                domainListFeats = []

                while (int(lparts[1]) > 0):
                    domainList.append(int(lparts[1]))
                    fts = lparts[3]
                    domainListFeats.append(fts)
                    line = file.readline()
                    if (line == ''):
                        break
                    else:
                        lparts = line.strip().split()

                tempDomain = Domain(StartPos, domainList.__len__(), domainList, domainListFeats)
                c = tempDomain.getUniFeatsCount()
                maxNoGapList.append(tempDomain)

            else:
                line = file.readline()

        maxNoGapList.sort(key=lambda x: (x.uniFeatsCount, x.size), reverse=True)

        tempList = []
        for i in range(0, maxNoGapList.__len__()):
            if maxNoGapList[i].getUniFeatsCount() >= minHMCount and maxNoGapList[i].getSize() >= minDomSize:
                tempList.append(maxNoGapList[i])

        result.write("Total Number of Domains = " + tempList.__len__().__str__() + "\n")
        result.write("StartPosition \t DomainSize(W) \t HMCount \t [ClassCombinations]\n")
        combo.write("Total Number of Domains = " + tempList.__len__().__str__() + "\n")
        combo.write("StartPosition \t DomainSize(W) \t HMCount \t [ClassCombinations]\n")
        for i in range(0, tempList.__len__()):
            result.write(tempList[i].getStartPos().__str__() + "\t" + tempList[i].getSize().__str__() + "\t\t" +
                         tempList[i].getUniFeatsCount().__str__() + "\t" + tempList[i].getfeats().__str__() + "\n")
            combo.write(tempList[i].getStartPos().__str__() + "\t" + tempList[i].getSize().__str__() + "\t\t" +
                        tempList[i].getUniFeatsCount().__str__() + "\t" + tempList[i].getfeatsSum() + "\t\t" + binToCat(
                key, tempList[i].getfeatsSum()) + "\n")

        file.close()
        result.close()
        combo.close()
# checked!
def binToCat(catKey, binStr):
    catStr = ""
    for i in range(0, binStr.__len__()):
        if binStr[i] == "1":
            catStr += catKey[i] + ","

    return catStr[0:catStr.__len__() - 1]

def main():
    wigFilesModifier()
    standarizeDataFiles()
    groupHMfiles()
    avgReplicas()
    combineHMFiles()
    hmList()
    peakCalling()
    #Apriori()
    maxDomainWithNoGaps()
    categoriesDef()
    listCategories()
    maxDomainWithNoGapsCat()