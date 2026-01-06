from io import open
import numpy as np
import const
import random

import utils


class GenAllData:
    @staticmethod
    def __convertBinStringToArray(bs):
        sz = len(bs)
        ar = np.zeros(sz, dtype=int)
        for i in range(sz):
            if bs[i] == "1":
                ar[i] = 1
        return ar
        
    def __convertNumStringToArray(self, bs):
        sz = len(bs)
        ar = np.zeros(sz, dtype=float)
        for i in range(sz):
                ar[i] = float(bs[i])
        return ar # xinya

    def paddingECEPFeatureToNumpyArray(self, features):
        numDrug = len(features)
        if numDrug == 1:
            return np.asarray(features)
        ar = np.zeros([numDrug, self.MAX_ATOMS, self.N_FEATURE])
        # for i,data in enumerate(features):
        #     for j,d in enumerate(data):
        #         ar[i][j][d] = 1
        for drugIdx in range(numDrug):
            drugData = features[drugIdx]
            nAtom, nFeature = drugData.shape

            for iAtom in range(nAtom):
                for iFeature in range(nFeature):
                    ar[drugIdx][iAtom][iFeature] = drugData[iAtom][iFeature]
        return ar

    def loadBio2RDFFeatures(self, path):
        fin = open(path)
        features = []
        while True:
            line = fin.readline()
            if line == "":
                break
            line = line.strip()
            line = line.split("|")[1]

            parts = line.split(",")
            arFeature = set()
            if len(line) > 0:
                for p in parts:
                    arFeature.add(p)

            features.append(arFeature)
        fin.close()
        return features


    def loadRDFkitFeatures(self, path):
        features = []
        with open(path) as fin:
            for line in fin:
                line = line.strip()
                if not line or "|" not in line:
                    continue  
                    
                feature_str = line.split("|")[1]    
                if feature_str:    
                    features.append([float(x) for x in feature_str.split(",")])
        return features  #xinya

    def loadECFPLiuData(self):
        ECFPFeatures = utils.load_obj(const.LIU_ECFP_PATH)
        BIO2RDFFeatures = self.loadBio2RDFFeatures(const.LIU_BIO2RDF_PATH)
        fADR = open(const.LIU_ADR_PATH)
        
        Chems = []
        ADRs = []
        while True:
            line = fADR.readline()
            if line == "":
                break
            line = line.strip()
            parts = line.split("|")
            chem = self.__convertBinStringToArray(parts[4])
            adr = self.__convertBinStringToArray(parts[3])

            Chems.append(chem)
            ADRs.append(adr)
        fADR.close()

        if len(ECFPFeatures) != len(Chems):
            print("Fatal error. Missmatched data")
            exit(-1)

        fin = open(const.LIU_INFO)
        self.N_DRUGS = int(fin.readline().split(":")[-1].strip())
        self.N_FEATURE = int(fin.readline().split(":")[-1].strip())
        self.MAX_ATOMS = int(fin.readline().split(":")[-1].strip())

        self.ECFPFeatures = ECFPFeatures
        self.Bio2RDFFeatures = BIO2RDFFeatures
        self.Chems = Chems
        self.ADRs = ADRs
        print("loading successfully")


    def loadECFPAEOLUSData(self):
        ECFPFeatures = utils.load_obj(const.AEOLUS_ECFP_PATH)
        BIO2RDFFeatures = self.loadBio2RDFFeatures(const.AEOLUS_BIO2RDF_PATH)
        fADR = open(const.AEOLUS_ADR_PATH)
        fChem = open(const.AEOLUS_CHEM_PATH)
        
        Chems = []
        ADRs = []
        while True:
            lineADR = fADR.readline()
            if lineADR == "":
                break
            lineADR = lineADR.strip()
            lineChem = fChem.readline().strip()

            adrString = lineADR.split("|")[-1].replace(",", "")

            # (1) raw feature (pubchem cids)
            # chemString = lineChem.split("|")[5].replace(",", "")[:881]
            
            # (2) ecfp 
            # chemString = lineChem.split("|")[6][:2048]
            # chemString = lineChem.split("|")[7][:2048]
            # chemString = lineChem.split("|")[8][:2048]
            
            # (3) maccs
            # chemString = lineChem.split("|")[9][:167]
            
            # (4) morgan fp
            chemString = lineChem.split("|")[10][:1024]
            
            # (5) Rdkit2D
            # chemString = lineChem.split("|")[12].split(",")
            # print(len(chemString), len(adrString))
            chem = self.__convertBinStringToArray(chemString)
            # chem = self.__convertNumStringToArray(chemString)
            adr = self.__convertBinStringToArray(adrString)

            Chems.append(chem)
            ADRs.append(adr)
        fADR.close()
        fChem.close()
        if len(ECFPFeatures) != len(Chems):
            print("Fatal error. Missmatched data")
            exit(-1)

        fin = open(const.AEOLUS_INFO)
        self.N_DRUGS = int(fin.readline().split(":")[-1].strip())
        self.N_FEATURE = int(fin.readline().split(":")[-1].strip())
        self.MAX_ATOMS = int(fin.readline().split(":")[-1].strip())

        self.ECFPFeatures = ECFPFeatures
        self.Bio2RDFFeatures = BIO2RDFFeatures
        self.Chems = Chems
        self.ADRs = ADRs
        print("loading successfully")

    def loadECFPSIDERData(self):
        ECFPFeatures = utils.load_obj(const.SIDER_ECFP_PATH)
        # BIO2RDFFeatures = self.loadBio2RDFFeatures(const.SIDER_BIO2RDF_PATH)
        BIO2RDFFeatures = self.loadRDFkitFeatures(const.SIDER_RDKIT_PATH)
        fADR = open(const.SIDER_ADR_PATH)
        fChem = open(const.SIDER_CHEM_PATH)
        
        Chems = []
        ADRs = []
        while True:
            lineADR = fADR.readline()
            lineChem = fChem.readline()
            if lineADR == "":
                break
            lineADR = lineADR.strip()
            lineChem = lineChem.strip()
            
            adrString = lineADR.split("|")[-1].replace(",", "")
            # chemString = lineChem.split("|")[3][:881]
            
            # (1) raw feature (pubchem cids)
            # chemString = lineChem.split("|")[4].replace(",", "")[:881]
            
            # (2) ecfp 
            chemString = lineChem.split("|")[4][:2048]
            # chemString = lineChem.split("|")[5][:2048]
            # chemString = lineChem.split("|")[6][:2048]
            
            # (3) maccs
            # chemString = lineChem.split("|")[7][:167]
            
            # (4) morgan fp
            # chemString = lineChem.split("|")[8][:1024]
            
            # (5) Rdkit2D
            # chemString = lineChem.split("|")[9].split(",")
            # print(len(chemString), len(adrString))
            
            
            # chemString = lineChem.split("|")[4][:2048] # xinya
            #print(len(chemString), len(adrString))
            chem = self.__convertBinStringToArray(chemString)
            adr = self.__convertBinStringToArray(adrString)

            Chems.append(chem)
            ADRs.append(adr)
        fADR.close()
        fChem.close()
        if len(ECFPFeatures) != len(Chems):
            print("Fatal error. Missmatched data")
            exit(-1)

        fin = open(const.SIDER_INFO)
        self.N_DRUGS = int(fin.readline().split(":")[-1].strip())
        self.N_FEATURE = int(fin.readline().split(":")[-1].strip())
        self.MAX_ATOMS = int(fin.readline().split(":")[-1].strip())

        self.ECFPFeatures = ECFPFeatures
        self.Bio2RDFFeatures = BIO2RDFFeatures
        self.Chems = Chems
        self.ADRs = ADRs # xinya
        print("loading successfully")

    
    def loadECFPTOXData(self):
        ECFPFeatures = utils.load_obj(const.TOX_ECFP_PATH)
        BIO2RDFFeatures = self.loadRDFkitFeatures(const.TOX_BIO2RDF_PATH)
        fADR = open(const.TOX_TOXCITY_PATH)
        fChem = open(const.TOX_CHEM_PATH)
        
        Chems = []
        ADRs = []
        while True:
            lineADR = fADR.readline()
            lineChem = fChem.readline()
            if lineADR == "":
                break
            lineADR = lineADR.strip()
            adrString = lineADR.replace(",", "")
            lineChem = lineChem.strip()
            chemString = lineChem.split("|")[2][:2048]
            #print(len(chemString), len(adrString))
            adr = self.__convertBinStringToArray(adrString)
            chem = self.__convertBinStringToArray(chemString)
            
            Chems.append(chem)
            ADRs.append(adr)
        fADR.close()
        fChem.close()
        if len(ECFPFeatures) != len(Chems):
            print("Fatal error. Missmatched data")
            exit(-1) 

        fin = open(const.TOX_INFO)
        self.N_DRUGS = int(fin.readline().split(":")[-1].strip())
        self.N_FEATURE = int(fin.readline().split(":")[-1].strip())
        self.MAX_ATOMS = int(fin.readline().split(":")[-1].strip())

        self.ECFPFeatures = ECFPFeatures
        self.Bio2RDFFeatures = BIO2RDFFeatures
        self.Chems = Chems
        self.ADRs = ADRs 
        print("loading successfully") # xinya


    def getTrainTestPathByIFold(self, ifold, root=const.KFOLD_FOLDER_EC_Liu):
        pTrainECFeature = "%s/%s_ec_%s" % (root, const.TRAIN_PREFIX_EC, ifold)
        pTrainChemFeature = "%s/%s_chem_%s" % (root, const.TRAIN_PREFIX_EC, ifold)
        pTrainADRs = "%s/%s_ADR_%s" % (root, const.TRAIN_PREFIX_EC, ifold)

        pTestECFeature = "%s/%s_ec_%s" % (root, const.TEST_PREFIX_EC, ifold)
        pTestChemFeature = "%s/%s_chem_%s" % (root, const.TEST_PREFIX_EC, ifold)
        pTestADRs = "%s/%s_ADR_%s" % (root, const.TEST_PREFIX_EC, ifold)

        pTrainBioRDFFeature = "%s/%s_biordf_%s" % (root, const.TRAIN_PREFIX_EC, ifold)
        pTestBioRDFFeature = "%s/%s_biordf_%s" % (root, const.TEST_PREFIX_EC, ifold)
        return pTrainECFeature, pTrainChemFeature, pTrainBioRDFFeature, pTrainADRs, pTestECFeature, pTestChemFeature,pTestBioRDFFeature, pTestADRs

    def exportKFold(self, root=const.KFOLD_FOLDER_EC_Liu):
        foldSize = self.N_DRUGS / const.KFOLD
        order = np.arange(0, self.N_DRUGS)
        random.seed(1)
        random.shuffle(order)

        for i in range(const.KFOLD):
            # pTrainECFeature, pTrainChemFeature, pTrainBioRDFFeature, pTrainADRs, pTestECFeature, pTestChemFeature,pTestBioRDFFeature, pTestADRs
            paths = self.getTrainTestPathByIFold(i, root)

            arTrain = []
            arTest = []

            start = i * foldSize
            end = (i + 1) * foldSize
            if i == const.KFOLD - 1:
                end = self.N_DRUGS
            for jj in range(self.N_DRUGS):
                ar = arTrain
                if start <= jj < end:
                    ar = arTest
                ix = order[jj]
                ar.append([self.ECFPFeatures[ix], self.Chems[ix], self.Bio2RDFFeatures[ix], self.ADRs[ix]])

            ars = [arTrain, arTest]

            for ii in range(2):
                for jj in range(4):
                    path = paths[ii * 4 + jj]
                    tmp = ars[ii]
                    data = []
                    for d in tmp:
                        data.append(d[jj])
                    if jj == 0 or jj == 2:
                        utils.save_obj(data, path)
                    else:
                        data = np.vstack(data)
                        np.savetxt(path, data)

    def loadFold(self, iFold, root=const.CURRENT_KFOLD):

        print("IFOLD: %s, FOLDER: %s" % (iFold, root))
        if root == const.KFOLD_FOLDER_EC_Liu:
            fin = open(const.LIU_INFO)
        elif root == const.KFOLD_FOLDER_EC_AEOLUS:
            fin = open(const.AEOLUS_INFO)
        elif root == const.KFOLD_FOLDER_EC_SIDER:
            fin = open(const.SIDER_INFO) # xinya
        elif root == const.KFOLD_FOLDER_EC_TOX:
            fin = open(const.TOX_INFO) # xinya
        else:
            print("data loading error")

        self.N_DRUGS = int(fin.readline().split(":")[-1].strip())
        self.N_FEATURE = int(fin.readline().split(":")[-1].strip())
        self.MAX_ATOMS = int(fin.readline().split(":")[-1].strip())

        paths = self.getTrainTestPathByIFold(iFold, root)
        data = []
        for i in range(2):
            for j in range(4):
                path = paths[i * 4 + j]
                print("attention: i = {}, j = {}, path = {}".format(i,j,path))
                if j == 0 or j == 2:
                    data.append(utils.load_obj(path))
                else:
                    data.append(np.loadtxt(path))

        self.N_ADRS = data[3].shape[1]

        print(self.N_DRUGS, self.N_ADRS)

        return data

    
def convertBioRDFSet2Array(sets):
    numDrugs = len(sets)
    if const.UPDATE:
        ar = np.array(sets).reshape(numDrugs, const.NUM_RDKIT2D_FEATURE) # xinya
    else: 
        ar = np.zeros((numDrugs, const.NUM_BIO2RDF_FEATURE))
        for i,vs in enumerate(sets):
            for v in vs:
                ar[i][int(v)] = 1
    return ar # xinya 

# def convertBioRDFSet2Array(sets):
#     numDrugs = len(sets)
#     if const.CURRENT_DATA == "TOX":
#         numFeatures = const.NUM_RDKIT2D_FEATURE
#     else:
#         numFeatures = const.NUM_BIO2RDF_FEATURE
#     if numFeatures == 6172: 
#         ar = np.zeros((numDrugs, const.NUM_BIO2RDF_FEATURE))
#         for i,vs in enumerate(sets):
#             for v in vs:
#                 ar[i][int(v)] = 1
#     elif numFeatures == 200:
#          ar = np.array(sets).reshape(numDrugs, numFeatures) # xinya
#     else:
#         print("x_bio data error!")
#     return ar # xinya 

# def convertBioRDFSet2Array(sets):
#     numDrugs = len(sets)
#     ar = np.zeros((numDrugs,const.NUM_BIO2RDF_FEATURE))
#     for i,vs in enumerate(sets):
#         for v in vs:
#             ar[i][int(v)] = 1
#     return ar


def genCombineData():
    dataLoader = GenAllData()
    for i in range(const.KFOLD):
        datas = dataLoader.loadFold(i)
        trainInp, trainKGInp, trainOut, testInp, testKGInp, testOut = datas[1], datas[2], datas[3], datas[5], datas[
            6], datas[7]
        trainInp2 = convertBioRDFSet2Array(trainKGInp)
        trainInp = np.concatenate([trainInp, trainInp2], axis=1)

        np.savetxt( "%s/trainInpAll_%s" % (const.RSCCA_DATA_DIR, i),trainInp)
        np.savetxt("%s/trainOutAll_%s" % (const.RSCCA_DATA_DIR, i), trainOut)


def genKFoldECFPLiu():
    data = GenAllData()
    data.loadECFPLiuData()
    data.exportKFold(const.KFOLD_FOLDER_EC_Liu)


def genKFoldECFPAEOLUS():
    data = GenAllData()
    data.loadECFPAEOLUSData()
    data.exportKFold(const.KFOLD_FOLDER_EC_AEOLUS)
    
def genKFoldECFPSIDER():
    data = GenAllData()
    data.loadECFPSIDERData()
    data.exportKFold(const.KFOLD_FOLDER_EC_SIDER) # xinya

def genKFoldECFPTOX():
    data = GenAllData()
    data.loadECFPTOXData()
    data.exportKFold(const.KFOLD_FOLDER_EC_TOX) # xinya
    
if __name__ == "__main__":    
    # Liu_Data With ECFP
    genKFoldECFPLiu()

    # AeolusData
    genKFoldECFPAEOLUS()

    # siderdata
    genKFoldECFPSIDER() # xinya

    # toxdata
    genKFoldECFPTOX() # xinya
    pass
