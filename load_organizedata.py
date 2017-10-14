
import csv
import numpy as np
import copy
from itertools import izip as zip, count  # izip for maximum efficiency


class loadorgdata:
    def __init__(self,csvname,assayname):
        self.csvname = csvname
        self.assayname = assayname

    def load_data(self):
        """ Load data from a given CSV file"""

        with open(self.csvname, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')

            for rowCt, row in enumerate(reader):
                if rowCt == 0:
                    self.alldata = []
                    self.colVals = row
                    self.indexAssaystart1 = row.index(self.assayname)

                else:
                    self.alldata.append(row)
                self.colCount = len(row)
                self.rowCount = rowCt

    def org_data(self):
        """Put data into numpy array and lists,
        Find index where assays start: This  will only work if
        diagnosis data comes before assay data."""
        self.load_data()
        # initilize arrays for diagnosis codes and genetic test codes
        self.diagCodes = np.zeros((self.rowCount, self.indexAssaystart1))
        # if there is a test performed or not
        genoCodes = np.zeros((self.rowCount, len(self.colVals) - self.indexAssaystart1))
        # categorize data
        genoCodesData = np.zeros((self.rowCount, len(self.colVals) - self.indexAssaystart1))

        # transpose list to access column data
        categoryData = map(list, zip(*self.alldata))

        # get unique elements of each column
        els = []
        for test in categoryData:
            els.append(list(set(test)))


        self.alldata_copy = copy.deepcopy(self.alldata)
        for rowNum, row in enumerate(self.alldata):
            for colNum, col in enumerate(row[0:self.indexAssaystart1 - 1]):
                if col == 'TRUE':
                    self.diagCodes[rowNum, colNum] = int(col == 'TRUE')

            colNum = self.indexAssaystart1
            for col in row[self.indexAssaystart1:]:
                if col == 'none':
                    genoCodes[rowNum, colNum - self.indexAssaystart1] = 0
                    self.alldata_copy[rowNum][colNum] = ''
                else:
                    genoCodes[rowNum, colNum - self.indexAssaystart1] = 1
                    for elnum, el in enumerate(els[colNum]):
                        if col == el:
                            genoCodesData[rowNum, colNum - self.indexAssaystart1] = elnum
                            #remove any data where the assay result is inconclusive
                            check = ['UND', 'N/A', 'INV']
                            if col in check:
                                genoCodesData[rowNum, colNum - self.indexAssaystart1] = 0
                                genoCodes[rowNum, colNum - self.indexAssaystart1] = 0
                                self.alldata_copy[rowNum][colNum] = ''
                try:
                    ind = els[colNum].index('none')
                    els[colNum][ind] = ['']
                except:
                    pass

                colNum += 1

        return rowNum, els, self.alldata_copy, genoCodes, genoCodesData, self.diagCodes,self.colVals

    def get_diags(self, diagThreshold):
        # looking at diagnosis codes
        commonCodes = []
        commonCodeinds = np.where(sum(self.diagCodes) > diagThreshold)[0]
        commonCodeinds_sorted = sorted(range(len(commonCodeinds)), key=sum(self.diagCodes)[commonCodeinds].__getitem__)
        for i in commonCodeinds_sorted:
            commonCodes.append(colVals[commonCodeinds[i]])

        # in reverse order with >50 diagnoses, so last is most observed in initial dataset

        commonDiags = ['Radiculopathy', 'Vitamin D deficiency', 'Hypercholesterolemia',
                       'Anxiety disorder', 'Fatigue', 'Low back pain','Mixed hyperlipidemia',
                       'Depressive disorder episode', 'Cocaine abuse', 'Long term opioid use',
                       'Opioid abuse with intoxication', 'Alcohol abuse', 'Hyperlipidemia NOS',
                       'Hypothyroidism','Opioid abuse','Esophageal reflux', 'Stimulant dependence',
                       'Diabetes', 'Anxiety', 'Chromosomal anomaly','Trichomoniasis',
                       'Sedative dependence', 'Opioid dependence', '', 'Candidiasis',
                       'STD screening', 'Hypertension','Vaginitis', 'Long-term drug therapy',
                       'Adverse effect of drugs or medication']
        groupNames = ['pain/opioid', 'std', 'psychiatric', 'alcohol', 'stimulant/cocaine', 'weight-related',
                      'back pain', 'other']

        # group indices for each diagnosis, in the same order as above
        diagGroups = [6, 7, 5, 2, 7, 6, 5, 2, 4, 0, 0, 3, 5, 5, 0, 7, 4, 5, 2, 7, 1, 0, 0, 7, 1, 1, 5, 1, 0, 0]
        codeColumn = [colVals.index(i) for i in commonCodes]
        diagDict = {'Diagnosis': commonDiags, 'Code': commonCodes, 'Group': diagGroups, 'GroupNames': groupNames,
                    'CodeColumn': codeColumn}
        return diagDict, commonCodeinds_sorted
