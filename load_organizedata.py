
import csv
import numpy as np
import pandas as pd
import copy
from itertools import izip as zip, count  # izip for maximum efficiency


class loadorgdata:
    def __init__(self,csvname,assayname):
        self.csvname = csvname
        self.assayname = assayname

    def load_data(self):
        """ Load data from a given CSV file"""
        self.alldata = pd.read_csv(csvname)
        self.colVals=list(alldata.columns.values)
        self.indexAssaystart1 = self.colVals.index(self.assayname)
        self.rowCount=len(newalldata)
        self.colCount=len(self.colVals)
        #
        # with open(self.csvname, 'rb') as csvfile:
        #     reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        #
        #     for rowCt, row in enumerate(reader):
        #         if rowCt == 0:
        #             self.alldata = []
        #
        #         else:
        #             self.alldata.append(row)
        #         self.colCount = len(row)
        #         self.rowCount = rowCt

    def org_data(self):
        """Put data into numpy array and lists,
        Find index where assays start: This  will only work if
        diagnosis data comes before assay data."""
        self.load_data()
        # initilize arrays for diagnosis codes and genetic test codes
        diagCodes = np.zeros((self.rowCount, self.indexAssaystart1))
        # if there is a test performed or not
        genoCodes = np.zeros((self.rowCount, len(self.colVals) - self.indexAssaystart1))
        # categorize data
        genoCodesData = np.zeros((self.rowCount, len(self.colVals) - self.indexAssaystart1))

        # get unique elements of each column
        els = []
        for test in self.colVals:
            els.append(list(set(self.alldata[test])))

        self.alldata_copy = copy.deepcopy(self.alldata)
        # running into problems here with pandas transformations to array
        # becomes problematic going into sklearn anyway
        for colNum, col in enumerate(row[0:self.indexAssaystart1 - 1]):
            if col == 'TRUE':
                diagCodes[0:, colNum] = int(col == 'TRUE')

        colNum = self.indexAssaystart1
        for col in colVals[colNum:]:
            for rowNum, row in enumerate(self.alldata):
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

        return rowNum, els, self.alldata_copy, genoCodes, genoCodesData, diagCodes,self.colVals
