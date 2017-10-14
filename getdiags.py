from load_organizedata import loadorgdata
csvname='genotype_report_080817.csv'
assayname='C__25986767_70'
lo=loadorgdata(csvname,assayname)

rowNum, els, alldata_copy, genoCodes, genoCodesData, diagCodes,colVals=lo.org_data()
