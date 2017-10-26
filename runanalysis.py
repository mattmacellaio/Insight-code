from load_organizedata import loadorgdata
from ML_analysis import ml_pipeline

#load and organize data
csvname='genotype_report_080817.csv'
assayname='C__25986767_70'
lo = loadorgdata(csvname,assayname)
genoCodes, genoCodesData, diagCodes,colVals = lo.org_data()
#selecting diagnoses with more than 50 patients
diagDict,commonCodeInds,commonCodeInds_sorted = lo.get_diags(50)
#using patients where more than 100 assays were performed
#this works out to a set of assays that had each at least 222 patients
categoryData,useTests = lo.choose_data(222,genoCodes,genoCodesData)
#one hot encode assay results
testResultData,testResultLabels = lo.onehot(categoryData,useTests)

#set up and run the ML pipeline
mlp = ml_pipeline(diagCodes,diagDict,commonCodeInds,commonCodeInds_sorted,testResultData,testResultLabels)
mlp.trainml_withfigs()



