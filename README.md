# Insight project: ForeCGAT
Prediction of medical diagnosis from genetic assays

## Dependencies, etc.
Written in Python 2

## Structure and use
Class loadorgdata(csvname, assayname)
1. load_data
-Called by org_data.
-Takes a CSV file in the current directory, organized with diagnoses and genetic assays as column headers and patients as rows. This is currently written to take a file organized such that all diagnosis columns come first. 
-assayname(name of the first assay column) is used to separate diagnoses from assays.
-Returns indexassaystart1 (index of first assay column).

2. org_data
-Separates genetic data from diagnosis data.
-Returns numpy arrays genoCodes (if a test was performed for each patient), genoCodesData (test results as categories), diagCodes (diagnosis codes), and colVals (column headers).

3. get_diags(diagThreshold)
-Selects diagnoses with more than diagThreshold patients. 
-Includes hard-coded diagnoses for the testing dataset with more than 50 patients, in ascending order (commonDiags). 
-Creates dictionary diagDict with commonDiags, commonCodes (the ICD10 codes for those diagnoses), diagGroups (diagnoses roughly grouped by type), groupNames, codeColumn (column number of each code relative to original structure). 
-Returns diagDict, commonCodeInds_sorted (column number of each code, sorted by number of diagnosed patients, ascending).

4. choose_data(ptThreshold, genoCodes, genoCodesData)
-Selects assays that had been performed on more than ptThreshold patients. 
-Returns categoryData (compressed version of genetic assay data), useTests (list of assay indices used, relative to genoCodesData array).

5. onehot(categoryData,useTests)
-One hot encodes genetic data as testResultData
-Generates label names from assay name and result as testResultLabels
-Returns testResultData, testResultLabels.







