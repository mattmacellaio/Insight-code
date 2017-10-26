# lasso
import numpy as np
import matplotlib as plt
from sklearn.linear_model import Lasso
from sklearn.metrics import roc_curve, auc, r2_score
import random

class ml_pipeline:
    def __init__(self,diagCodes,diagDict,commonCodeInds,commonCodeInds_sorted,testResultData,testResultLabels):
        self.diagCodes = diagCodes
        self.diagDict = diagDict
        self.commonCodeInds = commonCodeInds
        self.commonCodeInds_sorted = commonCodeInds_sorted
        self.testResultData = testResultData
        self.testResultLabels = testResultLabels

    def trainml_withfigs(self):
        # first checking for each of the major diagnoses
        diagCodes_use = self.diagCodes[usePtInds] #need usePtInds too
        falsepos = []
        truepos = []
        roc_auc = []
        codeTest = []
        predComps = []
        predCompsWeights = []
        r2_score_lasso = []
        self.coefs = []
        #is there a way to avoid using both CCI and CCI_S?
        for self.sortedInd in range(len(self.commonCodeinds_sorted)):
            self.coefs.append([])
            codeInd = self.commonCodeinds[self.commonCodeinds_sorted[self.sortedInd]]
            catPtDiagUse = np.zeros(len(diagCodes_use))
            catPtDiagUse[np.where(diagCodes_use.T[codeInd])[0]] = 1

            X = self.testResultData
            Y = catPtDiagUse

            trainfrac = 8
            diagInds = np.where(catPtDiagUse == 1)[0]
            if len(diagInds) < diagThreshold:
                print('not enough trials for ' + diagDict['Diagnosis'][self.sortedInd])
                continue
            numReps = 100
            # manually choosing training/testing samples
            for rep in range(numReps):
                nSampsD = len(diagInds)
                nondiagInds = np.where(catPtDiagUse == 0)[0]
                nSampsND = len(nondiagInds)
                randsampsD_OS = np.random.choice(range(0, nSampsD * trainfrac // 10), nSampsND)
                randsampsND = random.sample(range(0, nSampsND), nSampsND)
                # to balance, use nSampsD (number of samples with diagnosis) for selecting non-diagnosed pts too
                trainsamps = np.concatenate(([randsampsD_OS[n] for n in range(0, nSampsND * trainfrac // 10)],
                                             [randsampsND[n] for n in range(0, nSampsND * trainfrac // 10)]), axis=0)
                testsamps = np.concatenate(([randsampsD_OS[n] for n in range(nSampsND * trainfrac // 10, nSampsND)],
                                            [randsampsND[n] for n in range(nSampsND * trainfrac // 10, nSampsND)]), axis=0)

                X_train, X_test = X[trainsamps], X[testsamps]
                Y_train, Y_test = Y[trainsamps], Y[testsamps]
                Y_preds=run_lasso(X_train,Y_train,X_test)

                fpr, tpr, _ = roc_curve(Y_test, Y_preds)
                tpr_interp = np.interp(np.asarray(range(0, 10, 1)) / 10.0, fpr, tpr)
                falsepos.append(np.asarray(range(0, 10, 1)) / 10.0)
                truepos.append(tpr_interp)
                roc_auc.append(auc(np.asarray(range(0, 10, 1)) / 10.0, tpr_interp))

                codeTest.append(codeInd)
                r2_score_lasso.append(r2_score(Y_test, Y_preds))

            m = max(roc_auc[-numReps:])
            ind = (roc_auc[-numReps:].index(m))
            # only plot if the maximum AUC observed is >0.6
            if np.mean(roc_auc[-numReps:]) > 0.6:
                plotauc()

                # sorted indices of regression coefficients
                reginds_sorted = sorted(range(len(np.mean(coefs[-1], axis=0))),
                                        key=np.abs(np.mean(coefs[-1], axis=0)).__getitem__, reverse=True)

                # save weights and titles of components that provide predictive power
                predCompsWeights.append([np.mean(coefs[-1], axis=0)[i] for i in reginds_sorted])
                predComps.append([self.diagDict['Diagnosis'][self.sortedInd], [self.testResultLabels[i] for i in reginds_sorted]])
                if r2_score_lasso[-1] > 0.1:
                    print self.diagDict['Diagnosis'][self.sortedInd]
                    print("r^2 on test data : %f" % r2_score_lasso[-1])

                self.coefs[-1] = list(self.coefs[-1])
            else:
                print(self.diagDict['Diagnosis'][self.sortedInd] + ' AUC below 0.6')
            return predComps, predCompsWeights


    def run_lasso(self,X_train,Y_train,X_test):
        reg = Lasso(alpha=0.01, copy_X=True, fit_intercept=True, max_iter=1000, normalize=False,
                        positive=False, precompute=False, random_state=None, selection='cyclic', tol=0.0001,
                        warm_start=False)
        reg.fit(X_train, Y_train)
        Y_preds = reg.predict(X_test)
        self.coefs[-1].extend(reg.coef_)
        return Y_preds


    def plot_auc(self,falsepos,truepos,roc_auc):
        fig = plt.figure()
        cax = plt.subplot(1, 3, 1)
        plt.errorbar(np.mean(falsepos[-numReps:], axis=0), np.mean(truepos[-numReps:], axis=0),
                     np.std(truepos[-numReps:], axis=0),
                     color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % np.mean(roc_auc[-numReps:]))
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        xlo = plt.xlabel('False Positive Rate')
        ylo = plt.ylabel('True Positive Rate')
        titleobj = plt.title('ROC ' + diagDict['Diagnosis'][self.sortedInd])
        plt.legend(loc="lower right")
        coefs[-1] = np.reshape(np.asarray(self.coefs[-1]), (self.numReps, len(coefs[-1]) / numReps))
        cax = plt.subplot(1, 3, 3)
        for i in coefs[-1]:
            plt.plot(i, color='black')
        xlo = plt.xlabel('Feature')
        ylo = plt.ylabel('Weight')
        plt.ylim([-0.15, 0.15])

        titleobj = plt.title('Feature weights, n=' + str(nSampsD))
        for tick in cax.get_xticklabels():
            tick.set_rotation(45)

        fig.savefig(self.diagDict['Diagnosis'][self.sortedInd] + 'ROC_' + str(ptThreshold) + 'ptmin_oversampled.eps',
                    format='eps')
        plt.show()
