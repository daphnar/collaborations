import sklearn
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def get_pca_pvalue_from_SVC(PC1,PC2,Y,plot = True): #Y is case/control
    model = sklearn.svm.SVC(kernel='linear')
    X = np.column_stack((np.array(PC1), np.array(PC2)))
    model.fit(X, Y)
    PC2_classify = (-model.coef_[0][0] / model.coef_[0][1]) * xx - (model.intercept_[0]) / model.coef_[0][1]
    if plot:
        plt.scatter(PC1,PC2,c=Y)
        plt.plot(PC1,PC2_classify,'k')
        plt.show()
    prediction = [0 if PC2_classify[i]<pc2 else 1 for i,pc2 in enumerate(PC2)]
    statistic, pvalue = mannwhitneyu(prediction, Y)
    return pvalue
