import pandas as pd
from eranb import config
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
# from mne.stats import fdr_correction
# from extractData import extract
import seaborn as sns
from scipy.stats import mannwhitneyu
import os
from sklearn.decomposition import PCA
import xlrd
import sklearn
import numpy as np

def get_pca_pvalue_from_SVC(PC1,PC2,Y,plot = True): #Y is case/control
    model = sklearn.svm.SVC(kernel='linear')
    X = np.column_stack((np.array(PC1), np.array(PC2)))
    model.fit(X, Y)
    PC2_classify = (-model.coef_[0][0] / model.coef_[0][1]) * PC1 - (model.intercept_[0]) / model.coef_[0][1]
    if plot:
        #plt.scatter(PC1,PC2,c=Y)
        plt.plot(PC1,PC2_classify,'k')
        #plt.show()
    prediction = [0 if PC2_classify[i]<pc2 else 1 for i,pc2 in enumerate(PC2)]
    statistic, pvalue = mannwhitneyu(prediction, Y)
    return pvalue

def plot_PCA(X,Y=None, name="", impute = None,plot=True,format='pdf'):
    if impute is not None:
        X=X.fillna(X.mean())
    pca = PCA()
    allPCs = pca.fit_transform(X)
    all_pcoa_explained_var = pca.explained_variance_ratio_
    df_PCs=pd.DataFrame(allPCs, index=X.index)
    df_PCs.to_csv(os.path.join(config.analyses_path, 'PCA__%s.csv'%name))
    if plot:
        for pc_1 in range(5):
            for pc_2 in range(pc_1,5):
                if pc_1==pc_2:
                    continue
                plt.figure(facecolor=(1, 1, 1))
                if Y is not None:
                    colors = list(set(Y.dropna().values))
                    colors.sort()
                    cmap=plt.cm.get_cmap('Blues')#, len(colors))

                    plt.scatter(df_PCs[pc_1],df_PCs[pc_2],c=Y.values,#dataColors,
                                s=40,cmap=cmap)
                    ticksValues=list(set(Y.values))
                    plt.colorbar(ticks=ticksValues, label='Disease')
                    plt.clim(min(ticksValues)-0.5,max(ticksValues)+0.5)
                    print("PC%s-PC%s pvalue %s"%(pc_1+1,pc_2+1,get_pca_pvalue_from_SVC(df_PCs[pc_1].values, df_PCs[pc_2].values, Y.values.astype(int), plot=True)))
                else:
                    plt.scatter(df_PCs[pc_1].values, df_PCs[pc_2].values, marker='.', c='b', s=40,alpha=0.6)
                plt.xlabel('PC%s %.2f'%(pc_1+1,all_pcoa_explained_var[pc_1]),fontsize=14)
                plt.ylabel('PC%s %.2f'%(pc_2+1,all_pcoa_explained_var[pc_2]),fontsize=14)
                figname='PCA_%s_%s__%s.%s'%(pc_1+1,pc_2+1,name,format)
                plt.subplots_adjust(left=0.2)
                plt.savefig(os.path.join(config.analyses_path, figname), format=format)
                plt.close()
        #plt.show()
    return df_PCs, all_pcoa_explained_var


def getPhenotypes():
    xls = xlrd.open_workbook(config.phenotypes_df, on_demand=True)
    sheets = xls.sheet_names()
    for sheetname in sheets:
        pheno_df = pd.read_excel(config.phenotypes_df, sheet_name=sheetname).set_index('Animal ID')
    #pheno_df = pd.read_excel(config.phenotypes_df).set_index('Animal ID')
        for category ,pair in config.boolean_types.items():
            for key,value in pair.items():
                pheno_df.loc[pheno_df[category]==key,category] = value
        yield pheno_df.drop(config.y_category, axis=1), pheno_df[config.y_category], sheetname

def corrlatePhenotypesWithDisease(X,Y):
    pvalues=[]
    statistics = []
    for phenotype in X.columns.values:
            subset_phenotype_x = X[phenotype].dropna()
            statistic, pvalue = mannwhitneyu(subset_phenotype_x.values,Y.loc[subset_phenotype_x.index].values)
            statistics.append(statistic)
            pvalues.append(pvalue)

    qvalues=multipletests(pvalues,method='fdr_bh')[1]
    results=pd.DataFrame({'Phenotype':X.columns.values,
                          'statistics':statistics,
                          'pvalues':pvalues,
                          'qvalues':qvalues}).set_index('Phenotype')
    results.to_csv(os.path.join(config.analyses_path, 'Disease_associations.csv'))

if __name__=='__main__':
    for X,Y,name in getPhenotypes():
    #X,Y,name = getPhenotypes()
        plot_PCA(X,Y,name,impute='average')
    #corrlatePhenotypesWithDisease(X,Y)