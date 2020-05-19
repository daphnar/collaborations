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

def plot_PCA(X,Y=None, plot=True,format='pdf'):

    pca = PCA()
    allPCs = pca.fit_transform(X)
    all_pcoa_explained_var = pca.explained_variance_ratio_
    df_PCs=pd.DataFrame(allPCs, index=X.index)
    df_PCs.to_csv(os.path.join(config.analyses_path, 'PCA.csv'))
    if plot:
        plt.figure()
        if Y is not None:
            colors = list(set(Y.dropna().values))
            colors.sort()
            cmap=plt.cm.get_cmap('cubehelix', len(colors))

            plt.scatter(df_PCs[0],df_PCs[1],c=Y.values,#dataColors,
                        s=40,cmap=cmap)
            ticksValues=list(set(Y.values))
            plt.colorbar(ticks=ticksValues, label='Disease')
            plt.clim(min(ticksValues)-0.5,max(ticksValues)+0.5)
        else:
            plt.scatter(df_PCs[0].values, df_PCs[1].values, c='b', s=40)
        plt.xlabel('PCo1 %.2f'%all_pcoa_explained_var[0],fontsize=14)
        plt.ylabel('PCo2 %.2f'%all_pcoa_explained_var[1],fontsize=14)
        figname='PCA.%s'%(format)
        plt.savefig(os.path.join(config.analyses_path, figname), format=format)
        plt.close()
    return df_PCs, all_pcoa_explained_var


def getPhenotypes():
    pheno_df = pd.read_excel(config.phenotypes_df).set_index('Animal ID')
    for category ,pair in config.boolean_types.items():
        for key,value in pair.items():
            pheno_df[pheno_df[category]==key] = value
    return pheno_df.drop(config.y_category, axis=1), pheno_df[config.y_category]

def corrlatePhenotypesWithDisease(X,Y):
    pvalues=[]
    statistics = []
    for phenotype in X.columns.values:
            statistic, pvalue = mannwhitneyu(X.loc[Y == 0,phenotype],X.loc[Y == 1,phenotype])
            statistics.append(statistic)
            pvalues.append(pvalue)

    qvalues=multipletests(pvalues,method='fdr_bh')[1]
    results=pd.DataFrame({'Phenotype':X.columns.values,
                          'statistics':statistics,
                          'pvalues':pvalues,
                          'qvalues':qvalues}).set_index('Phenotype')
    results.to_csv(os.path.join(config.analyses_path, 'Disease_associations.csv'))

if __name__=='__main__':
    X,Y = getPhenotypes()
    plot_PCA(X,Y)
    corrlatePhenotypesWithDisease(X,Y)