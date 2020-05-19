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

def plot_PCA(X,Y=None, impute = None,plot=True,format='pdf'):
    if impute is not None:
        X=X.fillna(X.mean())
    pca = PCA()
    allPCs = pca.fit_transform(X)
    all_pcoa_explained_var = pca.explained_variance_ratio_
    df_PCs=pd.DataFrame(allPCs, index=X.index)
    df_PCs.to_csv(os.path.join(config.analyses_path, 'PCA.csv'))
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
                else:
                    plt.scatter(df_PCs[0].values, df_PCs[1].values, marker='.', c='b', s=40,alpha=0.6)
                plt.xlabel('PC%s %.2f'%(pc_1+1,all_pcoa_explained_var[pc_1]),fontsize=14)
                plt.ylabel('PC%s %.2f'%(pc_2+1,all_pcoa_explained_var[pc_2]),fontsize=14)
                figname='PCA_%s_%s.%s'%(pc_1+1,pc_2+1,format)
                plt.subplots_adjust(left=0.2)
                plt.savefig(os.path.join(config.analyses_path, figname), format=format)
        plt.show()
    return df_PCs, all_pcoa_explained_var


def getPhenotypes():
    pheno_df = pd.read_excel(config.phenotypes_df).set_index('Animal ID')
    for category ,pair in config.boolean_types.items():
        for key,value in pair.items():
            pheno_df.loc[pheno_df[category]==key,category] = value
    return pheno_df.drop(config.y_category, axis=1), pheno_df[config.y_category]

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
    X,Y = getPhenotypes()
    plot_PCA(X,Y,impute='average')
    corrlatePhenotypesWithDisease(X,Y)