from scipy.spatial.distance import squareform, pdist
import scipy.linalg as la
import warnings
def e_matrix(distance_matrix):
    """Compute E matrix from a distance matrix.
    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998."""
    return distance_matrix * distance_matrix / -2.0

def f_matrix(E_matrix):
    """Compute F matrix from E matrix.
    Centring step: for each element, the mean of the corresponding
    row and column are substracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
    row_means = E_matrix.mean(axis=1)#, keepdims=True)
    col_means = E_matrix.mean(axis=0)#, keepdims=True)
    matrix_mean = E_matrix.mean()
    return E_matrix - row_means - col_means + matrix_mean
def PCoA(D, suppress_warning=False, return_explained_variance=False):
    E_matrix = e_matrix(D)
    F_matrix = f_matrix(E_matrix)
    eigvals, eigvecs = la.eigh(F_matrix)
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0
    if (np.any(eigvals < 0) and not suppress_warning):
        warnings.warn(
            "The result contains negative eigenvalues."
            " Please compare their magnitude with the magnitude of some"
            " of the largest positive eigenvalues. If the negative ones"
            " are smaller, it's probably safe to ignore them, but if they"
            " are large in magnitude, the results won't be useful. See the"
            " Notes section for more details. The smallest eigenvalue is"
            " {0} and the largest is {1}.".format(eigvals.min(),
                                                  eigvals.max()),
            RuntimeWarning
            )


    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]
    num_positive = (eigvals >= 0).sum()
    coordinates = eigvecs * np.sqrt(eigvals)

    if return_explained_variance: return coordinates, eigvals/eigvals.sum()
    return coordinates

def BrayCurtis(X_orig, is_log_abundance=True, zero_min_value=True):
    if is_log_abundance: X = 10**X_orig
    else: X = X_orig.copy()
    if zero_min_value: X[X_orig==np.min(X_orig)]=0

    ####X = np.array([[1,2,3,0], [3,2,1,0], [0,0,0,1]])
    D = squareform(pdist(X, metric='braycurtis'))
    return D

def extractMPA(taxaLevel,bacteria=None,dataSource='PNP',onlySwab=True,threshold=1e-4,filterNonDetectableSpeciest=True,
               filterRatio=0.2,abundances='RA',scale='log'):

    dfmpa=pd.read_csv(ALSmpaPath,index_col=0)
    dfmpa=dfmpa.loc[dfmpa.index.map(lambda x: x.startswith('ALS') or x.startswith('CTRL'))]
    #dfmpa=dfmpa[dfmpa.index.get_level_values(0)==taxaLevel].reset_index().drop('TaxLevel',axis=1).set_index('Tax')
    #dfmpa=dfmpa.T
    dfmpa[dfmpa<=threshold]=np.nan
    dfmpa=dfmpa.dropna(how='all',axis=1)
    dfmpa=dfmpa.fillna(threshold)
    if dataSource=='PNP':
        dfmpa=dfmpa.groupby('RegNum').median()
    else:
        dfmpa.index.name='RegNum'
    if filterNonDetectableSpeciest and bacteria is None:
        appearInMpa=(dfmpa>threshold).sum()
        minPresent=dfmpa.shape[0]*filterRatio
        detectableSpecies=appearInMpa[appearInMpa>minPresent].index.values
        dfmpa=dfmpa[detectableSpecies]
    if bacteria is not None:
        dfmpa=dfmpa[bacteria]
    dfmpa=dfmpa.apply(np.log10)
    # dfmpa=dfmpa.T
    if abundances=='01':
        dfmpa = (dfmpa>-3.99999999).astype(np.int).astype(np.float)
    # if scale=='log':
    #     dfmpa = 10**dfmpa
#     if filterALLnonvariant:
#         #remove non-useful features
#         columns_common = dfmpa.apply(lambda s: s.value_counts().idxmax())
#         frac_common = (dfmpa==columns_common).mean(axis=0)
#         is_sparse_feature = (frac_common>0.97)
#         print 'Removing %d sparse features if not phylum'%(is_sparse_feature.sum())
#         print 'Removing ',dfmpa.loc[:,is_sparse_feature].columns.tolist(), 'if not phylum'
#         if taxaLevel!='p':
#             dfmpa = dfmpa.loc[:, ~is_sparse_feature]
#     dfmpa=dfmpa.T
    return dfmpa


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from kernel_utils import PCoA,BrayCurtis
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
ALSmpaPath='/Users/daphna/Downloads/MPA.df.species.csv'
# PNPmpaPath='/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/DFOut/MPA.dat'
threshold=1e-4

resolutionLevel='s'
abundances='RA'
scale='log'
filterNonDetectableSpeciest=False
filterRatio=0
usePNP=False
# pnpmpa=extractMPA(resolutionLevel,dataSource='PNP',threshold=threshold,
#                               filterNonDetectableSpeciest=filterNonDetectableSpeciest,filterRatio=filterRatio,
#                               abundances=abundances,scale=scale)
alsmpa=extractMPA(resolutionLevel,dataSource='ALS',threshold=threshold,
                              filterNonDetectableSpeciest=filterNonDetectableSpeciest,filterRatio=filterRatio,
                              abundances=abundances,scale=scale)
# alldata=pd.merge(pnpmpa.reset_index(),alsmpa.reset_index(),on=['RegNum']+list(set(pnpmpa.columns.tolist()).intersection(set(alsmpa.columns.tolist()))),how='outer')
# alldata=alldata.set_index('RegNum')
if not usePNP:
    alldata = alsmpa
alldata[alldata<-3.99999999]=np.nan
alldata.dropna(how='all',axis=1,inplace=True)
alldata.fillna(threshold,inplace=True)
allmpa=alldata
ALS = pd.Series(index=alsmpa.index,data=alsmpa.index.str.startswith('ALS')).astype(int)
BC=BrayCurtis(allmpa)
allPCs, all_pcoa_explained_var = PCoA(BC, return_explained_variance=True)
#pd.DataFrame(allPCs,index=allmpa.index).to_excel('/Users/daphna/Downloads/PCHuman.xlsx')
#allmpa.to_excel('/Users/daphna/Downloads/MPAHuman.xlsx')
allPCs=allPCs[:, :5]
allpcs_df=pd.DataFrame(allPCs[:,:2],index=allmpa.index)
ALS_BC=allpcs_df.loc[ALS[ALS==1].index]
C_BC=allpcs_df.loc[ALS[ALS==0].index]
plt.figure()
# if usePNP:
#         pnp_BC=allpcs_df.loc[pnpmpa.index]
#         plt.scatter(pnp_BC[0].values,pnp_BC[1].values,c='g',s=40)
plt.scatter(C_BC[0].values,C_BC[1].values,c='b',s=40)
plt.scatter(ALS_BC[0].values,ALS_BC[1].values,c='r',s=40)
plt.xlabel('PCo1 %.2f'%all_pcoa_explained_var[0])
plt.ylabel('PCo2 %.2f'%all_pcoa_explained_var[1])
plt.show()