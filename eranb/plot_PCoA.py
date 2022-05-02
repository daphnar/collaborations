import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

from scipy.spatial.distance import squareform, pdist
import scipy.linalg as la
import warnings

from skbio.stats.ordination import pcoa
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

def plot(source1=None,source2=None,indexName='index',
         threshold=1e-4,legend=None,figurename=None,
         method='PCoA',excelName=None,data=None,groups=None):
    if (source1 is not None) and (source2 is not None) and (data is None):
        if type(indexName) is str:
            alldata=pd.merge(source1.reset_index(),source2.reset_index(),on=[indexName]+list(set(source1.columns.tolist()).intersection(set(source2.columns.tolist()))),how='outer')
        else:
            alldata=pd.merge(source1.reset_index(),source2.reset_index(),on=indexName+list(set(source1.columns.tolist()).intersection(set(source2.columns.tolist()))),how='outer')
    else:
        alldata=pd.concat([data[group] for group in groups]).reset_index()
    alldata=alldata.set_index(indexName)
    shapebefore=alldata.shape[1]
    alldata.dropna(axis=1,inplace=True)
    shapeafter = alldata.shape[1]
    alldata[alldata<np.log10(threshold)+0.0000001]=np.nan
    alldata.dropna(how='all',axis=1,inplace=True)
    alldata.fillna(threshold,inplace=True)

    if method=='PCoA':
        BC = BrayCurtis(alldata)
        allPCs, all_pcoa_explained_var = PCoA(BC, return_explained_variance=True)
    if method=='PCA':
        pca=PCA()
        allPCs=pca.fit_transform(alldata)
        all_pcoa_explained_var=pca.explained_variance_ratio_
    if method == 'TSNE':
        allPCs = TSNE(n_components=2).fit_transform(alldata)
    allpcs_df=pd.DataFrame(allPCs[:,:2],index=alldata.index)
    # if excelName is not None:
    #     pd.DataFrame(allPCs, index=alldata.index).to_excel(excelName)
    plt.figure()
    if (source1 is not None) and (source2 is not None) and (data is None):
        source1_BC=allpcs_df.loc[source1.index]
        source2_BC=allpcs_df.loc[source2.index]
        plt.scatter(source1_BC[0].values,source1_BC[1].values,c='b',s=40)
        plt.scatter(source2_BC[0].values,source2_BC[1].values,c='r',s=40)
        ALSpc1=np.concatenate((source1_BC[0].values,source2_BC[0].values))
        ALSpc2=np.concatenate((source1_BC[1].values,source2_BC[1].values))
        ALS=np.concatenate((np.array([0]*len(source1_BC[0].values)),
                            np.array([1] * len(source2_BC[0].values))))
        print(spearmanr(ALSpc1, ALS))
        print(spearmanr(ALSpc2, ALS))
    else:
        colors=['r','g','b','k']
        for idx,group in enumerate(groups):
            source_BC = allpcs_df.loc[data[group].index]
            plt.scatter(source_BC[0].values, source_BC[1].values, c=colors[idx], s=40)
    if method == 'PCoA':
        plt.xlabel('PCo1 %.2f'%all_pcoa_explained_var[0])
        plt.ylabel('PCo2 %.2f'%all_pcoa_explained_var[1])
    if method == 'PCA':
        plt.xlabel('PC1 %.2f'%all_pcoa_explained_var[0])
        plt.ylabel('PC2 %.2f'%all_pcoa_explained_var[1])
    if method == 'TSNE':
        plt.xlabel('TSNE1')
        plt.ylabel('TSNE2')
    if legend is not None:
        plt.legend(legend)
    # if figurename is not None:
    #     plt.savefig(figurename,format='pdf')
    # else:
    plt.show()