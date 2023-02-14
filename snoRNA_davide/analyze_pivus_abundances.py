import pandas as pd
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from scipy.stats import ranksums, ttest_rel
from skbio.stats.ordination import pcoa
import os
import seaborn as sns
from statsmodels.stats.multitest import fdrcorrection

remote=False
local = not remote

def plot_PCoA(X,legend=None,n_pcs=5,plot_pc=None,ax=None,color_groups=None,s_special=1,titlename=None,
              plot_background=True,arrow=False,arrow_mapping=None,up_down_dict=None,plot_y=True,
              debug=False,output_pcos=None,output_pcos_proportion=None,X_is_pcos=False,color_set=None):
    if output_pcos is not None:
        if not os.path.exists(output_pcos):
            pcs = pcoa(X)
            df_PCs = pcs.samples
            df_PCs.index = X.index
            df_PCs.to_csv(output_pcos)
            all_pcoa_explained_var = pcs.proportion_explained
            pd.Series(all_pcoa_explained_var).to_csv(output_pcos_proportion,index=False)
        else:
            if X_is_pcos:
                df_PCs=X
                all_pcoa_explained_var = pd.read_csv(output_pcos_proportion, index_col=0).index
            else:
                df_PCs=pd.read_csv(output_pcos,index_col=0)
                all_pcoa_explained_var=pd.read_csv(output_pcos_proportion,index_col=0).index
    else:
        pcs = pcoa(X)
        df_PCs = pcs.samples
        df_PCs.index = X.index
        all_pcoa_explained_var = pcs.proportion_explained
    initial=df_PCs.columns.values[0]
    if not str(initial).startswith('PC'):
        df_PCs.columns=df_PCs.columns.map(lambda x: 'PC%s'%(int(x)+1))
    for pc_1 in range(n_pcs):
        for pc_2 in range(pc_1, n_pcs):
            if not debug:
                if (pc_1, pc_2) != plot_pc:
                    continue
                plt.sca(ax)
            else:
                if pc_1!= pc_2:
                    plt.figure()
                else:
                    continue
            if plot_background:
                plt.scatter(df_PCs['PC%s'%(pc_1+1)].values, df_PCs['PC%s'%(pc_2+1)].values, marker='.', c='gray', s=s_special, alpha=0.6)
            for i, group in enumerate(color_groups):
                if color_set is None:
                    plt.scatter(df_PCs.loc[group, 'PC%s'%(pc_1+1)].values, df_PCs.loc[group, 'PC%s'%(pc_2+1)].values, marker='.'
                                , s=s_special, alpha=0.6)
                else:
                    plt.scatter(df_PCs.loc[group, 'PC%s' % (pc_1 + 1)].values,
                                df_PCs.loc[group, 'PC%s' % (pc_2 + 1)].values, marker='.'
                                , s=s_special, alpha=0.6,color=color_set[i])
            plt.xlabel('PC%s %.2f' % (pc_1 + 1, all_pcoa_explained_var[pc_1]))
            if plot_y:
                plt.ylabel('PC%s %.2f' % (pc_2 + 1, all_pcoa_explained_var[pc_2]))
            else:
                plt.gca().axes.yaxis.set_visible(False)
            #ax.yaxis.labelpad = -1
            #ax.set_ylabel('PC%s %.2f' % (pc_2 + 1, all_pcoa_explained_var[pc_2]),labelpad=0)
            if legend is not None:
                plt.legend(legend,loc='lower right')
            # if arrow:
            #     count=0
            #     for age70 in arrow_mapping.index:
            #         age80 = arrow_mapping.loc[age70]
            #         if up_down_dict[age70]=='g':
            #             lines = sns.lineplot(y=[df_PCs.loc[age70,'PC%s'%(pc_2+1)], df_PCs.loc[age80,'PC%s'%(pc_2+1)]],
            #                              x=[df_PCs.loc[age70,'PC%s'%(pc_1+1)], df_PCs.loc[age80,'PC%s'%(pc_1+1)]], ax=ax, color='gray',
            #                              linewidth=0.5, alpha=0.5).get_lines()
            #             plt.scatter(y=[df_PCs.loc[age70, pc_2], df_PCs.loc[age80, pc_2]],
            #                                  x=[df_PCs.loc[age70, pc_1], df_PCs.loc[age80, pc_1]],
            #                                  color='white')
            #         # try:
            #         #     add_arrow(lines[-1], color=up_down_dict[age70], size=4) #'r'
            #         # except IndexError:
            #         #     add_arrow(lines[-1], color=up_down_dict[age70], size=4,direction='left') #'g'
            #             count+=1
            #     print(count)
            plt.title(titlename)


def add_pivus_annotations(ra_df,annotation_table):
    ra = pd.read_csv(ra_df, index_col=0)
    ra_copy = ra.copy()
    annotation_table_df = pd.read_csv(annotation_table, index_col=0,sep='\t')[['RNAseqID70','RNAseqID80']]
    age_dic = {}
    for column in annotation_table_df.columns:
        for sample_id in annotation_table_df[column]:
            age_dic[sample_id]=column[-2:]
    annotation_table_df=pd.Series(age_dic)
    ra_copy=ra_copy.loc[annotation_table_df.index[annotation_table_df.index.isin(ra_copy.index)]]
    ra_copy['age'] = ra_copy.index.map(lambda x: annotation_table_df[x])
    ra_copy=ra_copy.sort_values('age',ascending=True)
    return ra_copy

if remote:
    phenotypes_kallisto = '/oak/stanford/groups/pritch/users/daphna/snoRNA/data/pivus/PIVUS_MapGlobalID_RNAseqMetadata_merge_SeqSuccessOnly.txt'
    ra_df_path = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/pivus/pivus_clean_abundances.csv'
    annotated_ra_path='/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/pivus/pivus_clean_abundances_with_age.csv'
    ra_df = add_pivus_annotations(ra_df_path,phenotypes_kallisto)
    ra_df.to_csv(annotated_ra_path)

if local:
    annotated_ra_with_age_path='/Users/daphna/cluster2/users/daphna/snoRNA/pivus_clean_abundances_with_age.csv'
    annotated_ra_path='/Users/daphna/cluster2/users/daphna/snoRNA/pivus_abundances.csv'
    annotated_ra = pd.read_csv(annotated_ra_path,index_col=0)
    annotated_ra_with_age = pd.read_csv(annotated_ra_with_age_path,index_col=0)
    print(annotated_ra)
    print(annotated_ra.shape)

    min_abundance = 20 / 1500
    annotated_ra[annotated_ra < min_abundance] = min_abundance
    annotated_ra=annotated_ra.loc[:,annotated_ra.std()>0.01]
    age_groups = annotated_ra_with_age['age']
    #distmat = pd.DataFrame(squareform(pdist(annotated_ra.drop('age', 1), 'braycurtis'))).values
    distmat = pd.DataFrame(squareform(pdist(annotated_ra, 'braycurtis'))).values

    distmat_df = pd.DataFrame(index =annotated_ra.index,
                           columns = annotated_ra.index,
                           data=distmat)
    color_groups = [age_groups[age_groups==70].index.values,
                    age_groups[age_groups == 80].index.values]
    plot_PCoA(distmat_df,legend=['70','80'],
              color_groups=color_groups,debug=True,s_special=20)
    print(annotated_ra)
    res_pval = {}
    print(annotated_ra.shape)
    for snoRNA in annotated_ra.columns:
        statistic, pval = ttest_rel(annotated_ra.loc[age_groups[age_groups==70].index.values, snoRNA],
                                    annotated_ra.loc[age_groups[age_groups == 80].index.values, snoRNA])
        #pval = ranksums(annotated_ra.loc[age_groups[age_groups==70].index.values, snoRNA],
        #                annotated_ra.loc[age_groups[age_groups == 80].index.values, snoRNA])[1]
        res_pval[snoRNA]=[pval,1,annotated_ra.loc[age_groups[age_groups==70].index.values, snoRNA].median(),
                          annotated_ra.loc[age_groups[age_groups==80].index.values, snoRNA].median()]
    res_pval_df = pd.DataFrame(index = ['p-value','q-value','70 median RA','80 median RA'],data=res_pval).T
    res_pval_df['q-value']=fdrcorrection(res_pval_df['p-value'])[1]
    res_pval_df=res_pval_df.sort_values('q-value')
    print(res_pval_df)
    print(res_pval_df[res_pval_df['p-value']<0.05])
    print(res_pval_df[res_pval_df['p-value'] < 0.05].shape)
    plt.show()

