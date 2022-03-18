import pandas
import os
import numpy as np
# allmicepath='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/DFOut'
# miceByFigurePath='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/mice_extract_nameoffigure.xlsx'
# allMiceByFigure=pandas.read_excel(miceByFigurePath).set_index('uniting_key')
allHumanpath='/Users/daphna/Downloads'
humanByFigurePath='/Users/daphna/Downloads/HumanExtract.xlsx'
allHumanByFigure=pandas.read_excel(humanByFigurePath).set_index('uniting_key')
PNPK0math='/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017/KEGG_genes_Bacteria_protein_old_DF/KO_mat_by_KEGG_0.0001.df'
PNPmetadataPath='/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/DFOut/StoolMetadataDF.dat'
IsNexteraPath='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/HumanSamples/all_fds_barcode_df.dat'
WTvsSOD10PathMPA='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/newPipeline/DFOut3/MPA.df'
WTvsSOD10PathGenes='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/newPipeline/DFOut3/KEGGGene.df'
#WTvsSOD10PathGenes='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/newPipeline/KEGGGenesMappingFiles_HGF/KO_mat_by_KEGG_0.0001.df'
def extract(*args, **kwargs):
    if 'organism' not in kwargs:
        kwargs['organism']='mice'
    if 'figure' not in kwargs:
        kwargs['figure']='all'
    # if kwargs['organism']=='mice':
    #     organismPath=allmicepath
    #     organismByFigurePath=miceByFigurePath
    #     allOrganismByFigure=allMiceByFigure
    #     if 'threshold' not in kwargs:
    #         kwargs['threshold'] = 1e6
    #     if 'thresholdKegg' not in kwargs:
    #         kwargs['thresholdKegg'] = 20 / ((1e6) * 0.8)
    #     if 'thresholdMpa' not in kwargs:
    #         kwargs['thresholdMpa'] = 5e-4
    # elif kwargs['organism']=='human':
    if kwargs['organism'] == 'human':
        if 'cohort' not in kwargs:
            kwargs['cohort'] = 'ALS'
        organismPath=allHumanpath
        organismByFigurePath=humanByFigurePath
        if kwargs['cohort']=='ALS':
            allOrganismByFigure=allHumanByFigure
        elif kwargs['cohort']=='PNP':
            kwargs['KeepAll']=True
            samples=pandas.read_pickle(PNPK0math).index
            allOrganismByFigure=pandas.DataFrame({'uniting_key':samples.values,
                                                  'Genotype':['WT']*len(samples),
                                                  'NameOfFigure':[None]*len(samples),
                                                  'Timepoint':[None]*len(samples),
                                                  'Treatment':[None]*len(samples),
                                                  'Swab':[False]*len(samples),
                                                  'Nextera':[False]*len(samples)}).set_index('uniting_key')
            metadata = pandas.read_pickle(PNPmetadataPath).reset_index()[['FD', 'IsGenotek']].set_index('FD')
            allOrganismByFigure['FD'] = [fd_reg[0] for fd_reg in allOrganismByFigure.index.str.split('_')]
            allOrganismByFigure['UserID']=[fd_reg[1] for fd_reg in allOrganismByFigure.index.str.split('_')]

            allOrganismByFigure.loc[allOrganismByFigure['FD'].\
                                        isin(metadata[metadata['IsGenotek'] == 0].index), 'Swab']=True
            nextera = pandas.read_pickle(IsNexteraPath)
            nextera = nextera[nextera['Barcode'] == 'Nextera']
            allOrganismByFigure.loc[allOrganismByFigure['FD'].isin(nextera.index),'Nextera']=True

            allOrganismByFigure.drop('FD', inplace=True, axis=1)
            allOrganismByFigure.loc[(allOrganismByFigure['Nextera']==True)& \
                                (allOrganismByFigure['Swab'] == True),'NameOfFigure']='Human'
        if 'threshold' not in kwargs:
            if kwargs['organism']=='mice':
                kwargs['threshold'] = 1e6
            elif kwargs['organism']=='human':
                kwargs['threshold'] = 7*1e6
        if 'thresholdKegg' not in kwargs:
            #kwargs['thresholdKegg'] = 2 / (1e7) #20 / ((1e7) * 0.8)
            kwargs['thresholdKegg'] =20/((7*1e6)*0.8*0.55)
        if 'thresholdMpa' not in kwargs:
            kwargs['thresholdMpa'] = 5e-4

    if 'ratio' not in kwargs:
        kwargs['ratio']= 0.1#0.2
    if 'source' not in kwargs:
        kwargs['source']='mpa'

    if kwargs['source']=='kegg' or ('organism' in kwargs and kwargs['organism']=='human'):
        if 'diamond' not in kwargs:
            kwargs['diamond'] = False

        if kwargs['diamond']:
            #threshold = 0.1
            #allmicedf = pandas.read_pickle(
            #    '/net/mraid08/export/jafar/Microbiome/Analyses/PFP/Databases/AlsMice_DF_KOs_DIAMOND_threshold_%s_NO_Normalization_added_second_hit' % threshold)
            if 'IGCnKEGG' not in kwargs:
                kwargs['IGCnKEGG']=False
            if kwargs['organism']=='mice':
                # threshold=0.001
                # allmicedf =pandas.read_pickle(
                #     '/net/mraid08/export/jafar/Microbiome/Analyses/PFP/Databases/AlsMice_DF_KOs_DIAMOND_threshold_%s_NO_Normalization' % threshold)
                if kwargs['figure']=='WTvsSOD10':
                    allmicedf = pandas.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/newPipeline/KEGGGenesMappingFiles_HGF/KO_mat_by_KEGG_0.0001.df')
                    allmicedf.index=allmicedf.index.map(lambda x: x.split('~')[-1])
                    kwargs['thresholdKegg'] = 20 / ((5e6) * 0.8)
                    kwargs['KeepAll'] = True
                else:
                    allmicedf = pandas.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/tmp/Raw/KEGGGenesMapFiles/KO_mat_by_KEGG_0.0001.df')
                allmicedf.columns = allmicedf.columns.str.lstrip('ko:')
                allmicedf = allmicedf.drop('UNKNOWN', axis=1)
                allmicedf = allmicedf.div(allmicedf.sum(axis=1), axis=0)
            elif kwargs['organism'] == 'human':
                if kwargs['cohort'] == 'PNP':
                    allmicedf = pandas.read_pickle(PNPK0math)
                elif kwargs['cohort'] == 'ALS':
                    allmicedf = pandas.read_csv('/Users/daphna/Downloads/MPA.df.species.csv',index_col=0)
                    #allmicedf = pandas.read_pickle(os.path.join(allHumanpath, 'KEGGGene.df'))
                    # Before revision
                    # if kwargs['IGCnKEGG']:
                    #     allmicedf=pandas.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/PFP/Databases/AlsHuman_DF_KOs_DIAMOND_threshold_1e-07_NO_Normalization')
                    # else:
                    #     allmicedf=pandas.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/HumanSamples/tmp/Raw/KEGGGenesMapFiles/KO_mat_by_KEGG_0.0001.df')#KO_mat_by_KEGG_10.df')
                # allmicedf.columns=allmicedf.columns.str.lstrip('ko:')
                # allmicedf=allmicedf.drop('UNKNOWN',axis=1)
                # allmicedf=allmicedf.div(allmicedf.sum(axis=1),axis=0)
            else:
                if kwargs['figure'] == 'WTvsSOD10':
                    allmicedf = pandas.read_pickle(WTvsSOD10PathGenes)
                    kwargs['thresholdKegg'] = 20 / ((5e6) * 0.8)
                    kwargs['KeepAll'] = True
                else:
                    allmicedf = pandas.read_pickle(os.path.join(organismPath, 'KeggGeneSpidDF.dat'))
            kwargs['sourceThreshold'] = kwargs['thresholdKegg']
            if kwargs['source'] == 'mpa':
                if kwargs['figure'] == 'WTvsSOD10':
                    allmicedf = pandas.read_pickle(WTvsSOD10PathMPA)
                    kwargs['KeepAll'] = True
                    kwargs['thresholdMpa'] = 1e-4
                else:
                    allmicedf = pandas.read_csv('/Users/daphna/Downloads/MPA.df.species.csv',index_col=0)
                # if ('taxa' in kwargs.keys()):
                #     assert 'taxLevel' not in kwargs.keys(), \
                #         'taxa is mutual exclusive with taxLevel - t,s,g,f,o,c,p'
                #     allmicedf = allmicedf.loc[('s', kwargs['taxa'])]
                # else:
                #     if 'taxLevel' not in kwargs:
                #         kwargs['taxLevel'] = 's'
                #     allmicedf = allmicedf.T.loc[(kwargs['taxLevel'],)].T
                kwargs['sourceThreshold'] = kwargs['thresholdMpa']
            # if kwargs['source'] == 'count':
            #     allmicedf = pandas.read_pickle(os.path.join(organismPath, 'ReadCountDF.dat')).reset_index().set_index(
            #         'uniting_key')
            # if 'KeepAll' not in kwargs:
            #     countdf = pandas.read_pickle(os.path.join(organismPath, 'ReadCountDF.dat')).reset_index().set_index(
            #         'uniting_key')
            #     aboveThresholdSamples = countdf[countdf['PostSubSamp'] >= kwargs['threshold']].index
            #     allmicedf = allmicedf.loc[aboveThresholdSamples.intersection(allmicedf.index)]

            if kwargs['source'] != 'count':
                allmicedf[allmicedf < kwargs['sourceThreshold'] + 0.00000001] = np.nan
                allmicedf.dropna(how='all', axis=1, inplace=True)
                allmicedf.fillna(kwargs['sourceThreshold'], inplace=True)
    miceByFigure = allOrganismByFigure.copy()
    mice = []

    if kwargs['figure'] == 'WTvsSOD12':
        kwargs['treatment'] = 'PBS'
        kwargs['figure'] = 'AkkKEGG'
    if kwargs['figure'] == 'WTvsSOD1' or kwargs['figure'] == 'all':
        mice += miceByFigure[miceByFigure['NameOfFigure'] == 'WTvsSOD1'].index.values.tolist()
    if kwargs['figure'] == 'CoHousing' or kwargs['figure'] == 'all':
        mice += miceByFigure[miceByFigure['NameOfFigure'] == 'CoHousing'].index.values.tolist()
    if kwargs['figure'] == 'AkkKEGG' or kwargs['figure'] == 'all' or kwargs['figure'] == 'WTvsSOD12':
        mice += miceByFigure[miceByFigure['NameOfFigure'] == 'AkkKEGG'].index.values.tolist()
    if kwargs['figure'] == 'Human' or kwargs['figure'] == 'all':
        mice += miceByFigure[miceByFigure['NameOfFigure'] == 'Human'].index.values.tolist()
    if kwargs['figure'] == 'WTvsSOD10' or kwargs['figure'] == 'all':
        mice += miceByFigure[miceByFigure['NameOfFigure'] == 'WTvsSOD10'].index.values.tolist()

    miceByFigure = miceByFigure.loc[mice]
    if kwargs['organism'] == 'human':
        micedf = allmicedf.loc[list(set(miceByFigure.UserID.apply(lambda x: x.replace('Sample_','')).values))]  # miceByFigure.index]
    else:
        micedf = allmicedf.loc[miceByFigure.index]

    if kwargs['source'] != 'count':
        presence = (micedf > kwargs['sourceThreshold'] + 1e-7).astype(int).sum()
        presence = presence[presence > len(micedf.index) * kwargs['ratio']].index.values.tolist()
        micedf = micedf[presence]

    if 'timepoint' not in kwargs:
        kwargs['timepoint'] = None
    if kwargs['timepoint'] is not None:
        if type(kwargs['timepoint']) is int:
            miceByFigure = miceByFigure[miceByFigure['Timepoint'] == kwargs['timepoint']]
        elif kwargs['timepoint'] == 'D>0':
            miceByFigure = miceByFigure[miceByFigure['Timepoint'] > 0]
        elif kwargs['timepoint'] == 'sort':
            miceByFigure.sort_values('Timepoint', inplace=True, ascending=False)

    if kwargs['organism'] == 'human':
        micedf = micedf.loc[list(set(miceByFigure.UserID.values))]  # miceByFigure.index]
    else:
        micedf = micedf.loc[miceByFigure.index]

    if 'k0s' in kwargs:
        k0s = list(set(micedf.columns.values).intersection(set(kwargs['k0s'])))
        print
        "looking for: ", len(set(kwargs['k0s'])), "K0s, found ", len(k0s), " K0s"
        micedf = micedf[k0s]

    genotype = []
    if 'genotype' not in kwargs:
        kwargs['genotype'] = 'all'
    if kwargs['genotype'] == 'SOD1' or kwargs['genotype'] == 'all':
        genotype += miceByFigure[miceByFigure['Genotype'] == 'SOD1'].index.values.tolist()
    if kwargs['genotype'] == 'WT' or kwargs['genotype'] == 'all':
        genotype += miceByFigure[miceByFigure['Genotype'] == 'WT'].index.values.tolist()
    if kwargs['genotype'] == 'ALS' or kwargs['genotype'] == 'all':
        genotype += miceByFigure[miceByFigure['Genotype'] == 'ALS'].index.values.tolist()
    if kwargs['genotype'] == 'Control' or kwargs['genotype'] == 'all':
        genotype += miceByFigure[miceByFigure['Genotype'] == 'Control'].index.values.tolist()
    miceByFigure = miceByFigure.loc[genotype]

    treatment=[]
    if  kwargs['organism']=='human' or kwargs['figure']=='WTvsSOD1' \
            or kwargs['figure']=='WTvsSOD10':
        kwargs['treatment']=None
    if kwargs['treatment'] is None and \
            kwargs['organism']!='human' and \
            kwargs['figure']!='WTvsSOD1' and \
            kwargs['figure']!='WTvsSOD10':
        kwargs['treatment']='all'
    if kwargs['treatment'] == 'Akkermansia' or kwargs['treatment'] == 'all':
        treatment += miceByFigure[miceByFigure['Treatment'] == 'Akkermansia'].index.values.tolist()
    if kwargs['treatment'] == 'PBS' or kwargs['treatment'] == 'all':
        treatment += miceByFigure[miceByFigure['Treatment'] == 'PBS'].index.values.tolist()
    if kwargs['treatment'] == 'SH' or kwargs['treatment'] == 'all':
        treatment += miceByFigure[miceByFigure['Treatment'] == 'SH'].index.values.tolist()
    if kwargs['treatment'] == 'CH' or kwargs['treatment'] == 'all':
        treatment += miceByFigure[miceByFigure['Treatment'] == 'CH'].index.values.tolist()
    if kwargs['treatment'] is not None:
        miceByFigure = miceByFigure.loc[treatment]

    if kwargs['organism']=='human':
        micedf=micedf.loc[list(set(miceByFigure.UserID.values))]#miceByFigure.index]
    else:
        micedf=micedf.loc[miceByFigure.index]

    if 'apply' not in kwargs:
        kwargs['apply']=None
    if kwargs['apply'] is not None:
        micedf['UserID']=miceByFigure.loc[micedf.index,'UserID']
        micedf = micedf.groupby('UserID').apply(kwargs['apply'])
        if 'UserID' in micedf.columns:
            micedf.drop('UserID', inplace=True, axis=1)

    if kwargs['source']=='count':
        kwargs['scale']=None
    if 'scale' not in kwargs:
        kwargs['scale'] = 'log'
    if kwargs['scale'] == 'log':
        micedf = np.log10(micedf)
    if kwargs['source']!='count':
        micedf = micedf.dropna()

    if 'MostVariable' in kwargs:
        stds=micedf.std().sort_values(ascending=False)[:100].index.values
        micedf=micedf[stds]

    if 'MostAbundant' in kwargs:
        abun=micedf.median().sort_values(ascending=False)[:100].index.values
        micedf=micedf[abun]
    return micedf

if __name__=="__main__":
    # print extract(source='count',figure='WTvsSOD1').to_excel('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/ALSElinav/MiceRuns/mice_readCountsWTvsSOD1.xlsx')
    mice = extract(figure='Human', source='mpa', genotype='all',diamond=True,
                   organism='human')