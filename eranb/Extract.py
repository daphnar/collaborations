import pandas
import os
import numpy as np

def extract(genotype,sample_annotation,mpa_path,thresholdMpa = 5e-4,ratio = 0.1,taxLevel = 's'):
    samples_df=pandas.read_excel(sample_annotation).set_index('uniting_key')
    mpa_df = pandas.read_pickle(os.path.join(mpa_path, 'MPA.df'))
    mpa_df = mpa_df.T.loc[(taxLevel,)].T
    mpa_df[mpa_df < thresholdMpa + 1e-7] = np.nan
    mpa_df.dropna(how='all', axis=1, inplace=True)
    mpa_df.fillna(thresholdMpa, inplace=True)
    samples = samples_df[samples_df['NameOfFigure'] == 'Human'].index.values.tolist()
    samples_df = samples_df.loc[samples]
    mpa_df = mpa_df.loc[list(set(samples_df.UserID.values))]
    presence = (mpa_df > thresholdMpa + 1e-7).astype(int).sum()
    presence = presence[presence > len(mpa_df.index) * ratio].index.values.tolist()
    mpa_df = mpa_df[presence]
    mpa_df=mpa_df.loc[list(set(samples_df.UserID.values))]
    genotype_samples=[]
    if genotype == 'ALS':
        genotype_samples += samples_df[samples_df['Genotype'] == 'ALS'].index.values.tolist()
    if genotype == 'Control':
        genotype_samples += samples_df[samples_df['Genotype'] == 'Control'].index.values.tolist()
    mpa_df = mpa_df.loc[genotype_samples]
    mpa_df = np.log10(mpa_df)
    return mpa_df

