import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import pandas as pd
import glob,os
from eranb import config
from eranb.PCA_association import getPhenotypes

def run_ols_model(time,genotype,phenotype):
    data = pd.DataFrame({'time':time,'genotype':genotype,'phenotype':phenotype})
    data['timeMgenotype']=data['genotype'].mul(data['time'])
    stat = smf.ols("phenotype ~ time+genotype+timeMgenotype", data).fit()
    stats = stat.pvalues.drop(['time','genotype[T.1]','Intercept'])
    pvals_time = stat.pvalues['time']
    pvals_genotype = stat.pvalues['genotype[T.1]']
    pvals_timeMgenotype = min(stats)
    return pvals_time,pvals_genotype,pvals_timeMgenotype


def run_time_dependence(output_file):
    pvals_times=[]
    pvals_genotypes=[]
    pvals_timeMgenotypes=[]
    hypotheses=[]
    sheet = "Metadata_LumMFIav"
    for data, genotypes, sheetname in getPhenotypes():
        if sheetname!=sheet:
            continue

    run_groups = [(config.time_behavioral,config.behavioural_tests),
                  (config.time_blood,config.blood_tests)]
    for time,phenotypes in run_groups:
        for phenotype in phenotypes:
            hypotheses.append("%s~%s + %s + %s"%(phenotype,time,'genotype',time +'*genotype'))
            pvals_time,pvals_genotype,pvals_timeMgenotype=run_ols_model(data[time].values,genotypes.values,data[phenotype].values)
            pvals_times.append(pvals_time)
            pvals_genotypes.append(pvals_genotype)
            pvals_timeMgenotypes.append(pvals_timeMgenotype)
    qvals_times=multipletests(pvals_times,method='fdr_bh')[1]
    qvals_genotypes=multipletests(pvals_genotypes,method='fdr_bh')[1]
    qvals_timeMgenotypes=multipletests(pvals_timeMgenotypes,method='fdr_bh')[1]
    res = pd.DataFrame({'Test':hypotheses,
                  'P-value time':pvals_times,
                  'Q-value time': qvals_times,
                  'P-value genotype':pvals_genotypes,
                  'Q-value genotype': qvals_genotypes,
                  'P-value time*genotype':pvals_timeMgenotypes,
                  'Q-value time*genotype':qvals_timeMgenotypes}).set_index('Test')

    res.sort_values(by='Q-value time*genotype').to_csv(output_file)

if __name__=='__main__':
    output_file=os.path.join(config.analyses_path,'time_dependent_analysis.csv')
    run_time_dependence(output_file)