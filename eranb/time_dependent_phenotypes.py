import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import pandas as pd
import os
from eranb import config
from eranb.PCA_association import getPhenotypes
import matplotlib.pyplot as plt

def run_ols_model(time,genotype,phenotype):
    data = pd.DataFrame({'time':time,'genotype':genotype,'phenotype':phenotype}).dropna()
    data['timeMgenotype']=data['genotype'].mul(data['time'])
    data=data.astype(float)
    stat = smf.ols("phenotype ~ time+genotype+timeMgenotype-1", data).fit()
    pvals_time = stat.pvalues['time']
    pvals_genotype = stat.pvalues['genotype']
    pvals_timeMgenotype = stat.pvalues['timeMgenotype']
    return data,stat.params, pvals_time,pvals_genotype,pvals_timeMgenotype

def run_time_dependence(output_file):
    pvals_times=[]
    pvals_genotypes=[]
    pvals_timeMgenotypes=[]
    hypotheses=[]
    sheet = "Metadata_LumMFIav"
    for data, genotypes, sheetname in getPhenotypes():
        if sheetname==sheet:
            break
    data['genotype']=genotypes
    run_groups = [(config.time_behavioral,config.behavioural_tests),
                  (config.time_blood,config.blood_tests)]

    for time,phenotypes in run_groups:
        for phenotype in phenotypes:
            hypotheses.append("%s~%s + %s + %s"%(phenotype,time,'genotype',time +'*genotype'))
            model_data,coefs,pvals_time,pvals_genotype,pvals_timeMgenotype=run_ols_model(data[time].values,genotypes.values,data[phenotype].values)
            pvals_times.append(pvals_time)
            pvals_genotypes.append(pvals_genotype)
            pvals_timeMgenotypes.append(pvals_timeMgenotype)
            if pvals_timeMgenotype<0.05:
                plot_fit_ols(coefs,model_data,phenotype,pvals_timeMgenotype)
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

def plot_fit_ols(coefs,data,phenotype,pval):
    plt.figure()
    data['yfitted'] = data['genotype'].mul(coefs['genotype']).add(
        data['time'].mul(coefs['time'])).add( \
        data['timeMgenotype'].mul(coefs['timeMgenotype']))
    plt.plot(data[(data['genotype'] == 0)]['time'], \
             data[(data['genotype'] == 0)]['yfitted'], 'r')
    plt.scatter(data[(data['genotype'] == 0)]['time'], \
                data[(data['genotype'] == 0)]['phenotype'], c='r', alpha=0.25, s=50)
    plt.plot(data[(data['genotype'] == 1)]['time'], \
             data[(data['genotype'] == 1)]['yfitted'], 'g')
    plt.scatter(data[(data['genotype'] == 1)]['time'], \
                data[(data['genotype'] == 1)]['phenotype'], c='g', alpha=0.25, s=50)
    plt.legend(['WT','KO'])
    plt.xlabel('Time')
    plt.ylabel('Phenotype')
    plt.title(phenotype)
    plt.savefig(os.path.join(config.time_dependence_plots,phenotype.replace("/","$")+'pval_%.3f.pdf'%pval),format='pdf')
    plt.close()

if __name__=='__main__':
    output_file=os.path.join(config.analyses_path,'time_dependent_analysis.csv')
    run_time_dependence(output_file)