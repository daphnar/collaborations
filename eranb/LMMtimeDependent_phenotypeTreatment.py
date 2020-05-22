import pandas
import statsmodels.formula.api as smf
from mne.stats.multi_comp import fdr_correction
import glob,os

data_files='*.csv'
output_file='LMM_time_dependent_analysis.xlsx'
pvals_genotypeMtreatment=[]
pvals_timeMgenotypeMtreatment=[]
stats=[]
sheets = glob.glob(data_files)
for sheetname in sheets:
    data = pandas.read_csv(sheetname)
    data['genotypeMtreatment']=data['genotype2'].mul(data['treatment2'])
    data['timeMgenotype']=data['genotype2'].mul(data['Time'])
    data['timeMtreatment']=data['treatment2'].mul(data['Time'])
    data['timeMgenotypeMtreatment'] = data['timeMgenotype'].mul(data['treatment2'])
    if (data['Time']==data['timeMgenotype']).all():
        stat = smf.mixedlm("phenotype ~ Time+treatment2+timeMtreatment", data, groups=data["MouseID"]).fit()#treatment2
        pvals_genotypeMtreatment.append(stat.pvalues['treatment2'])
        pvals_timeMgenotypeMtreatment.append(stat.pvalues['timeMtreatment'])
    else:
        stat = smf.mixedlm("phenotype ~ Time+genotypeMtreatment+timeMgenotype+timeMtreatment+timeMgenotypeMtreatment", data, groups=data["MouseID"]).fit()
        pvals_genotypeMtreatment.append(stat.pvalues['genotypeMtreatment'])
        pvals_timeMgenotypeMtreatment.append(stat.pvalues['timeMgenotypeMtreatment'])
    stats.append(stat)

qvals_genotypeMtreatment=fdr_correction(pvals_genotypeMtreatment)[1]
qvals_timeMgenotypeMtreatment=fdr_correction(pvals_timeMgenotypeMtreatment)[1]

passing = [(os.path.basename(sheets[i]),qvals_genotypeMtreatment[i],qvals_timeMgenotypeMtreatment[i],stats[i]) for i in range(len(qvals_genotypeMtreatment)) if qvals_genotypeMtreatment[i]<0.05 or \
           qvals_timeMgenotypeMtreatment[i]<0.05]

pandas.DataFrame({'Test':[os.path.basename(sheets[i])[:-4] for i in range(len(sheets))],
                  'P-value Genotype*treatment':pvals_genotypeMtreatment,
                  'Q-value Genotype*treatment':qvals_genotypeMtreatment,
                  'P-value time*genotype*treatment':pvals_timeMgenotypeMtreatment,
                  'Q-value time*genotype*treatment':qvals_timeMgenotypeMtreatment}).set_index('Test').to_excel(output_file)