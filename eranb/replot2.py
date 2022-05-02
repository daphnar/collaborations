from eranb import plot_PCoA
from eranb.Extract import extract

sample_annotation='/Users/daphna/Downloads/reproduce/HumanExtract_revisions2.xlsx'
mpa_path='/Users/daphna/Downloads/reproduce/DFOut_revisions2'

genotype=('Control','ALS')
control = extract(genotype[0],sample_annotation,mpa_path)
als = extract(genotype[1],sample_annotation,mpa_path)
plot_PCoA.plot(control, als, legend=genotype, method='PCA')