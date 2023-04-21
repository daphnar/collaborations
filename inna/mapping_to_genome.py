#path_inna = '/Volumes/Favaro_HuBMAP/RNAseq_FFPE_bonemarrow'
#I created /Users/daphna/cluster2/users/daphna/inna with folders data and index_human.
#In index human there is the kallisto index found here:
#Kallisto index homo_sapiens download https://github.com/pachterlab/kallisto-transcriptome-indices/releases
#To download kallisto: (two lines)
#ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
#brew install kallisto

import pandas as pd
import os,glob

path_local = '/Users/daphna/cluster2/users/daphna/inna'
index_file = '/Users/daphna/cluster2/users/daphna/inna/index_human/transcriptome.idx'
fastq_files = glob.glob(os.path.join(path_local,'data','*_R1_001.fastq.gz'))
for pair_1 in fastq_files:
    pair_2 = pair_1.replace('_R1_','_R2_')
    basename = pair_1.split('_R1_')[0]
    kallisto_command = 'kallisto quant -i %s -o %s %s %s'%(index_file,basename,pair_1,pair_2)
    print(kallisto_command)
    os.system(kallisto_command)

for pair_1 in fastq_files:
    basename = pair_1.split('_R1_')[0]
    obj = pd.read_csv(os.path.join(basename,'abundance.tsv'), index_col=0, sep='\t')
    obj['length_eff_counts'] = (obj['est_counts'] / obj['eff_length'])
    obj['relative_expression'] = obj['length_eff_counts'] / obj['length_eff_counts'].sum()


