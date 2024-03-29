import pandas as pd
import os
import glob
def unite_kallisto_abundances(input_folder_re,output_ra):
    results = []
    for folder in glob.glob(input_folder_re):
        try:
            abundance = pd.read_csv(os.path.join(folder,'abundance.tsv'),sep='\t',index_col=0)[['eff_length','tpm']]
            abundance = abundance['tpm']/abundance['eff_length']
        except Exception:#FileNotFoundError: #EmptyDataError
            print(folder)
            continue
        ra = abundance/abundance.sum()
        ra.name = os.path.basename(folder)
        results.append(ra)
    pd.concat(results,axis=1).T.to_csv(output_ra)

def clean_abundances(input_abundances,output_abundances,min_count=20,sequencing_depth=750):
    obj = pd.read_csv(input_abundances, index_col=0)
    min_abundance=min_count / sequencing_depth
    obj[obj < min_abundance] = min_abundance
    #Throw columns with little information
    obj = obj.loc[:, obj.std() > 0.0001]
    obj.to_csv(output_abundances)

kallisto_index = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/kallisto/snoRNA'
kallisto_tool = '/oak/stanford/groups/pritch/users/daphna/tools/kallisto/kallisto'
sample_data_path = '/oak/stanford/groups/pritch/users/daphna/snoRNA/data/pivus/pivus_snoRNA_data_by_sample_readlen48_rarefaction750'
kallisto_output_folder = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/pivus/kallisto/'

chunks_file = '/oak/stanford/groups/pritch/users/daphna/snoRNA/data/pivus/samples.csv'
files = pd.read_csv(chunks_file, index_col=0)
for sample in files.index:
    sample_file=sample+'.fastq.gz'
    command = '%s quant -i %s -o %s --single -l 75 -s 0.01 %s'%(kallisto_tool,kallisto_index,
            os.path.join(kallisto_output_folder,sample), os.path.join(sample_data_path,sample_file))
    os.system(command)

input_folder_re = os.path.join(kallisto_output_folder,'PIVUS*')
output_ra = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/pivus/pivus_abundances.csv'
unite_kallisto_abundances(input_folder_re,output_ra)
output_norm = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/pivus/pivus_clean_abundances.csv'
clean_abundances(output_ra,output_norm)

