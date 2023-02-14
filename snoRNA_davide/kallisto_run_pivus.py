import pandas as pd
import os
kallisto_index = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/kallisto/snoRNA'
kallisto_tool = '/oak/stanford/groups/pritch/users/daphna/tools/kallisto/kallisto'
sample_data_path = '/oak/stanford/groups/pritch/users/daphna/snoRNA/pivus_snoRNA_data_by_sample_readlen48_rarefaction750'
kallisto_output_folder = '/oak/stanford/groups/pritch/users/daphna/snoRNA/analyses/kallisto/'

chunks_file = '/oak/stanford/groups/pritch/users/daphna/snoRNA/samples.csv'
files = pd.read_csv(chunks_file, index_col=0)
for sample in files.index:
    sample_file=sample+'.fastq.gz'
    command = '%s quant -i %s -o %s --single -l 75 -s 0.01 %s'%(kallisto_tool,kallisto_index,
            os.path.join(kallisto_output_folder,sample), os.path.join(sample_data_path,sample_file))
    os.system(command)