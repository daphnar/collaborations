import gzip
import numpy as np
from Bio import SeqIO
import os
import glob
import sys
import pandas as pd
import os

def trim_reads_to_equi_length(finput,foutput,readlen=75):
    with gzip.open(finput, "rt") as handle:
        try:
            new_records = [record[:readlen] for record in SeqIO.parse(handle, "fastq")]
        except Exception:
            print('%s file is bad' % finput)
            return 1
    with gzip.open(foutput, "wt") as handle:
        SeqIO.write(new_records, handle, 'fastq')
    return 1

def rarefaction(finput,foutput,keep_n=250000,throw_below=100000):
    with gzip.open(finput, "rt") as handle:
        try:
            records = list(SeqIO.parse(handle, "fastq"))
        except EOFError:
            print('%s file is bad'%finput)
            raise EOFError
    write=False
    if len(records)>keep_n:
        keep_reads = np.random.choice(len(records), size=keep_n, replace=False)
        keep_records = [records[i] for i in keep_reads]
        write=True
    elif len(records)>=throw_below:
        keep_records=records
        write = True
    if write:
        with gzip.open(foutput, "wt") as handle:
            SeqIO.write(keep_records,handle,'fastq')
    return 1

def unit_pivus_paired_end(sample,indir,outdir):
    return_val=0
    if os.path.exists('%s/%s_merge_R1.fastq.gz' % (indir,sample)) and os.path.exists('%s/%s_merge_R2.fastq.gz' %  (indir,sample)):
        os.system(
        'cat %s/%s_merge_R1.fastq.gz %s/%s_merge_R2.fastq.gz > %s/%s.fastq.gz' % (indir,sample, indir,sample,outdir,sample))
        return_val=1
    return return_val


if __name__=='__main__':
    chunk_count=0
    readlen = 75
    chunk_size=500
    chunk=0
    indir = '/oak/stanford/groups/pritch/users/daphna/snoRNA/pivus_snoRNA_data/'
    outdir_p = '/oak/stanford/groups/pritch/users/daphna/snoRNA/pivus_snoRNA_data_readlen48_rarefaction750'
    outdir_u = '/oak/stanford/groups/pritch/users/daphna/snoRNA/pivus_snoRNA_data_by_sample_readlen48_rarefaction750'
    chunks_file = '/oak/stanford/groups/pritch/users/daphna/snoRNA/samples.csv'
    files = pd.read_csv(chunks_file, index_col=0)
    for fin in files.iloc[chunk * chunk_size:(chunk + 1) * chunk_size, :].index:
        for p_end in [1,2]:
            finput = os.path.join(indir, fin)+'_merge_R%s.fastq.gz'%p_end
            foutput = os.path.join(outdir_p, os.path.basename(finput))
            print(finput,foutput)
            rarefaction(finput, foutput, keep_n=750, throw_below=0)
            trim_reads_to_equi_length(foutput,foutput,readlen=readlen)
        chunk_count+=unit_pivus_paired_end(fin,outdir_p, outdir_u)
