import gzip
import numpy as np
from Bio import SeqIO
import os
import glob
import sys
import pandas as pd
import os

def trim_reads_to_equi_length(finput,foutput,readlen=48):
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

def rarefaction_dataset(sample,indir,outdir):
    return_val=0
    if os.path.exists('%s/%s_merge_R1.fastq.gz' % (indir,sample)) and os.path.exists('%s/%s_merge_R2.fastq.gz' %  (indir,sample)):
        os.system(
        'cat %s/%s_1.fastq.gz %s/%s_merge_R1.fastq.gz > %s/%s.fastq.gz' % (indir,sample, indir,sample,outdir,sample))
        return_val=1
    else:
        if os.path.exists('%s/%s_merge_R1.fastq.gz' % (indir,sample)) and not os.path.exists('%s/%s_merge_R2.fastq.gz' % (indir,sample)):
            os.system('cat %s/%s_merge_R1.fastq.gz > %s/%s.fastq.gz' % ( indir,sample,outdir,sample))
            return_val=0
        elif os.path.exists('%s/%s_merge_R2.fastq.gz' % (indir,sample)) and not os.path.exists(
            '%s/%s_merge_R1.fastq.gz' % (indir,sample)):
            os.system('cat %s/%s_merge_R2.fastq.gz > %s/%s.fastq.gz'% ( indir,sample,outdir,sample))
            return_val = 0
    return return_val

if __name__=='__main__':
    chunk_count=0
    if False:
        indir='/home/users/daphna/data/GTEx/sequencing'
        chunks_file = os.path.join(indir,'lengths.csv')
        chunk=int(sys.argv[1])
        if len(sys.argv)>2:
            chunk_size = int(sys.argv[2])
        else:
            chunk_size = 50
        files = pd.read_csv(chunks_file,index_col=1)
        outdir='/home/users/daphna/data/GTEx/%s'
        for fin in files.iloc[chunk*chunk_size:(chunk+1)*chunk_size,:].index:
            finput = os.path.join(indir,fin)
            foutput=os.path.join(outdir,os.path.basename(finput))
            chunk_count+=rarefaction(finput, foutput)
    if True:
        indir = '/home/users/daphna/data/GTEx/sequencing_rarefaction_to_250K_throw_threshold_100K_ink'
        outdir = '/home/users/daphna/data/GTEx/sequencing_rarefaction_to_500K_throw_threshold_100K_ink'
        chunks_file = '/home/users/daphna/data/GTEx/sequencing_rarefaction_to_250K_throw_threshold_100K_ink/sample_names.csv'
        chunk = int(sys.argv[1])
        if len(sys.argv) > 2:
            chunk_size = int(sys.argv[2])
        else:
            chunk_size = 50
        files = pd.read_csv(chunks_file, index_col=0)
        for sample in files.iloc[chunk * chunk_size:(chunk + 1) * chunk_size, :].index:
            chunk_count+=unit_GTEx_paired_end(sample,indir, outdir)
    if False:
        indir = '/home/users/daphna/data/PIVUS/sequencing'
        outdir = '/home/users/daphna/data/PIVUS/sequencing_SE'
        chunks_file = '/home/users/daphna/data/PIVUS/samples.csv'
        chunk = int(sys.argv[1])
        if len(sys.argv) > 2:
            chunk_size = int(sys.argv[2])
        else:
            chunk_size = 50
        files = pd.read_csv(chunks_file, index_col=0)
        for sample in files.iloc[chunk * chunk_size:(chunk + 1) * chunk_size, :].index:
            chunk_count+=unit_PIVUS_paired_end(sample,indir, outdir)
    if False:
        # indir='/home/users/daphna/data/PIVUS/sequencing_SE'
        # chunks_file = '/home/users/daphna/data/PIVUS/samples.csv'
        # outdir = '/home/users/daphna/data/PIVUS/sequencing_rarefaction_to_124K/'

        # indir = '/home/users/daphna/data/DGN/sequencing_SE'
        # chunks_file = '/home/users/daphna/data/DGN/samples.csv'
        # # outdir = '/home/users/daphna/data/PIVUS/sequencing_rarefaction_to_250K_throw_threshold_100K/'
        # outdir = '/home/users/daphna/data/DGN/sequencing_rarefaction_to_441K/'
        # indir = '/home/users/daphna/data/TCGA/WGS/no_rarefaction'
        # chunks_file = '/home/users/daphna/data/TCGA/WGS/samples.csv'
        # # outdir = '/home/users/daphna/data/PIVUS/sequencing_rarefaction_to_250K_throw_threshold_100K/'
        # outdir = '/home/users/daphna/data/TCGA/WGS/readlen48_no_rarefaction/'
        # outdir = '/home/users/daphna/data/PIVUS/sequencing_rarefaction_to_250K_throw_threshold_100K/'
        chunk=int(sys.argv[1])
        if len(sys.argv)>2:
            chunk_size = int(sys.argv[2])
        else:
            chunk_size = 50
        if len(sys.argv)>3:
            readlen = int(sys.argv[3])
        else:
            readlen = 48
        if len(sys.argv)>4:
            chunks_file = sys.argv[4]
        else:
            chunks_file = '/home/users/daphna/data/TCGA/samples.csv'
        if len(sys.argv)>5:
            indir = sys.argv[5]
        else:
            indir = '/home/users/daphna/data/TCGA/sequencing_rarefaction_to_250K_throw_threshold_100K'
        if len(sys.argv)>6:
            outdir = sys.argv[6]
        else:
            outdir = '/home/users/daphna/data/TCGA/mrna_readlen48_rarefaction/'
        files = pd.read_csv(chunks_file,index_col=0)
        for fin in files.iloc[chunk*chunk_size:(chunk+1)*chunk_size,:].index:
            finput = os.path.join(indir,fin)
            if not finput.endswith('.fastq.gz'):
                finput+='.fastq.gz'
            foutput=os.path.join(outdir,os.path.basename(finput))
            if os.path.exists(finput):
                #chunk_count+=rarefaction(finput, foutput,keep_n=442000,throw_below=441000)
                #chunk_count += rarefaction(finput, foutput, keep_n=50000, throw_below=40000)
                chunk_count += rarefaction(finput, foutput, keep_n=250000, throw_below=100000)
            else:
                chunk_count+=1
            #chunk_count+=trim_reads_to_equi_length(finput,foutput,readlen=readlen)
    if chunk_count==chunk_size:
        os.system('touch %s'%os.path.join(outdir,'%s.readlen.done'%chunk))


finput2=finput+'_merge_R2.fastq.gz'
rarefaction(finput2,foutput2,keep_n=750,throw_below=0)