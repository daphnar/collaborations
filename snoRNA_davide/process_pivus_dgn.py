import os
import sys
from Bio import SeqIO,SeqRecord,Seq,bgzf
import pysam
from datetime import datetime
import subprocess
import pandas as pd

def run_shell_command(command,debug=True):
    if debug:
        print('Running: %s'%command)
        print('Starting: %s\n\n'%datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        output = subprocess.check_output(command, shell = True)
        print('End: %s\n\n'%datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    return output


def map_and_filter_sample_to_ribosomes(sample_path,output_directory,bowtie_executable,bowtie_index,samtools_executable,
                                       to_remove=True,is_zipped=False,max_keep_reads=10000000,fastq_output=None):
    mapped_sample=os.path.join(output_directory,os.path.basename(sample_path)+'.sam')
    sorted_mapped_sample=os.path.join(output_directory,os.path.basename(sample_path)+'.sorted.sam')
    #For fasta
    # command = '%s -x %s -f %s --threads 16 -S %s' % (
    #     bowtie_executable, bowtie_index, sample_path, mapped_sample)
    command = '%s --threads 2 -x %s -U %s -S %s' % (bowtie_executable, bowtie_index, sample_path, mapped_sample)
    run_shell_command(command)
    #No meed to sort - save time and space of not writing another file!
    # command = ('%s sort %s -o %s' % (samtools_executable, mapped_sample, sorted_mapped_sample))
    # run_shell_command(command)
    savereads=[]
    #mapiter = pysam.AlignmentFile(sorted_mapped_sample, 'r')
    mapiter = pysam.AlignmentFile(mapped_sample, 'r')
    for read in mapiter.fetch():
        if read.reference_name is None:
            continue
        savereads.append(read)
        if len(savereads)>=max_keep_reads:
            break
    fastq_records=[]
    for s in savereads:
        record = SeqRecord.SeqRecord(Seq.Seq(s.seq), id=s.qname, description='')
        record.letter_annotations["phred_quality"] = s.query_qualities
        fastq_records.append(record)
    if fastq_output is None:
        fastq_output = os.path.join(output_directory,os.path.basename(sample_path))
    # if is_zipped:
    #     fastq_output=fastq_output.replace('.gz','')
    # with open(fastq_output,'w') as f_h:
    if not fastq_output.endswith('.gz'):
        fastq_output = fastq_output+'.gz'
    with bgzf.BgzfWriter(fastq_output, "wb") as f_h:
        SeqIO.write(fastq_records,f_h,'fastq')
    if to_remove:
        #cmd = 'rm -rf %s %s'%(mapped_sample, sorted_mapped_sample)
        cmd = 'rm -rf %s' % mapped_sample
        run_shell_command(cmd)
    return fastq_output

def run_sample(sample_fastq_path, scratch_base,bowtie_executable,bowtie_index,samtools_executable,to_remove,dataset,rename_sample_id=None):
    output_ribo_filter_folder = os.path.join(scratch_base, 'daphna/rdna/analyses/ribo_filter/%s'%dataset)
    is_zipped=sample_fastq_path.endswith('.gz')
    if rename_sample_id is None:
        sample_id=os.path.basename(sample_fastq_path)
    else:
        sample_id = os.path.basename(sample_fastq_path)
        sample_id = os.path.join(os.path.dirname(sample_fastq_path),rename_sample_id+sample_id[sample_id.find('.'):])
    # if is_zipped:
    #     sample_id=sample_id.replace('.gz','')
    if not os.path.exists(output_ribo_filter_folder):
        os.makedirs(output_ribo_filter_folder)
    ribosome_reads = os.path.join(output_ribo_filter_folder, sample_id)
    if not os.path.exists(ribosome_reads):
        print(ribosome_reads)
        map_and_filter_sample_to_ribosomes(sample_fastq_path, output_ribo_filter_folder, bowtie_executable, bowtie_index,
            samtools_executable, to_remove,is_zipped,fastq_output=ribosome_reads)

if __name__=='__main__':
    sample_fastq_path = sys.argv[1]
    if len(sys.argv)>2:
        dataset = sys.argv[2]
    else:
        dataset = 'pivus'
    if len(sys.argv)>3:
        input_metadata = sys.argv[3]
        if input_metadata!='None':
            sample_id=sample_fastq_path
            sample_fastq_path=pd.read_csv(input_metadata,sep='\t',index_col=0,header=None).loc[sample_fastq_path,6]

    if len(sys.argv) > 4:
        to_remove = eval(sys.argv[4])
    else:
        to_remove = True
    if len(sys.argv)>5:
        bowtie_index = sys.argv[5]
    else:
        #/oak/stanford/groups/pritch/users
        bowtie_index = '/scratch/groups/pritch/daphna/rdna/analyses/' \
                       'URA/URA_DB_ribosomes_all_GIAB_T2T/bowtie_index/index'
    if len(sys.argv) > 6:
        scratch_base = sys.argv[6]
    else:
        scratch_base = '/scratch/groups/pritch'
    if len(sys.argv)>7:
        bowtie_executable = sys.argv[7]
    else:
        bowtie_executable = 'bowtie2'
    if len(sys.argv)>8:
        samtools_executable = sys.argv[8]
    else:
        samtools_executable = 'samtools'


    print(sys.argv)
    run_sample(sample_fastq_path, scratch_base, bowtie_executable, bowtie_index, samtools_executable, to_remove,dataset)


