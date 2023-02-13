#Data path: /oak/stanford/groups/pritch/users/daphna/snoRNA/pivus_snoRNA_data
import matplotlib.pyplot as plt
import glob
from Bio import SeqIO
import gzip

read_count = []
for sample in glob.glob('*R1*.fastq.gz'):
    with gzip.open(sample,'rt') as f_h:
        read_count.append(len(list(SeqIO.parse(f_h,'fastq'))))

plt.hist(read_count,30)
plt.xlabel('Number of snoRNA mapped reads')
plt.ylabel('Count')
plt.savefig('HistSnoRNAReadCount.pdf',format='pdf')