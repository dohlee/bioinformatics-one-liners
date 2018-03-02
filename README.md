# bioinformatics-one-liners
Useful bioinformatics one-line commands.

## samtools [(doc)](http://www.htslib.org/doc/samtools.html)
sort sam/bam files
```shell
samtools sort SRR1234567/SRR1234567.bam -o SRR1234567/SRR1234567.sorted.bam --threads 4
```

index sam/bam files (this creates SRR1234567.sorted.bam.bai file)
```shell
samtools index SRR1234567/SRR1234567.sorted.bam
```

count the number of reads in sam/bam files
```shell
samtools view -c SRR1234567/SRR1234567.bam
```

count the number of mapped reads in sam/bam files
```shell
samtools view -F 4 -c SRR1234567/SRR1234567.bam
```

count the number of unmapped reads in sam/bam files
```shell
samtools view -f 4 -c SRR1234567/SRR1234567.bam
```

show statistics of sam/bam files
```shell
samtools stats SRR1234567/SRR1234567.bam
```
## Dealing with fastq files
count the number of reads in fastq files
```shell
cat SRR1234567/SRR1234567.fastq | echo $((`wc -l` / 4))
```

count the number of reads in gzipped fastq files
```shell
zcat SRR1234567/SRR1234567.fastq.gz | echo $((`wc -l` / 4))
```

## cut
print only the second column of the tab-delimited file
```shell
cat file.txt | cut -f 2
```

convert tsv to csv
```shell
cat file.tsv | cut -f - --output-delimiter , > file.csv
```
## Linux basics
show the size of subdirectories
```shell
du -h
```

only show the size of current directory
```shell
du -sh
```
## Statistics
Multiple testing correction
```
from statsmodels.stats import multitest
reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, method='fdr_bh')
```
