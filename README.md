# bioinformatics-one-liners
Useful bioinformatics one-line commands.

## Indexing files
Indexing gtf files with tabix can be tricky. (from Book [Bioinformatics Data Skills](https://www.amazon.com/Bioinformatics-Data-Skills-Reproducible-Research/dp/1449367372))
```shell
(zgrep -v "^#" gencode.v32.annotation.gtf.gz | sort -k1,1 -k4,4n) | bgzip > gencode.v32.annotation.gtf.bgz \
&& tabix -p gff gencode.v32.annotation.gtf.bgz
```

## Fetch annotations
get CpG Island annotation in BED format ([ref](https://www.biostars.org/p/236141/), [BEDOPS](https://bedops.readthedocs.io/en/latest/) needed)
```shell
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
  | gunzip -c \
  | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, substr($0, index($0, $7)); }' \
  | sort-bed - \
  > cpgIslandExt.hg38.bed
```
## samtools [(doc)](http://www.htslib.org/doc/samtools.html)
sort sam/bam files
```shell
samtools sort SRR1234567/SRR1234567.bam -o SRR1234567/SRR1234567.sorted.bam --threads 4
```

index sam/bam files (this creates SRR1234567.sorted.bam.bai file)
```shell
samtools index SRR1234567/SRR1234567.sorted.bam
```

get alignments on chromosome 20 only
```shell
samtools view SRR1234567/SRR1234567.bam chr20 -b > SRR1234567/SRR1234567_chr20.bam
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

concatenate multiple fastq.gz files (just use cat!)
```shell
cat sample1.fastq.gz sample2.fastq.gz sample3.fastq.gz > concatenated.fastq.gz
```

fastq-dump best practice
```shell
fastq-dump --split-3 --skip-technical --gzip --readids --clip --read-filter pass <SRR->
```

use ascp to prefetch sra files ([aspera-connect](https://downloads.asperasoft.com/) required, and use absolute path for aspera)
```shell
prefetch --ascp-path \
'/path/to/home/.aspera/connect/bin/ascp|/path/to/home/.aspera/connect/etc/asperaweb_id_dsa.openssh' \
--ascp-options -l150M SRR1234567
```
## grep
print files in current directory which contain a pattern
```shell
grep -l <pattern> *
```
## sed
extract desired part of each line using substitution
```shell
cat file | sed -e 's/.*\([0-9]*\).*/\1/'
```
## join [(tutorial)](https://shapeshed.com/unix-join/)
join lines of two files on a common field in the first column
```shell
join file1.txt file2.txt > out.txt
```
## comm (mnemonics: `common`)
print only lines common to both files
```shell
comm -12 <(sort file1) <(sort file2)
```

print only lines unique to first file
```shell
comm -23 <(sort file1) <(sort file2)
```

print only lines unique to second file
```shell
comm -13 <(sort file1) <(sort file2)
```

print only lines common to three files
```shell
comm -12 <(sort file1) <(sort file2) | comm -12 <(sort file3) -
```
## tr (mnemonics: `translate`)
get complementary DNA sequence
```shell
echo ATGCTGTAGTC | rev | tr 'ACGT' 'TGCA'
```
## cut
print only the second column of the tab-delimited file
```shell
cat file.txt | cut -f 2
```

print only the first column of the comma-delimited file
```shell
cat file.txt | cut -d',' -f 1
```

print the first and second columns of the comma-delimited file
```shell
cat file.txt | cut -d',' -f 1,2
```

convert tsv to csv
```shell
cat file.tsv | cut -f - --output-delimiter , > file.csv
```

print unique values in second column
```shell
cat file.tsv | cut -f2 | uniq
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

print a file except first line
```shell
cat file.txt | tail -n +2
```

prepend a line to a file
```shell
sed '1i header' file.txt > prepended_file.txt
```

iterate over files in a directory
```shell
for f in *.sra; do fastq-dump $f; done
```

show from the a-th line to the b-th line of a file
```shell
head file -n $b | tail -n $(($b-$a+1))
```
## Statistics
Multiple testing correction
```python
from statsmodels.stats import multitest
reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(p_values, method='fdr_bh')
```
