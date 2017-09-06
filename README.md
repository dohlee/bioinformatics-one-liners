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
