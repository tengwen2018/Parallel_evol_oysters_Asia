**Phylogenetic analysis of southern-northern oyster species**

Whole genome sequencing (WGS) data of 21 *C. angulata* and 472 *C. gigas* were downloaded from the Sequence Read Archive (SRA) database under accession numbers PRJNA394055 (Li et al., 2018). WGS data of 261 *C. ariakensis* (69 southern samples and 192 northern samples) was downloaded from the SRA database under accession numbers PRJNA715058 (Wu et al., 2022).

**Data trimming**
```bash
fastp -i input_1.fastq.gz -I input_2.fastq.gz -o sample1_clean_1.fq.gz -O sample1_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
```
**Map reads to reference geneome**
```bash
bwa index ref.fa
bwa mem -M -t 4 ref.fa sample1_clean_1.fq.gz sample1_clean_2.fq.gz | samtools view -bS > sample1.bam
# mark duplications with samtools v1.7
samtools sort -@ 8 -n -o namesort.bam sample1.bam && \
rm -f sample1.bam && \
samtools fixmate -@ 8 -m namesort.bam fixmate.bam && \
rm -f namesort.bam && \
samtools sort -@ 8 -o positionsort.bam fixmate.bam && \
rm -f fixmate.bam && \
samtools markdup -@ 8 -r positionsort.bam sample1.dup.bam && \
rm -f positionsort.bam && \
samtools index -@ 8 sample1.dup.bam
# remove dup and low quality mapped reads
samtools view -@ 8 -h -F 0x100 -F 0x400 sample1.dup.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 8 -q 30 -bS > sample1.bam && \
samtools index sample1.bam && \
rm -f sample1.dup.bam
````
**Assemble the unmapped reads**
```bash
samtools fastq -@ 10 -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastqsample1.dup.bam
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq --careful --cov-cutoff auto -o spades_assembly -t 30
```

**Reference**

Li, L., Li, A., Song, K., Meng, J., Guo, X., Li, S., ... Zhang, G. (2018). Divergence and plasticity shape adaptive potential of the Pacific oyster. Nature Ecology & Evolution, 2, 1751-1760. https://doi.org/10.1038/s41559-018-0668-2

Wu, B., Chen, X., Yu, M. J., Ren, J. F., Hu, J., Shao, C. W., ... Liu, Z. H. (2022). Chromosome-level genome and population genomic analysis provide insights into the evolution and environmental adaptation of Jinjiang oyster Crassostrea ariakensis. Molecular Ecology Resources, 22, 1529-1544. https://doi.org/10.1111/1755-0998.13556
