**Phylogenetic analysis with COI sequences**

Whole genome sequencing (WGS) data of 21 *C. angulata* and 472 *C. gigas* were downloaded from the Sequence Read Archive (SRA) database under accession numbers PRJNA394055 (Li et al., 2018). WGS data of 261 *C. ariakensis* (69 southern samples and 192 northern samples) was downloaded from the SRA database under accession numbers PRJNA715058 (Wu et al., 2022).

**Data trimming**
```bash
for i in `cat sample.list`
do
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fq.gz -O ${i}_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
done
```
**Map reads to reference geneome**
```bash
bwa index ref.fa
for i in `cat sample.list`
do
bwa mem -M -t 4 ref.fa ${i}_clean_1.fq.gz ${i}_clean_2.fq.gz | samtools view -bS > $i.bam
# mark duplications with samtools v1.7
samtools sort -@ 8 -n -o namesort.bam $i.bam && \
rm -f $i.bam && \
samtools fixmate -@ 8 -m namesort.bam fixmate.bam && \
rm -f namesort.bam && \
samtools sort -@ 8 -o positionsort.bam fixmate.bam && \
rm -f fixmate.bam && \
samtools markdup -@ 8 -r positionsort.bam $i.dup.bam && \
rm -f positionsort.bam && \
samtools index -@ 8 $i.dup.bam
# remove dup and low quality mapped reads
samtools view -@ 8 -h -F 0x100 -F 0x400 $i.dup.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 8 -q 30 -bS > $i.bam && \
samtools index $i.bam && \
rm -f $i.dup.bam
done
````
**Assemble the unmapped reads**
```bash
for i in `cat sample.list`
do
mkdir $i
cd $i
samtools fastq -@ 10 -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq $i.dup.bam
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq --careful --cov-cutoff auto -o spades_assembly -t 30
cd ../
done
```
**Extract COI sequences**
```bash
for i in `cat sample.list`
do
cd $i
makeblastdb -in scaffolds.fasta -dbtype nucl
blastn -query FJ841964.1.fa -db scaffolds.fasta -out $i.blt -outfmt 6
head -n 1 $i.blt > $i.head1.blt
coi=`sed -n -e 1p $i.head1.blt | awk '{print $2}'`
start=`sed -n -e 1p $i.head1.blt | awk '{print $9}'`
end=`sed -n -e 1p $i.head1.blt | awk '{print $10}'`
samtools faidx scaffolds.fasta && samtools faidx scaffolds.fasta ${coi}:${start}-${end} > ${i}_coi.dna.fa
python cds2prot.py ${i}_coi.dna.fa ${i}_coi.prot.fa ${i}_coi.nucl.fa
cd ../
done
```
**Multiple alignment with MUSCLE v3.8.31 (Edgar & Soc, 2004).**
```bash
muscle -in merge.fa -out merge.afa
```
****
convertFasta2Phylip.sh merge.afa > merge.phy
raxmlHPC-PTHREADS -f a -# 100 -m GTRGAMMA -p 12345 -x 12345 -s merge.phy -n merge.tree -T 30 


**Reference**

Edgar, R. C. (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5, 1-19. https://doi.org/10.1186/1471-2105-5-113

Li, L., Li, A., Song, K., Meng, J., Guo, X., Li, S., ... Zhang, G. (2018). Divergence and plasticity shape adaptive potential of the Pacific oyster. Nature Ecology & Evolution, 2, 1751-1760. https://doi.org/10.1038/s41559-018-0668-2

Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. Current protocols in bioinformatics, 70, e102. https://doi.org/10.1002/cpbi.102

Wu, B., Chen, X., Yu, M. J., Ren, J. F., Hu, J., Shao, C. W., ... Liu, Z. H. (2022). Chromosome-level genome and population genomic analysis provide insights into the evolution and environmental adaptation of Jinjiang oyster Crassostrea ariakensis. Molecular Ecology Resources, 22, 1529-1544. https://doi.org/10.1111/1755-0998.13556
