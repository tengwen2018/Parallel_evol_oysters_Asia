**Genome assembly of Fujian oyster (*Crassostrea angulata*), also named Portuguese oyster.**

Teng Wen 14/07/2023

**1. Genome size estimation with Illumina 150bp PE reads**

Data trimming (All Illumina 150bp PE reads used in the following analysis were trimmed with this criteria)
```bash
fastp -i input_1.fastq.gz -I input_2.fastq.gz -o output_clean_1.fq.gz -O output_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
```
Genome size and complexity were estimated with Jellyfish (Marçais & Kingsford, 2011) and GenomeScop (Vurture et al., 2017).
```bash
for i in 17 19 21 23 25 27 29 31 33 35 37 39 41
do
jellyfish count -m $i -s 1000M -t 10 -C *.fq.gz -o kmer$i.reads.jf
jellyfish histo kmer$i.reads.jf > kmer$i.reads.histo
Rscript ~/tools/genomescope-master/genomescope.R kmer$i.reads.histo $i 150 genomescope_kmer$i
done
```
**2. De novo assembly**

Primary genome was assembled by hifiasm v0.16 with default parameters (Cheng et al., 2021).
```bash
hifiasm -o output.asm -t32 input1.hifi.fastq.gz input2.hifi.fastq.gz
```
Haploid sub-assembly was reconstructed by HaploMerger2 (Huang et al., 2017).
```bash
windowmasker \
    -checkdup true \
    -mk_counts \
    -in input.fasta \
    -out masking-library.ustat \
    -mem 200000

windowmasker \
  -ustat masking-library.ustat \
  -in input.fasta \
  -out contigs_wm.fa \
  -outfmt fasta \
  -dust true

gzip contigs_wm.fa

sh ./hm.batchA1.initiation_and_all_lastz contigs_wm
sh ./hm.batchA2.chainNet_and_netToMaf contigs_wm
sh ./hm.batchA3.misjoin_processing contigs_wm

sh ./hm.batchA1.initiation_and_all_lastz contigs_wm_A 
sh ./hm.batchA2.chainNet_and_netToMaf contigs_wm_A 
sh ./hm.batchA3.misjoin_processing contigs_wm_A

sh ./hm.batchA1.initiation_and_all_lastz contigs_wm_A_A
sh ./hm.batchA2.chainNet_and_netToMaf contigs_wm_A_A
sh ./hm.batchA3.misjoin_processing contigs_wm_A_A

sh hm.batchB1.initiation_and_all_lastz contigs_wm_A_A_A
sh hm.batchB2.chainNet_and_netToMaf contigs_wm_A_A_A
sh hm.batchB3.haplomerger contigs_wm_A_A_A
sh hm.batchB4.refine_unpaired_sequences contigs_wm_A_A_A

sh hm.batchB5.merge_paired_and_unpaired_sequences contigs_wm_A_A_A
```
The chromosome‐level assembly was performed by Juicer v1.6 with default parameters (Durand et al., 2016).
```bash
bwa index contigs.fa
python ~/tools/juicer/misc/generate_site_positions.py DpnII contigs contigs.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' contigs_DpnII.txt > contigs.size
~/tools/juicer/scripts/juicer.sh \
   -g asm.contigs \
   -s DpnII \
   -z reference/contigs.fa \
   -y reference/contigs_DpnII.txt \
   -p reference/contigs.size \
   -D /public/home/sundongf/tools/juicer \
   -t 20
cd ./i15k_r0
~/tools/3d-dna/run-asm-pipeline.sh -i 15000 -r 0 ../reference/contigs.fa ../aligned/merged_nodups.txt
```
Genome assembly completeness was accessed by BUSCO v5.1.3 (Simão et al., 2015) using the “mollusca_odb10” BUSCO gene set collection (Waterhouse et al., 2018).

```bash
busco --offline -f -c 10 -i input.fa -l mollusca_odb10 -o output_busco -m genome
```
Here is what the outpu looked like:
```
# BUSCO version is: 5.1.2 
# The lineage dataset is: mollusca_odb10 (Creation date: 2020-08-05, number of genomes: 7, number of BUSCOs: 5295)
# Summarized benchmarking in BUSCO notation for file /yourpath/crassostrea_angulata.fna
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

	***** Results: *****

	C:99.0%[S:98.1%,D:0.9%],F:0.4%,M:0.6%,n:5295	   
	5242	Complete BUSCOs (C)			   
	5193	Complete and single-copy BUSCOs (S)	   
	49	Complete and duplicated BUSCOs (D)	   
	21	Fragmented BUSCOs (F)			   
	32	Missing BUSCOs (M)			   
	5295	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.3
	metaeuk: 4.a0f584d
```

**Reference**

Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature Methods, 18, 170-175. https://doi.org/10.1038/s41592-020-01056-5

Durand, N. C., Robinson, J. T., Shamim, M. S., Machol, I., Mesirov, J. P., Lander, E. S., & Aiden, E. L. (2016). Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom. Cell systems, 3, 99-101. https://doi.org/10.1016/j.cels.2015.07.012

Huang, S., Kang, M., & Xu, A. (2017). HaploMerger2: rebuilding both haploid sub-assemblies from high-heterozygosity diploid genome assembly. Bioinformatics, 33, 2577-2579. https://doi.org/10.1093/bioinformatics/btx220

Marçais, G., & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27, 764-770. https://doi.org/10.1093/bioinformatics/btr011

Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31, 3210-3212. https://doi.org/10.1093/bioinformatics/btv351

Vurture, G. W., Sedlazeck, F. J., Nattestad, M., Underwood, C. J., Fang, H., Gurtowski, J., & Schatz, M. C. (2017). GenomeScope: fast reference-free genome profiling from short reads. Bioinformatics, 33, 2202-2204. https://doi.org/10.1093/bioinformatics/btx153

Waterhouse, R. M., Seppey, M., Simão, F. A., Manni, M., Ioannidis, P., Klioutchnikov, G., ... Zdobnov, E. M. (2018). BUSCO applications from quality assessments to gene prediction and phylogenomics. Molecular Biology and Evolution, 35, 543-548. https://doi.org/10.1093/molbev/msx319
