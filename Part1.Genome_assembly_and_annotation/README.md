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
Genome assembly completeness was accessed by BUSCO v5.1.3 (Simão et al., 2015) using the "mollusca_odb10" BUSCO gene set collection (Waterhouse et al., 2018)(https://busco-data.ezlab.org/v5/data/lineages/mollusca_odb10.2020-08-05.tar.gz).

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
Transposable elements de novo prediction and homology searching were performed by RepeatModeler v2.0.3 and RepeatMasker v4.1.1 (https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.1.tar.gz) pipeline (Smit AFA).
```bash
BuildDatabase -name ref input.fna
RepeatModeler -database ref -pa 20 -LTRStruct
RepeatMasker -e ncbi -lib ref-families.fa -gff -dir 00-RepeatMask_RepeatModeler -pa 15 -a input.fna 
```
```
==================================================
file name: input.fna                  
sequences:           412
total length:  624377014 bp  (624245514 bp excl N/X-runs)
GC level:         33.58 %
bases masked:  310143980 bp ( 49.67 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements        60645     36303191 bp    5.81 %
   SINEs:             2420       361749 bp    0.06 %
   Penelope           9215      3440732 bp    0.55 %
   LINEs:            36167     17074402 bp    2.73 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex        927       359090 bp    0.06 %
     R1/LOA/Jockey     189        53400 bp    0.01 %
     R2/R4/NeSL        180       144925 bp    0.02 %
     RTE/Bov-B       12984      5346668 bp    0.86 %
     L1/CIN4            72        11223 bp    0.00 %
   LTR elements:     22058     18867040 bp    3.02 %
     BEL/Pao          1959      2281547 bp    0.37 %
     Ty1/Copia         225       195687 bp    0.03 %
     Gypsy/DIRS1     15108     13152263 bp    2.11 %
       Retroviral        0            0 bp    0.00 %

DNA transposons     365920    118872436 bp   19.04 %
   hobo-Activator     9890      3711019 bp    0.59 %
   Tc1-IS630-Pogo    37403     10920771 bp    1.75 %
   En-Spm                0            0 bp    0.00 %
   MuDR-IS905            0            0 bp    0.00 %
   PiggyBac           5280      1318534 bp    0.21 %
   Tourist/Harbinger 15668      3870853 bp    0.62 %
   Other (Mirage,      131        36952 bp    0.01 %
    P-element, Transib)

Rolling-circles     129838     82442325 bp   13.20 %

Unclassified:       234562     64598689 bp   10.35 %

Total interspersed repeats:   219774316 bp   35.20 %


Small RNA:             741       152705 bp    0.02 %

Satellites:           1654       356882 bp    0.06 %
Simple repeats:     134493      6637749 bp    1.06 %
Low complexity:      16771       784708 bp    0.13 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element
                                                      

RepeatMasker version 4.1.1 , default mode
                                        
run with rmblastn version 2.10.0+
The query was compared to classified sequences in "merged.fa"
```
Alignment assembly was first completed using transcriptome data by PASA pipeline v2.5.1 (Haas et al., 2003). The genome was also annotated by MAKER v3.01.03 (Cantarel et al., 2008).
```bash

```

**Reference**

Cantarel, B. L., Korf, I., Robb, S. M., Parra, G., Ross, E., Moore, B., ... Yandell, M. (2008). MAKER: an easy-to-use annotation pipeline designed for emerging model organism genomes. Genome research, 18, 188-196. https://doi.org/10.1101/gr.6743907

Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature Methods, 18, 170-175. https://doi.org/10.1038/s41592-020-01056-5

Durand, N. C., Robinson, J. T., Shamim, M. S., Machol, I., Mesirov, J. P., Lander, E. S., & Aiden, E. L. (2016). Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom. Cell systems, 3, 99-101. https://doi.org/10.1016/j.cels.2015.07.012

Haas, B. J., Delcher, A. L., Mount, S. M., Wortman, J. R., Smith Jr, R. K., Hannick, L. I., ... Town, C. D. (2003). Improving the Arabidopsis genome annotation using maximal transcript alignment assemblies. Nucleic Acids Research, 31, 5654-5666. https://doi.org/10.1093/nar/gkg770

Huang, S., Kang, M., & Xu, A. (2017). HaploMerger2: rebuilding both haploid sub-assemblies from high-heterozygosity diploid genome assembly. Bioinformatics, 33, 2577-2579. https://doi.org/10.1093/bioinformatics/btx220

Marçais, G., & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27, 764-770. https://doi.org/10.1093/bioinformatics/btr011

Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31, 3210-3212. https://doi.org/10.1093/bioinformatics/btv351

Smit AFA, H. R. RepeatMasker. http://www.repeatmasker.org/RepeatModeler/. Accessed 23 June 2022

Vurture, G. W., Sedlazeck, F. J., Nattestad, M., Underwood, C. J., Fang, H., Gurtowski, J., & Schatz, M. C. (2017). GenomeScope: fast reference-free genome profiling from short reads. Bioinformatics, 33, 2202-2204. https://doi.org/10.1093/bioinformatics/btx153

Waterhouse, R. M., Seppey, M., Simão, F. A., Manni, M., Ioannidis, P., Klioutchnikov, G., ... Zdobnov, E. M. (2018). BUSCO applications from quality assessments to gene prediction and phylogenomics. Molecular Biology and Evolution, 35, 543-548. https://doi.org/10.1093/molbev/msx319
