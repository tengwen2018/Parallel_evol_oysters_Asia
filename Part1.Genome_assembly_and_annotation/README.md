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

**Primary genome was assembled by hifiasm v0.16 with default parameters (Cheng et al., 2021).**
```bash
hifiasm -o output.asm -t32 input1.hifi.fastq.gz input2.hifi.fastq.gz
```
**Haploid sub-assembly was reconstructed by HaploMerger2 (Huang et al., 2017).**
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
**The chromosome‐level assembly** was performed by Juicer v1.6 with default parameters (Durand et al., 2016).
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
**Genome assembly completeness was accessed by BUSCO v5.1.3 (Simão et al., 2015) using the "mollusca_odb10" BUSCO gene set collection (Waterhouse et al., 2018)(https://busco-data.ezlab.org/v5/data/lineages/mollusca_odb10.2020-08-05.tar.gz).**

```bash
busco --offline -f -c 10 -i ref.fasta -l mollusca_odb10 -o output_busco -m genome
```
Here is what the outpu looked like:
```
# BUSCO version is: 5.1.2 
# The lineage dataset is: mollusca_odb10 (Creation date: 2020-08-05, number of genomes: 7, number of BUSCOs: 5295)
# Summarized benchmarking in BUSCO notation for file crassostrea_angulata.fna
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
**Transposable elements de novo prediction and homology searching were performed by RepeatModeler v2.0.3 and RepeatMasker v4.1.1 (https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.1.tar.gz) pipeline (Smit AFA).**
```bash
BuildDatabase -name ref ref.fasta
RepeatModeler -database ref -pa 20 -LTRStruct
RepeatMasker -e ncbi -lib ref-families.fa -gff -dir 00-RepeatMask_RepeatModeler -pa 15 -a ref.fasta
```
Here is what the outpu looked like:
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
**Alignment assembly was first completed using transcriptome data by PASA pipeline v2.5.1 (Haas et al., 2003). The genome was also annotated by MAKER v3.01.03 (Cantarel et al., 2008).**
```bash
# Data trimming of RNA HiFi reads
lima pt1mix10.ccs.bam primers.fasta pt1mix10.flnc.bam --isoseq --peek-guess
# Convert bam file to fastq file
bam2fastq -o pt1mix10.flnc pt1mix10.flnc.bam
# Mapping reads to our genome reference
minimap2 -ax splice -t 40 ref.fasta pt1mix10.flnc.fastq.gz | samtools view -h -F 4 | samtools sort -@ 10  > pt1mix10.flnc.mapped.bam
# Transcriptome assembly
Trinity --genome_guided_bam pt1mix10.flnc.mapped.bam --max_memory 1000G --genome_guided_max_intron 10000 --CPU 50
# Annotation with PASA
~/tools/PASApipeline/Launch_PASA_pipeline.pl -c ~/tools/PASApipeline/pasa_conf/pt1mix10_withdatatrimming.config -C -R -g ref.fasta -t Trinity-GG.rename.fasta.clean -T -u Trinity-GG.rename.fasta --ALIGNERS blat,gmap --CPU 50
```
**Genome Annotation with MAKER**

Step1. Initial MAKER Analysis
```
~/tools/RepeatMasker/util/rmOutToGFF3.pl ref.unmasked.fa.out > ref.unmasked.fa.gff3
grep -v -e "Satellite" -e ")n" -e "-rich" ref.unmasked.fa.gff3 > ref.unmasked.fa.complex.gff3
cat ref.unmasked.fa.complex.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > ref.unmasked.fa.complex.reformat.gff3

# rnd1_maker_opts.ctl
#-----Genome (these are always required)
genome==ref.fasta
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=Trinity-GG.rename.fasta
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=GCF_902806645.1_cgigas_uk_roslin_v1_protein.faa,GCF_002022765.2_C_virginica-3.0_protein.faa,S_glomerata_protein.faa,P_fucata_protein.faa,GCF_002113885.1_M_yessoensis_protein.faa
in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=simple #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=ref.unmasked.fasta.complex.reformat.gff3 #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
```
# maker_bopts.ctl
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'
use_rapsearch=0 #use rapsearch instead of blastx, 1 = yes, 0 = no

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
```
```
# maker_exe.ctl
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=~/tools/rmblast-2.11.0/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=~/tools/rmblast-2.11.0/bin/blastn #location of NCBI+ blastn executable
blastx=~/tools/rmblast-2.11.0/bin/blastx #location of NCBI+ blastx executable
tblastx=~/tools/rmblast-2.11.0/bin/tblastx #location of NCBI+ tblastx executable
formatdb=~/tools/bin/formatdb #location of NCBI formatdb executable
blastall=~/tools/bin/blastall #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=~/tools/RepeatMasker/RepeatMasker #location of RepeatMasker executable
exonerate=~/tools/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=~/tools/snap/snap #location of snap executable
gmhmme3=~/tools/gmes_linux_64/gmhmme3 #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=~/tools/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
evm= #location of EvidenceModeler executable
tRNAscan-SE=~/tools/bin/tRNAscan-SE #location of trnascan executable
snoscan=~/tools/bin/snoscan #location of snoscan executable

#-----Other Algorithms
probuild=~/tools/gmes_linux_64/probuild #location of probuild executable (required for genemark)
```
```
mpiexec -n 26 maker -base rnd1 rnd1_maker_opts.ctl maker_bopts.ctl maker_exe.ctl | tee round1_run_maker.log
cd rnd1.maker.output
gff3_merge -s -d *_master_datastore_index.log > rnd1.all.maker.gff
fasta_merge -d *_master_datastore_index.log
gff3_merge -n -s -d *_master_datastore_index.log > rnd1.all.maker.noseq.gff
# transcript alignments
awk '{ if ($2 ~ /est2genome/) print $0 }' rnd1.all.maker.noseq.gff > rnd1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 ~ /protein2genome/) print $0 }' rnd1.all.maker.noseq.gff > rnd1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' rnd1.all.maker.noseq.gff > rnd1.all.maker.repeats.gff
```
Step2. Training Gene Prediction Software
```
# SNAP
mkdir -p snap/round1
cd snap/round1
# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../../rnd1.maker.output/rnd1_master_datastore_index.log
rename 's/genome/rnd1.zff.length50_aed0.25/g' *
# gather some stats and validate
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# assembly the HMM
hmm-assembler.pl rnd1.zff.length50_aed0.25 params > rnd1.zff.length50_aed0.25.hmm

# AUGUSTUS
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ../../rnd1.maker.output/rnd1.all.maker.noseq.gff | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools getfasta -fi ../../ref.fasta -bed - -fo rnd1.all.maker.transcripts1000.fasta
~/tools/busco/scripts/run_BUSCO.py -i rnd1.all.maker.transcripts1000.fasta -f -o rnd1_maker -l ./busco_downloads/lineages/mollusca_odb10 -m genome -c 48 --long -sp human -z --augustus_parameters='--progress=true'

cd snap/round1/run_rnd1_maker/augustus_output/retraining_parameters
rename 's/Boa_constrictor/snap_round1/g' *
sed -i 's/Boa_constrictor/snap_round1/g' snap_round1_parameters.cfg
sed -i 's/Boa_constrictor/snap_round1/g' snap_round1_parameters.cfg.orig1
cp -f * ~/tools/augustus/config/species/snap_round1/
```
Step3. MAKER With Ab Initio Gene Predictors
```
# rnd2_maker_opts.ctl
#-----Genome (these are always required)
genome=ref.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=rnd1.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=rnd1.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=rnd1.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=snap/round1/rnd1.zff.length50_aed0.25.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=snap_round1 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=20000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
```
mpiexec -n 26 maker -base rnd2 rnd2_maker_opts.ctl maker_bopts.ctl maker_exe.ctl

cd rnd2.maker.output
gff3_merge -s -d *_master_datastore_index.log > rnd2.all.maker.gff
fasta_merge -d *_master_datastore_index.log
# GFF w/o the sequences
gff3_merge -n -s -d *_master_datastore_index.log > rnd2.all.maker.noseq.gff
# transcript alignments
awk '{ if ($2 ~ /est2genome/) print $0 }' rnd2.all.maker.noseq.gff > rnd2.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 ~ /protein2genome/) print $0 }' rnd2.all.maker.noseq.gff > rnd2.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' rnd2.all.maker.noseq.gff > rnd2.all.maker.repeats.gff
```
Repat step2 and step3 three times.


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
