**SNP calling and filtering with VCFtools v0.1.16 (Danecek et al., 2011).**
```bash
bcftools mpileup --threads 10 -a AD,DP,SP -Ou -f ref.fa *.dup.bam | bcftools call --threads 10 -f GQ,GP -mO z -o oyster.vcf.gz
VCF_IN=oyster.vcf.gz
VCF_OUT=oyster.filtered.vcf.gz
MAF=0.1
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT
```

**Principal components analysis (PCA) with PLINK v1.9(Purcell et al., 2007).**
```bash
# perform linkage pruning - i.e. identify prune sites
plink --vcf oyster.filtered.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out oyster
# prune and create pca
plink --vcf oyster.filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract oyster.prune.in \
--threads 10 \
--make-bed --pca --out oyster
```
**Individual ancestries estimation with ADMIXTURE (Alexander et al., 2009).**
```bash
FILE=oyster.filtered
# Generate the input file in plink format
plink --vcf $FILE.vcf.gz --make-bed --out $FILE --allow-extra-chr
# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim
for i in {1..5}
do
admixture --cv $FILE.bed $i > log${i}.out
sh admixture_cangcgig.$i.sh
done
```
**Estimation of per-site *F*<sub>ST</sub> and sliding-window *F*<sub>ST</sub>.**
```bash
# filter individuals with missingness above 20%
vcftools --gzvcf oyster.filtered.vcf.gz --missing-indv
sed '1d' out.imiss | awk '{if($5>0.2)print $1}' > lowDP.indv
vcftools --gzvcf oyster.filtered.vcf.gz --remove lowDP.indv --recode --stdout | gzip -c > filterlowDP.vcf.gz
# per site fst
time vcftools --gzvcf filterlowDP.vcf.gz \
--weir-fst-pop south \
--weir-fst-pop north \
--out ./south.vs.north
# sliding windows of 2kb
parseVCF.py -i filterlowDP.vcf.gz | bgzip > filterlowDP.geno.gz
popgenWindows.py -g filterlowDP.geno.gz -o filterlowDP.2k.Fst.Dxy.pi.csv.gz \
   -f phased -w 2000 -m 10 -s 2000 \
   -p south -p north \
   --popsFile pop.info
```


**Reference**

Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome research, 19, 1655-1664. https://doi.org/10.1101/gr.094052.109

Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., ... Sherry, S. T. (2011). The variant call format and VCFtools. Bioinformatics, 27, 2156-2158. https://doi.org/10.1093/bioinformatics/btr330

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., ... Sham, P.C. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81, 559-575. https://doi.org/10.1086/519795
