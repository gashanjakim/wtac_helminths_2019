# Helminths 2019 - Genetic Diversity module

## Software needed
samtools-1.6
bcftools-1.9
bwa 0.7.17-r1188



## Working dir
```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity
```

## Get and prepare some raw data
Extracting mtDNA reads already mapped to the HCON_V4 genome from my populaiton genetics analysis. 
``` shell
# extact reads on mtDNA 
#--- new bams: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/MAPPING/MERGED_BAMS
#--- GS bams: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/MAPPING

for i in $( cd MERGED_BAMS ; ls -1 *.merged.bam) ; do \
	samtools view -f3 -b MERGED_BAMS/${i} hcontortus_chr_mtDNA_arrow_pilon | \
    samtools sort -n -o ${i}.sorted - ; \
    samtools fastq -1 ${i%%.merged.bam}_1.fq.gz -2 ${i%%.merged.bam}_2.fq.gz ${i}.sorted ; \
    rm ${i}.sorted; \
    done

for i in $( cd MAPPING ; ls -1 *.merged.bam) ; do \
	samtools view -f3 -b MAPPING/${i} hcontortus_chr_mtDNA_arrow_pilon | \
    samtools sort -n -o ${i}.sorted - ; \
    samtools fastq -1 ${i%%.merged.bam}_1.fq.gz -2 ${i%%.merged.bam}_2.fq.gz ${i}.sorted ; \
    rm ${i}.sorted; \
    done

# found that downsampling to 5000 reads is sufficient for mapping and SNP calling, while at the same time minimising the file size footprint. Using seqtk to downsample

for i in *.gz; do  seqtk sample -s100  ${I} 5000 > ${i%%.fq.gz}.fastq; gzip ${i%%.fq.gz}.fastq; done &
``` 
    
    
## mapping to the mtDNA genome
Will perfomr this in two parts for the module, first with a single sample, and second looping over all samples


## mapping - single sample
``` shell
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity/MAPPING

#--- get reference - need to provide this
cp ../../../../../hc/GENOME/REF/hcontortus_chr_mtDNA_arrow_pilon.fa .

#--- get reads
ln -s ../RAW_READS/ZAI_ZAI_OA_014_1.fastq.gz
ln -s ../RAW_READS/ZAI_ZAI_OA_014_2.fastq.gz


#--- mapping 

bwa index hcontortus_chr_mtDNA_arrow_pilon.fa

bwa mem hcontortus_chr_mtDNA_arrow_pilon.fa ZAI_ZAI_OA_014_1.fastq.gz ZAI_ZAI_OA_014_2.fastq.gz > ZAI_ZAI_OA_014.tmp.sam

samtools view -q 15 -b -o ZAI_ZAI_OA_014.tmp.bam ZAI_ZAI_OA_014.tmp.sam

samtools sort ZAI_ZAI_OA_014.tmp.bam -o ZAI_ZAI_OA_014.sorted.bam

samtools index ZAI_ZAI_OA_014.sorted.bam 

#--- this can be viewed in artemis
```


## SNP calling - single sample
```
bcftools-1.9 mpileup -Ou -f hcontortus_chr_mtDNA_arrow_pilon.fa ZAI_ZAI_OA_014.sorted.bam  | bcftools-1.9 call -v -c --ploidy 1 -Ob --skip-variants indels > ZAI_ZAI_OA_014.bcf

bcftools index ZAI_ZAI_OA_014.bcf


bcftools-1.9 view ZAI_ZAI_OA_014.bcf -Oz > ZAI_ZAI_OA_014.vcf.gz
tabix -p vcf ZAI_ZAI_OA_014.vcf.gz
#--- this can be viewed in artemis
```


## Mapping - multiple samples
``` 

for i in $( cd ../RAW_READS ; ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//g'); do
# map reads
bwa mem hcontortus_chr_mtDNA_arrow_pilon.fa ../RAW_READS/${i}_1.fastq.gz ../RAW_READS/${i}_2.fastq.gz > ${i}.tmp.sam ;
# convert sam to bam
samtools view -q 15 -b -o ${i}.tmp.bam ${i}.tmp.sam ;
#sort reads
samtools sort ${i}.tmp.bam -o ${i}.sorted.bam ; 
# index the bam
samtools index ${i}.sorted.bam ;
# remove files you dont need
rm *tmp*;
done


## SNP calling - multiple samples
```
# make a list of bam files that you want to include in the analysis
ls -1 *.sorted.bam > bam.list
# call SNPs from the sampels listed in the bam list, to generate a multi-sample bcf
bcftools-1.9 mpileup -Ou --annotate FORMAT/DP -f hcontortus_chr_mtDNA_arrow_pilon.fa --bam-list bam.list | bcftools-1.9 call -v -c --ploidy 1 -Ob --skip-variants indels  > all_samples.bcf
# index the multi-sample bcf
bcftools index all_samples.bcf
# convert the bcf to a vcf
bcftools-1.9 view all_samples.bcf -Oz > all_samples.vcf.gz
# index the vcf
tabix -p vcf all_samples.vcf.gz

```


vcftools-0.1.14 --gzvcf all_samples.vcf.gz --maf 0.05


## Analysis
``` R
R-3.5.0
library(adegenet)
library(vcfR)

vcf_file <- "all.vcf.gz"
dna_file <- "hcontortus_chr_mtDNA_arrow_pilon.fa"

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")

chrom <- create.chromR(name="hcontortus_mtDNA", vcf=vcf, seq=dna, verbose=TRUE)



```
