# Helminths 2019 - genetic Diverity module

```shell
# Working dir
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity
```

``` shell
# GET AND PREPARE RAW DATA

# extact reads on mtDNA 
#--- new bams: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/MAPPING/MERGED_BAMS
#--- GS bams: /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/POPULATION_DIVERSITY/GS_ORIGINAL/MAPPING
```


```
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

for i in *.gz; do  seqtk sample -s100  ${I} 5000 > ${I%%.fq.gz}.fastq; gzip ${I%%.fq.gz}.fastq; done &
``` 
    
    




## MAPPING
``` shell
cd /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity/MAPPING

# get reference
cp ../../../../../hc/GENOME/REF/hcontortus_chr_mtDNA_arrow_pilon.fa .


ln -s /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity/RAW_READS/CH_SWI_003_1.fq.gz
ln -s /nfs/users/nfs_s/sd21/lustre118_link/WTAC/HELMINTHS_2019/data/Module_GeneticDiversity/RAW_READS/CH_SWI_003_2.fq.gz

#--- subsample 10000 reads for comparison

seqtk sample -s100  CH_SWI_003_1.fq.gz 5000 > sub1.fq
seqtk sample -s100  CH_SWI_003_2.fq.gz 1000 > sub2.fq


#--- mapping 

bwa index hcontortus_chr_mtDNA_arrow_pilon.fa

bwa mem hcontortus_chr_mtDNA_arrow_pilon.fa CH_SWI_003_1.fq.gz CH_SWI_003_2.fq.gz > all.sam
bwa mem hcontortus_chr_mtDNA_arrow_pilon.fa sub1.fq sub2.fq > sub.sam


samtools view -q 15 -b -o all.bam all.sam
samtools view -q 15 -b -o sub.bam sub.sam


samtools sort –o all2.bam all.bam
samtools sort –o sub2.bam sub.bam

samtools index all2.bam
samtools index sub2.bam
#--- this can be viewed in artemis
```







## SNP calling
```
bcftools-1.9 mpileup -Ou -f hcontortus_chr_mtDNA_arrow_pilon.fa all2.bam | bcftools-1.9 call -v -c --ploidy 1 -Ob --skip-variants indels > all.bcf

bcftools index all.bcf


bcftools-1.9 view all.bcf -Oz > all.vcf.gz
tabix -p vcf all.vcf.gz
#--- this can be viewed in artemis
```






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
