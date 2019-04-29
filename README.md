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
```shell
bcftools-1.9 mpileup -Ou -f hcontortus_chr_mtDNA_arrow_pilon.fa ZAI_ZAI_OA_014.sorted.bam  | bcftools-1.9 call -v -c --ploidy 1 -Ob --skip-variants indels > ZAI_ZAI_OA_014.bcf

bcftools index ZAI_ZAI_OA_014.bcf


bcftools-1.9 view ZAI_ZAI_OA_014.bcf -Oz > ZAI_ZAI_OA_014.vcf.gz
tabix -p vcf ZAI_ZAI_OA_014.vcf.gz
#--- this can be viewed in artemis
```


## Mapping - multiple samples
``` shell

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
```

## SNP calling - multiple samples
```shell
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

```shell
#filter SNPs by MAF
vcftools-0.1.14 --gzvcf all_samples.vcf.gz --maf 0.05 --min-alleles 2 --max-alleles 2 --recode --out all_samples.filtered
```

## Analysis
``` R
R-3.5.0
library(adegenet)
library(vcfR)
library(ggplot2)
library(patchwork)
library(poppr)
library(ape)
library(RColorBrewer)

library(dplyr)
library(ggrepel)
library(ggtree)
library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")
library(reshape2)

# prepare your data
#--- specify your input files
vcf_file <- "all_samples.filtered.recode.vcf"
dna_file <- "hcontortus_chr_mtDNA_arrow_pilon.fa"
metadata<-read.table("../sample_metadata.txt",header=T)

#--- read the input files into R
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")

#--- convert your vcf file into a dataframe that R can interpret
vcf.gl <- vcfR2genlight(vcf)
pop(vcf.gl) <- metadata$country
ploidy(vcf.gl) <- 1


#--- lets setup some colours to use throughout
cols<-colorRampPalette(brewer.pal(8, "Set1"))(18)


# PCA analysis
#--- perform a PCA
vcf.pca <- glPca(vcf.gl, nf = 10)
vcf.pca.scores <- as.data.frame(vcf.pca$scores)
vcf.pca.scores$country <- metadata$country

#--- Lets calculate the amount of variance each principal component describes. We will use this in the plots
var_frac <- vcf.pca$eig/sum(vcf.pca$eig)*100
PC1_variance <- formatC(head(vcf.pca$eig)[1]/sum(vcf.pca$eig)*100)
PC2_variance <- formatC(head(vcf.pca$eig)[2]/sum(vcf.pca$eig)*100)
PC3_variance <- formatC(head(vcf.pca$eig)[3]/sum(vcf.pca$eig)*100)
PC4_variance <- formatC(head(vcf.pca$eig)[4]/sum(vcf.pca$eig)*100)

#--- plot eigenvectors
barplot(100*vcf.pca$eig/sum(vcf.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)




# make a basic plot of the first two principal components
pca_plot_pc12 <- ggplot(vcf.pca.scores,aes(PC1,PC2,col=country))+geom_point()
pca_plot_pc12
# lets add some labels
pca_plot_pc12 <- pca_plot_pc12 + labs(x=paste0("PC1 variance = ",PC1_variance,"%"),y=paste0("PC2 variance = ",PC2_variance,"%"))
pca_plot_pc12
# lets add some colours per country
pca_plot_pc12 <- pca_plot_pc12 + geom_point(aes(col=country)) + stat_ellipse(level = 0.95, size = 1)
pca_plot_pc12

# plot principal components 3 and 4, and we will compare them to 1 and 2
pca_plot_pc34 <- ggplot(vcf.pca.scores,aes(PC3,PC4,col=country))+geom_point()+ stat_ellipse(level = 0.95, size = 1)+ labs(x=paste0("PC3 variance = ",PC3_variance,"%"),y=paste0("PC4 variance = ",PC4_variance,"%"))
pca_plot_pc34

# calculate the mean value of the principal components for each country. We can use this to make some labels for our plots
means = vcf.pca.scores %>% group_by(country) %>% summarize(meanPC1 = mean(PC1),meanPC2 = mean(PC2),meanPC3 = mean(PC3),meanPC4 = mean(PC4))

# label plot with mean values labels for each country
pca_plot_pc12 <- pca_plot_pc12 + geom_label_repel(data=means,aes(means$meanPC1,means$meanPC2,col=means$country,label=means$country))
pca_plot_pc12

# lets look at the large cluster in more detail
pca_plot_pc12_zoom <- pca_plot_pc12 + xlim(5,10) + ylim(-1,0.5) + geom_label_repel(data=means,aes(means$meanPC1,means$meanPC2,col=means$country,label=means$country))		
pca_plot_pc12_zoom

#Questions
#--- 1.



# generate a tree based on genetic distances between samples

#--- generate a pairwise distance matrix from the mtDNA vcf, which we will use to make a tree
tree_data <- aboot(vcf.gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

#--- plot the tree
#plot.phylo(tree_data, cex = 0.4, font = 2, adj = 0, tip.color =  cols[pop(vcf.gl)])
#nodelabels(tree_data$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.4,font = 3, xpd = TRUE)

ggtree(tree_data)+ geom_tiplab(size=2,color=cols[pop(vcf.gl)])+xlim(-0.1, 0.3)+geom_nodelab(size=2,nudge_x=-0.006,nudge_y=1)+theme_tree2(legend.position='centre')
#--- add the node labels


#--- add a legend





#--- quantitative differentiation
#--- allele frequencies
myDiff_pops <- genetic_diff(vcf,pops=vcf.gl@pop)
AF_data <- myDiff_pops[,c(1:20)]
AF_data <- melt(AF_data)
colnames(AF_data) <- c("CHROM","POS","country","value")
AF_data$country <- gsub("Hs_","",AF_data$country)

#--- extract the latitude and longitude for each country from the metadata file
coords <- data.frame(metadata$country,metadata$latitude,metadata$longitude)
coords <- unique(coords)
colnames(coords) <- c("country","latitude","longitude")

#--- join the allele frequency data and the latitude/longitude data together
AF_data_coords <- dplyr::left_join(AF_data, coords, by = "country")


#--- calculate Gst and Gst'
myDiff_pairwise <- pairwise_genetic_diff(vcf,pops=vcf.gl@pop)
colMeans(myDiff[,c(4:ncol(myDiff))], na.rm = TRUE)
as.matrix(colMeans(myDiff[,c(4:ncol(myDiff))], na.rm = TRUE))


# make a map of the sampling locations, and plot some data on it

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(dplyr)
library(ggrepel)
library(plotrix)

#mapWorld <- borders("world", colour="gray85", fill="gray85")
#ggplot() + mapWorld + geom_point(aes(metadata$longitude, metadata$latitude,col=metadata$country))+scale_colour_manual(values=cols)+theme_bw()+labs(x="Longitude",y="Latitude")

# select a SNP of interest based on its position
AF_SNP_coords <- AF_data_coords[AF_data_coords$POS=="10076",]
AF_SNP_coords <- AF_data_coords[1175,]

par(fg = "black")
map("world",col="grey85",fill=TRUE, border=FALSE)
map.axes()
points(metadata$longitude, metadata$latitude, cex=1,pch=20,col=cols[pop(vcf.gl)])
for (i in 1:nrow(AF_SNP_coords)){
   add.pie(z=c(AF_SNP_coords$value[i],1-AF_SNP_coords$value[i]),x=AF_SNP_coords$longitude[i]+10,y=AF_SNP_coords$latitude[i],radius=5,col=c(alpha("orange", 0.5), alpha("blue", 0.5)), labels="")
}
legend( x="left", legend=AF_SNP_coords$country,col=cols[as.factor(AF_SNP_coords$country)], lwd="1", lty=0, pch=20,box.lwd = 0,cex = 0.9)
```
