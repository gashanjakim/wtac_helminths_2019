# hcontortus_qtl

This primary focus of this project on QTL mapping ivermectin response in Haemonchus contortus.

The origin of the material analysed is the MHco3/MHco18 genetic cross performed by the BUG Consortium. The analysis of these data have been described previously in different ways, including
- genetic map (Doyle et al. 2018 GBE)
- XQTL genomics (Doyle, Laing et al. *in prep*)
- XQTL transcriptomics (Laing et al. *in prep*)

The main QTL experiment here is a dose response larval development experiment, in which two groups of larvae were collected:
- L3 that developed normally at high doses of ivermectin
- L1/L2 that were not developing well on low doses of ivermectin
We originally performed a poolseq experiment, where pools of larvae (n=200 per pool) were sequenced and compared, revealing a consistent chromosome 5 signature of selection as the XQTL experiments. The rationale for this new QTL experiment was that genotypes from individual worms, rather than allele frequencies from pooled worms, may provide additional support for variants in the region under selection.

The added benefit of sequencing individual larvae is that there a lot more than can be done with the data, including (and this forms my to-do list)
- determine LD and explore recombination
- perform phasing
- sex-specific analyses, which can be determined from X chromosome coverage and heterozygosity
- check for variation in ploidy
- genetic relatedness
- more accurate CNV analyses
- perform GWAS, taking into account genetic relatedness.

Further, this is going to be an excellent resource for some benchmarking, including
- SNP calling
- genome graphs


We also prepared sequencing libraries from larvae obtained from US farms that have been phenotypically tested using drenchrite assays by Ray Kaplan. These are the same farms for which we have sequenced pools and analysed in the XQTL paper. Library prep revealed these didn't work that well - the libraries were quite weak especially compared with the QTL samples - and so we have taken a different pooling approached to try and maximise samples. There were more samples prepared than sequenced, however, becasue of the library strength we have put these on hold for now. For those samples that were sequecned, these probably wont work that well.


## Sequencing data  
### Sequencing data - QTL
- susceptible (267 samples)
     - 36342_3 - 94 samples (SUS)
     - 36808_1 - 79 samples (SUS) - QC PENDING
     - 35990_1 - 94 samples (SUS)
- resistance (235 samples)
     - 36342_1 - 75 samples (RES)
     - 36342_2 - 82 samples (RES)
     - 36342_4 - 78 samples (RES)


### Sequencing data - Farm
- 35887_2 - 96 samples (FARM)




```bash
# get reference
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/01_REFERENCE

cp ../../GENOME/REF/HAEM_V4_final.chr.fa .

# make a index and a dict file
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```
[↥ **Back to top**](#top)



```bash
# get the raw data
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW

for i in 36342_3 36808_1 35990_1 36342_1 36342_2 36342_4 35887_2; do
     pf data --type lane --id ${i} --symlink ./ --rename --filetype fastq;
     done

# extract lane/barcode and sample IDs to create some metadata
>sample_lanes.list
for i in 36342_3 36808_1 35990_1 36342_1 36342_2 36342_4 35887_2; do
     pf supplementary --type lane --id ${i} | grep -v "Sample" | sed 's/\#/_/g' | awk '{print $3,$6}' OFS="\t" >> sample_lanes.list;
     done

```
[↥ **Back to top**](#top)



```bash
# mapping
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/03_MAPPING

ln -s ../02_RAW/sample_lanes.list

screen
while read LANE NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter \
     ${NAME} \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW/${LANE}_1.fastq.gz \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW/${LANE}_2.fastq.gz;
     done < sample_lanes.list &


```
