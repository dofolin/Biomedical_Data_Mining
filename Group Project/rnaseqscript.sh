#112-1 Summer Biomedical Data Mining_2023_08_08_RNA-seq
## conda environment
[STEP1] create personal environment
conda create --name rnas
[STEP2] activate your own environment
conda activate rnas

## get the sample file

[STEP1] cp course file 
cp -r /home/r11631063/rnaseq .
[STEP2] make a folder to save sample file
cd rnaseq
mkdir data
[STEP3] unzip tar file(hcc1395 fastq data)
tar -xvf practical.tar -C data/
[STEP4] remove tar file (optional)
rm practical.tar


## install FastQC
[STEP1] install FastQC
conda install -c bioconda fastqc
[STEP2] create FastQC output folder
mkdir fastqc_raw
[STEP3] analyze with FastQC
fastqc -t 2 data/*.fastq.gz -o fastqc_raw/

## install Trimmomatic
[STEP1] install Trimmomatic
conda install -c bioconda trimmomatic
[STEP2] make sure there is shell script used and addapter fasta
 trimmomatic.sh 
 TruSeq3-PE.fa 

#!/bin/bash
conda activate rnas

cd ./data/
gunzip *.fastq.gz
for f1 in *r1.fastq
do
        f2=${f1%%r1.fastq}"r2.fastq"
        trimmomatic PE -threads 2 $f1 $f2 ${f1%%.fastq}"_clean_paired.fastq" ${f1%%.fastq}"_unpaired.fastq" ${f2%%.fastq}"_clean_paired.fastq" ${f2%%.fastq}"_unpaired.fastq" ILLUMINACLIP:../TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done
mkdir clean_fastq
mv *_clean_paired.fastq clean_fastq

[STEP3]  run trimmomatic script
bash trimmomatic.sh

## rerun fastQC
fastqc -t 2 ./data/clean_fastq/*_clean_paired.fastq -o ./data/clean_fastq/

## install HISAT2
[STEP1] check  folder
make sure there is a folder call hisat
[STEP2] get into hisat2
cd hisat
[STEP3] install
conda install -c bioconda hisat2
[STEP4] reference genome:
        GRCh38_latest_genomic.fna
#[STEP5] and [STEP6] takes a long time we've done it for you        
[STEP5] build library   (or download from https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz)
        hisat2-build.sh 

#!/bin/sh
conda activate rnas

hisat2-build -p 2 GRCh38_latest_genomic.fna genome

[STEP6] run hisat2-build 
bash hisat2-build.sh
[STEP7] check file existence
    make sure there is genome.(1~8).ht2 in hisat file
[STEP8]get back to rnaseq
    cd ..
[STEP9] prepare to run hisat2 
       hisat2.sh 

#!/bin/sh

conda activate rnas

mkdir hisat2_log
for f1 in ./data/clean_fastq/*r1_clean_paired.fastq
do
        f2=${f1%%r1_clean_paired.fastq}"r2_clean_paired.fastq"
        hisat2 -p 2 --dta -x hisat/genome -1 $f1 -2 $f2 -S ${f1:19:-21}align.sam 2>hisat2_log/${f1:20:-21}align.log
done

[STEP11] run hisat2 
    bash hisat2.sh
[STEP11]
    check if there is hcc1395_normal_rep1_align.sam 
## run StringTie
[STEP1] install StringTie and samtools
conda install -c bioconda stringtie
conda install -c bioconda samtools
[STEP2] convert sam file to bam
    samtobam.sh 

#!/bin/sh

conda activate rnas

for i in *.sam
do
        i=${i%_align.sam*}
        samtools sort -@ 2 -o ${i}.bam ${i}_align.sam
done
[STEP3] run samtobam.sh
    bash samtobam.sh
    (if there is no bam file ,try conda deactivate then bash again)

[STEP4] check gtf file
    make sure there is GRCh38_latest_genomic.gff

[STEP5] StringTie
    stringtie.sh 

#!/bin/sh

conda activate rnas

for i in *.bam
do
        i=${i%.bam*}
        stringtie -p 2 -G GRCh38_latest_genomic.gff -o ${i}.gff ${i}.bam
done
[STEP6] run stringtie
    bash stringtie.sh


[STEP7] gff file list
    vim gff_list.txt 
    copy paste below in to gff_list.txt

hcc1395_normal_rep1.gff
hcc1395_normal_rep2.gff
hcc1395_normal_rep3.gff
hcc1395_tumor_rep1.gff
hcc1395_tumor_rep2.gff
hcc1395_tumor_rep3.gff

[STEP8] merge all gff files
stringtie --merge -p 2 -G GRCh38_latest_genomic.gff -o stringtie_merged.gff  gff_list.txt

[STEP10] transcript assembly and quantification with StringTie
conda install -c bioconda gffcompare
gffcompare -r GRCh38_latest_genomic.gff -G -o merged stringtie_merged.gff

[STEP11] restringtie
    restringtie.sh 

#!/bin/sh

conda activate rnas
mkdir final_gff
for i in *.bam 
do 
        i=${i%.bam*}
        stringtie -e -B -p 2 -G stringtie_merged.gff -o final_gff/${i}/${i}.gff  ${i}.bam 
done
[STEP12] run restringtie
    bash restringtie.sh

## Differential expression analyze with DESeq2 in R
[STEP1] get into final_gff folder
    cd final_gff
[STEP2] prepare file list
    vim final_gff_list.txt 
    copy paste below in to gff_list.txt

normal_rep1 ./hcc1395_normal_rep1/hcc1395_normal_rep1.gff
normal_rep2 ./hcc1395_normal_rep2/hcc1395_normal_rep2.gff
normal_rep3 ./hcc1395_normal_rep3/hcc1395_normal_rep3.gff
tumor_rep1 ./hcc1395_tumor_rep1/hcc1395_tumor_rep1.gff
tumor_rep2 ./hcc1395_tumor_rep2/hcc1395_tumor_rep2.gff
tumor_rep3 ./hcc1395_tumor_rep3/hcc1395_tumor_rep3.gff

[STEP3] move and check prepDE.py
    mv ../prepDE.py .

[STEP4] prepare for DESeq2
python prepDE.py -i final_gff_list.txt -g gene_results.csv -t transcript_results.csv
[STEP5] install R & DESeq2
conda deactivate
conda create --name deseq2
conda activate deseq2
conda install -c conda-forge r-base
conda install -c bioconda bioconductor-deseq2
[STEP6] move R script to final_gff
    mv ../deseq2.R .

library(DESeq2)
database <- read.table(file = "transcript_results.csv", sep = ",", header = TRUE, row.names = 1)
database <- round(as.matrix(database))

condition <- factor(c("normal", "normal", "normal", "tumor", "tumor", "tumor"),
  levels = c("normal", "tumor"))
coldata <- data.frame(row.names = colnames(database), condition)
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)
res <- results(dds)

res <- res[order(res$padj),]
diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene <- row.names(diff_gene)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
write.csv(resdata,file = "gene_diff.csv",row.names = FALSE)

[STEP7] run R script
R < deseq2.R --no-save

[STEP8] see the results
          
#draw volcanoplot
[STEP1] download your csv
    gene_diff.csv
[STEP2] make minus log(padj)
    open the csv let N=-LOG10(G2)
    Rowname= minus_logpadj
[STEP3]
    go to https://huygens.science.uva.nl/VolcaNoseR2/
