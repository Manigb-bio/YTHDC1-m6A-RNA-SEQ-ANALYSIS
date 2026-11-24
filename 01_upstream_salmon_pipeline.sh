#!/bin/bash
# Upstream Bulk RNA-Seq Pipeline: Raw SRA -> FastQC/MultiQC -> Salmon Quantification
# Dataset: NCBI GEO GSE74397 / SRA SRP065329 (Mus musculus mESCs)

conda activate rnaseq
mkdir -p ~/ythdc1_work/{raw_fastq,qc_reports,salmon_index,salmon_quant}

# 1. Download subsampled raw paired-end FASTQ reads
cd ~/ythdc1_work/raw_fastq
fastq-dump -X 500000 --split-files --gzip SRR2830633

# 2. Quality Control (FastQC + MultiQC)
fastqc ~/ythdc1_work/raw_fastq/*.fastq.gz --outdir ~/ythdc1_work/qc_reports
cd ~/ythdc1_work/qc_reports && multiqc .

# 3. Download GENCODE Mouse Transcriptome & Build Salmon Index
cd ~/ythdc1_work/salmon_index
rsync -av --progress rsync://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.pc_transcripts.fa.gz .
salmon index -t gencode.vM35.pc_transcripts.fa.gz -i mouse_index -p 4 --gencode

# 4. Quantify Transcript Expression with Salmon
cd ~/ythdc1_work
salmon quant -i salmon_index/mouse_index -l A \
  -1 raw_fastq/SRR2830633_1.fastq.gz \
  -2 raw_fastq/SRR2830633_2.fastq.gz \
  -p 4 --validateMappings \
  -o salmon_quant/SRR2830633_quant

# 5. Copy Outputs to D: Drive & Clean Up Temporary Linux Storage
cp -r ~/ythdc1_work/qc_reports/* /mnt/d/Ythdc1_SRA_Pipeline/qc_reports/
cp -r ~/ythdc1_work/salmon_quant/* /mnt/d/Ythdc1_SRA_Pipeline/salmon_quant/
rm -rf ~/ythdc1_work
