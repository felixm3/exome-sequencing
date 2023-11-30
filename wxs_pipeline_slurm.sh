#!/bin/bash
#
#SBATCH --job-name=wxs_pipeline_slurm.sh
#SBATCH --output=/home/fmbuga/tools/slurm_scripts/slurm_out_err/wxs_pipeline_%j.out
#SBATCH --error=/home/fmbuga/tools/slurm_scripts/slurm_out_err/wxs_pipeline_%j.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G

datenow=$(date)
echo -e "START TIME: $datenow\n"

srun hostname
echo ""

start=$(date +%s)
echo -e "hostname: $HOSTNAME\n"

# SRA Accession ID of WXS data
sra_accession=$1
echo -e "SRA Accession ID: $sra_accession\n"

# activate conda environment with wxs variant-calling tools
eval "$(conda shell.bash hook)"
conda activate wxs-pipeline

# log info on conda env and package versions
command="conda info"
echo -e "$command\n"
eval $command

command="conda list"
echo -e "\n$command\n"
eval $command

# make directory to store output files
current_date=$(date +"%Y_%m_%d_%H_%M_%S")
echo -e "\n$current_date\n"
mkdir $current_date
cd $current_date
echo -e "Current working directory: $(pwd)\n"

## prefetch + fasterq-dump SRA file
command="prefetch $sra_accession"
echo $command
eval $command
echo -e "\n prefetch SRA accession COMPLETE.\n"

command="fasterq-dump $sra_accession \
            --outdir 00_fastq_raw"
echo $command
eval $command
echo -e "\n fasterq-dump SRA accession COMPLETE.\n"

## QC: FastQC
input_folder="00_fastq_raw"
output_folder="01_fastq_raw_FastQC"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Use wildcards to match the FASTQ files
forward_fastq="$input_folder/*_1.fastq"
reverse_fastq="$input_folder/*_2.fastq"
echo $forward_fastq
echo $reverse_fastq

# Run FastQC on the matched FASTQ files
fastqc \
    --outdir "$output_folder" \
    --threads 2 \
    $forward_fastq $reverse_fastq

echo -e "\nFastQC analysis completed. Results are in $output_folder\n"

## Trimming: TrimGalore
output_folder="02_fastq_trimmed"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

trim_galore \
    --fastqc \
    --output_dir "$output_folder" \
    --cores 8 \
    --paired \
    $forward_fastq $reverse_fastq

echo -e "\nTrimGalore trimming and post-trim FastQC completed. Results are in $output_folder\n"

## Mapping: bwa mem
input_folder="02_fastq_trimmed"
output_folder="03_bam_raw"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Use wildcards to match the FASTQ files
forward_fastq="$input_folder/*_1.fq"
reverse_fastq="$input_folder/*_2.fq"
echo $forward_fastq
echo $reverse_fastq

out_name=$(basename $forward_fastq _1_val_1.fq)

bwa mem \
    -t 28 \
    -R "@RG\tID:$out_name\tSM:$out_name" \
    ~/tools/GIAB/GRCh38/refs/GRCh38_GIABv3.fasta \
    $forward_fastq $reverse_fastq | \
        samtools view -uh | \
            samtools sort \
                -o "$output_folder/$out_name".bam \
                -@ 28 \
                --write-index \
                --verbosity 9 \
                -m 2G

echo -e "\nmapping with bwa mem completed. Results are in $output_folder\n"

## Mark duplicates: GATK MarkDuplicatesSpark
input_bam="$output_folder/$out_name".bam
output_folder="04_bam_markdup"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Extract the file name without extension
file_name=$(basename "$input_bam" .bam)

echo -e "\nProcessing $input_bam\n"

# Run GATK SortSam
gatk SortSam \
    -I "$input_bam" \
    -O "$output_folder/${file_name}_temp.bam" \
    -SO queryname

# Run GATK MarkDuplicatesSpark
gatk MarkDuplicatesSpark \
    -I "$output_folder/${file_name}_temp.bam" \
    -O "$output_folder/${file_name}_markdup.bam" \
    -M "$output_folder/${file_name}_markdup_metrics.txt" \
    --conf "spark.local.dir=$output_folder"

rm "$output_folder/${file_name}_temp.bam"

echo -e "\nmarking duplicates with gatk MarkDuplicatesSpark completed. Results are in $output_folder\n"

## Call variants: bcftools mpileup --> bcftools call
input_bam="$output_folder/${file_name}_markdup.bam"
output_folder="05_vcf_raw"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Extract the file name without extension
file_name=$(basename "$input_bam" .bam)

echo -e "\nProcessing $input_bam\n"

# Run bcftools mpileup --> bcftools call
min_MQ=0 
bcftools mpileup \
    -Ou \
    --min-MQ $min_MQ \
    --config illumina \
    -f ~/tools/GIAB/GRCh38/refs/GRCh38_GIABv3.fasta \
    $input_bam |\
    bcftools call \
        --ploidy GRCh38 \
        -m -v \
        -P 0.99 \
        -o "$output_folder/$file_name"_bcftools_P_099_minMQ_"$min_MQ".vcf

echo -e "\nvariant calling with bcftools mpileup --> bcftools call completed. Results are in $output_folder\n"

## Intersect with vendor exome regions BED and GIAB hi-confidence BED
input_vcf="$output_folder/$file_name"_bcftools_P_099_minMQ_"$min_MQ".vcf
output_folder="06_vcf_isec"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Extract the file name without extension
file_name=$(basename $input_vcf .vcf)

echo -e "\nProcessing $input_vcf\n"

# bedtools intersect with WXS regions
bedtools intersect \
    -a $input_vcf \
    -b ~/tools/Agilent/SureSelect/SureSelect_Human_All_Exon_V4/Agilent_S04380110_hs_hg19/S04380110_Covered_hg19toHg38.bed \
    -header > "$output_folder/$file_name"_temp.vcf

# bedtools intersect with GIAB
bedtools intersect \
    -a "$output_folder/$file_name"_temp.vcf \
    -b ~/tools/GIAB/GRCh38/HG001_NA12878/v4_2_1/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    -header > "$output_folder/$file_name"_isec.vcf

# create GIAB WXS vcf for comparison
bedtools intersect \
    -a ~/tools/GIAB/GRCh38/HG001_NA12878/v4_2_1/HG001_GRCh38_v4_2_1.vcf \
    -b ~/tools/Agilent/SureSelect/SureSelect_Human_All_Exon_V4/Agilent_S04380110_hs_hg19/S04380110_Covered_hg19toHg38.bed \
    -header > "$output_folder/"GIAB_WXS_isec.vcf

rm "$output_folder/$file_name"_temp.vcf

echo -e "\n VCF intersection with bedtools intersect complete. Results are in $output_folder \n"

echo -e "\n WORKFLOW COMPLETE \n"

##################

echo ""
end=$(date +%s)
echo -e "END TIME: $(date)\n"

runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"


