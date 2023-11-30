# Whole Exome Sequencing (WXS/WES) data analysis pipeline

I wrote this SLURM bash script to automate a workflow for Whole Exome Sequencing (WXS) data analysis on a high-performance computing (HPC) cluster. It activates a conda environment that I created with the bioinformatics tools required to performs several tasks:

- **Prefetching and Downloading**: Downloads raw sequence data in FASTQ format from NCBI's Sequence Read Archive (SRA) database given your desired run accession.
- **Quality Control (QC)**:
  - Runs FastQC to assess the quality of the raw FASTQ files.
  - Trims the raw FASTQ files using TrimGalore and generates post-trim FastQC reports.
- **Read Alignment**:
  - Aligns trimmed reads to a human reference genome (GRCh38) using BWA-MEM.
- **Mark Duplicates**: Uses GATK's MarkDuplicatesSpark to mark PCR duplicates in the aligned BAM files.
- **Variant Calling**:
  - Uses bcftools mpileup and bcftools call to call variants from the aligned and marked BAM files. The variant calling parameters were selected to prioritize heightened sensitivity, even though this might result in increased false positives. The intention is to subsequently employ machine learning models to filter out these false positives.
- **Variant Subsetting**:
  - Intersects the called variants with vendor exome regions BED and Genome-in-a-Bottle (GIAB) high-confidence BED files using bedtools intersect to subset variants to those within specific regions. The intersection with the GIAB high-confidence BED file is to enable determination of false positive and false negative variant calls for the purpose of building machine learning models to mitigate the artifacts.
- **Output**: Generates intermediate files for each step placing them in their respective folders.
<img width="443" alt="image" src="https://github.com/felixm3/exome-sequencing/assets/43228120/64ade30c-31a4-46dd-817d-c94ac0100a84">


**Input Files Required**:
- Human reference genome (GRCh38)
- Vendor exome regions BED file
- GIAB high-confidence BED file

**Required Bioinformatics Tools**:
- BWA
- FastQC
- TrimGalore
- Samtools
- GATK
- bcftools
- bedtools

**Outputs**:
- FASTQ quality reports (FastQC)
- Trimmed FASTQ files
- Aligned BAM files
- Marked duplicate BAM files
- Variant calling results in VCF format
- VCF files subsetted to specified regions



The script logs system information, information on the various processing steps, information on the versions of the various tools used, and execution time for the full run. I used the output VCF files to create machine learning models of sequencing artifacts using Bioconductor and R packages (VariantAnnotation, caret, etc).
