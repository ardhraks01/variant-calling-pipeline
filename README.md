# Variant Calling Pipeline

A reproducible variant-calling workflow on simulated SARS-CoV-2 reads: quality control, trimming, alignment, variant calling, filtering and variant summary, run in WSL2 Ubuntu with a conda environment.

## Workflow
1. **Quality control and trimming:** fastp
2. **Alignment to the reference:** BWA
3. **Sorting and indexing:** samtools
4. **Variant calling and filtering:** bcftools (filter: QUAL >= 30 and DP >= 10)
5. **Summary and plots:** Python (biopython, pandas, matplotlib, pysam)

## Results
On simulated reads (wgsim, 1% mutation rate), **327 variants** were retained after filtering: **281 SNPs** and **46 indels**.
The filtered VCF is in `work/`, with plots of read quality and variant density.

## Installation
git clone https://github.com/YOURUSERNAME/variant-calling-pipeline.git
cd variant-calling-pipeline
conda env create -f environment.yml
conda activate varpipe

## Usage
Input reads are in `data/` and the reference is in `ref/`. Outputs are written to `work/`.

## Data
See `data/README.md` for the source of the test data.

## Limitations and future work
- Small practice dataset, not a clinical analysis
- Future: variant annotation (for example SnpEff) and a Snakemake or Nextflow version

## Author
Ardhra K S | MSc Genomic Science | [LinkedIn](https://www.linkedin.com/in/ardhra-k-s-3b4358275/) | ardhraks01@gmail.com

## Tool versions
- fastp: fastp 1.3.7
- bwa: [main] unrecognized command '--version'
- samtools: samtools 1.24
- bcftools: bcftools 1.24
