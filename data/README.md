# Data

The pipeline runs on **simulated** Illumina reads generated from a real reference genome.

## Reference (from nf-core test-datasets, SARS-CoV-2)

Source: https://github.com/nf-core/test-datasets/tree/modules/data/genomics/sarscov2/genome

    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta -O ref/genome.fa
    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.gff3 -O ref/genome.gff3

## Simulated reads (wgsim)

30,000 paired-end reads (2 x 150 bp), 1% mutation rate, random seed 42, haplotype mode (-h, one haploid genome like a virus):

    wgsim -h -S 42 -N 30000 -1 150 -2 150 -r 0.01 ref/genome.fa data/sim_R1.fq data/sim_R2.fq
    gzip data/sim_R1.fq data/sim_R2.fq

Because the reads are simulated, the variants called are simulated mutations, not real biological variants.

## Not used by the pipeline

sample_R1.fastq.gz and sample_R2.fastq.gz are the small real nf-core test reads (test_1/test_2.fastq.gz). They were downloaded but are not used in the pipeline run.
