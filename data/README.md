# Data

All test data comes from the nf-core test-datasets repository (SARS-CoV-2, `modules` branch):
https://github.com/nf-core/test-datasets/tree/modules/data/genomics/sarscov2

- Reference genome: genome/genome.fasta (saved as ref/genome.fa)
- Gene annotation: genome/genome.gff3
- Reads: illumina/fastq/test_1.fastq.gz and test_2.fastq.gz (saved as data/sample_R1.fastq.gz and data/sample_R2.fastq.gz)

Download:
wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta -O ref/genome.fa
wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.gff3 -O ref/genome.gff3
wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz -O data/sample_R1.fastq.gz
wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz -O data/sample_R2.fastq.gz
