import subprocess, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam

WORK = Path("work"); WORK.mkdir(exist_ok=True)
REF  = "ref/genome.fa"
GFF  = "ref/genome.gff3"
R1, R2 = "data/sim_R1.fq.gz", "data/sim_R2.fq.gz"
THREADS = 4

def run(cmd):
    print(">>", cmd)
    subprocess.run(cmd, shell=True, check=True)

def read_fastq_stats(path, n=50000):
    lengths, quals = [], []
    with gzip.open(path, "rt") as fh:
        for i, rec in enumerate(SeqIO.parse(fh, "fastq")):
            if i >= n: break
            lengths.append(len(rec.seq))
            quals.append(np.mean(rec.letter_annotations["phred_quality"]))
    return np.array(lengths), np.array(quals)

# 1. QC
def qc_plot(r1, label):
    lengths, quals = read_fastq_stats(r1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(lengths, bins=30); ax[0].set_title(f"{label}: read length")
    ax[1].hist(quals, bins=30);   ax[1].set_title(f"{label}: mean Phred quality")
    plt.tight_layout(); plt.savefig(WORK / f"qc_{label}.png", dpi=150); plt.close()
    print(f"{label}: reads={len(lengths)}, mean len={lengths.mean():.1f}, "
          f"mean Q={quals.mean():.1f}, %reads Q>=30: {(quals>=30).mean()*100:.1f}")

qc_plot(R1, "raw")
pass  # fastqc skipped

# 2. Trimming
T1, T2 = WORK/"trim_R1.fq.gz", WORK/"trim_R2.fq.gz"
run(f"fastp -i {R1} -I {R2} -o {T1} -O {T2} "
    f"--detect_adapter_for_pe -q 15 -l 30 "
    f"-w {THREADS} -j {WORK}/fastp.json -h {WORK}/fastp.html")
qc_plot(str(T1), "trimmed")

# 3. Alignment
run(f"bwa index {REF}")
BAM = WORK / "sample.sorted.bam"
run(f"bwa mem -t {THREADS} -R '@RG\\tID:s1\\tSM:sample1\\tPL:ILLUMINA' "
    f"{REF} {T1} {T2} | samtools sort -@ {THREADS} -o {BAM}")
run(f"samtools index {BAM}")
run(f"samtools flagstat {BAM} > {WORK}/flagstat.txt")
run(f"samtools coverage {BAM} > {WORK}/coverage.txt")

# 4. Variant calling
run(f"samtools faidx {REF}")
VCF = WORK / "variants.vcf.gz"
run(f"bcftools mpileup -f {REF} -Ou {BAM} | bcftools call -mv -Oz -o {VCF}")
FILT = WORK / "variants.filtered.vcf.gz"
run(f"bcftools filter -i 'QUAL>=30 && DP>=10' -Oz -o {FILT} {VCF}")
run(f"bcftools index {FILT}")

# 5. Annotation
def load_genes(gff):
    cols = ["chrom","src","type","start","end","score","strand","phase","attr"]
    df = pd.read_csv(gff, sep="\t", comment="#", names=cols)
    df = df[df["type"].isin(["gene", "CDS"])].copy()
    df["name"] = df["attr"].str.extract(r"(?:Name|gene|ID)=([^;]+)")
    return df[["chrom","start","end","strand","name"]]

genes = load_genes(GFF)

def annotate(chrom, pos):
    hit = genes[(genes.chrom == chrom) & (genes.start <= pos) & (genes.end >= pos)]
    return hit.iloc[0]["name"] if len(hit) else "intergenic"

# 6. Variant analysis
rows = []
for rec in pysam.VariantFile(str(FILT)):
    alt = rec.alts[0]
    is_snp = len(rec.ref) == 1 and len(alt) == 1
    rows.append({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt,
                 "qual": rec.qual, "depth": rec.info.get("DP"),
                 "type": "SNP" if is_snp else "INDEL",
                 "gene": annotate(rec.chrom, rec.pos)})

df = pd.DataFrame(rows)
if df.empty:
    print("\nNo variants passed the filters (QUAL>=30, DP>=10).")
else:
    df.to_csv(WORK / "annotated_variants.csv", index=False)
    print("\nTotal variants:", len(df))
    print(df["type"].value_counts())
    print("\nVariants per gene (top 10):")
    print(df["gene"].value_counts().head(10))

    purines, pyrimidines = {"A","G"}, {"C","T"}
    snps = df[df["type"] == "SNP"]
    ti = sum(((r in purines and a in purines) or (r in pyrimidines and a in pyrimidines))
             for r, a in zip(snps.ref, snps.alt))
    tv = len(snps) - ti
    print(f"\nTi/Tv ratio: {ti/tv:.2f}" if tv else "\nNo transversions")

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.hist(df["pos"], bins=100)
    ax.set_xlabel("Genomic position"); ax.set_ylabel("Variant count")
    ax.set_title("Variant density")
    plt.tight_layout(); plt.savefig(WORK / "variant_density.png", dpi=150)
print("\nDone. Outputs are in the work/ folder.")
