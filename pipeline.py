import subprocess, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam

# ---------- Settings ----------
WORK = Path("work"); WORK.mkdir(exist_ok=True)
for f in WORK.iterdir():              # clear old outputs so results are never stale
    if f.is_file():
        f.unlink()
REF = "ref/genome.fa"
GFF = "ref/genome.gff3"
R1, R2 = "data/sim_R1.fq.gz", "data/sim_R2.fq.gz"
THREADS = 4
MIN_QUAL, MIN_DEPTH = 30, 10

def run(cmd):
    print(">>", cmd)
    subprocess.run(cmd, shell=True, check=True)

# ---------- 1-2. Import + Quality control ----------
def read_fastq_stats(path, n=50000):
    lengths, quals = [], []
    with gzip.open(path, "rt") as fh:
        for i, rec in enumerate(SeqIO.parse(fh, "fastq")):
            if i >= n: break
            lengths.append(len(rec.seq))
            quals.append(np.mean(rec.letter_annotations["phred_quality"]))
    return np.array(lengths), np.array(quals)

def qc_plot(path, label):
    lengths, quals = read_fastq_stats(path)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(lengths, bins=30); ax[0].set_title(f"{label}: read length")
    ax[1].hist(quals, bins=30);   ax[1].set_title(f"{label}: mean Phred quality")
    plt.tight_layout(); plt.savefig(WORK / f"qc_{label}.png", dpi=150); plt.close()
    print(f"{label}: reads={len(lengths)}, mean len={lengths.mean():.1f}, "
          f"mean Q={quals.mean():.1f}, %reads Q>=30: {(quals>=30).mean()*100:.1f}")

qc_plot(R1, "raw")

# ---------- 3. Trimming ----------
T1, T2 = WORK/"trim_R1.fq.gz", WORK/"trim_R2.fq.gz"
run(f"fastp -i {R1} -I {R2} -o {T1} -O {T2} "
    f"--detect_adapter_for_pe -q 15 -l 30 "
    f"-w {THREADS} -j {WORK}/fastp.json -h {WORK}/fastp.html")
qc_plot(str(T1), "trimmed")

# ---------- 4. Alignment ----------
run(f"bwa index {REF}")
BAM = WORK / "sample.sorted.bam"
run(f"bwa mem -t {THREADS} -R '@RG\\tID:s1\\tSM:sample1\\tPL:ILLUMINA' "
    f"{REF} {T1} {T2} | samtools sort -@ {THREADS} -o {BAM}")
run(f"samtools index {BAM}")
run(f"samtools flagstat {BAM} > {WORK}/flagstat.txt")
run(f"samtools coverage {BAM} > {WORK}/coverage.txt")

# ---------- 5. Variant calling ----------
run(f"samtools faidx {REF}")
VCF = WORK / "variants.vcf.gz"
run(f"bcftools mpileup -f {REF} -Ou {BAM} | bcftools call -mv --ploidy 1 -Oz -o {VCF}")
FILT = WORK / "variants.filtered.vcf.gz"
run(f"bcftools filter -i 'QUAL>={MIN_QUAL} && DP>={MIN_DEPTH}' -Oz -o {FILT} {VCF}")
run(f"bcftools index {FILT}")

# ---------- 6. Annotation ----------
def parse_attr(s):
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def pick_name(attr):
    d = parse_attr(attr)
    for key in ("gene", "Name", "gene_name", "locus_tag", "ID"):
        if key in d:
            return d[key]
    return "unnamed_feature"

def load_features(gff):
    cols = ["chrom","src","type","start","end","score","strand","phase","attr"]
    df = pd.read_csv(gff, sep="\t", comment="#", names=cols)
    df = df[df["type"] == "CDS"].copy()
    df["name"] = df["attr"].apply(pick_name)
    df["rank"] = (df["type"] != "gene").astype(int)   # prefer gene over CDS
    return df.sort_values("rank")[["chrom","start","end","name"]]

features = load_features(GFF)

def annotate(chrom, pos):
    hit = features[(features.chrom == chrom) & (features.start <= pos) & (features.end >= pos)]
    return hit.iloc[0]["name"] if len(hit) else "intergenic"

# ---------- 7. Variant analysis ----------
rows = []
for rec in pysam.VariantFile(str(FILT)):
    if rec.alts is None:
        continue
    alt = rec.alts[0]
    is_snp = len(rec.ref) == 1 and len(alt) == 1
    rows.append({"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt,
                 "qual": rec.qual, "depth": rec.info.get("DP"),
                 "type": "SNP" if is_snp else "INDEL",
                 "gene": annotate(rec.chrom, rec.pos)})

df = pd.DataFrame(rows)
if df.empty:
    print(f"\nNo variants passed the filters (QUAL>={MIN_QUAL}, DP>={MIN_DEPTH}).")
else:
    df.to_csv(WORK / "annotated_variants.csv", index=False)
    print("\nTotal variants:", len(df))
    print(df["type"].value_counts())
    print("\nVariants per gene:")
    print(df["gene"].value_counts())

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

print("\n===== ALIGNMENT SUMMARY =====")
print(open(WORK / "coverage.txt").read())
print("\nDone. Outputs are in the work/ folder.")
