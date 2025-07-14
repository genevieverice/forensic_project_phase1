import pysam
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.Primary.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.Primary.bam"
}

colors = {
    "enzyme_treated": "cornflowerblue",
    "sheared": "green",
    "control": "orange"
}

mt_chromosomes = ["chrM", "MT"]
read_lengths = {}

for label, path in bam_files.items():
    read_lengths[label] = []
    try:
        bam = pysam.AlignmentFile(path, "rb")
        for chrom in mt_chromosomes:
            if chrom in bam.references:
                for read in bam.fetch(chrom):
                    if not read.is_unmapped and read.query_length is not None:
                        read_lengths[label].append(read.query_length)
                break
        bam.close()
    except Exception as e:
        print(f"Error with {label}: {e}")

plt.figure(figsize=(12, 6))

for label, lengths in read_lengths.items():
    sns.histplot(
        lengths, bins=100, kde=False, stat="density", element="step", fill=True,
        label=label, color=colors[label], alpha=0.6
    )

plt.axvline(16569, color='gray', linestyle='--', label='Full-length mtDNA (16,569 bp)')
plt.title("mtDNA Read Length Distribution")
plt.xlabel("Read Length (bp)")
plt.ylabel("Density")
plt.xlim(0, 17500)
plt.legend()
plt.tight_layout()
plt.show()
