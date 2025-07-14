import pysam
import matplotlib.pyplot as plt

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.Primary.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.Primary.bam"
}

mt_chromosomes = ["chrM", "MT"]

plt.figure(figsize=(14, 6))
colors = {"enzyme_treated": "blue", "sheared": "orange", "control": "green"}

for label, path in bam_files.items():
    bamfile = pysam.AlignmentFile(path, "rb")
    mt_ref = None

    for chrom in mt_chromosomes:
        if chrom in bamfile.references:
            mt_ref = chrom
            break

    if mt_ref is None:
        print(f"{label}: chrM/MT not found.")
        bamfile.close()
        continue
    chr_len = bamfile.get_reference_length(mt_ref)
    coverage = [0] * chr_len

    for pileupcolumn in bamfile.pileup(mt_ref, truncate=True):
        pos = pileupcolumn.reference_pos
        coverage[pos] = pileupcolumn.nsegments

    bamfile.close()
 
    plt.plot(range(1, chr_len + 1), coverage, label=label, color=colors[label], alpha=0.8)

plt.title("Coverage Across chrM (mtDNA)")
plt.xlabel("Position on mtDNA (bp)")
plt.ylabel("Coverage Depth")
plt.legend()
plt.tight_layout()
plt.xlim(0,17500)
plt.grid(True)
plt.show()
