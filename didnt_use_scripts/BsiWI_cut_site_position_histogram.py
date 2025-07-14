import pysam
import matplotlib.pyplot as plt
import numpy as np

bam = pysam.AlignmentFile("/Users/genevieverice/Documents/RP/Enzymes.hg38.bam", "rb")

starts = []
ends = []

for read in bam.fetch("chrM"):
    if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
        starts.append(read.reference_start)
        ends.append(read.reference_end)

bam.close()

cut_site = 8997
buffer = 10
bins = np.arange(cut_site - buffer, cut_site + buffer + 1)

start_counts, _ = np.histogram(starts, bins=bins)
end_counts, _ = np.histogram(ends, bins=bins)

plt.figure(figsize=(10, 5))
plt.bar(bins[:-1], start_counts, width=1, color="blue", alpha=0.6, label="Start Positions")
plt.bar(bins[:-1], end_counts, width=1, color="red", alpha=0.6, label="End Positions")
plt.axvline(cut_site, color="gray", linestyle="--", label="BsiWI Cut Site (8997)")
plt.xlabel("Position on chrM")
plt.ylabel("Read Count")
plt.title("Read Start and End Positions Near BsiWI Cut Site")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
