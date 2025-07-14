!pip install pysam

import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore

bam_files = {
    "Sample 1 (A)": "/content/sample_data/24D5647.chrM.bam",
    "Sample 3 (A)": "/content/sample_data/24D5655.chrM.bam",
    "Sample 5 (B)": "/content/sample_data/24D5649.chrM.bam"
}

reference = "chrM"
genome_length = 16569
bin_size = 50
num_bins = genome_length // bin_size

def get_binned_coverage(bam_path, ref, bin_size, genome_length):
  bamfile = pysam.AlignmentFile(bam_path, "rb")
  binned_coverage = np.zeros(genome_length // bin_size)

  for pileupcolumn in bamfile.pileup(ref, truncate = True):
    bin_index = pileupcolumn.reference_pos // bin_size
    if bin_index < len(binned_coverage):
      binned_coverage[bin_index] += pileupcolumn.nsegments

  bamfile.close()
  return binned_coverage

coverage_matrix =[]
labels = []

for name, path in bam_files.items():
  cov = get_binned_coverage(path, reference, bin_size, genome_length)
  coverage_matrix.append(cov)
  labels.append(name)

coverage_array = np.array(coverage_matrix)

zscore_coverage = zscore(coverage_array, axis=1)

plt.figure(figsize=(14, 5))
sns.heatmap(
    zscore_coverage,
    cmap="coolwarm",  # diverging colormap is useful for z-scores
    yticklabels=labels,
    cbar_kws={'label': 'Z-score of Read Depth'}
)
plt.title("Z-score Normalized Read Depth Heatmap Across mtDNA")
plt.xlabel("Genome Position (bin index)")
plt.ylabel("Sample")
plt.xticks(
    ticks=np.linspace(0, zscore_coverage.shape[1], 6),
    labels=[f"{int(x)}bp" for x in np.linspace(0, genome_length, 6)]
)

bin_count = normalized_coverage_array.shape[1]
xtick_locs = np.linspace(0, bin_count, 6)
xtick_labels = [f"{int(x)}bp" for x in np.linspace(0, 16569, 6)]

plt.tight_layout()
plt.show()
