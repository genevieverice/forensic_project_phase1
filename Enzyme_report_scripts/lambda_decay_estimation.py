import pysam
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

samples = {
    "Sample 1 (A)": "/content/sample_data/24D5647.chrM.bam",
    "Sample 3 (A)": "/content/sample_data/24D5655.chrM.bam",
    "Sample 5 (B)": "/content/sample_data/24D5649.chrM.bam"
}

min_length = 50
max_length = 2000
bin_size = 10
bins = np.arange(min_length, max_length + bin_size, bin_size)
centers = (bins[:-1] + bins[1:]) / 2

def decay_model(L, N0, lamb):
    return N0 * np.exp(-lamb * L)

lambda_values = {}

for sample_name, path in samples.items():
    lengths = []

    bam = pysam.AlignmentFile(path, "rb")
    for read in bam:
        if not read.is_unmapped and read.is_proper_pair == False:
            length = read.query_length
            if min_length <= length <= max_length:
                lengths.append(length)
    bam.close()

    hist, _ = np.histogram(lengths, bins=bins)

    try:
        popt, _ = curve_fit(decay_model, centers, hist, p0=(max(hist), 0.001), maxfev=10000)
        N0_fit, lambda_fit = popt
    except RuntimeError:
        N0_fit, lambda_fit = np.nan, np.nan
        print(f"Fit failed for {sample_name}")

    lambda_values[sample_name] = lambda_fit

    plt.figure(figsize=(8, 4))
    plt.bar(centers, hist, width=bin_size, alpha=0.6, label="Observed")
    if not np.isnan(lambda_fit):
        plt.plot(centers, decay_model(centers, *popt), 'r--', label=f"Fit: λ = {lambda_fit:.4f}")
    plt.title(f"Fragment Length Decay – {sample_name}")
    plt.xlabel("Fragment Length (bp)")
    plt.ylabel("Read Count")
    plt.legend()
    plt.tight_layout()
    plt.show()

plt.figure(figsize=(7, 5))
plt.bar(lambda_values.keys(), lambda_values.values(), color='mediumseagreen')
plt.ylabel("Estimated λ (Decay Constant)")
plt.title("DNA Degradation Score per Sample")
plt.xticks(rotation=15)
plt.tight_layout()
plt.show()

for sample, lam in lambda_values.items():
    print(f"{sample}: λ = {lam:.5f}")
