import pysam
import matplotlib.pyplot as plt
import seaborn as sns

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.bam"
}

for label, bam_path in bam_files.items():
    print(f"Processing {label}...")

    starts = []
    ends = []

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_name in ["MT", "chrM"]:
                starts.append(read.reference_start)
                ends.append(read.reference_end)

        bam.close()

        plt.figure(figsize=(12, 5))
        sns.kdeplot(starts, label="Start positions", color="blue", bw_adjust=0.2)
        sns.kdeplot(ends, label="End positions", color="red", bw_adjust=0.2)

        plt.title(f"mtDNA Read Start and End Position Density: {label}")
        plt.xlabel("Position on chrM")
        plt.ylabel("Density")
        plt.legend()
        plt.grid(True)
        plt.xlim(0, 16569)
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Error processing {label}: {e}")
