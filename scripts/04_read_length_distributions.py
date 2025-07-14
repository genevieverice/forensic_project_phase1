import matplotlib.pyplot as plt
import pysam

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.Primary.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.Primary.bam"
}

colors = {
    "enzyme treated": "cornflowerblue",  
    "sheared": "orange",        
    "control": "green"         
}

mt_read_lengths = {}

# Extract read lengths from chrM
for label, bam_path in bam_files.items():
    mt_read_lengths[label] = []
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.reference_name in ["MT", "chrM"]:
                if read.query_length:
                    mt_read_lengths[label].append(read.query_length)
        bam.close()
    except Exception as e:
        print(f"Error processing {label}: {e}")

# Plotting raw counts
for label, lengths in mt_read_lengths.items():
    plt.figure(figsize=(10, 5))
    plt.hist(
        lengths, bins=100, alpha=0.7,
        color=colors.get(label), histtype='stepfilled'
    )
    plt.axvline(16569, color='gray', linestyle='--', label='Full-length mtDNA (16,569 bp)')
    plt.title(f"mtDNA Read Length Distribution (including alignments): {label}")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Count")
    plt.xlim(0, 17500)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
