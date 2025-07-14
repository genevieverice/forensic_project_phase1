import pysam

bam_path = "/Users/genevieverice/Documents/RP/Enzymes.hg38.bam"

enzymes = {
    "BsiWI (GTAC)": 8997,
    "PmeI (blunt)": 10414,
    "SnaBI (blunt)": 10734
}
buffer = 10

enzyme_start_counts = {name: 0 for name in enzymes}
enzyme_end_counts = {name: 0 for name in enzymes}

bam = pysam.AlignmentFile(bam_path, "rb")

for read in bam.fetch("chrM"):  
    if read.is_unmapped:
        continue
    
    start = read.reference_start
    end = read.reference_end

    for enzyme, site in enzymes.items():
        if abs(start - site) <= buffer:
            enzyme_start_counts[enzyme] += 1
        if abs(end - site) <= buffer:
            enzyme_end_counts[enzyme] += 1

bam.close()

print("Read starts near cut sites (10 bp):")
for enzyme, count in enzyme_start_counts.items():
    print(f"  {enzyme}: {count} reads")

print("\nRead ends near cut sites (10 bp):")
for enzyme, count in enzyme_end_counts.items():
    print(f"  {enzyme}: {count} reads")
