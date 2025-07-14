import pysam

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.Primary.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.Primary.bam"}

mt_chromosomes = ["chrM", "MT"]

for label, bam_path in bam_files.items():
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    total_reads = sum(1 for read in bamfile.fetch(until_eof=True) if not read.is_unmapped)
    
    bamfile.reset()
    mt_count = 0
    mt_found = False

    for chrom in mt_chromosomes:
        if chrom in bamfile.references:
            mt_found = True
            mt_count = bamfile.count(contig=chrom)
            break

    bamfile.close()

    if mt_found:
        percent_mt = (mt_count / total_reads * 100) if total_reads > 0 else 0
        print(f"{label}:\n"
              f"  Total mapped reads: {total_reads}\n"
              f"  Mitochondrial reads ({chrom}): {mt_count}\n"
              f"  Percentage mtDNA: {percent_mt:.2f}%\n")
    else:
        print(f"{label}: mtDNA chromosome not found in BAM references.\n")
