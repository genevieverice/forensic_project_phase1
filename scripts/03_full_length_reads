import pysam

bam_files = {
    "enzyme_treated": "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam",
    "sheared": "/Users/genevieverice/Documents/RP/gTUBE.hg38.Primary.bam",
    "control": "/Users/genevieverice/Documents/RP/Untreated.hg38.Primary.bam"
}

mt_chromosomes = ["chrM", "MT"]
MIN_LENGTH = 16000  

for label, path in bam_files.items():
    bamfile = pysam.AlignmentFile(path, "rb")
    mt_found = False

    for chrom in mt_chromosomes:
        if chrom in bamfile.references:
            mt_found = True
            total_mt_reads = 0
            full_length_reads = 0

            for read in bamfile.fetch(chrom):
                if not read.is_unmapped and read.query_length:
                    total_mt_reads += 1
                    if read.query_length >= MIN_LENGTH:
                        full_length_reads += 1

            percent = (full_length_reads / total_mt_reads * 100) if total_mt_reads > 0 else 0
            print(f"{label}:\n"
                  f"  Total mtDNA reads: {total_mt_reads}\n"
                  f"  Full-length reads (≥{MIN_LENGTH} bp): {full_length_reads}\n"
                  f"  Percent full-length: {percent:.2f}%\n")
            break

    if not mt_found:
        print(f"{label}: chrM/MT not found in BAM references.\n")
