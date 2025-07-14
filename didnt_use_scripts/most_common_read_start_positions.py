import pysam
from collections import Counter

bam = pysam.AlignmentFile("/Users/genevieverice/Documents/RP/Enzymes.hg38.bam", "rb")
starts = [read.reference_start for read in bam if read.reference_name in ["MT", "chrM"]]
bam.close()

Counter(starts).most_common(10)
