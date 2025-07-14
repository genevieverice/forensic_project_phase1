import pysam
import numpy as np
import pandas as pd

bam_path = "/Users/genevieverice/Documents/RP/Enzymes.hg38.Primary.bam"
mt_chromosomes = ["chrM", "MT"]
cut_sites = {
    "BsiWI": (8897, 9001),
    "PmeI": (10416, 10417),
    "SnaBI": (10736, 10737),
}
window_size = 100

bamfile = pysam.AlignmentFile(bam_path, "rb")
mt_ref = next((chrom for chrom in mt_chromosomes if chrom in bamfile.references), None)
if mt_ref is None:
    raise ValueError("Mitochondrial chromosome not found")

chr_len = bamfile.get_reference_length(mt_ref)
coverage = [0] * chr_len
for pileupcolumn in bamfile.pileup(mt_ref, truncate=True):
    pos = pileupcolumn.reference_pos
    coverage[pos] = pileupcolumn.nsegments
bamfile.close()

drop_results = {}
for enzyme, (cut_before, cut_after) in cut_sites.items():
    before_region = coverage[max(cut_before - window_size,0):cut_before]
    after_region = coverage[cut_after:cut_after +window_size]
    
    avg_before = np.mean(before_region)
    avg_after = np.mean(after_region)
    drop = avg_before - avg_after
    percent_drop = (drop / avg_before *100) if avg_before !=0 else 0
    
    drop_results[enzyme] = {
        "Cut before": cut_before,
        "Cut after": cut_after,
        "avg before": avg_before,
        "avg after": avg_after,
        "coverage drop": drop,
        "% drop": percent_drop
    }


df = pd.DataFrame.from_dict(drop_results, orient="index")
print(df)
