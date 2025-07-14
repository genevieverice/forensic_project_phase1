!pip install biopython

from Bio import SeqIO
import re
import pandas as pd

fasta_file = "/content/snp_targets_sequences(+-200bp).fa"
motif = "CACGTC" 
output_file = "/content/motif_positionsBmgBI_snps.csv"

results = []

for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = str(record.seq).upper()
    snp_id = record.id
    snp_position = len(sequence) // 2  

    for match in re.finditer(motif, sequence):
        motif_start = match.start()
        distance_to_snp = motif_start - snp_position

        results.append({
            "SNP_ID": snp_id,
            "SNP_Position": snp_position,
            "Motif_Start": motif_start,
            "Distance_Motif_to_SNP": distance_to_snp
        })

df = pd.DataFrame(results)
df.to_csv(output_file, index=False)

df.head()
