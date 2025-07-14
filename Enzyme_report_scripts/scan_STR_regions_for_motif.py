from Bio import SeqIO

fasta_file = '/content/ForenseqSTR_target_sequences(+-200bp).fa'
motif = 'GTTTAAAC'
flank_size = 200
motif_len = len(motif)

results = []

for record in SeqIO.parse(fasta_file, "fasta"):
    seq_id = record.id
    sequence = str(record.seq).upper()

    for i in range(len(sequence) - motif_len + 1):
        if sequence[i:i+motif_len] == motif:
            # Determine region
            if i < flank_size:
                region = '5_flank'
            elif i >= len(sequence) - flank_size:
                region = '3_flank'
            else:
                region = 'STR_region'

            results.append({
                'STR_ID': seq_id,
                'motif_position': i,
                'region': region,
                'relative_to_STR': i - flank_size
            })

import pandas as pd
df = pd.DataFrame(results)
df.head()
