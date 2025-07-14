from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

fasta_file = '/Users/genevieverice/hg38.fa'   		
bed_file = '/Users/genevieverice/horvathmarkers.bed' 	
output_fasta = '/Users/genevieverice/horvath_target_sequences(+-200bp).fa'  
								
flank_size = 200	# 200 bases up and downstream

genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

bed = pd.read_csv(bed_file, sep='\t', header=None, names=["Chromosome", "Start", "End", "CpG_ID"])  

records = []
for _, row in bed.iterrows():
    chrom = row['Chromosome']
    start = int(row['Start']) - flank_size
    end = int(row['End']) + flank_size
    name = row['CpG_ID']

    start = max(0, start)

    if chrom in genome:
        seq = genome[chrom].seq[start:end]        # Extract the sequence
        header = f"{name}|{chrom}:{start}-{end}"
        record = SeqRecord(seq, id=header, description="")
        records.append(record)
    else:
        print(f"Warning: Chromosome {chrom} not found in genome.")

SeqIO.write(records, output_fasta, "fasta")

print(f"{len(records)} sequences saved with ±200 bp flanks to {output_fasta}")
