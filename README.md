# forensic_project_phase1
This repository contains analysis scripts and documentation for a proof-of-concept experiment assessing mitochondrial DNA (mtDNA) degradation via enzyme digestion and nanopore sequencing.

Mitochondrial DNA (mtDNA) is circular and commonly used in forensic analyses. Estimating time since deposition (TsD) of blood is challenging. This study explores enzyme-based linearization methods (BsiWI, PmeI, SnaBI) combined with long-read sequencing to track fragment patterns over time.

Reads aligned against hg38 using Minimap2 and sorted to BAM.
- Only primary alignments retained:

To avoid overrepresentation of reads that span the circular mitochondrial genome junction (position 16,569 to 1), only primary alignments were retained for downstream analyses. During linear alignment, such junction-spanning reads are split into two alignments: one for the longer fragment (assigned as the primary alignment) and one for the shorter fragment (marked as secondary or supplementary).

To prevent double-counting of these reads, secondary (0x100) and supplementary (0x800) alignments were excluded using the samtools command. 

Because all three restriction enzyme cut sites were located in the second half of the mitochondrial genome, primary alignments of junction-spanning reads generally originated in the first half (starting near position 1 and ending at the enzyme-specific cut site).

All analyses that follow have been done on Primary alignments only. 
