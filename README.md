# forensic_project_phase1
This repository contains analysis scripts and documentation for a proof-of-concept experiment assessing mitochondrial DNA (mtDNA) degradation via enzyme digestion and nanopore sequencing.

Mitochondrial DNA (mtDNA) is circular and commonly used in forensic analyses. Estimating time since deposition (TsD) of blood is challenging. This study explores enzyme-based linearization methods (BsiWI, PmeI, SnaBI) combined with long-read sequencing to track fragment patterns over time.

Reads aligned against hg38 using Minimap2 and sorted to BAM.
- Only primary alignments retained:

To avoid overrepresentation of reads that span the circular mitochondrial genome junction (position 16,569 to 1), only primary alignments were retained for downstream analyses. During linear alignment, such junction-spanning reads are split into two alignments: one for the longer fragment (assigned as the primary alignment) and one for the shorter fragment (marked as secondary or supplementary).

To prevent double-counting of these reads, secondary (0x100) and supplementary (0x800) alignments were excluded using the samtools command. 

Because all three restriction enzyme cut sites were located in the second half of the mitochondrial genome, primary alignments of junction-spanning reads generally originated in the first half (starting near position 1 and ending at the enzyme-specific cut site).

All analyses that follow have been done on Primary alignments only. 

All bioinformatic analyses were conducted using Python 3. Data handling and visualization were performed with a combination of open-source packages: pysam for BAM file parsing, numpy and pandas for numerical computations and data organization, and matplotlib and seaborn for generating plots.

Fragment length distributions were generated using pysam and matplotlib. Aligned reads were extracted from each BAM file (enzyme-treated, sheared, and untreated control) and only reads mapped to the mitochondrial chromosome ("chrM") were included. Read lengths were retrieved using the query_length attribute of each read. Histograms were plotted to visualize the distribution of mtDNA fragment sizes for each sample. A reference line at 16,569 bp was added to indicate full-length mitochondrial genomes. Because all reads had been filtered to include only primary alignments, the resulting distributions reflect representative molecule lengths while minimizing overcounting.

Coverage depth across the mitochondrial genome was calculated using pysam.pileup, which provides per-base read depth. For each sample, coverage was calculated across the full length of chrM. Plots were generated using matplotlib, with drop-offs in coverage observed at enzyme-specific cut sites. These drop-offs were used to evaluate the relative efficiency of each enzyme. For each cut site, the average coverage within a 100 bp window upstream ("before") and downstream ("after") of the cut was calculated using numpy. The absolute difference in coverage before and after the cut site was used as an indication for enzyme cutting efficiency.

