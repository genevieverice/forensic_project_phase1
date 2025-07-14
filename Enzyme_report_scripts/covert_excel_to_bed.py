import pandas as pd

excel_file = '/content/AllTargetRegionsWithoutFlanksCopy.xlsx'  
sheet_name = 'CpGs Horvath'                                  
output_bed = '/content/horvathmarkers.bed'                    

df = pd.read_excel(excel_file, sheet_name=sheet_name)

print(df.columns)
df.head()

bed_df = df[['Chromosome', 'StartPosition_hg38', 'EndPosition_hg38', 'CpG_ID']]

bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
