

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and Upregulation in scRNA-seq >
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.com/NinaGilshteyn](https://g>


import pandas as pd
import numpy as np 



df = pd.read_csv('cd56NKcellsdf_len8385.csv', header = 0, index_col = 0) # genes are columns and cells are rows 
df = df.drop('labels', axis = 1)

df = df.T  # cells are columns and genes are rows
Total_RNA_per_Cell = df.sum()

#print(df)
#print(Total_RNA_per_Cell)

Total_RNA_per_Cell.to_csv("Total_RNA_per_Cell.csv")

