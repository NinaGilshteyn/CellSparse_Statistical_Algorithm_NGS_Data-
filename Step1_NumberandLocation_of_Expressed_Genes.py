

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and Upregulation in scRNA-seq D>
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.com/NinaGilshteyn](https://gi>


import pandas as pd
import numpy as np


df = pd.read_csv('cd56NKcellsdf_len8385.csv', header = 0, index_col =0)
print('df pretranspose', df.head(3))

df = df.drop('labels', axis = 1)
print('df pretranspose', df.head(3))



df = df.T  #transpose df so that cells are the columns and genes are the rows 
print("df post transpose", df.shape, df.head(3))

def find_non_zero_indices(matrix):
    num_rows = len(matrix)
    num_columns = len(matrix[0])
    #print(num_rows)
    #print(num_columns)

    non_zero_indices_per_column = []

    for column in range(num_columns):
        non_zero_indices = []
        for row in range(num_rows):
            if matrix[row][column] != 0:
                non_zero_indices.append(row)  ## row value in column that does not equal zero 
        non_zero_indices_per_column.append(non_zero_indices)

    return non_zero_indices_per_column

# Get the non-zero indices for each column
matrix = df.values
results = find_non_zero_indices(matrix)
#print("nonzero row indices for each columns", results)

#now we have the information about the data, we can count how many non-zero row values are in each column 
def count_sublist_elements(data):  
    counts = []  # To store the counts of elements in each sublist
    for sublist in data:
        sublist_count = len(sublist)  # Count the elements in the current sublist
        counts.append(sublist_count)
    return counts

element_counts = count_sublist_elements(results)
#print("Element counts:", element_counts)

### writing everything to disk 
element_counts =pd.DataFrame(element_counts)
#print("element_counts", element_counts)
element_counts.to_csv("AmountExpressedGenes.csv")

results = pd.DataFrame(results)
results.to_csv('Coordinates_of_Expressed_Genes.csv')


