import pandas as pd 
import numpy as np
import seaborn as sns


df = pd.read_csv('Full.csv', header = 0, index_col =0)

#df = pd.DataFrame({'gene1':[0,0,0], 'gene2': [0,0,0], 'gene3':[2,2,2], 'gene4': [0,0,0]})

print('df', df.head(3))
#df = df.T  #transpose df so that cells are the columns and genes are the rows
#print("df post transpose", df.shape, df.head(3))


# get sums of each column to see if any are zero 


sumz = df.sum() 

sumz = np.array(sumz)

zero = np.sum(sumz == 0)

print(zero)
