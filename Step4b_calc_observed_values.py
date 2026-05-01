import pandas as pd 
import numpy as np 
import seaborn as sns 




df =  pd.read_csv('Full.csv', header = 0, index_col = 0)
#df = pd.DataFrame({'gene1':[1,1,1], 'gene2':[0,0,0], 'gene3': [2,3,7]})
print('df', df.head(3))

# find genes that are expressed in the data set 
net_transcripts = pd.DataFrame(df.sum())
net_transcripts.columns = ['net_transcripts']
non_zeros = net_transcripts[net_transcripts['net_transcripts'] != 0] #subset data to netcounts largers than zero 
print('non_zeros', non_zeros)

## convert non_zeros variable into to correct formatting to filter later 
genes_of_interest = non_zeros.index
del non_zeros
genes_of_interest = pd.DataFrame(genes_of_interest)
genes_of_interest.columns = ['genes']
print('genes of interest', genes_of_interest)

# format original data for filtering step
df = df.T
allgenes = pd.DataFrame(df.index)
allgenes.columns = ['allgenes']
allgenes.index = df.index

df = pd.concat([allgenes,df], axis = 1)



#filter data for only cols whose sums are non zero
df = df[df['allgenes'].isin(genes_of_interest['genes'])]


## reformat data to calculate the weight

df = df.drop('allgenes', axis = 1)
df = df.T

print('df subset')
print(df)

## calculate weights on this data set 

means = df.mean()
medians = df.median()
weights = means - medians 

print(weights)
weights.to_csv('obs_weight.csv')



