import pandas as pd
import numpy as np
import seaborn as sns 
import ast


weights_df = pd.read_csv('correct_null_weights.csv', header = 0, index_col = 0)
## we had to repeat it in shell and I talked myself into makeing a seperate script for t

print('null weights')
print(weights_df.head(3))#
print('null weights')
#print(weights_df.head(3))

weights_df.drop('0', axis = 1)



def unlist(nested_list):
    return [item for sublist in nested_list for item in sublist] # each item is sublist is put to the list

print('null weights')
print(weights_df.head(3))


# 1st should have p val of 0.5  and second zero, and third should be 1 since they are all larger)
obs = pd.read_csv('obs_weight.csv', header = 0, index_col = None)
obs.columns = ['geneName', 'weight']

A = obs['geneName'].values.tolist()
B = obs['weight'].values.tolist()

print(A)
print(B)


keys = list(obs.index)
print(keys)


ziped = zip(keys,A,B)
dictt = {keys:[A,B] for keys,A,B in ziped}
print('dict')
print(dictt)




# i need to covert instead the obs into a dictionary where there keys correspond to the column index of the null df. Then the item can be the count and the gene name so 
# then I can can extract the item that is the obs value for the pval calculation and the gene name to to then link it and I can just cvqrgb

pvals = {}

sim_total = len(weights_df)
for i in weights_df.columns:
    print('i', i)
    null_values = weights_df[i]
    #null_values = weights_df
    print('null values')
    print(null_values) 
    gene_index = int(i)
    obsarr = dictt.get(gene_index)
    print('obs_ars', obsarr)
    obs_val= obsarr[1]
    print('obs_val', obs_val)        
    num1 = np.sum(null_values >= obs_val)
    num2 = np.sum(null_values <= -obs_val)
    p_val = (num1 + num2) / (2*sim_total) 
    print('n1:', num1, 'n2:', num2, 'p', p_val)
    #print('dict key is gene index')
    pvals[obsarr[0]] = ['observed weight:', obs_val, 'two-sided pval', p_val]
   

print(pvals)
print(pvals)
pvals = pd.DataFrame(pvals)
print(pvals) 

print('len pvals in dataframe', len(pvals))
pvals.to_csv('newnewpvals_for_gene_index_as_column_index.csv')
