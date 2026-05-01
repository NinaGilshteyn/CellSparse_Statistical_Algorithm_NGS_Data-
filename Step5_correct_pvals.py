
#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and Upregulation in scRNA-seq >
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.com/NinaGilshteyn](https://g>


import pandas as pd 
import seaborn as sns 
import numpy as np 



df = pd.read_csv('newnewpvals_for_gene_index_as_column_index.csv', header = 0, index_col =0)

print(df.head(4))


val = df.iloc[3,:]
val = pd.DataFrame(val)
print(val.head(5))
val.columns = ['pvals']

print(type(val['pvals']))
print(val['pvals'])

print(val['pvals'].dtypes)

val['pvals'] = pd.to_numeric(val['pvals'])

print('length of pvals', len(val['pvals']))

a = val[val['pvals'] == 0.0]
print('sum of zeros', a)


b = val[val['pvals'] != 0.0]

binz = int(len(b)**0.5)
print('nonzeropval', len(b), 'binz', binz)
c = sns.displot(b, bins =binz)
c.savefig('histogram_of_nonzerop_vals.png') 


#d = b[b['pvals'] > 1]
#print(d)


import statsmodels.stats.multitest as smm
arr = smm.multipletests(val['pvals'], alpha = 0.01, method = 'fdr_bh')
print(arr)


A = arr[0]
print(A)

z = np.sum(A == True)
print('num sig after correct', z)

B = arr[1]

g = sns.displot(B, bins = 22) 
g.savefig('corrected_pvals.png')


Names = pd.DataFrame(df.columns)

A = pd.DataFrame(A)
B = pd.DataFrame(B)
C = pd.read_csv('obs_weight.csv', header = 0, index_col = 0)
print(C.head(4))
C.index = B.index
B = pd.concat([Names,A,B,C], axis=1)
B.columns = ['gene_name','significant', 'corrected_pval', 'observed weight']
B = B.sort_values(by = 'corrected_pval', ascending = True)
print(B.head(4))
B.to_csv('fulldf_corrected_pval.csv')


D = B[B['significant'] == True]
D = D.sort_values(by = 'observed weight', ascending = False)
print(D.head(4))
D.to_csv('significant_genes.csv')

 
D = B[B['significant'] == False]
D = D.sort_values(by = 'observed weight', ascending = False)
print(D.head(4))
D.to_csv('not_significant_genes.csv')



