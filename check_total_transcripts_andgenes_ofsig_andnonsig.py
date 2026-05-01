import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 

small_sig = pd.read_csv('redosmall_significantweight.csv', header = 0, index_col =0)

non_sig = pd.read_csv('not_significant_genes.csv', header = 0, index_col = 0)

#print('small_sig')
#print(small_sig.head(3))
#print('non_sig')
#print(non_sig.head(3))


#fulldf = pd.read_csv('Full.csv', header = 0, index_col = 0)
#fulldf = fulldf.drop('Cell Type', axis = 1)
#fulldf = fulldf.T

#print('full df')
#print(fulldf.head(3))


# get gene subsets 

#sig_gene = small_sig['gene_name']
#non_sig_gene = non_sig['gene_name']

#sig_data = fulldf.loc[sig_gene]
#print('sig_gene')
#print(sig_gene.head(3))
#print('sig_data')
#print(sig_data.head(3))

#sig_data.to_csv('sigdata.csv')

#non_sig_data = fulldf.loc[non_sig_gene]
#print('sig_gene')
#print(sig_gene.head(3))
#print('sig_data')
#print(sig_data.head(3))

#non_sig_data.to_csv('non_sigdata.csv')

##### read in the files now

nsdf = pd.read_csv('non_sigdata.csv', header =0, index_col = 0)
sdf = pd.read_csv('sigdata.csv', header = 0, index_col = 0) 

print('nsdf')
print(nsdf.head(3))
print('sdf')
print(sdf.head(3))

# put sums from each group on 1 df 
#A = sns.displot(nsw, bins = binz, color = 'firebrick', edgecolor = 'black', log_scale = (False, True))
#plt.axvline(med, color = 'darkorange', label = f'Median:{med}')





ns_transcripts = np.sum(nsdf, axis = 0)
ns_activegenes = np.sum(nsdf!=0, axis = 0)
plt.scatter(ns_activegenes, ns_transcripts)
plt.show

s_transcripts = np.sum(sdf, axis = 0)
s_activegenes = np.sum(sdf!=0, axis = 0)



#### here for the histograms 
ns_transcripts = np.sum(nsdf, axis = 0)
print('ns transcripts')
print(ns_transcripts)

Binz = int(len(ns_transcripts)**0.5)
A = sns.displot(ns_transcripts, bins = Binz, color = 'firebrick', edgecolor = 'black')
mean = np.mean(ns_transcripts)
plt.axvline(mean, color = 'dodgerblue', label = f'Mean:{mean}')
med = np.median(ns_transcripts)
plt.axvline(med, color = 'darkorange', label = f'Median:{med}')
plt.title('Net_transcripts')
plt.suptitle('Hydra_nonsig_cellsparse')
plt.legend()
A.savefig('non_significant_netranscripts.png') 

s_transcripts = np.sum(sdf, axis =0)
Binz = int(len(s_transcripts)**0.5)
A = sns.displot(s_transcripts, bins = Binz, color = 'darkgreen', edgecolor = 'black')
mean = np.mean(s_transcripts)
plt.axvline(mean, color = 'dodgerblue', label = f'Mean:{mean}')
med = np.median(s_transcripts)
plt.axvline(med, color = 'darkorange', label = f'Median:{med}')
plt.title('Net_transcripts')
plt.suptitle('Hydra_sig_cellsparse')
plt.legend()

A.savefig('significant_netranscripts.png')




#put number of genes in group on 1 df 






ns_transcripts = np.sum(nsdf!=0, axis = 0)
print('ns transcripts')
print(ns_transcripts)

Binz = int(len(ns_transcripts)**0.5)
A = sns.displot(ns_transcripts, bins = Binz, color = 'firebrick', edgecolor = 'black')
mean = np.mean(ns_transcripts)
plt.axvline(mean, color = 'dodgerblue', label = f'Mean:{mean}')
med = np.median(ns_transcripts)
plt.axvline(med, color = 'darkorange', label = f'Median:{med}')
plt.title('Net_activegenes')
plt.suptitle('Hydra_nonsig_cellsparse')
plt.legend()
A.savefig('non_significant_netactivegenes.png')

s_transcripts = np.sum(sdf!=0, axis =0)
Binz = int(len(s_transcripts)**0.5)
A = sns.displot(s_transcripts, bins = Binz, color = 'darkgreen', edgecolor = 'black')
mean = np.mean(s_transcripts)
plt.axvline(mean, color = 'dodgerblue', label = f'Mean:{mean}')
med = np.median(s_transcripts)
plt.axvline(med, color = 'darkorange', label = f'Median:{med}')
plt.title('Net_Active genes')
plt.suptitle('Hydra_sig_cellsparse')
plt.legend()

A.savefig('significant_netactivegenes.png')





# repeast his with output steps from earlier in the algorith m

