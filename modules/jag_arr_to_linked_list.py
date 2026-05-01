
# CellSparse: A Statistical Algorithm for Detecting Upregulated Genes in scRNA-seq >

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and>
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.c>



import pandas as pd
from collections import defaultdict

def reformat_jagged_array_into_linked_list(coordinates, null_counts):
    #print('genes as position:')
    #print(coordinates)
    #print('pseduo count for gene at position:')
    #print(null_counts)
    # Initialize an empty dictionary to hold the result

# Initialize result as defaultdict
    result = defaultdict(list)

# Iterate through the coordinates and null_counts together
    for coord_list, null_list in zip(coordinates, null_counts):  # iterate through sublists
        for coord, null_count in zip(coord_list, null_list):     # iterate through items within sublist
            result[coord].append(null_count)  # Directly append to the list

    #print('result before loops below', result)
    # Convert any single values to lists for uniformity, genes with only 1 pseudo count value wouldn't have the value in list format because it is just being put as = in braces and list isn't forced  
    result = {k: v if isinstance(v, list) else [v] for k, v in result.items()}
    #print('output of reformatting jagged array')
    #print(result)
    #keys_ofresult = result.keys()
    #keys_ofresult = pd.DataFrame(keys_ofresult)
    #print('amount of keys', keys_ofresult)
    return result



#### dummy to test 

# Input lists of lists
#genes = [[1, 2, 3], [1, 3], [1,2,3], [54]]
#pseduo_counts = [[4, 6, 4], [5, 7], [8,8,8], [2]]
#print('genes')
#print(genes)
#print('pseudo_counts')
#print(pseduo_counts)

#A = reformat_jagged_array_into_linked_list(genes, pseduo_counts)


# practice extracting Keys to see amount of unique keys 

#keyz = A.keys()

#keyz = pd.DataFrame(keyz)

#print('keys')

#print(keyz)

######


















