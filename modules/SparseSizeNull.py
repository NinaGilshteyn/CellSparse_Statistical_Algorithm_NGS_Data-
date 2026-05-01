

# CellSparse: A Statistical Algorithm for Detecting Upregulated Genes in scRNA-seq >

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and>
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.c>


import numpy as np
import random
def sparse_size_null(joint):
    null_list = []
    for i in joint.items():  # _ is another iterator instead of i but I am not sure if it can be replace with another value. but basically it is going through each row in this 2D object and is taking out the results and sumz. I wonder if this could be used as a dictkeys instead
        
        Array_Containing_Stats = i[1]
        #print('what enters sparse null')
        #print(Array_Containing_Stats)     
        Number_of_Activated_Genes = int(Array_Containing_Stats[1])
        Number_of_Transcripts = int(Array_Containing_Stats[0]) 
        #print('Number of active genes:', Number_of_Activated_Genes, 'number of transcripts:', Number_of_Transcripts)
        pre_null_list = [None]*1  # this will be filled each iteration then appended to the meta null_list
        array = np.zeros((1, Number_of_Activated_Genes))
        index_list = [l for l in range(Number_of_Activated_Genes)]
        for _ in range(Number_of_Transcripts):  
            #index_0list = [l for l in range(Number_of_Activated_Genes)]            
            random_index = random.choice(index_list)  
            #print("random index", random_index)
            array[0, random_index] += 1  # placing transcript into gene type, randomly. 
            #pre_null_list[0] = array # storing current distribution into memory so the next count is added to what is rememebered
            #print("pre_null_list", pre_null_list)
        #loop is done for this gene - meaning we have the pseudo-counts so now we qill append it to the main list
        #print('array is appended to nulllist', array)
        #null_list.append(pre_null_list) 
        null_list.append(array)

    # now null list is weired format 
    #print('null list in weird format', null_list)  
    null_list = [item[0].tolist() for item in null_list] # I think this extracts arrays to release nested list
    #print('what is returned by sparse null', null_list)
    return null_list


## dummydata 




