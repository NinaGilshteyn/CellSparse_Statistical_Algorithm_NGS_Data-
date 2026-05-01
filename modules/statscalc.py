

# CellSparse: A Statistical Algorithm for Detecting Upregulated Genes in scRNA-seq >

#**Author:** Nina Gilshteyn, M.S.
#**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology
#**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and>
#**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.c>




import pandas as pd
import numpy as np



def calc_stat(number_of_cells, dictionary):# cell number is the number of cells in the sample
    Means = []
    Medians = []
   
    #print(number_of_cells)
    #print(dictionary)
    value = dictionary[1]
    length_list_of_genes =  len(value)
    #print('Pseudo_counts')
    #print(value)
    cell_number = number_of_cells 
    length_50 = cell_number/2
    if len(value) < length_50:
        med = 0
        #print('median is zero')
        mean = np.sum(value)/cell_number
        Means.append(mean)
        Medians.append(med)    
    else:
        #cell_number =  5
        empty_cells = (cell_number - int(len(value)))
        #print('number empty_cells')
        zeroz = np.zeros(empty_cells)
        zeroz = list(zeroz)
        #print('zeros', zeroz)
        #value = np.array(value)
        population = value + zeroz
        population = np.array(population)
        population = list(population)
        med = np.median(population) # this will automatically sort the values and find the middle value or average the middle values for a list with an even amount of elements.  
        #print('population', population)  
        #print('median', med)
        Medians.append(med)
        #print('null counts inside stats function before mean calc:', value)
        mean = np.mean(population)
        Means.append(mean)
        
    stats_pair = [Means,Medians]    
    return stats_pair

     

##### dummy
#result = {1: [1,4, 5, 8], 2: [2,6, 8, 8], 3: [3,4, 7], 4: [4,2]}
#a = calc_stat(result,6)
#print('TADA', a)
