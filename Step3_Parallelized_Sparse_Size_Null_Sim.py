import pandas as pd
import numpy as np
import multiprocessing
import random
from functools import partial
import sys
import os
import gzip
import csv
import time
from itertools import chain
import argparse 

from jag_arr_to_linked_list import reformat_jagged_array_into_linked_list
from statscalc import calc_stat
from compressedwriter import compressed_write
from SparseSizeNull import sparse_size_null


def divide_iteration_across_cpus(chunkz, num_cpus, iterator, coordinates_of_Genes, numberOfCells):
    ### the loop above diving the objects into chunks that are processing it into the parallelizing tool 
  
    pool = multiprocessing.Pool(num_cpus)
    resultss = pool.map(sparse_size_null, chunkz)  #results for 1 simulation
    pool.close()
    pool.join()
    
    #### ok so basically each process adds another bracket or nests each sublist so we need to extract the values 
    #result = []
    #for sub_nest_list in resultss:
    #    for element in sub_nest_list:
    #        result.append(element)
    # Because each process add another bracket to the chunk of data handled by the process.  so we need to flatten it and this code apparently does it faster 
    ## apparently this function is faster for larger data sets because it avoids creating intermediate lists during the flattening 
    result = list(chain.from_iterable(resultss))
    #print('original form of result')
    #print(resultss)
    #print('final format of result from sparsenull')
    #print(result)  
    del resultss


    from multiprocessing import Pool
    from functools import partial
    
    if __name__ == "__main__":

      ## turning stacked pseudo_counts into a dictionary 
        pseduo_counts_linked_list = reformat_jagged_array_into_linked_list(coordinates_of_Genes, result)
        number_of_cells_in_whole_data  = numberOfCells
        #print('number of cells within divide iteration', number_of_cells_in_whole_data)
        #print('coordinates', coordinates_of_Genes)
        pseduo_counts_linked_list = dict(sorted(pseduo_counts_linked_list.items()))  # this is so that the indicides are done in order so the output stats are in order of the genes are in order also 
        #print('input to pool', pseduo_counts_linked_list)
        #print('len of dict', len(pseduo_counts_linked_list))
        partial_function = partial(calc_stat, number_of_cells_in_whole_data)
        with Pool() as pool:
            results = pool.map(partial_function, pseduo_counts_linked_list.items())

   
    return results     


# this function creates an interator that tracks of the amount of simulations 
# for each gene and so within 1 of the iterations it is going to feed into the
# paralellizing machine called divide_iteration_across_cpus


def divide_and_save_iterations(joints, num_cpus, num_iterations, coordinates, numberofCells, simulation_points_for_disk_writing):
    null_stats = []
    disk_writing_points_set = set(simulation_points_for_disk_writing)  # Use a set for faster membership checks
    parser = argparse.ArgumentParser(description = 'input number simulation')
    parser.add_argument('number', help = 'number of simulations')
    args = parser.parse_args()

    for citerator in range(num_iterations):
        # Perform the simulation for the current iteration
        stats_for_1_simulation = divide_iteration_across_cpus(joints, num_cpus, citerator, coordinates, numberofCells)
        null_stats.append(stats_for_1_simulation)
        #parser = argparse.ArgumentParser(description = 'input number simulation')
        #parser.add_argument('number', help = 'number of simulations')
        #args = parser.parse_args()


        # If the iteration is in disk_writing_points, write to disk
        if citerator in disk_writing_points_set:
            statsdf = pd.DataFrame(null_stats)
            #print(statsdf)
            jobnumber = args.number
            statsdf.to_csv(f'jobnumber_{jobnumber}_simulationnumber_{citerator}_null_stats.csv', header=False, index=False)
            null_stats = []  # Reset the null_stats for the next batch
      

# Here we are reading in the files from previous steps and reset index so we can join the files together                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
Number_Activated_Genes = pd.read_csv("AmountExpressedGenes.csv", index_col=0, header=0)  #this is the amount of expressed genes or the length of the array we are depositing transcripts into
Number_Transcripts = pd.read_csv("Total_RNA_per_Cell.csv", index_col=0, header=0)  # is the total amount of RNA in each cell
Number_Transcripts = Number_Transcripts.reset_index(drop=True)
Number_Activated_Genes = Number_Activated_Genes.reset_index(drop=True) 
coordinates = 'Coordinates_of_Expressed_Genes.csv'
coordinatesdf = pd.read_csv(coordinates, header = 0, index_col = 0)
df_as_list = coordinatesdf.values.tolist()
del coordinates
del coordinatesdf
# Remove NaN values from the list
corr_list = [[value for value in row if not pd.isna(value)] for row in df_as_list] # goes through each sublist and 
corr_list = [[int(value) for value in sublist] for sublist in corr_list]
del df_as_list 

## this is so that is executes in order or not automatically 
if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description = 'input number simulation')
    #parser.add_argument('number', help = 'number of simulations')
    #args = parser.parse_args()
 

    num_cpus = 22  # Number of CPUs to use to get a single null value for each gene 
    num_iterations = 200  # equivalent to the number of null values for each gene 
    joint_df = pd.concat([Number_Transcripts, Number_Activated_Genes], axis=1)
    joint_df.columns = ['Number_Transcripts', 'Number_Activated_Genes']
    #print('joint_df', joint_df)
    print('num_iterations:', num_iterations)
    number_of_cells = len(joint_df)     

    chunk_size = len(joint_df) // num_cpus   # this calculates the amount of cells each CPU should recieve
    print('num cpus', num_cpus)

    chunks = [joint_df[i:i+chunk_size] for i in range(0, len(joint_df), chunk_size)]   
    ##  i changes by increments of chunk size and puts a comma to slice the data so that the slices can be fed into the divide_and_save_iterations function
    list_of_dicts = [chunk.T.to_dict(orient='list') for chunk in chunks]
    del chunks
    when_to_write = []
    breaks = 4
    for i in range((breaks - 1),(num_iterations -1), breaks): #breaks -1 is the timepoint I want to 1st write and iterations are numbered starting with 0 so I need to deduct 
        when_to_write.append(i)  
    g = int(num_iterations - 1) # the last time point is excluded in the loop so I need to manually add it because the code only right when there is a match
    when_to_write.append(g)
    #print(list_of_dicts) 
    divide_and_save_iterations(list_of_dicts, num_cpus, num_iterations, corr_list, number_of_cells, when_to_write) 




