
import pandas as pd
import numpy as np
import seaborn as sns
import ast



df = pd.read_csv('total_combined.csv', header = None, index_col = None)
print(df.head(4))

def unlist(nested_list):
    return [item for sublist in nested_list for item in sublist] # each item is sublist is put to the list

def process_element(x):
        # Check if x is a string, and convert it to a list if necessary
    if isinstance(x, str):
        x = ast.literal_eval(x)  # Safely evaluate the string into a list
        # Ensure x is a valid nested list
    if isinstance(x, list):
        MeanandMedian = unlist(x)  # Flatten the nested list

        Mean = MeanandMedian[0]
        Median = MeanandMedian[1]
        weight = Mean - Median
        return weight


# Apply the process_element function to each element of the DataFrame, this is will go through all the sublists in each coloumn and that the difference in the mean and
# median. creaing a new df that replaces each nested list with the weight or difference between two lists
weights_df = df.applymap(process_element)
weights_df.to_csv('correct_null_weights.csv')

