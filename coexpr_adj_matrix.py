# Comparing Adjacency Matrices

import numpy as np
import pandas as pd


__all__ = [
    "corr_input",
    "create_corr_matrix",
    "corr_to_adj",
    "flatten_matrix",
    "sub_matrices",
    "add_matrices",
    "retrieve_occurences",
    "combine_dict",
    "freq_dict",
    "add_freq_occurences",
    "compute_perc",
    "comp_matrices",
    "pairwise_comp",
]
 
 
def corr_input (file):
    """
    Finds the correlation of the protein
    against itself in a bi-condition comparison
    -------------
    Parameters:
    file: csv of correlation between same
    protein of two groups
    -------------
    Returns: corr (int): the correlation values
    """

    corr = pd.read_csv(file)
    corr = corr[['Corr']]

    return corr
def create_corr_matrix (data_frame, file):
    """
    Creates a correlation matrix for a dataframe, for when the
    correlation is the protein against itself, the correlation from the
    file inserts the diagonal.
    -------------
    Parameters:
    dataframe: array of integers representing correlation
    file: csv file including correlation for protein against itself
    -------------
    Return: corr_matrix (array): an array where protein correlation
    values are compared to itself and other protein correlation
    """

    corr = corr_input(file)
    corr_matrix = data_frame.T.corr(method = 'pearson')
    np.fill_diagonal(corr_matrix.values, corr)

    return corr_matrix
def corr_to_adj (corr_matrix, threshold):
  """
  Creates a correlation matrix based on the given dataframe
  using Pearson Correlation coeffecient.
  -------------
  Parameters:
  corr_matrix: (array) array of integers representing correlation
  threshold: (int) threshold to pass/not pass rhymicity in both conditions
  ------------
  Return: adj_matrix (array): array of 1's where corr value is > then the
  threshold, and 0's where corr value is < then the threshold
  """

  adj_matrix = corr_matrix.map(lambda x: 1 if abs(x) > threshold else 0)
  return adj_matrix

def flatten_matrix (matrix):
    """
    Flattens a matrix into a numpy list
    -------------
    Parameters:
    matrix1: (array) adjacency matrix for one group condition
    -------------
    Returns: (array) a numpy 1D array of the upper triangle elements
    """

    return matrix[np.triu_indices_from(matrix, k=1)]

def sub_matrices(adj_matrix1, adj_matrix2):
    """
    Subtracts an adjacency matrix from another.
    -------------
    Parameters:
    adj_matrix1: (array) adjacency matrix for one group condition
    adj_matrix2: (array) adjacency matrix for one group condition
    ------------
    Returns: array: the resulting adjacency matrix after
    subtracting one from another
    """

    return (adj_matrix1 - adj_matrix2).to_numpy()
def add_matrices(adj_matrix1, adj_matrix2):
    """
    Adds an adjacency matrix to another.
    -------------
    Parameters:
    adj_matrix1: (array) adjacency matrix for one group condition
    adj_matrix2: (array) adjacency matrix for one group condition
    -------------
    Returns: array: the resulting adjacency matrix after
    adding one to another
    """
    return (adj_matrix1 + adj_matrix2).to_numpy()
def retrieve_occurences (flat_matrix):
    """
    Makes a dictionary of the occurences of values within a list.
    -------------
    Parameters:
    flat_matrix: (array) adjacency matrix for one group condition
    ------------
    Returns: dictionary: key = int value, value = count of how often key
    appears in list
    """

    unique, counts = np.unique(flat_matrix, return_counts=True)
    return dict(zip(unique, counts))
def combine_dict (dict1, dict2):
    """
    Combines two dictionaries.
    ------------
    Parameters:
    dict1 :(dictionary) 1st dictionary
    dict2 :(dictionary) 2nd dictionary
    ------------
    Returns: dictionary: a combined dictionary of dict1 and dict2
    """
    return dict1 | dict2
def freq_dict (freq_dict, num1, num2):
    """
    Makes a dictionary of occurences for two specific
    value of the integers given. Num1 and Num2 should not
    be the same integer.
    ------------
    Parameters:
    freq_dict :(dictionary) adjacency matrix for one group condition
    num1 :(int) occurence of number given
    num2 :(int) occurence of number given
    ------------
    Returns: dictionary: key = int value, value = count of how often key
    appears in list
    """

    dict = {int(k): int(v) for k, v in freq_dict.items()}
    dict = {k: v for k, v in freq_dict.items() if k in [num1,num2]}

    return dict

def add_freq_occurences(freq_dict):
    """Adds the amount of occurences for each key value within the
    dictionary.
    -------------
    Parameters:
        freq_dict (dcitionary): adjacency matrix for one group condition
    -------------
    Returns: total: (int) count of how many values in a key value pair are in the
    dictionary

    """
    total = 0
    for key in freq_dict:
        total += freq_dict[key]

    return total

def compute_perc(count,total):
    """
    Makes a dictionary of the occurences of values within a list
    -----
    Parameters:
    matrix1: (array) adjacency matrix for one group condition
    -----
    Returns: freq: (dictionary) key = int value, value = count of how often key
    appears in list
    """

    freq = (count/total)*100

    return freq
def comp_matrices(adj_matrix1, adj_matrix2):
    """
    Compares two adjacency matrices. Analyzes how proteins
    are correlated from one group to another.
    --------------
    Parameters:
    adj_matrix1: (array) adjacency matrix for one group condition
    adj_matrix2: (array) adjacency matrix for a different group condition
    --------------
    Prints: Outcome of comparing matricy (adding and subtracting) in frequency int
    """

    #### creating dictionaries and finding counts
    # 1 = 1-->0, -1 = 0-->1
    # 2 = 1-->1, 0 = 0--> 0

    sub = sub_matrices(adj_matrix1, adj_matrix2)
    sub_flat = flatten_matrix(sub)
    sub_dict = retrieve_occurences(sub_flat)
    sub_dict = freq_dict(sub_dict,1,-1)

    add = add_matrices(adj_matrix1, adj_matrix2)
    add_flat = flatten_matrix(add)
    add_dict = retrieve_occurences(add_flat)
    add_dict = freq_dict(add_dict,0,2)


    comb_dict = combine_dict(sub_dict, add_dict)
    total = add_freq_occurences(comb_dict)

    occur_2 = comb_dict[2]
    occur_0 = comb_dict[0]
    occur_1 = comb_dict[1]
    occur_m1 = comb_dict[-1]

    freq_2  = compute_perc(occur_2, total)
    freq_0  = compute_perc(occur_0, total)
    freq_1  = compute_perc(occur_1, total)
    freq_m1 = compute_perc(occur_m1, total)

    print("Frequency of 1-->1:", f"{freq_2:.2f}")
    print("Frequency of 0-->0:", f"{freq_0:.2f}")
    print("Frequency of 0-->1:", f"{freq_1:.2f}")
    print("Frequency of 1-->0:", f"{freq_m1:.2f}")


def pairwise_comp(control, experimental, file1, file2, intersect_data, threshold=0.8):

    file1 = pd.read_excel(file1)
    file2 = pd.read_excel(file2)

    df1 = file1[[col for col in file1 if col.startswith('Fitted')]] #type: ignore
    df2 = file2[[col for col in file2 if col.startswith('Fitted')]] #type: ignore

    matrix1 = create_corr_matrix(df1,intersect_data) # Control vs. Experimental w/Intersect
    matrix2 = create_corr_matrix(df2,intersect_data) # Control vs. Experimental w/Intersect

    adjacency_matrix1 = corr_to_adj(matrix1, threshold) #Control
    adjacency_matrix2 = corr_to_adj(matrix2, threshold) #Experimental

    print(control,"VS", experimental)
    comp_matrices(adjacency_matrix2, adjacency_matrix1)

