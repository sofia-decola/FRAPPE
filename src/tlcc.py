"""
Time Lagged Cross Correlation (TLCC)
"""
import pandas as pd
import numpy as np
from scipy import stats
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import math

__all__ = [
    "tlcc",
    "tlcc_datasets",
    "clean_protein_data",
    "pears_corr",
    "spman_corr",
    "wrap_lag",
    "plot_tlcc_histogram",
    "get_indiv_protein_tlcc",
]
 
"""
Time Lagged Cross Correlation: A set of time series data analysis tools that compares correlation at different time lags for two datasets of proteomics time series data,
with the two data sets representing proteomics data in different experimental groups

Input Data Set Required Setup: 
- any .csv file of data that includes a column of protein names (default label is 'Gene Name') as well as columns for each time point.
- the protein_col name must be consistent across datasets

This colelction of functions can run TLCC on individual proteins or whole datasets. The goal is to analyze the optimal time lag distribution across two datasets, 
one from a control condition and one from an experimental condition
"""


def tlcc(protein, cols, exper_path, contr_path, corr=0, lag=1):
    """
    Time Lagged Cross Correlation (TLCC) Function - comparing the correlation of two proteins over time lags
    ----------------
    Parameters: 
    protein: (string) name of protein to run TLCC on
    cols: (array of strings) columns to be dropped, any column from the data set that isn't the protein name or the time points
    exper_path: (string) path to a .csv file with time series data for the experimental condition, does not need to include '.csv'
    contr_path: (string) path to a .csv file with time series data for the control condition, does not need to include '.csv'
    corr: (int) indicate whether the code uses Spearman (corr=0) or Pearson (corr=1) correlation in the TLCC calculations
    default value corr=0, meaning Spearman correlation is used
    lag: (int) number of hours between timepoints in the data sets
    default value lag=1, meaning time points are 1 hour apart
    ----------------
    Return: array of protein name, correlation of exp group data to control group data as is, the max correlation once the time lags are introduced,
    the associated time lag that caused the max correlation (in hours), and the p value for the max corr statistic
    """
    exp_df = clean_protein_data(protein, exper_path, cols)
    control_df = clean_protein_data(protein, contr_path, cols)
    exp_y = exp_df.iloc[0].values
    contr_y = control_df.iloc[0].values
    if corr == 0:
        output = spman_corr(contr_y, exp_y, protein, lag)
    elif corr == 1:
        output = pears_corr(contr_y, exp_y, protein, lag)
    return output
    

def tlcc_datasets(exp_data, control_data, cols, output_csv, corr=0, lag=1, protein_col='Gene Name'):
   """
    Runs TLCC across two datasets of proteins in different experimental groups, a control against an experimental group. 
    Finds common proteins between data set and runs TLCC on each protein, comparing the experimental time 
    --------------
    Parameters
    exp_data: (string) path to time series .csv with data from the experimental group condition, need to include '.csv'
    control_data: (string) path to time series .csv with data from the control group condition, need to include '.csv'
    cols: (array of strings) columns to be dropped, any column from the data set that isn't the protein name or the time points
    output_csv: (string) name of the desired output file, does not require the .csv ending
    corr: (int) indicate whether the code uses Spearman (corr=0) or Pearson (corr=1) correlation in the TLCC calculations
    default value corr=0, meaning Spearman correlation is used
    lag: (int) number of hours between timepoints in the data sets
    default value lag = 1, meaning time points are 1 hour apart
    protein_col: (string) name of the column with all the protein names in exp_data and control_data
    --------------
    Result: output .csv where each row contains the protein name, correlation, max correlation, and p value for max correlation 
    for every protein in both datasets
    """
   exp_df = pd.read_csv(exp_data)
   control_df = pd.read_csv(control_data)
   common_proteins = set(exp_df[protein_col]).intersection(set(control_df[protein_col]))
   results_df = pd.DataFrame(columns=["Protein Name", "Corr", "Max Corr", "Time Lag For Max Corr", "P-Value Max Corr"])
   print("Running...")
   with alive_bar(len(common_proteins)) as bar: 
        for protein in common_proteins:
            results_df.loc[len(results_df.index)] = (tlcc(protein, cols, exp_data, control_data, corr, lag))
            bar()
   results_df.to_csv(f"{output_csv}.csv")
   print("Done!")
   


def clean_protein_data(protein, path, cols, protein_col='Gene Name'):
    """
    Cleans the dataframe to get the necessary values for Time Lagged Cross Correlation, removing excess columns, helper function for TLCC
    Only gets time points for one protein, drops all other columns based on specified parameter
    --------------
    Parameters: 
    protein: (string) name of the protein being acquired 
    path: (string) path to data file
    cols: (array of strings) columns to be dropped, any column from the data set that isn't the gene name or the time points
    --------------
    Returns: 
    row of data from a pandas DataFrame, containing time points for TLCC analysis
    """
    df = pd.read_csv(path)
    data_row = (df.loc[df[protein_col] == protein])
    data_row = data_row.drop(columns=cols) 
    data_row = data_row.loc[data_row[protein_col] == protein].drop(protein_col, axis=1)
    return data_row


def pears_corr(contr, exp, protein, lag=1):
    """
    Helper function for TLCC, calculates Pearson correlation (using SciPy stats) for experimental condition versus control condition time points, 
    rolling the experimental time points, returning the max correlation of the time points and the time lag (in hours) associated with the max correlation
    -------------
    Parameters:
    contr: (string) path to a .csv file with time series data for the control condition 
    exp: (string) path to a .csv file with time series data for the experimental condition
    protein: (string) name of the protein running the correlation value on
    lag: (int) number of hours between timepoints in the data sets
    default value lag=1, meaning time points are 1 hour apart
    -------------
    Returns: array of protein name, correlation of exp group data to control group data as is, the max correlation once the time lags are introduced,
    the optimal time lag that caused the max correlation (in hours), and the p value for the max corr statistic
    """
    nan = False
    if np.std(contr) == 0.0 or np.std(exp) == 0.0:
        corr = np.nan
        nan = math.isnan(corr)
    else:
        corr  = stats.pearsonr(contr, exp)
    max_corr = corr
    optimal_lag = 0
    if (nan == False):
        for i in range(1, len(exp) - 1):       
            exp = np.roll(exp, -1)  
            roll_corr = stats.pearsonr(contr, exp)
            if roll_corr.statistic >= max_corr.statistic: # type: ignore
                max_corr = roll_corr
                optimal_lag = i

    if nan:
        return [protein, corr, max_corr, optimal_lag * lag, math.nan] 
    else:
        return [protein, corr.statistic, max_corr.statistic, optimal_lag * lag, max_corr.pvalue] # type: ignore


def spman_corr(contr, exp, protein, lag=1):
    """
    Helper function for TLCC, calculates Spearman correlation (using SciPy stats) for experimental condition versus control condition time points, 
    rolling the experimental time points, returning the max correlation of the time points and the time lag (in hours) associated with the max correlation
    -------------
    Parameters:
    contr: (string) path to a .csv file with time series data for the control condition, does not need to include '.csv' 
    exp: (string) path to a .csv file with time series data for the experimental condition, does not need to include '.csv'
    protein: (string) name of the protein running the correlation value on
    lag: (int) number of hours between timepoints in the data sets
    default value lag=1, meaning time points are 1 hour apart
    -------------
    Returns: array of protein name, correlation of exp group data to control group data as is, the max correlation once the time lags are introduced,
    the optimal time lag that caused the max correlation (in hours), and the p value for the max corr statistic
    """
    
    nan = False
    if np.std(contr) == 0.0 or np.std(exp) == 0.0:
        corr = np.nan
        nan = math.isnan(corr)
    else:
        corr  = stats.spearmanr(contr, exp)
    max_corr = corr
    optimal_lag = 0
    if (nan == False):
        for i in range(1, len(exp) - 1):       
            exp = np.roll(exp, -1)  
            roll_corr = stats.spearmanr(contr, exp)
            if roll_corr.statistic >= max_corr.statistic: # type: ignore
                max_corr = roll_corr
                optimal_lag = i
    if nan:
        return [protein, corr, max_corr, optimal_lag * lag, math.nan] 
    else:
        return [protein, corr.statistic, max_corr.statistic, optimal_lag * lag, max_corr.pvalue] # type: ignore

def wrap_lag(x):
    return x if x <= 12 else x - 24


def plot_tlcc_histogram(data_path, contr_condition="Control", exp_condition="Experimental", lag=1, male_path=None, female_path=None, color='skyblue'):
    """
    Plots the histogram of the distribution of optimal time lags between two data sets, either for one data set or two
    --------------
    Parameters:
    data_path: (string) path to TLCC output file
    contr_condition: (string) indicating the control condition (for plot title)
    default value = "Control"
    exp_condition: (string) indicating the experimental condition (for plot title)
    default value = "Experimental"
    lag: (int) number of hours between timepoints in the data sets
    default value lag = 1, meaning time points are 1 hour apart
    -------------- 
    """
    if male_path is None or female_path is None:
        df = pd.read_csv(data_path)
        wrapped_lags = df["Time Lag For Max Corr"].apply(wrap_lag)

        bins = np.arange(-12.5, 13.5, lag)
        bin_centers = np.arange(-9, 12 + lag, lag)

        bins = np.concatenate([
                bin_centers - lag / 2,
                [bin_centers[-1] + lag / 2 + 1e-6]  # ensure +12 is included
            ])

        __, ax = plt.subplots()
        ax.hist(
            x=wrapped_lags,
            bins=bins, #type: ignore
            rwidth=0.85,
            color=color,
            edgecolor='black'
        )
        ax.set_xlabel("Time Lag for Max Correlation (in hours)")
        ax.set_ylabel("Number of Proteins")
        ax.set_title(f"Time-Lagged Cross-Correlation: {contr_condition} vs. {exp_condition}")
        
        ax.set_xticks(bin_centers)
        ax.axvline(0, color='black', linewidth=1)  # visually centers zero

        plt.show()
        return

    male_df = pd.read_csv(male_path)
    female_df = pd.read_csv(female_path)
    male_lags = male_df["Time Lag For Max Corr"].apply(wrap_lag)
    female_lags = female_df["Time Lag For Max Corr"].apply(wrap_lag)

    bin_centers = np.arange(-9, 12 + lag, lag)

    bins = np.concatenate([
        bin_centers - lag / 2,
        [bin_centers[-1] + lag / 2 + 1e-6]  # ensure +12 is included
    ])


    __, ax = plt.subplots()

    ax.hist(
        [male_lags, female_lags],
        bins=bins, #type: ignore
        rwidth=0.85,
        label=["Male", "Female"],
        histtype='bar',
        edgecolor='black',
        color=["black", "white"]
    )

    ax.set_xlabel("Time Lag for Max Correlation (in hours)")
    ax.set_ylabel("Number of Proteins")
    ax.set_title(f"Time-Lagged Cross-Correlation: {contr_condition} vs. {exp_condition}")
    ax.set_xticks(bin_centers)
    ax.axvline(0, color='black', linewidth=1)
    ax.legend()
    plt.show()

def get_indiv_protein_tlcc(contr_df, exp_df, protein):
    contr = contr_df.loc[contr_df["Protein Name"] == protein]
    exp = exp_df.loc[exp_df["Protein Name"] == protein]
    contr = contr.set_index('Protein Name')
    exp = exp.set_index('Protein Name')
    print(contr)
    print(exp)
   
