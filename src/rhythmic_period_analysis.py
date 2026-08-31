# Rhythmic Period Analysis (RPA)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import anderson

__all__ = ["rpa", "plot_period_delta"]

def rpa(data1, data2, condition1, condition2, upper_bound, lower_bound, plot_color="skyblue", output="output", protein_col="Gene Name", period_col="Period"):
    """
    Generalized function to take in two datasets and compare the change in period in proteins as a result of the presence of experimental stimuli
    --------------
    Parameters:
    data1: (string) first data set inputted as .csv file path
    data2: (string) second data set inputted as .csv file path
    condition1: (string) a string delineating the condition correlating to data1 (i.e. "Control Group)
    condition2: (string) a string delineating the condition correlating to data2 (i.e. "Experimental Group)
    upper_bound: (int) the upper bound of hours for a normal period. In circadian biology the period range is 18-26 hours,
    so 26 hours would be the upper bound
    lower_bound: (int) the lower bound of hours for a normal period. 18 hours in circiadian biology
    output: (string) the string for the output file of the rhythmic period analysis, do not include ".csv"
    protein_col: (string) name of the column with all the protein names in data1 and data2, case sensitive.
    period_col: (string) name of the column with period length values in data1 and data2, case sensitive.
    --------------
    Function will save a csv file of the filtered dataframe, and produce a histogram of delta period when run
    """
    df1 = pd.read_csv(data1)
    df1 = df1[[protein_col, period_col]]
    df2 = pd.read_csv(data2)
    df2 = df2[[protein_col, period_col]]
    df3 = df1.merge(df2, left_on=protein_col, right_on=protein_col, suffixes=(f"_{condition1}", f"_{condition2}"), how="inner")
    df3[f"{condition1} Within Bound"] = ((df3[f"{period_col}_{condition1}"] >= lower_bound) & (df3[f"{period_col}_{condition1}"] <= upper_bound))
    df3[f"{condition2} Within Bound"] = ((df3[f"{period_col}_{condition2}"] >= lower_bound) & (df3[f"{period_col}_{condition2}"] <= upper_bound))
    print(df3)
    filtered_df = df3[(df3[f"{condition1} Within Bound"]) | (df3[f"{condition2} Within Bound"])].copy()
    filtered_df["Change in Period"] = (
        filtered_df[f"{period_col}_{condition2}"] -
        filtered_df[f"{period_col}_{condition1}"]
    )    
    filtered_df.to_csv(f"{output}.csv")
    print(filtered_df.head())
    print(filtered_df.dtypes)
    print(filtered_df["Change in Period"].isna().sum())
    plot_period_delta(filtered_df, condition1, condition2, plot_color)


# Plotting Code
def plot_period_delta(filtered_df, condition1, condition2, plot_color="skyblue"):
    """
    Plotting function to show the distribution of changes in length of period across an inputted proteome data set
    Automatically called when rpa or circadian rpa are called
    --------------
    Parameters:
    filtered_df: (DataFrame) created in rpa/circadian_rpa functions, contains data on changes in period
    condition1: (string) string denotion of control condition
    condition2: (string) string denotion of experimental condition
    --------------
    Produces histogram of changes in period across proteome
    """
    
    fig, ax = plt.subplots(figsize=(11, 6))

    min_val = filtered_df["Change in Period"].min()
    max_val = filtered_df["Change in Period"].max()
    bins = np.arange(np.floor(min_val/3)*3, np.ceil(max_val/3)*3 + 3, 3)
    ax.hist(
        filtered_df["Change in Period"],
        bins=bins,
        color=plot_color,
        edgecolor='black',
        align='mid'
    )
    anderson_val = anderson(filtered_df["Change in Period"])
    print(anderson_val)
    mean_val = filtered_df["Change in Period"].mean()

    # Anderson-Darling normality/skew check at the 5% and 1% significance levels
    sig_levels = list(anderson_val.significance_level) #type: ignore
    ad_stat = anderson_val.statistic  #type: ignore
    idx_5 = sig_levels.index(5.0)
    idx_1 = sig_levels.index(1.0)
    crit_5 = anderson_val.critical_values[idx_5]  #type: ignore
    crit_1 = anderson_val.critical_values[idx_1]  #type: ignore

    verdict_5 = "reject" if ad_stat > crit_5 else "fail to reject"
    verdict_1 = "reject" if ad_stat > crit_1 else "fail to reject"

    ad_label = (
        f"Anderson-Darling stat={ad_stat:.2f}\n"
        f"5% crit={crit_5:.2f} ({verdict_5} normality)\n"
        f"1% crit={crit_1:.2f} ({verdict_1} normality)"
    )

    ax.set_title(f"Distribution of Change in Period ({condition2} - {condition1})")
    ax.set_xlabel("Change of Period (hours)")
    ax.set_ylabel("Number of Proteins")
    ax.axvline(0, color='red', linestyle='--', label='No Change')
    ax.axvline(mean_val, color='black', linestyle=':', label=f'Mean ({mean_val:.2f})')
    ax.plot([], [], ' ', label=ad_label)  # invisible handle just to add text to legend
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, handlelength=1.5, borderpad=0.6)
    ax.grid(axis='y', linestyle='--', alpha=0.6)

    # Explicitly set axes position so the histogram takes up most of the figure width,
    # leaving a narrow strip on the right just for the compact legend
    fig.subplots_adjust(left=0.08, right=0.82, top=0.9, bottom=0.12)
    plt.show()
