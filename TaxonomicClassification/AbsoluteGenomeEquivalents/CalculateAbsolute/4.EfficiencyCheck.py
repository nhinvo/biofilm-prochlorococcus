"""
Purpose: to see if there are any patterns in 
DNA efficiency (if any factor impacted efficiency). 

Conclusion: filter size = biggest influence on extraction efficiency
"""
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path 

def distribution_plots(df, fname):
    """
    """
    # 3 subplots 
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    # violin plot 
    axs[0].violinplot(df['Efficiency'])
    axs[0].set_title("Violin Plot of \nStandard Recovery Efficiency")

    # box plot 
    axs[1].boxplot(df['Efficiency'])
    axs[1].set_title("Boxplot of \nStandard Recovery Efficiency")

    # histogram plot 
    axs[2].hist(df['Efficiency'], bins=15, color='skyblue', edgecolor='black')
    axs[2].set_title("Histogram of Efficiency")
    axs[2].set_ylabel("Frequency")
    axs[2].set_xlabel("Efficiency")

    plt.tight_layout()

    plt.savefig(f"data/EfficiencyPlots/Distrubution/{fname}.png")
    plt.close()

def distribution_plotting(df):
    """
    Plot: box plot, violin plot, histogram of distribution of efficiencies. 
    """
    Path('data/EfficiencyPlots/Distrubution').mkdir(exist_ok=True, parents=True)
    # plot all efficiencies
    distribution_plots(df, "AllEfficiencies")

    # Filter 1: remove the 1 outlier point 
    df = df[df['Efficiency'] < 2e-7] 
    distribution_plots(df, "Filter1")

    # Filter 2: remove the 4 outlier points 
    df = df[df['Efficiency'] < 6e-8] 
    distribution_plots(df, "Filter2")

    # Filter 3: remove the 2 outlier points 
    df = df[df['Efficiency'] < 4.5e-8] 
    distribution_plots(df, "Filter3")

def categorical_correlation_plot(df, outlier_df, category_col):
    """
    """
    # 2 subplots 
    fig, axs = plt.subplots(2, 1, figsize=(6, 8))

    sns.boxplot(x=category_col, y='Efficiency', data=df, ax=axs[0])
    sns.stripplot(x=category_col, y='Efficiency', data=outlier_df, 
        ax=axs[0], color='red', size=6, jitter=True, marker='o')  

    sns.violinplot(x=category_col, y='Efficiency', data=df, ax=axs[1])
    sns.stripplot(x=category_col, y='Efficiency', data=outlier_df, 
        ax=axs[1], color='red', size=6, jitter=True, marker='o')  

    # plt.title(f"Box Plot of Efficiency by {category_col}")

    plt.tight_layout()

    plt.savefig(f"data/EfficiencyPlots/CategoryCorr/{category_col}.png")
    plt.close()

def categorical_correlation_plotting(df):
    """
    """
    Path('data/EfficiencyPlots/CategoryCorr').mkdir(exist_ok=True, parents=True)
    
    # remove the 1 outlier that heavily skews plot
    df = df[df['Absolute Efficiency Z-scores'] < 7].copy()

    # obtain whether or not samples are combined
    df['Combined'] = np.where(df['Notes'].str.contains('combined', case=False, na=False), 'Combined', 'Individual')

    # obtain df of outliers
    outlier_df = df[df['Absolute Efficiency Z-scores'] > 0.5].copy()

    categorical_correlation_plot(df, outlier_df, 'date extracted')
    categorical_correlation_plot(df, outlier_df, 'Filter size (um)')
    categorical_correlation_plot(df, outlier_df, 'Depth')
    categorical_correlation_plot(df, outlier_df, 'Combined')

def continuous_correlation_plot(df, continuous_col):
    """
    """
    df[continuous_col] = df[continuous_col].astype(float)

    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Efficiency', y=continuous_col, hue='color', palette={'Outlier': 'red', 'Non-Outlier': 'blue'}, data=df, s=100, alpha=0.7)

    plt.title("Scatter Plot")
    plt.xlabel("Efficiency")
    plt.ylabel(continuous_col)
    plt.grid(True)

    plt.savefig(f"data/EfficiencyPlots/ContinuousCorr/{continuous_col}.png")
    plt.close()
    

def continuous_correlation_plotting(df):
    """
    """
    Path('data/EfficiencyPlots/ContinuousCorr').mkdir(exist_ok=True, parents=True)

    # remove the 1 outlier that heavily skews plot
    df = df[df['Absolute Efficiency Z-scores'] < 7].copy()

    # create col for outlier & color 
    df['outlier'] = df['Absolute Efficiency Z-scores'] > 0.5 
    df['color'] = np.where(df['outlier'], 'Outlier', 'Non-Outlier')
    
    continuous_correlation_plot(df, 'Standard Added Vol (uL)')

def main():
    # import table 
    df = pd.read_table('data/EfficiencyMetadata/StandardData.tsv')
    # df['Efficiency'] = df['Standard Added Molecules'] / df['Standard Recovered Molecules']

    # distribution plots
    # distribution_plotting(df)

    # categorical correlation plots
    categorical_correlation_plotting(df)

    # continuous correlation plots
    continuous_correlation_plotting(df)




main()