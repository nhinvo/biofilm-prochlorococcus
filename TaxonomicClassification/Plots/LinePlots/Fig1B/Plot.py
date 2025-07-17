"""
Purpose: to plot of HL and LL particle/free-living fractional data across depth.  
"""
import pandas as pd 
import numpy as np 
from pathlib import Path 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

PLOT_DATA = '../../../StandardizeData/data/AllNormalizedCounts_PlotFormat.tsv'

# PLOT DETAILS / SPECIFICATIONS #
# colors
FL_color = '#e66101'
PB_color = '#b2abd2'
HL_color = '#bbede8'
LL_color = '#edebbb'

# lines 
std_shading_alpha = 0.15  # alpha of std error shading 
linewidth=2.5  # width of line 
elinewidth=1.5  # width of error bar

# subplot 
axis_title_fsize = 12  # axis (y/x) titles 
subplot_title_pos = 1.01  # position of axes title
axes_title_fsize = 10  # subplot title (i.e. HL/LL)


def import_plot_data(fpath):
    """
    Split data into sub-groups for plotting. 
    """
    df = pd.read_table(fpath)

    # obtain HL_FL rows 
    HL_FL = df[(df['classification'] == 'HL') & (df['cell_state'] == 'Free-living')].copy()

    # obtain HL_PB rows 
    HL_PB = df[(df['classification'] == 'HL') & (df['cell_state'] == 'Particle-bound')].copy()

    # obtain LL_FL rows 
    LL_FL = df[(df['classification'] == 'LL') & (df['cell_state'] == 'Free-living')].copy()

    # obtain LL_PB rows 
    LL_PB = df[(df['classification'] == 'LL') & (df['cell_state'] == 'Particle-bound')].copy()

    return HL_FL, HL_PB, LL_FL, LL_PB

def main():
    # import plot data 
    HL_FL, HL_PB, LL_FL, LL_PB = import_plot_data(PLOT_DATA)

    # init figure and subplots 
    fig, (HL_ax, LL_ax) = plt.subplots(1, 2, figsize=(5, 6), sharex=True)    
    

    #######################
    ### 1. plot HL axes ###
    ## 1a. HL Free-living ## 
    # i. obtain mean/std + plot line and error bar 
    HL_FL = HL_FL.groupby('binned_depth')['percentage'].agg(['mean', 'std']).reset_index()
    HL_ax.errorbar(
        HL_FL['mean'], HL_FL['binned_depth'], xerr=HL_FL['std'], 
        fmt='-o', color=FL_color, elinewidth=elinewidth, linewidth=linewidth
    )

    # ii. shading the standard deviation
    HL_ax.fill_betweenx(
        HL_FL['binned_depth'], HL_FL['mean'] - HL_FL['std'], HL_FL['mean'] + HL_FL['std'], 
        color=FL_color, alpha=std_shading_alpha
    )

    ## 1b. HL Particle-bound ## 
    # i. obtain mean/std + plot line and error bar 
    HL_PB = HL_PB.groupby('binned_depth')['percentage'].agg(['mean', 'std']).reset_index()
    HL_ax.errorbar(
        HL_PB['mean'], HL_PB['binned_depth'], xerr=HL_PB['std'], 
        fmt='-o', color=PB_color, elinewidth=elinewidth, linewidth=linewidth
    )

    # ii. shading the standard deviation
    HL_ax.fill_betweenx(
        HL_PB['binned_depth'], HL_PB['mean'] - HL_PB['std'], HL_PB['mean'] + HL_PB['std'], 
        color=PB_color, alpha=std_shading_alpha
    )

    ## 1c. Other HL axes specifications ## 
    HL_ax.tick_params(axis='y', which='both', left=True, right=True, direction='in')
    HL_ax.set_title('High light-adapted clades', size=axes_title_fsize, y=subplot_title_pos, backgroundcolor=HL_color)


    #######################
    ### 2. plot LL axes ###
    ## 2a. HL Free-living ## 
    # i. obtain mean/std + plot line and error bar
    LL_FL = LL_FL.groupby('binned_depth')['percentage'].agg(['mean', 'std']).reset_index()
    LL_ax.errorbar(
        LL_FL['mean'], LL_FL['binned_depth'], xerr=LL_FL['std'], 
        fmt='-o', color=FL_color, elinewidth=elinewidth, linewidth=linewidth
    )

    # ii. shading the standard deviation
    LL_ax.fill_betweenx(
        LL_FL['binned_depth'], LL_FL['mean'] - LL_FL['std'], LL_FL['mean'] + LL_FL['std'], 
        color=FL_color, alpha=std_shading_alpha
    )

    ## 2b. LL Particle-bound ## 
    # i. obtain mean/std + plot line and error bar 
    LL_PB = LL_PB.groupby('binned_depth')['percentage'].agg(['mean', 'std']).reset_index()
    LL_ax.errorbar(
        LL_PB['mean'], LL_PB['binned_depth'], xerr=LL_PB['std'], 
        fmt='-o', color=PB_color, elinewidth=elinewidth, linewidth=linewidth
    )

    # ii. shading the standard deviation
    LL_ax.fill_betweenx(
        LL_PB['binned_depth'], LL_PB['mean'] - LL_PB['std'], LL_PB['mean'] + LL_PB['std'], 
        color=PB_color, alpha=std_shading_alpha
    )

    ## 2c. Other LL axes specifications ## 
    LL_ax.tick_params(axis='y', which='both', labelleft=False, left=True, right=True, direction='in')
    LL_ax.set_title('Low light-adapted clades', size=axes_title_fsize, y=subplot_title_pos, backgroundcolor=LL_color)


    #######################
    ### 3. Plot Details ### 
    # axis edits 
    fig.supxlabel('Fraction of classified $Prochlorococcus$', size=axis_title_fsize, ha='center')  # shared x-axis label 
    fig.supylabel('Depth (m)', size=axis_title_fsize)

    # x-axis edits for both subplots 
    for ax in [HL_ax, LL_ax]:
        # set axis limits  
        ax.set_xlim((0, 100))
        ax.set_ylim((-15, 1015))  # extra space above/below 0-1000

        # x-axis: major ticks every 20 | minor ticks every 10 
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(10))

        # x-axis ticks
        ax.tick_params(axis='x', bottom=True, direction='in', which='both')

        # y-axis: make major and minor ticks same size
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.tick_params(axis='y', which='both', length=3.5, width=1)
    
    # invert y so 0m depth is on top
    HL_ax.invert_yaxis()   
    LL_ax.invert_yaxis()  

    plt.tight_layout()

    # space between subplots
    plt.subplots_adjust(wspace=0.15)  

    plt.savefig('Fig1B.png')
    plt.savefig("Fig1B.svg")
    plt.close()

main()