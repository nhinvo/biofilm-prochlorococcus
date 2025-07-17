"""
"""
import pandas as pd 
from pathlib import Path 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

PLOT_DATA = "../../../StandardizeData/data/AllSummaryReadCount.tsv"

# PLOT DETAILS / SPECIFICATIONS #
# colors
FL_color = '#e66101'
PB_color = '#b2abd2'

# lines 
linewidth=2.5  # width of line 
elinewidth=1.5  # width of error bar


def import_plot_data(fpath):
    """
    """
    df = pd.read_table(fpath)
    
    # filter for read count only (out of percent and reads) & depth 
    df = df[(df['summary_type'] == 'reads')].copy()

    # obtain total read count 
    df['classified_reads'] = df['other_genus'] + df['Prochlorococcus'] + df['Synechococcus'] 

    # obtain pro percent (out of all classified reads )
    df['pro_percent'] = (df['Prochlorococcus'] / df['classified_reads']) * 100 

    FL = df[df['cell_state'] == 'Free-living'].copy()
    FL = FL.groupby('binned_depth')['pro_percent'].agg(['mean', 'std']).reset_index()

    # obtain HL_PB rows 
    PB = df[df['cell_state'] == 'Particle-bound'].copy()
    # PB = df[df['cell_state'] == 'Particle-bound']
    PB = PB.groupby('binned_depth')['pro_percent'].agg(['mean', 'std']).reset_index()

    return FL, PB


def main():
    # import plot data 
    FL, PB = import_plot_data(PLOT_DATA)

    # init figure  
    fig, ax = plt.subplots(1, 1, figsize=(3, 6))

    # plot the mean lines
    ax.plot(FL['mean'], FL['binned_depth'], label='Free-living', marker='o', color=FL_color)
    ax.plot(PB['mean'], PB['binned_depth'], label='Particle-bound', marker='o', color=PB_color)

    # plot error bars on line 
    ax.errorbar(
        FL['mean'], FL['binned_depth'], xerr=FL['std'], label='Free-living', fmt='-o', color=FL_color,
        elinewidth=elinewidth, linewidth=linewidth
    )
    ax.errorbar(
        PB['mean'], PB['binned_depth'], xerr=PB['std'], label='Particle-bound', fmt='-o', color=PB_color,
        elinewidth=elinewidth, linewidth=linewidth
    )

    # shade the standard deviation
    ax.fill_betweenx(FL['binned_depth'], FL['mean'] - FL['std'], FL['mean'] + FL['std'], color=FL_color, alpha=0.2)
    ax.fill_betweenx(PB['binned_depth'], PB['mean'] - PB['std'], PB['mean'] + PB['std'], color=PB_color, alpha=0.2)

    # titles
    ax.set_title('Total $Prochlorococcus$', size=12)
    ax.set_ylabel('Depth (m)', size=12)
    ax.set_xlabel('Fraction of microbial community', size=12, x=0.35)

    # tick edits 
    ax.tick_params(axis='x', bottom=True, top=True, direction='in', which='both')
    ax.xaxis.set_major_locator(MultipleLocator(10))

    ax.yaxis.set_major_locator(MultipleLocator(200))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.tick_params(axis='y', which='both', length=3.5, width=1, left=True, right=True, direction='in')

    ax.set_ylim((-15, 1015))  # extra space above/below 0-1000
    
    # invert y so 0m depth is on top
    ax.invert_yaxis()   

    plt.tight_layout()

    plt.savefig(f'Fig1A.png')
    plt.savefig(f'Fig1A.svg')
    plt.close()
    


main()