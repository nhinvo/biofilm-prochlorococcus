"""
Purpose: to plot of HL and LL particle/free-living fractional data across depth.  
"""
import pandas as pd 
import numpy as np 
from pathlib import Path 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
   

def main():
    Path('data').mkdir(exist_ok=True)

    # import table 
    fpath = '../CalculateAbsolute/data/AbsoluteGE/CladeAbosluteGenomeEquivalent.tsv'
    df = pd.read_table(fpath)

    # obtain samples up to 150m depths
    df = df[df['depth'] <= 150]

    PB = df[df['cell_state'] == 'Particle-bound'].copy()
    PB = PB.groupby('depth')['log10 Absolute Genome Equivalent'].agg(['mean', 'std']).reset_index()
    print(PB)

    FL = df[df['cell_state'] == 'Free-living'].copy()
    FL = FL.groupby('depth')['log10 Absolute Genome Equivalent'].agg(['mean', 'std']).reset_index()
    print(FL)


    # plotting 
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    # colors
    FL_color = '#d73027'
    PB_color = '#4575b4'
    
    # ax.errorbar(FL['mean'], FL['depth'], xerr=FL['std'], label='Free-living', fmt='-o', color=FL_color)
    # ax.errorbar(PB['mean'], PB['depth'], xerr=PB['std'], label='Particle-bound', fmt='-o', color=PB_color)

    # Plot the mean lines
    ax.plot(FL['mean'], FL['depth'], label='Free-living', marker='o', color=FL_color)
    ax.plot(PB['mean'], PB['depth'], label='Particle-bound', marker='o', color=PB_color)

    # shade the standard deviation
    ax.fill_betweenx(FL['depth'], FL['mean'] - FL['std'], FL['mean'] + FL['std'], color=FL_color, alpha=0.2)
    ax.fill_betweenx(PB['depth'], PB['mean'] - PB['std'], PB['mean'] + PB['std'], color=PB_color, alpha=0.2)

    ax.set_title('Free-living vs. Particle-bound', size=14)
    ax.set_ylabel('Depth (m)', size=12)
    ax.set_xlabel('Log10(Absolute Genome Equivalent)', size=12)

    # y-axis edits (set ticks where data points are)
    ax.set_yticks(FL['depth'].values)  
    
    # invert y so 0m depth is on top
    ax.invert_yaxis()   

    # add legend to right of subplot 
    custom_lines = [Line2D([0], [0], color=FL_color, lw=10),
                    Line2D([0], [0], color=PB_color, lw=10)]
    legend_labels = ['Free-living', 'Particle-bound']
    ax.legend(custom_lines, legend_labels, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1), frameon=False)

    plt.tight_layout()

    plt.savefig(f'data/LinePlot.png')
    plt.savefig(f'data/LinePlot.svg')
    plt.close()
    


main()