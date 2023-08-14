import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def double_axis_plot(df,figsize,outname):
    # Set "Pathways" as the index of the dataframe
    df_pathways = df.set_index('Pathways')

    # Sort the dataframe in descending order based on Log2FE
    df_pathways_sorted = df_pathways.sort_values(by='Log2FE', ascending=False)
    # Create a new figure and a subplot
    fig, ax1 = plt.subplots(figsize=figsize)

    # Plot Log2FE as a horizontal bar plot
    # Invert the y-axis to have the pathway with the highest Log2FE at the top
    ax1.plot(df_pathways_sorted['P-value'][::-1], df_pathways_sorted.index[::-1], color='orange', marker='o')
    ax1.set_xlabel('P-value', color='orange', fontweight='bold',fontsize=18)
    ax1.tick_params('x', colors='orange', direction='out', labelsize = 18)
    # Custom P-value ticks from 0 to 0.05
    ax1.set_xticks(np.arange(0, 0.1, 0.02))
    # Make all fonts bold
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        label.set_fontweight('bold')
    for label in ax1.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(18)
        # Create a second x-axis for the same y-axis
    ax2 = ax1.twiny()

    # Plot P-value as a line plot on the second x-axis
    # Invert the y-axis to have the pathway with the highest Log2FE at the top

    ax2.barh(df_pathways_sorted.index[::-1], df_pathways_sorted['Log2FE'][::-1], color='navy', alpha=0.7)
    ax2.set_xlabel('Log2FE', color='navy', fontweight='bold', fontsize=18)
    ax2.tick_params('x', colors='navy', direction='out',labelsize = 18)
    ax2.tick_params('y', labelsize=20)  # Increase y-axis font size


    # Make all fonts bold for the second axis
    for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        label.set_fontweight('bold')

    # Remove grid
    ax1.grid(False)
    ax2.grid(False)

    # Bring the line plot to front
    ax1.set_zorder(ax2.get_zorder() + 1)  # Make the ax1 overlay on ax2
    ax1.patch.set_visible(False)

    # Show the plot
    plt.tight_layout()
    plt.savefig(outname,dpi=300)

datdir = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig6/'
# Load the original dataframe
cs_go = pd.read_excel("{}/kun_reformat.xlsx".format(datdir), sheet_name="CS_GO")
double_axis_plot(cs_go.head(5),(13,8),'{}/CoreSig_GO.png'.format(datdir))

# Load the original dataframe
cs_kegg = pd.read_excel("{}/kun_reformat.xlsx".format(datdir), sheet_name="CS_KEGG")
double_axis_plot(cs_kegg.head(5),(10,7),'{}/CoreSig_KEGG.png'.format(datdir))

# Load the original dataframe
cs_reactome = pd.read_excel("{}/kun_reformat.xlsx".format(datdir), sheet_name="CS_REACTOME")
double_axis_plot(cs_reactome.head(5),(10,7),'{}/CoreSig_REACTOME.png'.format(datdir))