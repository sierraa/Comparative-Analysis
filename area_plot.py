# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 10:09:57 2015

@author: Sierra Anderson

Plot area plot. Provides methods to generate either static or interactive plots.
"""

# General imports
from random import random 
import sys

# specific imports that must be pre-installed
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt

MAX_DATA_POINTS = 20 # Maximum number of attributes that can be plotted on the area plot. Anymore and the plot is not legible and very slow to create.

WIDTH = 0.8 # Width of the bars

# Colors for area plot:
# blue, yellow, red, green, magenta, aqua, orange, purple, lime green,
# hot pink, cyan, dark red, dark blue, peach, gray, dark green, lavendar

colors = ["#0000ff", "#ffff00", "#ff0000", "#00ff00", "#ff0066", "#99ff99", "#ff9900", 
         "#660066", "#99ff00", "#ff0099", "#99ffff", "#990000", "#000066", "#ff9966",
         "#c0c0c0", "#006600", "#cc99ff"]

# Helper methods

def __find_most_abundant(df, row_label):
    """ Finds the most abundant attribute in this row of the DataFrame.
    
    Requires:
        DataFrame entries are non-negative numbers.
    
    Args:
        df: pandas DataFrame
        row_label: label of the row to search for the most abundant attribute in.

    Returns:
        Label of the most abundant column value in this row_label.
    """
    sample = df.loc[row_label]
    maximum = 0
    max_name = ""
    for k in sample.index:
        if sample[k] >= maximum:
            max_name = k
            maximum = sample[k]
    
    return max_name

def __sort_by_most_abundant(df):
    """ Sorts the DataFrame rows by the most abundant column.
    
    Args:
        df: pandas DataFrame.

    Effects:
        df rows are now sorted by the DataFrame's most abundant column.
    """
    column_label = __find_most_abundant(df, df.index[0])
    df.sort(columns=column_label, axis=0, inplace=True)

def __generate_colors(n):
    """ Generate a list of n colors from a predefined color list.
    
    Requires:
        n >= 0
        
    Args:
        n: number of colors to generate.
    
    Returns:
        a list of at least n colors
    """
    assert n >= 0
    
    if len(colors) >= n:
        return colors

    for i in range(n - len(colors)):
        colors.append((random(), random(), random()))

    return colors
    
def __sort_for_plot(df):
    """ Sorts the columns in roughly reverse order such that
    the first column in the DataFrame has the largest numerical
    entries. 

    Args:
        df: pandas DataFrame.

    Returns:
        A list of the new order of the columns. 
    """
    series = df.iloc[0]
    series.sort(ascending=False)
    return list(series.index)

def __get_xticks(profile):
    """ Get xticks for this plot.
    
    Args:
        profile: metagenomic profile instance.

    Effects:
        Plots separating line on current axes. 

    Returns:
        List of values to put ticks at.

    Note:
        Not compatible with mpld3's interactive plot. 
    """
    ticks = list()
    ticks.append(0)
    running = 0

    for cls in list(profile.references.keys()):
        running += len(profile.references[cls])*WIDTH
        ticks.append(running)
        plt.axvline(x=running, color='black')

    for i in range(len(ticks) - 1): # center the labels
        ticks[i] = (ticks[i] + ticks[i + 1]) / 2

    return ticks

def __plot_bars(profile, colors, interactive=False):
    """ Plot bars for area plot.

    Args:
        profile: metagenomic_profile instance 
        colors: list of colors to use to plot bars
        interactive: 
    
    Effects:
        Plots bars on current axes.

    Returns:
        legend labels, and list of labels, ids for use with interactive plotting
    """
    
    if interactive:
        import mpld3
    
    w = 0 # x-coordinate to plot each new bar on    
    prev = dict()
    plt.clf()
    
    lgd_labels = dict() #(keys, values) = (labels, handles) for plotting legend
    
    # interactive specific data structures
    interactive_labels = list()
    ids = []
    
    last_class = list(profile.references.keys())[-1]

    for cls in list(profile.references.keys()):
        class_df = profile.abundance_data.loc[profile.references[cls]]
        __sort_by_most_abundant(class_df)

        # change order of columns so most abundant attribute is plotted first
        l = list(class_df.columns)[::-1]
        class_df = class_df[l]

        last_sample = class_df.index[-1]

        for sample in class_df.index:
            for i in range(len(class_df.columns)):
                attr = class_df.columns[i]
                if i == 0:
                    prev[sample] = 0
                bars = plt.bar(w, profile.abundance_data.loc[sample, attr], linewidth=0, 
                               bottom=prev[sample], color=colors[i])
                if interactive:
                    ids.extend([mpld3.utils.get_id(bar) for bar in bars])
                prev[sample] += class_df.loc[sample, attr]
                if attr not in list(lgd_labels.keys()):
                    lgd_labels[attr] = mpatches.Patch(color=colors[i], label=attr)
                #interactive_labels.append(attr + " ("+cls+", " + sample + ")")
                interactive_labels.append(sample)
            
            # plot separating line 
            if sample == last_sample and cls != last_class:
                w += WIDTH
                plt.bar(w, 1, linewidth=0, color='black')
            w += WIDTH
            
    return lgd_labels, interactive_labels, ids

    
# Public methods

def area_plot(profile, output_dir):
    """ Create a png area plot of the attributes on this data. Assumes that
    the data is already normalized. 

    Args:
        profile: metagenomic profile instance.
        output_dir: output directory where png will be saved.
    
    Raises:
        ValueError: if the number of attributes is greater than MAX_DATA_POINTS.
        
    Returns:
        path to output image
    """
    
    plt.clf()
    
    df = profile.abundance_data
    if len(df.columns) > MAX_DATA_POINTS:
        raise ValueError("Too many attributes to generate area plot.")
    
    col_label = __sort_for_plot(df)
    df.sort(columns=col_label, axis=0, inplace=True)
    colors = __generate_colors(len(df.columns))
    
    lgd_labels = __plot_bars(profile, colors)[0]

    ticks = __get_xticks(profile)

    plt.xticks(ticks, list(profile.references.keys()))

    plt.xlim(0, len(profile.abundance_data.index)*WIDTH)
    plt.ylim(0, 1)

    plt.xlabel("Samples")
    plt.ylabel("Relative Abundance")

    lgd = plt.legend(title="Features", handles=list(lgd_labels.values()),
                     loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2,
                     fontsize=8)
    
    fname = output_dir + "/" + "area_plot.png"
    plt.savefig(fname, bbox_extra_artists=(lgd,),
                bbox_inches='tight', dpi=(400), figsize=(24,24))
                
    return fname
    
def area_plot_interactive(profile, output_dir):
    """ Create an interactive area plot of the attributes on this data.
    Assumes that the data is already normalized. 

    Args:
        profile: a metagenomic profile instance.
        output_dir: output directory to save html file. 
    
    Returns:
        path to output, path to legend
    
    Raises:
        ValueError: if the number of attributes is greater than MAX_DATA_POINTS.
    """
    try:
        import mpld3
        import plugins
    except ImportError:
        print("Error: Could not import mpld3. Please install mpld3 or set 'interactive_plots' option to 'false.'")
        sys.exit(0)
        
    plt.clf() 
    plt.figure(figsize=(12,12))
    
    df = profile.abundance_data
    if len(df.columns) > MAX_DATA_POINTS:
        raise ValueError("Too many attributes to generate area plot.")
    
    # Do some sorting
    col_label = __sort_for_plot(df)
    df.sort(columns=col_label, axis=0, inplace=True)
    colors = __generate_colors(len(df.columns))
    
    # Main plotting
    lgd_labels, interactive_labels, ids = __plot_bars(profile, colors, interactive=True)
    
    plt.xticks([])    
    
    plt.xlim(0, len(profile.abundance_data.index)*WIDTH + WIDTH*(len(list(profile.references.keys())) - 1)) # last part of sum to account for separating line
    plt.ylim(0, 1)
    
    # Padding for x and y labels 
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', pad=15)
    
    plt.ylabel("Relative Abundance", fontsize=16)
    
    # add plugins
    fig = plt.gcf()

    # plugin for labels on bar plot
    mpld3.plugins.connect(fig, plugins.BarLabelToolTip(ids, interactive_labels))
    
    # plugin to unhide toolbar        
    mpld3.plugins.connect(fig, plugins.TweakToolbar()) 
    
    fname = output_dir + "/" + "area_plot_interactive.html"
    mpld3.save_html(plt.gcf(), fname)
    
    lgd_fig = plt.figure(figsize=(10,8))

    # save legend as separate feature
    lf = lgd_fig.legend(title="Features", handles=list(lgd_labels.values()),
                        labels=list(lgd_labels.keys()), loc="center", ncol=2, fontsize=12)
    
    lf.get_frame().set_linewidth(0.0) # remove box around legend
    
    lgd_fig.canvas.draw()
        
    lgd_fname = output_dir + "/" + "legend.png"
    lgd_fig.savefig(lgd_fname, bbox_inches=lf.get_window_extent().transformed(lgd_fig.dpi_scale_trans.inverted()))
                
    return fname, lgd_fname