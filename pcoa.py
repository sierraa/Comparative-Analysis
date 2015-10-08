# -*- coding: utf-8 -*-

"""
@author: Sierra Anderson

Plot principal coordinate analysis (PCoA) or principle component analysis (PCA). 
Provides methods to generate either static or interactive plots.
"""
# General imports
import warnings
import os
import math

# specific imports that must be pre-installed
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

# complex numbers must be cast to real in order to plot
warnings.simplefilter("ignore", np.ComplexWarning)

# colors for plot markers (in order):
# blue, yellow, red, green, magenta, sea green, orange, lime green,
# hot pink, cyan, dark red, dark blue, peach, gray, dark green, lavendar
colors = ["#0000ff", "#ff0000", "#ffff00", "#00ff00", "#ff0066", "#99ff99", 
          "#ff9900", "#660066", "#99ff00", "#ff0099", "#99ffff", "#990000", 
          "#000066", "#ff9966", "#c0c0c0", "#006600", "#cc99ff"]

# plot marker shapes
markers = ['o', 'D', 'v', 'd', '<', 'h', '+', 's', '>', '|', 'p', 'H', '.',
           'x', '*', '^', ',', '_']

# Helper methods

def __check_input(output_dir, num_of_loadings=0):
    """ Checks input for errors.
    
    Args:
        output_dir (str): Directory to save the output in. 
        num_of_loadings (int, default=0): Number of loadings to plot on the PCA.
        
    Raises:
        IOError: if output directory is not valid.
        ValueError: if num_of_loadings < 0 or num_of_loadings >= 5
    """
    if num_of_loadings < 0:
        raise ValueError("Negative number of loadings.")
    elif num_of_loadings >= 5:
        raise ValueError("Maximum number of loadings is 5.")
    elif not os.path.isdir(output_dir):
        raise IOError("Output directory '" + output_dir + "'not found.")

def __partition_abundance_data(profile):
    """Partition data by sample class.
    
    Args:
        profile (metagenomic_profile): profile to be partitioned. .
        
    Returns:
        New partitioned DataFrame.
    """
    abundances_to_class = dict()

    for k in list(profile.references.keys()):
        abundances_to_class[k] = profile.abundance_data.loc[profile.references[k]]
    
    return pd.concat(list(abundances_to_class.values()))

def __get_eig_pairs(data, dist_type):
    """ Computes eigenvalues and eigenvectors for this matrix. Calculated
    using methods described in Numerical Ecology (pp 391-443, Legendre 1998).
    
    Args:
        data (pandas.DataFrame): a sorted (by sample class) matrix 
            containing abundance data.
        dist_type (str) : distance metric to use (Euclidean for PCA)
    
    Returns:
        List of pairs of eigenvalues and eigenvectors.
    """
    dist_matrix = pairwise_distances(data, metric=dist_type)
    
    # 9.20
    A_matrix = np.linalg.matrix_power(dist_matrix, 2) / -2
    n = int(A_matrix.shape[0])
    a_mean = np.mean(A_matrix)
    
    # 9.21
    ctr_matrix = [[0 for i in range(n)] for j in range(n)]
    
    for i in range(n):
        for j in range(n):
            s = A_matrix[i][j] - np.mean(A_matrix[i][:]) - np.mean(A_matrix[:][j]) + a_mean
            ctr_matrix[i][j] = s
    
    eig_val, eig_vec = np.linalg.eig(ctr_matrix)
    eig_pairs = [(eig_val[i], eig_vec[:,i]) for i in range(len(eig_val))]
            
    return eig_pairs

def __get_loadings(rotation):
    """ Get loadings for this matrix.     
    
    Args:
        rotatation (pandas.DataFrame): matrix of features x principal components
        
    Returns: 
        List of all loading factors in the form: 
        (feature name, loading vector norm, loading vector) sorted by vector norm.

    """
    all_loadings = list()
    for i in rotation.index:
        all_loadings.append((i, float(np.linalg.norm(rotation.loc[i])), list(rotation.loc[i])))
        
    all_loadings = sorted(all_loadings, reverse=True, key= lambda triple: triple[1])
    
    return all_loadings
    
def __plot_markers(profile, PC1, PC2, msize=5, a=0.9):
    """Plot markers on the current plot.
    
    Args:
        PC1 (list or seq): 1st PC/PCo
        PC2 (list or seq): 2nd PC/PCo
        msize (int, default=5): size of plot markers. 
        a (int or float, default=0.9): transparency of markers.
    """
    
    marker_index = color_index = i = prev = 0

    for k in list(profile.references.keys()):
        n = len(profile.references[k]) # number of samples in this class

        plt.plot(PC1[prev:n+prev], PC2[prev:n+prev], markers[marker_index], 
                     markersize=5, alpha=0.9, color=colors[color_index], label=k)

        marker_index = (marker_index + 1) % len(markers) # enable wrap around
        color_index = (color_index + 1) % len(colors)
        prev += n
        i += 1

def __plot_markers_interactive(profile, PC1, PC2, a=0.9):
    """Plot markers on the current plot.
    
    Args:
        PC1 (list or seq): 1st PC/PCo
        PC2 (list or seq): 2nd PC/PCo
        a (int or float, default=0.9): transparency of markers.
        
    Returns:
        Dictionary of legend labels mapping labels to color patches. 
    """
    import mpld3 # Provides interactive graphs 
    
    fig, ax = plt.subplots()
    
    fig.set_size_inches(12, 12)    
    
    marker_index = color_index = prev = i = 0    
    
    lgd_labels = dict() # for creating the legend    
    
    for k in list(profile.references.keys()):
        n = len(profile.references[k])
        
        scatter = ax.scatter(PC1[prev:n+prev], PC2[prev:n+prev], marker=markers[marker_index],
                             alpha=0.9, color=colors[color_index], label=k)
        
        lgd_labels[k] = mpatches.Patch(color=colors[color_index], label=k)
                     
        marker_index = (marker_index + 1) % len(markers) # enable wrap around
        color_index = (color_index + 1) % len(colors)
        prev += n
        i += 1
        
        labels = profile.references[k]
        
        tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels)   
        mpld3.plugins.connect(fig, tooltip)
        
    return lgd_labels
        
def __scale_loading(x_range, y_range, x_cor, y_cor):
    """ Scale this loading vector to fit on the plot.
    
    Args:
        x_range: range of the x-axis. 
        y_range: range of the y-axis. 
        x_cor: x-coordinate of the loading.
        y_cor: y-coordinate of the loading.
    
    Returns:
        Two floats x, y that are resized coordinates for the vector 
    """
    magnitude = math.sqrt(x_cor**2 + y_cor**2)
    
    normalized = [x_cor /magnitude, y_cor / magnitude] # get the unit vector
    
    scale_factor = min(x_range, y_range) / 2
    
    return normalized[0]*scale_factor, normalized[1]*scale_factor
        
def __plot_loadings(loadings, num_of_loadings, sz=6, interactive=False):
    """Plot loadings for this plot.
    
    Args:
        loadings (list[ ]):
        num_of_loadings (int): number of loadings to be plotted. 
        sz (int, default=6): Font size for annotations.
        interactive (bool, default=False): If True arrow dimensions are altered for better 
            visibility on interactive plots.
    """
    
    x_range = abs(plt.xlim()[0]) + abs(plt.xlim()[1])
    y_range = abs(plt.ylim()[0]) + abs(plt.ylim()[1])
    
    for i in range(num_of_loadings):
        x_cor = loadings[i][2][0] 
        y_cor = loadings[i][2][1] 
        
        vec_x, vec_y = __scale_loading(x_range, y_range, x_cor, y_cor)        
        
        if not interactive:
            plt.arrow(0, 0, vec_x, vec_y, head_width=0.000625, 
                      head_length=0.000625, width=0.0001, color='black')
            name = loadings[i][0].replace(" ", "\n")
        else:
            plt.plot([0, vec_x], [0, vec_y], color='black')
            name = loadings[i][0]
            
        plt.annotate(xy=(vec_x, vec_y), s=name, size=sz, 
                   bbox=dict(facecolor='w', edgecolor='none', alpha=0.75))

def __create_legend(lgd_labels, output_dir):
    
    lgd_fig = plt.figure(figsize=(6,4))
    lf = lgd_fig.legend(handles=list(lgd_labels.values()), labels=list(lgd_labels.keys()),
                        loc="center", ncol=2, fontsize=12)
                        
    lf.get_frame().set_linewidth(0.0)
    lgd_fig.canvas.draw()
    
    lgd_fig.savefig(output_dir + "/legend.png", bbox_inches=lf.get_window_extent().transformed(lgd_fig.dpi_scale_trans.inverted()))
        
    lgd_fname = output_dir + "/" + "legend.png"
    lgd_fig.savefig(lgd_fname, bbox_inches=lf.get_window_extent().transformed(lgd_fig.dpi_scale_trans.inverted()))
    
    return lgd_fname
    
# General methods 

def pca_plot(profile, output_dir, filename="pca.png", num_of_loadings=3):
    """Generate PCA plot. A PCA is a PCoA with a Euclidean distance metric. Non-interactive. 
    
    Args:
        profile (metagenomic_profile): profile instance containing data. 
        output_dir (str): path to directory to save output
        filename (str, default="pca.png"): name of png file to be saved. 
        num_of_loadings (int, default=3, max=5, min=0): number of loadings to display.
        
    Returns:
        Path to output.
    """
    __check_input(output_dir, num_of_loadings)
    
    df = __partition_abundance_data(profile)

    my_pca = PCA(n_components=2)
    df_new = my_pca.fit_transform(df).T

    # matrix of variable loadings
    rotation = pd.DataFrame(my_pca.components_.T, index=df.columns, columns=[1, 2])
    
    if num_of_loadings > 0:
        loadings = __get_loadings(rotation)
    
    PC1 = df_new[:][0]
    PC2 = df_new[:][1]
    PC1_variance, PC2_variance = my_pca.explained_variance_ratio_[0]*100, my_pca.explained_variance_ratio_[1]*100

    # Begin plotting

    # Clear any current figures from the plot
    plt.clf()

    ax = plt.subplot(111)
    
    __plot_markers(profile, PC1, PC2) # Main plotting 

    # Plot the loadings 
    if num_of_loadings > 0:
        __plot_loadings(loadings, num_of_loadings)
    
    plt.xlabel("PC1" + " (" + str(round(PC1_variance, 2)) + "%)")
    plt.ylabel("PC2" + " (" + str(round(PC2_variance, 2)) + "%)")
    
    # Set up legend 
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.1), numpoints=1) 

    # Rescale axes to fit loading vectors
    ax = plt.gca()
    ax.relim()
    ax.autoscale_view()
    
    fname = output_dir + "/" + "pca.png"    
    
    plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=(200))
    
    return fname
    
def pcoa_plot(profile, output_dir, dist_type):
    """Generate PCoA plot. Non-interactive. Saves png file "pcoa_(dist_type).png."
    
    Args:
        profile (metagenomic profile):  profile instance containing abundance data. 
        output_dir (str): directory to save output. 
        dist_type (str): distance metric to use for PCoA. 
        
    Returns:
        Path to output.
    """
    __check_input(output_dir)    
    
    df = __partition_abundance_data(profile)
    eig_pairs = __get_eig_pairs(df, dist_type)
    
    eig_pairs.sort()
    eig_pairs.reverse()    
            
    PCo1 = eig_pairs[0][1]
    PCo2 = eig_pairs[1][1]
    
    # Begin plotting

    # Clear any current figures from the plot
    plt.clf()

    ax = plt.subplot(111)
    
    __plot_markers(profile, PCo1, PCo2) # Main plotting
    
    plt.xlabel("PCo1")
    plt.ylabel("PCo2")
    
    # Set up legend 
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.1), numpoints=1)
    
    fname = output_dir + "/" + "pcoa_" + dist_type + ".png"    
    
    plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=(200))
    
    plt.clf()
    
    return fname

def pca_plot_interactive(profile, output_dir, num_of_loadings=3):
    """Generate interactive PCA plot. Saves html file "pca.html."
    
    Args:
        profile (metagenomic_profile): profile containing abundance data to be plotted. 
        output_dir (str): path to directory to save output.
        num_of_loadings (int, default=3, max=5, min=0): number of loadings to plot.
        
    Returns:
        Path to output file.
    """
    
    try:
        import mpld3 # Provides interactive graphs 
        import plugins # Custom mpld3 plugins
    except ImportError:
        print("Could not import mpld3. Please install mpld3 or set 'interactive_plots' option to 'false.'")
        return
        
    __check_input(output_dir, num_of_loadings)
    
    df = __partition_abundance_data(profile)
    
    my_pca = PCA(n_components=2)
    df_new = my_pca.fit_transform(df).T

    # matrix of variable loadings
    rotation = pd.DataFrame(my_pca.components_.T, index=df.columns, columns=[1, 2])
    
    if num_of_loadings > 0:
        loadings = __get_loadings(rotation)

    PC1 = df_new[:][0]
    PC2 = df_new[:][1]
    
    PC1_variance, PC2_variance = my_pca.explained_variance_ratio_[0]*100, my_pca.explained_variance_ratio_[1]*100
    
    # Begin plotting

    # Clear any current figures from the plot
    plt.clf()
    
    lgd_labels = __plot_markers_interactive(profile, PC1, PC2, a=0.7) # Main plotting 
    
    # Plot the loadings 
    if num_of_loadings > 0:
        __plot_loadings(loadings, num_of_loadings, sz=10, interactive=True)
    
    # Padding for x and y labels 
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', pad=15)
    
    plt.xlabel("PC1" + " (" + str(round(PC1_variance, 2)) + "%)", fontsize=16)
    plt.ylabel("PC2" + " (" + str(round(PC2_variance, 2)) + "%)", fontsize=16)
    
    fname = output_dir + "/" + "pca.html"    
    
    fig = plt.gcf()
    
    mpld3.plugins.connect(fig, plugins.TweakToolbar())
    mpld3.save_html(fig, fname)
        
    lgd_fname = __create_legend(lgd_labels, output_dir)
    
    return fname, lgd_fname

def pcoa_plot_interactive(profile, output_dir, dist_type):
    """Generate interactive PCoA plot.
    
    Args:
        profile (metagenomic_profile): Profile instance containing data. 
        output_dir (str): path to directory to save output
        dist_type (str): distance metric to use in PCoA.
        
    Returns:
        Path to output file. 
    """
    try:
        import mpld3 # Provides interactive graphs 
        import plugins # Custom mpld3 plugins    
    except ImportError:
        print("Could not import mpld3. Please install mpld3 or set 'interactive_plots' option to 'false.'")
        return
        
    __check_input(output_dir)
    
    df = __partition_abundance_data(profile)
    eig_pairs = __get_eig_pairs(df, dist_type)
    
    eig_pairs.sort()
    eig_pairs.reverse()
        
    PCo1 = eig_pairs[0][1]
    PCo2 = eig_pairs[1][1]
    
    # Begin plotting

    # Clear any current figures from the plot
    plt.clf()
    
    lgd_labels = __plot_markers_interactive(profile, PCo1, PCo2, a=0.7) # Main plotting 

    # Padding for x and y labels    
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', pad=15)    
    
    plt.xlabel("PCo1", fontsize=16)
    plt.ylabel("PCo2", fontsize=16)
    
    fname = output_dir + "/" + "pcoa_" + dist_type + ".html"
    
    fig = plt.gcf()
    
    mpld3.plugins.connect(fig, plugins.TweakToolbar())
    mpld3.save_html(fig, fname)
    
    # Create legend 
    lgd_fname = __create_legend(lgd_labels, output_dir)
    
    return fname, lgd_fname