# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 10:46:38 2015

@author: Sierra Anderson

Perform enrichment tests for this data and write results to a .tab file.
"""

# General imports
import math
import operator

# Specific imports that must be pre-installed 
from scipy import stats
from numpy import nanmean
import pandas as pd

# Helper methods

def __round_sig(x, n=4):
    """ Rounds x to n significant digits.
    
    Args: 
        x: value to be rounded 
        n: number of significant digits to round to (default=4).
        
    Returns:
        x rounded to n significant digits
    """
    return round(x, n - int(math.floor(math.log10(x))) - 1) 

def __write_to_file(output_dir, p_values, nans, fname):
    """ Write the p-values to file. 
    
    Args:
        output_dir: path to directory to save output.
        p_values: list of tuple'd p_values in the form (p, attribute name, class it is enriched in)
        nans: list of tuple'd NaN values in the form (NaN, 'attribute name)
        fname: filename to save output
    """
    fname = output_dir + "/" + fname
    
    f = open(fname, 'w')
    f.write('name\tp-val\tenrinched in\n')
    p_values.sort()
    
    for tp in p_values:
        pval = ("%.12f" % __round_sig(tp[0])).rstrip('0')
        attr_name = str(tp[1])
        enriched_in = str(tp[2])
        f.write(attr_name + "\t" + pval + "\t" + enriched_in + "\n")

    for n in nans:
        attr_name = str(n[1])
        f.write(attr_name + "\tn/a\n")

    f.close()

def __split_data(profile):
    """ Split and return the abundance data into a list of DataFrames based on
    class and return them. 
    
    Args:
        profile: a metagenomic_profile instance 
        
    Returns:
        list of DataFrames
    """
    dataframes = list()
    key_list = list(profile.references.keys())
    for key in key_list:
        df = profile.abundance_data.loc[profile.references[key]]
        dataframes.append(df)
        
    return dataframes

def __get_means(df1, df2):
    """ Returns dictionary mapping attribute to a tuple with the mean
    for that attribute from the first argument and the mean for
    that attribute from the second argument.
    
    Requires: 
        Both DataFrames have the same number of columns.
        
    Args:
        df1, df2: pandas DataFrames
        
    Returns:
        Dictionary with attributes as keys and a tuple containing the 
        mean for this attribute in df1 and the mean in df2
    """
    means = dict()
    
    for attr in df1.columns:
        df1_mean = nanmean(df1[attr])
        df2_mean = nanmean(df2[attr])
        means[attr] = (df1_mean, df2_mean)
    
    return means

def __bonferonni_correction(pvalues):
    """ Performs a Bonferroni correction on the data.
    
    Args:
        pvalues: list of p-values to be corrected
        n: number of hypotheses (i.e. attributes)
        
    Returns:
        Bonferroni corrected list of p-values.
    """
    result = list()
    n = len(pvalues)
    
    for p in pvalues:
        p_adjust = p*n if p*n < 1.0 else 1.0
        result.append(p_adjust)
        
    return result
    
def __fdr_correction(pvalues, FDR=0.1):
    """ Performs Benjamini-Hochberg procedure to adjust p-values.
    
    Args:
        pvalues: dictionary of p-values to be corrected, in the form {label1:p-value1, ...}
        FDR: the false discovery rate, a float x such that 0 < x < 1 (default=0.1)
        
    Returns:
        Benjamini-Hochberg adjusted p-values with corresponding attribute in a list of tuples.
    """
    sorted_values = sorted(pvalues.items(), key=operator.itemgetter(1)) # tuple representation of dict
    result = list()
    n = len(pvalues)
    
    for i in range(len(sorted_values)):
        lbl = sorted_values[i][0]
        p = sorted_values[i][1]
        result.append((lbl, p * (float(n) / (i + 1))))
    
    return result 
    
def __enrichment(df1, df2, label1, label2, output_dir, enrichment_type, correction=None):
    """ Helper method to perform the enrichment.
    
    Args:
        df1: The dataframe containing the first class of samples (i.e. control).
        df2: The dataframe containing the second class of samples (i.e. case).
        label1 (str): Label for df1.
        label2 (str): Label for df2. 
        output_dir: directory output is saved to. 
        enrichment_type: type of enrichment test to perform ('ranksums' or 'ttest').
        correction: type of correction to be performed. Options: "bonferroni", "fdr-0.1", 
            "fdr-O.05", "fdr-0.01"
            
    Effects:
        Writes out to the output_dir a .tab file containing the p-values. 
    Returns:
        Filename (str).
    """
    
    means = __get_means(df1, df2)
    pvals = list()
    nans = list()
    for attr in df1.columns:
        directionality = "n/a" # which class is the attribute enriched in
        if means[attr][0] > means[attr][1]:
            directionality = label1
        elif means[attr][1] > means[attr][0]:
            directionality = label2
        p = enrichment_type(df1[attr], df2[attr])[1]
        if math.isnan(p):
            nans.append((p, attr))
        else:
            pvals.append((p, attr, directionality))
            
    if correction == "bonferroni":
        ps = __bonferonni_correction([x[0] for x in pvals])
        pvals = [(ps[i], pvals[i][1], pvals[i][2]) for i in range(len(ps))]        
    elif correction.split("-")[0] == "fdr":
        rate = float(correction.split("-")[1])
        directions = dict([(pvals[i][1], pvals[i][2]) for i in range(len(pvals))]) # preserve directionality
        ps = dict([(pvals[i][1], pvals[i][0]) for i in range(len(pvals))])
        corrected = __fdr_correction(ps, FDR=rate)
        pvals = [(corrected[i][1], corrected[i][0], directions[corrected[i][0]]) for i in range(len(corrected))]
    __write_to_file(output_dir, pvals, nans, enrichment_type.__name__ + ".tab")
    return output_dir + "/" + enrichment_type.__name__ + ".tab"


def __pairwise(profile, output_dir, enrichment_type):
    """ Performs pairwise comparisons. Save results to file. 
    
    Args:
        profile: metagenomic profile instance.
        output_dir: directory output is saved to. 
        enrichment_type: type of enrichment test to perform ('ranksums' or 'ttest'). 
    
    Returns:
        a list of output files  
    
    """
    output_files = list()
    reference_keys = list(profile.references.keys())
    dataframes = __split_data(profile)
    for i in range(len(reference_keys) - 1):
        for j in range(i + 1, len(reference_keys)):
            df1 = dataframes[i]
            df2 = dataframes[j]
            out = __enrichment(df1, df2, reference_keys[i], reference_keys[j], output_dir, stats.ttest_ind)
            output_files.append(out)
            
    return output_files
    
# Public methods 

def ranksums(profile, output_dir, correction=None):
    """Perform the Wilcoxon Rank-Sum test. Save results to file.
    
    Args:
        profile: metagenomic profile instance.
        output_dir: directory output is saved to.
        correction: type of correction to be performed. Options: "bonferroni", "fdr-0.1", 
            "fdr-O.05", "fdr-0.01"
        
    Returns:
        Path to output.
    """
    if len(list(profile.references.keys())) > 2:
        outfiles = __pairwise(profile, output_dir, stats.ranksums)
        dfs = list()
        for fname in outfiles:
            df = pd.DataFrame.from_csv(fname, sep="\t")
            df['Comparison'] = pd.Series([fname.split(".")[0] for i in range(len(df.index))], index=df.index)
            dfs.append(df)
            
        dfs = pd.concat(dfs)
        fname = output_dir + "/enrichment_master.tab"
        dfs.to_csv(fname, sep="\t")
        return fname
    else:
        dfs = __split_data(profile)
        df1, df2 = dfs[0], dfs[1]
        ref_keys = list(profile.references.keys())
        return __enrichment(df1, df2, ref_keys[0], ref_keys[1], output_dir, stats.ranksums, correction=correction)
    

def ttest(profile, output_dir, correction=None):
    """Perform the Student's t-test. Save results to file.
    
    Args:
        profile: metagenomic profile instance.
        output_dir: directory output is saved to.
        correction: type of correction to be performed. Options: "bonferroni", "fdr-0.1", 
            "fdr-O.05", "fdr-0.01"
            
    Returns:
        Path to output.
    """
    if len(list(profile.references.keys())) > 2: # then do a pairwise comparison 
        outfiles = __pairwise(profile, output_dir, stats.ttest_ind)
        dfs = list()        
        for fname in outfiles:
            df = pd.DataFrame.from_csv(fname, sep="\t")
            df['Comparison'] = pd.Series([fname.split(".")[0] for i in range(len(df.index))], index=df.index)
            dfs.append(df)
            
        dfs = pd.concat(dfs)
        fname = output_dir + "/enrichment_master.tab"
        dfs.to_csv(fname, sep="\t")
        return fname
    else:
        dfs = __split_data(profile)
        df1, df2 = dfs[0], dfs[1]
        ref_keys = list(profile.references.keys())
        return __enrichment(df1, df2, ref_keys[0], ref_keys[1], output_dir, stats.ttest_ind, correction=correction)