# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:50:08 2015

@author: Sierra Anderson

Normalizes the data according to a user specified technique.
Updates the abundance data for the metagenomic profile instance
passed into the file and writes normalized data to file.
"""

import pandas as pd

# Helper methods

def __normalize_dataframe(df):
    """ Normalizes the pandas DataFrame df relatively.
    
    Args:
        df (pandas.DataFrame): data to be normalized
    """
    return df.div(df.sum(axis=1), axis=0)
            
# Public methods 

def relative_normalization(profile, output_dir):
    """Perform a relative normalization on this data. Write normalized data to file
    and update profile's abundance data to be normalized. 
    
    Args:
        profile (metagenomic_profile): profile containing the abundance data
        output_dir (str): output directory where the tab file will be saved.
        
    Effects: 
        profile's abundance data is now normalized. 
    """
    profile.abundance_data = __normalize_dataframe(profile.abundance_data)
    profile.to_file_abundance_data("normalized_abundance_data.tab", output_dir)

def musicc_normalization(profile, in_file, output_dir, musicc_inter=True, input_format='tab', output_format='tab', 
                         musicc_intra='use_generic', compute_scores=False, verbose=False):    
    """Perform a MUSiCC normalization on this profile's abundance data. Write normalized
    data to file and update profile's abundance data to be normalized.
    
    Requires:
        MUSiCC module can be imported and abundance features are KOs. 
    
    Args:
        profile (metagenomic_profile): profile containing the abundance data
        in_file (str): path to original abundance data file 
        output_dir (str): output directory where normalized data will be saved
    
    Effects: 
        profile's abundance data is now normalized.
    
    Note:
        For more details see "MUSiCC: A marker genes based framework for 
        metagenomic normalization and accurate profiling of gene abundances in the 
        microbiome." Ohad Manor and Elhanan Borenstein. Genome Biology.
    """
    
    from musicc.core import correct_and_normalize
    
    musicc_args = {'musicc_inter':musicc_inter, 'input_format':input_format, 'output_format':output_format,
               'musicc_intra':musicc_intra, 'compute_scores':compute_scores, 'verbose':verbose}    
    
    musicc_args['input_file'] = in_file 
    musicc_args['output_file'] = output_dir + "//" + "musicc_normalized_abundance.tab"
    correct_and_normalize(musicc_args)
    
    profile.set_abundance_data(pd.DataFrame.from_csv(musicc_args['output_file'], sep='\t'))