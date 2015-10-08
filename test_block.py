# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:04:38 2015

@author: Sierra Anderson

"""

# Internal imports 
import test_runner as tr

class test_block:
    """ Represents a block of tests to be run on this data.
    
    Attributes:
        params (dictionary): parameters for this block of tests 
    """
    
    def __init__(self, parameters):
        """ Create a new instance of a test_block with specified parameters.
        
        Args:
            parameters: dictionary containing the parameters for this test block
        """
        self.params = parameters
        
    def set_general_parameters(self, general_parameters):
        """ Set the general global parameters for this test_block.
        
        Args:
            general_parameters: a dictionary of general parameters for this dataset.
        
        """
        self.gen_params = general_parameters
        
    def set_metagenomic_profile(self, mgprofile):
        """ Set the metagenomic profile containing the abundance and metadata for
        for this test_block.
        
        Args:
            mgprofile: metagenomic_profile instance containing data 
            
        """
        self.metagenomic_profile = mgprofile
    
    def get_type(self):
        """ Return the type of test (PCOA, PCA, Area plot, Enrichment)
        """
        return self.params["test_type"]
        
    def get_name(self):
        """ Return the name of this test. 
        """
        return self.params["test_name"]
        
    def run(self):
        """ Run tests specified in this test block. 
        
        Returns:
            list of results
        """
        return tr.test_runner(self).run_tests()
        