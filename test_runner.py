# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:14:35 2015

@author: Sierra Anderson

The class that performs the tests for the comparative analysis. 
"""

# internal imports
import pcoa
import enrichment 
import normalization
import area_plot
from result import png_result, html_result, table_result

# General imports 
import os 

class test_runner:    
    """ A class containing functions to perform the tests for the comparative analysis. 
    
    Attributes:
        block: test_block instance containing instructions for this run. 
        new_dir: directory where the output from this test run will be saved. 
        result: list of result class instances containing the results for this test block. 
    """
    
    def __init__(self, test_block):
        """ Create a new test_runner instance. 
        
        Args:
            test_block: test_block instance containing instructions for this run. 
        """
        self.block = test_block
        os.chdir(self.block.gen_params["output_directory"])
        self.new_dir = self.block.get_name()  
        os.mkdir(self.new_dir)
        self.result = None
    
    # Helper methods 
    
    def __generate_about(self):
        """ Return string formatted with information about this test for HTML page.
        """
        result = ""
        
        if self.block.get_type() == "pca":
            loadings = self.block.params["number_of_loadings"]
            result = "PCA with " + (str(loadings) if loadings > 0 else "no") + " loadings shown.\n"
        elif self.block.get_type() == "pcoa":
            dist_metric = self.block.params["distance_metric"]
            result = "PCoA with " + dist_metric + " distance metric shown.\n"
        elif self.block.get_type() == "enrichment":
            correction = self.block.params["correction"]
            test = "student's t-test" if self.block.params["test"] == "ttest" else "Wilcoxon ranksums test"
            result = "Enrichment was performed using " + test + ".\n"
            if correction == "bonferroni":
                result += "P-values adjusted using the Bonferroni correction.\n" 
            elif correction.split("-")[0] == "fdr":
                result += "P-values adjusted using the Benjamini-Hochberg method using a false discovery rate = " 
                result += correction.split("-")[1] + ".\n"
        return result
    
    def __plot_static(self):
        """ Creates static plots. 
        """
        mgprofile = self.block.metagenomic_profile
        
        if self.block.get_type() == "area_plot":
            display_name = "Area Plot" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            area_img = area_plot.area_plot(mgprofile, self.new_dir)
            self.result = png_result(area_img, self.__generate_about(), display_name, self.block.get_name())        
        
        elif self.block.get_type() == "pca":
            loadings = int(self.block.params["number_of_loadings"])
            display_name = "PCA" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            pca_img = pcoa.pca_plot(mgprofile, self.new_dir, num_of_loadings=loadings)
            self.result = png_result(pca_img, self.__generate_about(), 
                                           display_name, self.block.get_name())            
        
        elif self.block.get_type() == "pcoa":
            dist = self.block.params["distance_metric"]
            display_name = "PCoA: " + dist.capitalize() if self.block.get_name() == self.block.get_type() else self.block.get_name()
            pcoa_img = pcoa.pcoa_plot(mgprofile, self.new_dir, dist_type=dist)
            self.result = png_result(pcoa_img, self.__generate_about(), 
                                           display_name, self.block.get_name())
            
            
    def __plot_dynamic(self):
        """ Creates interactive plots. 
        """
        mgprofile = self.block.metagenomic_profile
        
        if self.block.get_type() == "pca":
            loadings = int(self.block.params["number_of_loadings"])
            display_name = "PCA" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            pca_html, lgd_png = pcoa.pca_plot_interactive(mgprofile, self.new_dir, num_of_loadings=loadings)
            self.result = html_result(pca_html, self.__generate_about(), 
                                            display_name, self.block.get_name(), lgd=lgd_png)
        
        elif self.block.get_type() == "pcoa":
            dist = self.block.params["distance_metric"]
            display_name = "PCoA: " + dist.capitalize() if self.block.get_name() == self.block.get_type() else self.block.get_name()
            pcoa_html, lgd_png = pcoa.pcoa_plot_interactive(mgprofile, self.new_dir, dist_type=dist)
            self.result = html_result(pcoa_html, self.__generate_about(), 
                                            display_name, self.block.get_name(), lgd=lgd_png)
        
        elif self.block.get_type() == "area_plot":
            x_label = ""
            for cls in list(self.block.metagenomic_profile.references.keys()):
                x_label += cls.ljust(80) # add spaces between labels
            display_name = "Area Plot" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            area_html, lgd_png = area_plot.area_plot_interactive(mgprofile, self.new_dir)
            self.result = html_result(area_html, self.__generate_about(), display_name, self.block.get_name(), lgd=lgd_png, 
                                            x_lbl=x_label)
    
    def __perform_enrichment(self):
        """ Performs a user-specified enrichment test on this metagenomic profile.
        """
        mgprofile = self.block.metagenomic_profile
        correction_method = self.block.params["correction"]
        if self.block.params["test"] == "ranksums":
            display_name = "Wilcoxon rank-sum test" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            enrich_table = enrichment.ranksums(mgprofile, self.new_dir, correction=correction_method)
            self.result = table_result(enrich_table, self.__generate_about(), 
                                             display_name, self.block.get_name())
        else:
            display_name = "Student's t-test" if self.block.get_name() == self.block.get_type() else self.block.get_name()
            enrich_table = enrichment.ttest(mgprofile, self.new_dir, correction=correction_method)
            self.result = table_result(enrich_table, self.__generate_about(), 
                                             display_name, self.block.get_name())
    
    def __perform_normalization(self, normalization_type, musicc_intra='use_generic'):
        """ Normalizes the data in the metagenomic profile of this test block.
        """
        mgprofile = self.block.metagenomic_profile
        if normalization_type == "relative":
            normalization.relative_normalization(mgprofile, self.new_dir)
        elif normalization_type == "musicc":
            normalization.musicc_normalization(mgprofile, self.block.gen_params['abundance_data'], 
                                               self.new_dir, musicc_intra=musicc_intra)
    
    # Public methods 
      
    def run_tests(self):
        """ Run tests on this test block according to the specifications in the parameters file.
        
        Returns:
            Result instance for this test run. 
        """
        if "normalization" in self.block.params: # Normalization needs to be performed 1st
            self.__perform_normalization(self.block.params["normalization"])
        if "interactive_plots" in list(self.block.gen_params.keys()) and self.block.gen_params["interactive_plots"]:
            if "static_plots" in list(self.block.gen_params.keys()) and self.block.gen_params["static_plots"]:
                self.__plot_static()
            self.__plot_dynamic()
        else: 
            self.__plot_static()
        if self.block.get_type() == "enrichment":
            self.__perform_enrichment()
            
        return self.result
        
    def get_result(self):
        """ Return the result for this test_runner instance.  
        """
        return self.result