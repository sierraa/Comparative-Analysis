# -*- coding: utf-8 -*-
"""
Created on Wed Oct 07 15:21:55 2015

@author: Sierra Anderson

Test that everything needed to run the comparative analysis is installed and
all the internal files are in the correct directory. 

"""

def main():
    modules = ["matplotlib", "sklearn", "numpy", "scipy", "pandas"]
    internal_files = ["area_plot", "check_parameters", "comparative_analysis", 
                      "pcoa", "enrichment", "test_block", "test_runner", "normalization",
                      "metagenomic_profile", "result", "generate_html"]
    success = True
    
    for m in modules:
        try:        
            __import__(m)
        except ImportError:
            print(("Could not import " + m + ".\n"))
            success = False
            
    for m in internal_files:
        try:
            __import__(m)
        except ImportError:
            print(("Could not find " + m + ".py in directory.\n"))
            success = False
            
    if success:
        print("Test successful. Ready to run.")
    
    
if __name__ == '__main__':
    main()