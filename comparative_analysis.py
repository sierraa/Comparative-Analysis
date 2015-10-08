# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:02:13 2015

@author: Sierra Anderson

Main module for comparative analysis.
"""

# General imports 
from sys import argv
script, filename = argv
import webbrowser
import os
import shutil 

# Internal imports 
from check_parameters import parse
import metagenomic_profile as mgp
import generate_html 
import normalization 

def main():
    """ Main script for running comparative analysis from the command line.
    """    
    
    tests, genparams = parse(filename)      
    
    # dictionary to hold tests and test results
    test_results = dict() 
    
    abundance_sep = genparams['abundance_sep']
    metadata_sep = genparams['metadata_sep']
    abundance_data_path = genparams["abundance_data"]
    metadata_path = genparams["sample_metadata"]
    metadata_hdr = genparams["metadata_header"]
    output_dir = genparams["output_directory"]
    lbl = genparams["metadata_label"]
    class_names = genparams["class_names"]
    
    # give each test block a copy of the profile
    
    run_order = list()
    
    master_mp = mgp.metagenomic_profile(abundance_data_path, metadata_path, a_sep=abundance_sep, 
                                     m_sep=metadata_sep, metadata_header=metadata_hdr, metadata_label=lbl,
                                     class_names=class_names)
                                     
    if genparams["normalization"] == "relative":
        normalization.relative_normalization(master_mp, output_dir)
    elif genparams["normalization"] == "musicc":
        normalization.musicc_normalization(master_mp, abundance_data_path, output_dir)        
    
    for tb in tests:
         
        tb_class_names = tb.params["class_names"] if "class_names" in list(tb.params.keys()) else genparams["class_names"] 
        if "metadata_label" in tb.params:
            lbl = tb.params["metadata_label"]
            
        mp = mgp.metagenomic_profile(abundance_data_path, metadata_path, a_sep=abundance_sep, 
                                     m_sep=metadata_sep, metadata_header=metadata_hdr, metadata_label=lbl,
                                     class_names=tb_class_names)
                                     
        if "feature_metadata" in list(genparams.keys()):
            mp.add_feature_metadata(genparams["feature_metadata"])
        tb.set_metagenomic_profile(mp)
        
        # work around for plug-in issues
        if "area_plot" in tb.params:
            run_order.insert(0, tb)
        else:
            run_order.append(tb)
    
    # Copy parameters file to new directory.
    shutil.copyfile(filename, output_dir + "/" + filename)    
    
    os.chdir(output_dir)
    
    # Call run on each test.
    for tb in run_order:
        result = tb.run()
        test_results[tb] = result
    try:
        if genparams["to_html"][0] == "t":    
            results_web_page = generate_html.create_page(test_results, output_dir, filename, order=tests)
            try:            
                if genparams["open_page"][0] == "t":
                    webbrowser.open(results_web_page)
            except IndexError:
                print("Warning: HTML page could not be opened. Value for 'open_page' missing in parameters file.")
    except IndexError:
        print("Warning: HTML page could not be created. Value for 'to_html' missing in parameters file.")
        
    # Done
    print("Tests complete.")
    print(("Output saved at " + genparams["output_directory"]))

if __name__ == "__main__":
    main()