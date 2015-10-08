# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:24:01 2015

@author: Sierra Anderson

Reads and parses the parameter file. 
"""

# General imports
import os
import sys
import test_block
import string 
import datetime

supported_distance_metrics = ["cityblock", "cosine", "euclidean", "braycurtis",
                              "canberra", "chebyshev", "correlation", "dice",
                              "kulsinki", "mahalanobis", "matching", "minkowski",
                              "rogerstanimoto", "seuclidean", "sokalmichener", 
                              "sokalsneath", "sqeuclidean"]
                              
supported_test_types = ["pcoa", "pca", "area_plot", "enrichment"]

supported_enrichment_tests = ["ttest", "ranksums"]

# Helper methods

def __string_to_dict(s, int_keys=True):
    """Create a dictionary from a string in the form {key1:val1, key2:val2, ...}
    
    Args:
        s: string to be turned into a dictionary
        int_keys (default=True): if True, assumes keys are integer values
    
    Return: 
        dictionary from s 
        
    Raises:
        ValueError: if string is formatted uncorrectly
    """
    d = dict()
    try:    
        l = s.split(',')
        for x in l:
            x = x.strip()
            x = x.strip("{")
            x = x.strip("}")
            xl = x.split(":")
            if int_keys:
                d[int(xl[0])] = xl[1] 
            else:
                d[xl[0]] = xl[1]        
        return d
    except: 
        raise ValueError("Invalid input. Dictionary could not be parsed.")

def __check_general_parameters(gen_params):
    """Checks if the file names are valid and the necessary parameter options
    are included.
    
    Raises:
        NameError: if sample_metadata or abundance_data files are not included
        or if these files or the feature_metadata file is not found at the path specified 
        or if one of the other expected parameters
        is missing from the file (output_directory, metadata_header, open_page, to_html,
        interactive_plots, abundance_sep, metadata_sep). 
    
    """
    if "sample_metadata" not in gen_params:
        print("Error: Metadata file not specified in parameters. Please specify metadata file.")
        sys.exit(0)
    elif not os.path.isfile(gen_params["sample_metadata"]):
        print(("Error: file " + str(gen_params["sample_metadata"]) + " could not be found. Please check that path is correct."))
        sys.exit(0)
    if "abundance_data" not in gen_params:
        print("Error: Abundance file not specified in parameters. Please specify abundance file.")
        sys.exit(0)
    elif not os.path.isfile(gen_params["abundance_data"]):
        print(("Error: file " + str(gen_params["abundance_data"]) + " could not be found. Please check that path is correct."))
        sys.exit(0)
    if "feature_metadata" in gen_params and not os.path.isfile(gen_params["feature_metadata"]):
        print(("Warning: file " + str(gen_params["feature_metadata"]) + " could not be found."))
        
    expected = ["metadata_header", "output_directory", "open_page", "to_html",
                "interactive_plots", "abundance_sep", "metadata_sep", 
                "normalization", "class_names", "metadata_label"]
    
    for e in expected:
        if e not in gen_params: 
            print(("Error: Expected "  + e + " in parameters. Please include " + e + " in parameters file."))
            sys.exit(0)

def __check_test(test_block):
    """ Checks that the test block is created correctly 
    
    Args:
        test_block
    Effects:
        Prints a message to the console if a test cannot be performed.
    Returns:
        True/False
    """
    test_type = test_block.params["test_type"]
    if test_type == "pcoa":
        if "distance_metric" not in test_block.params:
            print("Warning: Could not create PCoA plot. Distance metric not specified.")
            return False
        elif test_block.params["distance_metric"] not in supported_distance_metrics:
            print(("Warning: Could not create PCoA plot. Distance metric '" + test_block.params["distance_metric"] + "' not supported."))
            return False
    if test_type == "pca":
        if "number_of_loadings" not in test_block.params:
            test_block.params["number_of_loadings"] == 0
            return True
        elif test_block.params["number_of_loadings"] < 0 or test_block.params["number_of_loadings"] > 5:
            print(("Warning: Could not create PCA plot. Invalid number of loadings: " + str(test_block.params["number_of_loadings"])))
            return False
    if test_type == "enrichment":
        if "test" not in test_block.params:
            print("Warning: Could not perform enrichment test. Enrichment test not specified.")
            return False
        elif test_block.params["test"] not in supported_enrichment_tests:
            print(("Warning: Could not perform enrichment test. Enrichment test '" + test_block.params["test"] + "' not supported."))
            return False
        if "multiple_hypothesis_correction" not in test_block.params:
            test_block.params["multiple_hypothesis_correction"] = None
    if test_type not in supported_test_types:
        print(("Warning: Unknown test type '" + test_type + "'")) 
        return False
    
    return True

def __create_output_dir():
    """ Return a string representing the output directory for the results.
    """
    today = datetime.datetime.now()
    dirname = "comparative analysis results "
    dirname += "{0}-{1}-{2}-{3}.{4}.{5}".format(today.month, today.day, today.year, today.hour, today.minute, today.second)
    os.mkdir(dirname)    
    return dirname

def __name_tests(result):
    """ Give unique names to tests.
    
    Args:
        result (list): 
    """
    num_of_tests = dict.fromkeys(supported_test_types, 0)
    for tblock in result:
        test_type = tblock.params["test_type"]
        num_of_tests[test_type] += 1 
        suffix = ""
        if num_of_tests[test_type] > 1:
            suffix = "_" + str(num_of_tests[test_type])
        tblock.params["test_name"] = test_type + suffix

# Public methods

def parse(filename):
    """Check parameters passed to the script for correctness.
    
    Args:
        filename: name of the file containing the parameters.

    Raises:
        NameError: if filenames for data are not included or not found.
    
    Returns:
        List containing test_block instances to be run.
    """
    if not os.path.isfile(filename):
        print("Parameters file not found.")
        sys.exit(0)
        
    filename = open(filename, 'rU')
    
    result = list()
    general_parameters = dict()    
    tests = False # flag for when tests in the file are reached 
    line_number = 0 # keep track of line number in file 
    
    for line in filename:
        try:
            if line not in string.whitespace and line[0] != "#" and not tests:
                line = line.rstrip()
                line = line.split("=")
                line[0] = line[0].lower()
                if line[0] == "abundance_sep" or line[0] == "metadata_sep":
                    general_parameters[line[0]] = ',' if line[1].rstrip() == 'csv' else '\t' # default to tab
                elif line[0] == "metadata_header":
                    general_parameters[line[0]] = line[1].rstrip() if line[1].rstrip() == "true" else None
                elif line[0] == "interactive_plots":
                    general_parameters[line[0]] = True if line[1].rstrip().lower()[0] == "t" else False
                elif line[0] == "output_directory":
                    if line[1].rstrip() == "current":
                        general_parameters[line[0]] = os.getcwd()
                    elif os.path.isdir(line[1].rstrip()):
                        general_parameters[line[0]] = line[1].rstrip()
                    else:
                        print(("Error: Could not find directory '" + line[1].rstrip() + ".' Please check that directory is correct."))
                        sys.exit(0)
                    general_parameters[line[0]] += "/" + __create_output_dir()
                elif line[0] == "metadata_label" and line[1].rstrip().lower() == "n/a":
                    general_parameters[line[0]] = None
                elif line[0] != "title" and line[0] != "class_names" and line[0] != "test_type":
                    try:
                        general_parameters[line[0]] = line[1].rstrip().lower()
                    except IndexError:
                        print(("Error: Invalid formatting in parameters file at line " + str(line_number) + ": '" + line[0] + "'"))
                        sys.exit(0)
                elif line[0] == "class_names":
                    try: 
                        general_parameters[line[0]] = __string_to_dict(line[1])
                    except ValueError:
                        try: 
                            general_parameters[line[0]] = __string_to_dict(line[1],int_keys=False)
                        except ValueError:
                            print(("Warning: Class names improperly formated (line " + str(line_number) + ")."))
                            general_parameters[line[0]] = None
                elif line[0] != "test_type":
                    general_parameters[line[0]] = line[1].rstrip()
                
                elif line[0] == "test_type": # Handle the case of the first test in the file
                    test_params = dict()
                    test_params["test_type"] = line[1].lower().rstrip()
                    tests = True
            elif line not in string.whitespace and line[0] != "#" and tests:
                # test parameters parsed here 
                line = line.split("=")
                line[0] = line[0].lower()
            
                if line[0] == "test_type":
                    # add previous test and create a new one
                    tblock = test_block.test_block(test_params)
                    result.append(tblock)            
                    test_params = dict()
                    test_params["test_type"] = line[1].lower().rstrip()
                elif line[0] == "number_of_loadings":
                    try:
                        test_params[line[0]] = int(line[1])
                    except ValueError:
                        print("Warning: Must specify integer value for 'number_of_loadings' in PCA test.")
                        test_params[line[0]] = 0
                elif line[0] == "metadata_label" and line[1] == "n/a":
                    test_params[line[0]] = None
                elif line[0] == "class_names":
                    try: 
                        test_params[line[0]] = __string_to_dict(line[1])
                    except ValueError:
                        try:
                            test_params[line[0]] = __string_to_dict(line[1],int_keys=False)
                        except ValueError:
                            print(("Warning: Class names improperly formated (line " + str(line_number) + ")."))
                else:
                    test_params[line[0]] = line[1].lower().rstrip()
                
            line_number += 1
        except IndexError:
            print(("Invalid formatting in parameters file at line " + str(line_number) + ": '" + line[0] + "'"))
            sys.exit(0)
            
   # After EOF is reached
    try:
        tblock = test_block.test_block(test_params)
        result.append(tblock)
    except UnboundLocalError: 
        print("Error: Empty parameters file. Please check that the correct file was specified.")
        sys.exit(0)
        
    __name_tests(result)    
    
    __check_general_parameters(general_parameters)

    final_result = list()

    for tblock in result:
        if __check_test(tblock):            
            tblock.set_general_parameters(general_parameters)
            final_result.append(tblock)
        
    return final_result, general_parameters