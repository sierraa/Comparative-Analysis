# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:27:34 2015

@author: Sierra Anderson

Class instance representing metagenomic profile containing all data.

"""

import sys

import pandas as pd

class metagenomic_profile(object):
    """ Represents a metagenomic profile containing all data from the experiment.
    
    Attributes:
        abundance_data (pandas.DataFrame): n x m matrix holding data across n samples and m attributes.
        metadata (pandas.DataFrame): n x m matrix holding metadata across n samples and m attributes.
        references (dict[str, list[str]]): maps class types to a list of sample labels in that class.
    """ 
    
    def __init__(self, abundance_data_path, metadata_path, a_sep='\t', m_sep='\t',
                 metadata_header=False, metadata_label=None, class_names=None, 
                 filter_rules=None, filter_labels=None):
        """Create new instance of a metagenomic profile
        
        Args:
            abundance_data_path (str): Path to file containing abundance data
            metadata_path (str): Path to file containing the metadata
            metadata_header (bool): Does the metadata file contain a header? 
                Defaults to False.
            metadata_label (str): Use metadata column with this label for tests. 
                Defaults to None. If None, uses first label in metadata DataFrame.
            class_names (dict[int, str] or dict[str, str]): 
                dictionary mapping metadata values to string names
            a_sep (str): separating character in the abundance data file (ex: '\t', ',').
                Defaults to '\t'.
            m_sep (str): separating character in the metadata file. (ex: '\t', ',').
                Defaults to '\t'.
            filter_rules: list of rules to filter out by, where a rule is a tuple in the 
                form (label, operator, value)
            filter_labels: list of labels to filter out
        """
        self.abundance_data = pd.DataFrame.from_csv(path=abundance_data_path, sep=a_sep)

        if not metadata_header:
            self.metadata = pd.DataFrame.from_csv(path=metadata_path, sep=m_sep, header=None)
        else:
            self.metadata = pd.DataFrame.from_csv(path=metadata_path, sep=m_sep)

        self.__check_abundance_data_shape()
        
        if filter_labels != None:    
            self.__filter_samples_by_name(filter_labels)
        
        if filter_rules != None:
            for rule in filter_rules:
                label, op, val = rule
                self.__filter_samples_by_rule(op, val, label)
                
        if metadata_label != None:
            try: 
                self.metadata.sort(columns=metadata_label, inplace=True)
                self.metadata = self.metadata[metadata_label]
            except KeyError:
                print(("Error: No metadata label '" + metadata_label + "' found. Please check 'metadata_label' option in parameters file."))
                sys.exit(0)
        else:
            self.metadata.sort(columns=self.metadata.columns[0], inplace=True)
            self.metadata = self.metadata[self.metadata.columns[0]]
        
        # keys = class labels, values = labels of class members len(value) = class size
        self.references = dict()
        
        i = 1
        ref_type = self.metadata[0]
        
        if class_names == None:
            self.references["class_" + str(ref_type)] = list()
            self.references["class_" + str(ref_type)].append(self.metadata.index[0])
        else:
            self.references[class_names[ref_type]] = list()
            self.references[class_names[ref_type]].append(self.metadata.index[0])
        while i < len(self.metadata):
            if (self.metadata[i] == ref_type):
                if class_names == None:
                    self.references["class_" + str(ref_type)].append(self.metadata.index[i])
                else:
                    self.references[class_names[ref_type]].append(self.metadata.index[i])
                i += 1
            else:
                ref_type = self.metadata[i]
                if class_names == None:
                    self.references["class_" + str(ref_type)] = list()
                    self.references["class_" + str(ref_type)].append(self.metadata.index[i])
                else:
                    self.references[class_names[ref_type]] = list()
                    self.references[class_names[ref_type]].append(self.metadata.index[i])
                i += 1
                
        self.num_of_classes = len(self.references)
        self.total_sample_count = 0
        
        for k in list(self.references.keys()):
            self.total_sample_count += len(self.references[k])
        
    def __check_abundance_data_shape(self):
        """ Test modules expect the abundance data to have samples as rows
        and attributes as columns. Transforms abundance data if this is not
        the case. 
        """
        if self.abundance_data.index[0] not in self.metadata.index:
            self.abundance_data = self.abundance_data.T
    
    def __filter_samples_by_name(self, names):
        """ Filter out samples according to label.
        
        Args:
            names: list of sample labels to filter out 
        """
        self.abundance_data.drop(names)
        self.metadata.drop(names)
        
    def __filter_samples_by_rule(self, op, value, label=None):
        """ Filter out samples according to rule. 
        
        Args:
            op: String representing the operator filtering samples ("=", "!=", ">", "<")
            label: class label to reference in metadata, defaults to the first column if set to None 
            value: value to compare against using operator 
        """
        
        # iterate over metadata, checking rule        
        for sample in self.metadata.index:
            lbl = self.metadata.columns[0] if label == None else label
            this_val = self.metadata[lbl][sample]
            boolean = False
            if op == "=":
                boolean = this_val == value
            elif op == "!=":
                boolean = this_val != value
            elif op == ">":
                boolean = this_val > value 
            elif op == "<":
                boolean = this_val < value
            if boolean:
                self.abundance_data.drop(sample)
                self.metadata.drop(sample)
        
    def add_feature_metadata(self, path, sp='\t'):
        """ Add a third dataframe to this profile, holding feature metadata (as opposed to 
        sample metadata).
        
        Args:
            path (string): path to file containing feature metadata
            sp (string): delimiting character in file. Defaults to '\t'.
        """
        self.feature_metadata = pd.DataFrame.from_csv(path, sep=sp)
    
    def to_file_abundance_data(self, filename, output_dir, separator="\t"):
        """ Write abundance data to file.
        
        Args:
            filename: name of the output file
            output_dir: directory where the output file should be saved
            separator (default="\t"): separating character for the out file 
        """
        if output_dir != "current":
            filename = output_dir + "//" + filename
        
        self.abundance_data.to_csv(path_or_buf=filename, sep=separator)
    
    def set_abundance_data(self, dataframe):
        """ Sets a new abundance data matrix for this metagenomic profile.
        
        Args:
            dataframe (pandas.DataFrame): new abundance data matrix
        """
        self.abundance_data = dataframe
        self.__check_abundance_data_shape()
        
    def set_metadata(self, dataframe):
        """ Sets a new metadata matrix for this metagenomic profile.
        
        Args:
            dataframe (pandas.DataFrame): new metadata matrix
        """
        self.metadata = dataframe