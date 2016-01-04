# ------------------EXPERIMENT METADATA------------------ 

# Data to be displayed on html page.
Title=T2D Dataset
Name=One Person
Collaborator=Another Person
Date=2015
Sequence_type=Illumina

# ------------------GENERAL PARAMETERS------------------

# Path to the abundance data
Abundance_data=brite_subcat_normalized.tab

# Path to the sample metadata, i.e. sample labels vs class
Sample_metadata=Age_metadata.tab

# Add path to feature metadata here if it exists
# Feature_metadata=Not_a_real_file.tab

# Does the sample metadata contain a header for the columns?
# Options: true/false
Metadata_header=true

# Separating character for the abundance data
# Options: tab/csv
abundance_sep=tab

# Separating character for the metadata
# Options: tab/csv
metadata_sep=tab

# Path to the output directory where results should be saved.
# NOTE: Each test will be put in a separate folder within a unique
# folder for this run. 
# Options: 'current' for the current working directory or custom file path
output_directory=current

# -------------- DEFAULT PARAMETERS --------------

# Defined here as defaults, include these parameters under the
# 'test_type' field if an override is desired

# How should the data be normalized?
# Options: None, MUSiCC, relative
normalization=relative

# Class names for the indicators in the metadata column being used.
class_names={0:Control, 1:T2D}

# Label of the column to use in the metadata to define classes
# if not specified, first column is used by default.
class_label=n/a

# ------------------PLOTTING OPTIONS------------------

# Create a html page displaying the results of the analysis?
# Options: true/false
To_html=True

# If to_html=True, open html page in default browser when tests are complete?
# Options: true/false
Open_page=True

# Create interactive plots? 
# Options: true/false
# If false, plots will be png images
Interactive_plots=True

#If interactive_plots=True, generate .png plots as well? 
Static_plots=True

# --------------------TESTS----------------------------

# The type of test to be performed
# Options: area_plot, pcoa, pca, enrichment
test_type=area_plot
# Unique display name for this test (optional)
test_name=Area Plot

test_type=pcoa
# Distance metric to be used in plotting PCoA
distance_metric=Cosine
test_name=Filtered PCoA
# Filter out samples by sample labels (optional)
filter_labels=DLM007,CON-050,T2D-097

test_type=pcoa
distance_metric=chebyshev
test_name=Chebyshev PCoA

test_type=pca
test_name=Filtered PCA
# Filter out samples by rule (optional)
# Operators: 'is','isnot','<','>'
# Separate multiple rules using a comma
filter_rules=AGE isnot 45

test_type=pca
# Number of loadings to be plotted in PCA
# Options: an integer between 0 and 5 (inclusive)
Number_of_loadings=3
test_name=PCA with loadings

test_type=enrichment
# Type of enrichment test to be performed
# Options: ttest (Student's t-test), ranksums (Wilcoxon Ranksums Test)
test=ttest
test_name=FDR-corrected t-test
# Type of correction to be used on p-values (optional)
# Options: bonferroni, FDR-0.1, FDR-0.05, FDR-0.01
correction=FDR-0.1

test_type=enrichment
# Type of enrichment test to be performed
# Options: ttest (Student's t-test), ranksums (Wilcoxon Ranksums Test)
test=ranksums
test_name=Bonferroni-corrected ranksums test
# Type of correction to be used on p-values (optional)
# Options: bonferroni, FDR-0.1, FDR-0.05, FDR-0.01
correction=bonferroni
