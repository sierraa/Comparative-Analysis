# Comparative Analysis

Comparative analysis is a tool developed to run tests on and generate plots for metagenomic abundance data.

# Getting started

The Github package includes everything needed, including a sample parameters file and sample data. Python 3.x is reccommended, but the package is backwards compatible with Python 2.x. 

## Dependencies 

Required:
<ul> 
<li><a href="http://pandas.pydata.org/">pandas</a> (v0.16.0)</li>
<li><a href="http://matplotlib.org/">matplotlib</a> (v1.4.0)</li>
<li><a href="http://scikit-learn.org/stable/">sklearn</a> (v0.15.2)</li>
<li><a href="http://www.numpy.org/">numpy </a>(v1.9.2)</li>
<li><a href="http://www.scipy.org/">scipy </a>(v0.14.0)</li>
</ul>

Optional:
<ul> 
<li><a href="http://mpld3.github.io/">mpld3</a> (v0.2)
    <ul><li>For interactive plotting</li></ul>
</li>

<li><a href="https://github.com/omanor/MUSiCC">MUSiCC</a></li>
    <ul><li>For MUSiCC normalization</li></ul>
</ul>


## Testing installation

To test that all the packages needed are installed and all the internal files are in the correct place,
run the following command:

        $ python test_ca.py

If everything is set up correctly, it should output the following:

        Test successful. Ready to run. 

Note that the installation of optional packages such as mpld3 are not tested by this script. 

## Files

Two files must be specified in the parameters file. Samples of these are included, using data from
<a href="http://www.nature.com/nature/journal/v490/n7418/abs/nature11450.html">Qin et al. 2012</a>.

### Abundance data

The sample abundance data file has been mapped from the KO-level to the <a href="http://www.genome.jp/kegg/brite.html">BRITE</a>
subcatergory of metabolism. User specified data is expected to be labeled. Either samples as rows and
attributes as columns or vice versa is acceptable. Comma or tab seperated data is acceptable. 

### Metadata

If there is only one column, a header for the metadata is not required, but this fact must be noted via:

        metadata_header=false


# Running the analysis

The analysis is run directly from the command line:

        $ python comparative_analysis.py parameters.sh 

## Parameters file

A default parameters file "parameters.sh" is contained in the Github package. Use this as a template 
to set parameters for your own data. 

## Overriding parameters

A subset of parameters are overridable, meaning they can be defined globally and be applied to all
tests, or they can be specified for individual tests. These can be found in the default parameters 
section of the sample parameters file. 

### Normalization 

For more information about normalization, see the normalization section under "Statistical functions."

### Class names

In the sample parameters file, the following line appears:

        class_names={0:Control, 1:T2D}

This creates a mapping between 0, 1 in the metadata file to the labels "Control" and "T2D." These labels
will now be used in plots and tables generated.
This can be overridden to change which classes are being considered in the metadata file. For instance if we
had a dataset with three classes of samples, we might change it to read

        class_names={0:Control, 2:Pre-Diabetic}

in order to compare a new class to the control class. 

If left blank or not included, names will automatically be generated for plots based on markers in metadata file. 

### Class label

In many datasets, there will be several dimensions of the metadata to consider. Use this option
to specify which dimension in the metadata file should be used to compare samples. The sample dataset
uses only case/control, so

        class_label=n/a

suffices. If multiple metadata columns are included and "n/a" is still specified, the first column
is chosen by default. The first column is also chosen in the case that the label is left blank or not included.

# Output options 

## Output directory

Each time the analysis is run, a new timestamped folder containing a copy of the parameters file is created. 
To save this folder in the current working directory, the parameters file should contain

        output_directory=current

in the general parameters. To specify another location, simply change the line to the desired file path.

## HTML results

The plots and generated files can be displayed in an HTML page titled "results.html" by setting

        to_html=true

Otherwise the results will be saved to file only. If 

        open_page=true

the default browser will attempt to open "results.html" upon the test completion. 

## Static vs dynamic plotting

Interactive plots saved as html files can be created if mpld3 is installed. These plots allow the user
to identify individual samples in the plots, zoom and move the data around the axes. This can be turned on or off via

        interactive_plots=true

in the plotting options section of the parameters file. When set to false, any plots will be saved as .png files.

## Naming tests 

Custom names can be specified for each test by including the following keyword after the test type has been specified:

        test_name=My Test

Test names will be used for file naming and HTML display. Test names must be unique. Note that test names must not include 
characters forbidden in directory names. For instance, on a Windows machine, the characters 

        < > : " / \ | * 

are disallowed. Similarly, "=" has a special meaning in the parameters file and is also reserved. 

# Filtering out samples

If desired, samples can be filtered out by rule or by sample label before a test is performed. Filtering is test-specific,
and must be included under the "test_type=..." keyword. The two filtering methods can be used in tandem if needed. 

## Filtering by label

Filter out specific samples by their label. 

        filter_labels=Sample1,Sample2,...

## Filtering by rule 

Filter out samples using simple logical rules. 

        filter_rules=Rule1,Rule2,...

A "rule" is of the form:

        Attribute Op Value

Where "Attribute" is a column in the metadata matrix, "Op" is one of ">", "<", "is", or "isnot", and Value is some string or number
that each sample will be compared against. Note that filtering by rule is only possible in datasets with multidimensional metadata. 

Example:

        filter_labels=AGE is 18

# Statistical functions

Here are all the statistical functions currently supported. 

## Normalization

As mentioned previously, normalization can be defined globally for all tests in the 
overridable parameters section, or it may defined for a specific test by defining it 
below the test defintion. 

Example:

        test_type=pca
        normalization=relative

There are two options for normalizing abundance data. If the data is pre-normalized or the user wishes
to leave it as it, specify

        normalization=none

### Relative Normalization

Relative abundances are calculated for each sample by summing the total abundances
and then dividing each value by this sum, such that the total of the abundances for the sample
add up to 1. To enable relative normalization, specify

        normalization=relative

in the overrideable parameters section of the parameters file or under a specific test.

### MUSiCC Normalization

Abundance data can also be normalized via <a href="https://github.com/omanor/MUSiCC">MUSiCC</a>. 
More information about MUSiCC can be found <a href="http://www.genomebiology.com/2015/16/1/53">here</a>.
To enable MUSiCC normalization, specify

        normalization=musicc
    
in the overrideable parameters section of the parameters file or under a specific test. Note
that MUSiCC must be installed prior. 

## PCA

<b> Keyword </b>

        test_type=pca

<b> Options </b>

A principle component analysis can be generated without loadings or with up to five loadings by
including the following line under the test definition:

        number_of_loadings=3

PCA is supported for an arbitrary number of sample classes. 

## PCoA

<b> Keyword </b>

        test_type=pcoa

<b> Options </b>

A principle coordinate analysis can be generated using several different distance metrics
by including the following line under the test definition:

        distance_metric=cosine

Current supported distance metrics include: "cityblock", "cosine", "euclidean", "braycurtis", 
"canberra", "chebyshev", "correlation", "dice", "kulsinki", "mahalanobis", "matching", "minkowski",
"rogerstanimoto", "seuclidean", "sokalmichener", "sokalsneath", "sqeuclidean". 
See the <a href="http://docs.scipy.org/doc/scipy/reference/spatial.distance.html">scipy documentation</a> for more information on these metrics. 

PCoA is supported for an arbitrary number of classes. 

## Enrichment

<b> Keyword </b>

        test_type=enrichment

<b> Options </b>

Either a student's t-test or a wilcoxon ranksums test can be performed by including either of the following
keywords under the test definition:

        test=ttest      /*OR*/      test=ranksums

The resulting p-values from an enrichment test can be adjusted by including the following optional keyword:

        correction=bonferroni

Current supported p-value adjustments include Bonferroni correction ("bonferroni") and Benjamini-Hochberg correction
with a false discovery rate of 0.1, 0.05, or 0.01 ("fdr-0.1", "fdr-0.05", and "fdr-0.01" respectively). 

## Area plot

<b> Keyword </b>

        test_type=area_plot

An area plot can only be created from data with a limited number of attributes, the max being 20. 
Note that generating an area plot on data with a large number of attributes will add
significantly to the running time of the analysis. 
