# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:36:15 2015

@author: Sierra Anderson

Generate HTML file displaying comparative analysis results.
"""

import datetime
import string 

# Helper methods

# old colors = D4DCE1, ADBBC3

def __get_css(color_1="#CCFF66", color_2="#99CC00"):
    """ Return the CSS for this page.
    
    Args:
        color_1 (default="#CCFF66"): string representing 1st color in this color scheme.
        color_2 (default="#99CC00"): string representing 2nd color in this color scheme.
        
    Returns:
        String holding CSS for the HTML page (stored internally in the .html file).
    """
    contents = '''<style>
    body {
        background-color:''' + color_1 + ''';
        font-family: Verdana, sans-serif;        
    }
    
    h1 {
        /* main heading for the page */
        padding-top: 100px;
    }
    
    h2 {
        /* heading for each test block */
        padding-top: 100px;
    }
    
    table.center {
        border: 5px solid ''' + color_2 +''';
    }
    
    td {
        font-size: 95%;
        padding: 5px;
        border: 1px solid black;
    }
    
    th {
      font-weight: bold;  
    }
    
    img.center {
        display: block;
        margin-left: 3cm;
        border: 5px solid ''' + color_2 + ''';
    }
    
    /* make mpld3 toolbar buttons more visible */
    
    image.mpld3-resetbutton {
        opacity: 0.9 !important;
    }
    
    image.mpld3-zoombutton {
        opacity: 0.9 !important;
    }
    
    image.mpld3-boxzoombutton {
        opacity: 0.9 !important;
    }
    
    image.active {
        opacity: 0.2 !important;
    }
    
    image.pressed {
        opacity: 0.2 !important;
    }
        
    /* for navigation bar links */
    
    a.nav:link {
        padding: 10px;
        font-weight: bold;
        color: black;
        text-decoration: none;
        background-color: ''' + color_2 + ''';
    }
    
    a.nav:visited {
        color: black;
        font-weight: bold;
        text-decoration: none;
    }
    
    a.nav:hover {
        background-color: ''' + color_1 + ''';    
    }
    
    /* general links */    
    
    a:link {
        font-weight: bold;
        text-decoration: none;
        color: black;
    }
    
    a:visited {
        color: black;
        font-weight: bold;
    }
    
    table.links {
        border-spacing: 10px;
        background-color: ''' + color_1 + ''';
        font-color: black;
    }
    
    p.about {
        margin-left: 3cm;
    }
    
    table.center {
        margin-left: 1cm;
        border: 5px solid ''' + color_2 + ''';
        background-color: white;
    }
    
    th.links {
        display: inline-block;
        padding: 12px;
    }
    
    div.x_label {
        margin-left: 5cm;
        font-size: 16px; 
        margin-top: -75px;
        margin-bottom: 50px;
    }
    
    .fixed-nav-bar {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        z-index: 9999;
        width: 100%;
    }
    
    div.interactive {
        display: block;
        width: 50%;
        font-size: 90%;
    }
    </style>
    '''
    return contents 
    
def __create_header(test_blocks, results):
    """ Creates the navigation bar for the html results page.
    
    Args:
        test_blocks: list of test block instances (defines the order).
        results: mapping of test blocks to their results.
        
    Returns:
        String representing the HTML for the navigation bar.
    """
    contents = '<nav class="fixed-nav-bar"><p><table style="width:"700" class="links">'
    contents += '<th class="links"><a href="#top" class="nav">^</a></th>' # back to top button
    
    for test in test_blocks:    
        contents += '<th class="links"><a href="#' + test.get_name() + '" class="nav">' + results[test].get_result_name() + '</a></th>'
    
    contents += '</table></p></nav><div><a name="top"></a></div>'
    return contents
    
def __format_experiment_metadata(gen_params):  
    """ Formats the experiment metadata
    
    Args:
        gen_params: a dictionary of general parameters 
        
    Returns:
        HTML string of formatted experiment metadata.
    """
    contents = "<b>Title:</b> " + gen_params["title"] + "<br>"
    contents += "<b>Name:</b> " + string.capwords(gen_params["name"]) + "<br>"
    contents += "<b>Collaborator(s):</b> " + string.capwords(gen_params["collaborator"]) + "<br>"
    contents += "<b>Date:</b> " + string.capwords(gen_params["date"]) + "<br>"
    contents += "<b>Sequence Type:</b> " + string.capwords(gen_params["sequence_type"]) + "<br>"
    
    current_date = datetime.datetime.now()    
    contents += "<b>Date generated:</b> {0}-{1}-{2} {3}:".format(current_date.month, current_date.day,
            current_date.year, current_date.hour) + "%02d" % current_date.minute + "<br>"
    
    return contents
    
def __get_parameters_info(fname):   
    """ Generates a HTML string containing info about the parameters file.
    
    Args:
        fname (str): path to the parameters file
    """
    return 'Parameters used to generate results can be found <a href=' + fname + ' target="_blank">here.</a><br><br>'
    
# Public methods    
    
def create_page(results, output_dir, parameters_file, order=None):
    """ Create full HTML/CSS for this page. Save result to file. 
    
    Args:
        results: Dictionary mapping test blocks to a list of results for that test block
        output_dir: path to directory to save this HTML page. 
        parameters_file (str): path to parameters file
        order (list[ ]): Order to add tests to page. Defaults to None in which case it orders
        them alphabetically.
        
    Returns:
        path to output HTML file
    """
    
    if order == None:
        order = sorted(list(results.keys()))
        
    contents = '<!DOCTYPE html><html><head>' + __get_css()
    contents += '<title>Results</title></head><body>' + __create_header(order, results)
    
    # Experimenting with table stuff 
    contents += '<table style="width:100%" id="main_table"><tr>'    
    
    contents += '<h1>Comparative Analysis Results</h1></tr><tr>' + __get_parameters_info(parameters_file)
    contents += __format_experiment_metadata(list(results.keys())[0].gen_params) + '</tr>'
    
    for test in order:
        contents += '<tr>'
        contents += '<a name="' + test.get_name() + '"></a><h2>' + results[test].get_result_name() + '</h2>'
        contents += results[test].to_html()
        contents += '</tr>'
        
    contents += '</table></body></html>'
    
    fname = output_dir + "/" + "results.html"    
    
    out = open(fname, 'w')
    out.write(contents)
    out.close()
    
    return fname