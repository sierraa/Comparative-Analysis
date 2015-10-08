# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:02:54 2015

@author: Sierra

Module containing several classes representing types of results. A result
can take the form of a table, a static image, or an interactive plot. The final
comparative analysis html page is composed from several instances of result class. 

"""

class abstract_result(object):
    """ Abstract class representing a result from a test. 
    
    Attributes:
        output: path to raw output for this result. 
        meta: information about this result. 
        result_name: specific name of this result. 
        test_name: type of test that generated this result.
    """
    def __init__(self, path, meta, result_name, test_name):
        """Create a new result to display.
        
        Args:
            path: path to file containing raw output for this result. 
            meta: information about this result.
            result_name: specfic name of this result. 
            test_name: type of test that generated this result. 
        """        
        self.output = path
        self.meta = meta
        self.result_name = result_name
        self.test_name = test_name 
        
    def get_output(self):
        """ Return the path to the output for this result. 
        """
        return self.output
        
    def get_result_name(self):
        """ Return name of this result. 
        """
        return self.result_name
        
    def get_meta(self):
        """ Return metadata from this result. 
        """
        return self.meta
        
    def get_result_id(self):
        """ Used for internal linking.
        
        Returns:
            String for internal linking.
        """
        return self.test_name + "_"
        
class png_result(abstract_result):
    """ A result consisting of a png image. 
    """
    def to_html(self):
        """ Returns html formatting for this result. 
        """
        contents = '<div><a name="' + self.get_result_id() + '"></a></div>' # internal linking
        contents += '<p><div class="image"><img src="' + self.get_output() + '" width="900" class="center"></div></p>' 
        contents += '<p class="about">' + self.get_meta() + '</p><br>'

        return contents
        
class html_result(abstract_result):
    """ A result containing an html file, typically an interactive plot. 
    
    Attributes:
        legend: path to png legend for this interactive plot. 
        x_label: string representing the x_label for the interactive plot. 
    """

    def __init__(self, path, meta, result_name, test_name, lgd=None, x_lbl=None):
        """ Create a new html_result to display.
        
        Args:
            lgd (default=None): path to png legend for this interactive plot.
            x_lbl (default=None): string representing the x_label for the interactive plot. 
        """
        abstract_result.__init__(self, path, meta, result_name, test_name)
        self.legend = lgd
        self.x_label = x_lbl
        
    def to_html(self):
        """ Returns html formatting for this result. 
        """
        contents = '<a name="' + self.get_result_id() + '"></a>' # for internal linking
        
        contents += '<div class="interactive">'
        
        for line in open(self.get_output(), 'r'):
            contents += line
            
        contents += '</div>'
        
        if self.x_label != None:
            contents += '<div class="x_label">' + self.x_label.replace(" ", "&nbsp;") + '</div>'
        
        if self.legend != None:
            contents += '<div class="image"><img src="' + self.legend + '" class="center"></div><br>'
            
        contents += '<p class="about">' + self.get_meta() + '</p><br>'
        
        return contents
    
class table_result(abstract_result):
    """ A result containing a table.
    """ 
    
    def to_html(self):
        """ Returns html formatting for this result. 
        """
        contents = '<div><a name="' + self.get_result_id() + '"></a>' # for internal linking
        contents += '<div id="enrich"><table style="width:"900" class="center"><tbody>' 
        
        i = 0        
        for line in open(self.get_output(), 'r'):
            line = line.split('\t')
            contents += '<tr>'            
            if i == 0:
                for word in line:
                    contents += '<th>' + word.capitalize() + '</th>'
            else:
                for word in line:
                    contents += '<td>' + word + '</td>'
            contents += '</tr>'
            
            i += 1
            
        contents += '<tbody></table></div>'
        contents += '<p class="about">' + self.get_meta() + '</p><br>'
        return contents