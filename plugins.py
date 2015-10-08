# -*- coding: utf-8 -*-
"""
Created on Thu Jul 02 13:55:03 2015

@author: Sierra Anderson 

Utitlity module containing several mpld3 plugin classes.

"""
import mpld3.plugins as plugins

class TweakToolbar(plugins.PluginBase):
    """Plugin for changing toolbar.
    """

    JAVASCRIPT = """
    mpld3.register_plugin("tweaktoolbar", TweakToolbar);
    TweakToolbar.prototype = Object.create(mpld3.Plugin.prototype);
    TweakToolbar.prototype.constructor = TweakToolbar;
    function TweakToolbar(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    TweakToolbar.prototype.draw = function(){
      // the toolbar svg doesn't exist
      // yet, so first draw it
      this.fig.toolbar.draw();
      
      // show toolbar
      this.fig.toolbar.buttonsobj.transition(750).attr("y", 0);
      
      // remove event triggers
      this.fig.canvas
        .on("mouseenter", null)
        .on("mouseleave", null)
        .on("touchenter", null)
        .on("touchstart", null);
      
      // then remove the draw function,
      // so that it is not called again
      this.fig.toolbar.draw = function() {}
    }
    """
    def __init__(self):
        self.dict_ = {"type": "tweaktoolbar"}
        
class BarLabelToolTip(plugins.PluginBase):
    """ Plugin for labeling bars in area plot. 
    """
    
    JAVASCRIPT = """
    mpld3.register_plugin("barlabeltoolTip", BarLabelToolTip);
    BarLabelToolTip.prototype = Object.create(mpld3.Plugin.prototype);
    BarLabelToolTip.prototype.constructor = BarLabelToolTip;
    BarLabelToolTip.prototype.requiredProps = ["ids","labels"];
    BarLabelToolTip.prototype.defaultProps = {
        hoffset: 0,
        voffset: 10,
        location: 'mouse'
    };
    function BarLabelToolTip(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    BarLabelToolTip.prototype.draw = function(){
        var svg = d3.select("#" + this.fig.figid);
        var objs = svg.selectAll(".mpld3-path");
        var loc = this.props.location;
        var labels = this.props.labels
        
        test = this.fig.canvas.append("text")
            .text("hello world")
            .style("font-size", 72)
            .style("opacity", 0.5)
            .style("text-anchor", "middle")
            .attr("x", this.fig.width / 2)
            .attr("y", this.fig.height / 2)
            .style("visibility", "hidden");
        
        function mousemove(d) {
            if (loc === "mouse") {
                var pos = d3.mouse(this.fig.canvas.node())
                this.x = pos[0] + this.props.hoffset;
                this.y = pos[1] - this.props.voffset;
            }

            test
                .attr("x", this.x)
                .attr("y", this.y);
        };

        function mouseout(d) {
            test.style("visibility", "hidden")
        };
            
        this.props.ids.forEach(function(id, i) {
            
            
            var obj = mpld3.get_element(id);
            
            function mouseover(d) {
                test.style("visibility", "visible")
                    .style("font-size", 16)
                    .style("opacity", 1)
                    .text(labels[i])
            };
            
            obj.elements().on("mouseover", mouseover.bind(this))
                          
        });
            
       objs.on("mousemove", mousemove.bind(this)) 
           .on("mouseout", mouseout.bind(this));     
        
    }       
    """
    def __init__(self, ids, labels=None, location="mouse"):

        self.dict_ = {"type": "barlabeltoolTip",
                      "ids": ids,
                      "labels": labels,
                      "location": location}