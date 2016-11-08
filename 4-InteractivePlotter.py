# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 10:48:04 2016

@author: Will Kew
will.kew@gmail.com

    Copyright Will Kew, 2016

    This file is part of FTMS Visualisation (also known as i-van Krevelen).

    FTMS Visualisation is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FTMS Visualisation is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FTMS Visualisation.  If not, see <http://www.gnu.org/licenses/>.

This script, like 2-StaticPlotter, reads in a pre-formatted CSV to produce scatter plots.
However, this one is designed to work on single files, and produces an HTML file with interactive graphics in it, allowing for better understanding of the data.
This functionality is novel and detailed in a recent letter (in preparation). 

This script relies on the bokeh package
http://bokeh.pydata.org/en/latest/
https://github.com/bokeh/bokeh
"""
from __future__ import print_function    # Python 2 compatibility
from __future__ import absolute_import   # Python 2 compatibility

import os, sys
#import numpy as np
#import pandas as pd
#import matplotlib 
#from matplotlib import cm
from bokeh.embed import file_html
from bokeh.util.browser import view
from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, OpenURL, TapTool, layouts, CustomJS, PrintfTickFormatter, ColorBar, LinearColorMapper, FixedTicker
from bokeh.models.widgets import Panel, Tabs, Button, DataTable, TableColumn, NumberFormatter#, Dropdown
#from bokeh.palettes import Viridis10, Inferno8
import bokeh.palettes as bp
from bokeh.resources import CDN #,JSResources,INLINE
from bokeh.io import reset_output
from jinja2 import Template 

from collections import OrderedDict

"""
# We import also the FTMSVizProcessingModule which contains a few useful functions.
# here we define where the scripts are stored. 
# Make sure to change this to where you have saved these scripts.
"""
try: #test if running in ipython
    __IPYTHON__
except NameError: #if not running in ipython....
    import FTMSVizProcessingModule as FTPM
    path  = os.getcwd()+"data\\" #example data location
else: #if running in ipython
    #homepath
    #scriptlocation = "C:\\Users\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\\FTMS\\DataProcessingScripts"
    #OfficeDesktopPath
    scriptlocation = "F:\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\FTMS\\DataProcessingScripts"
    sys.path.append(scriptlocation)
    import FTMSVizProcessingModule as FTPM
    #OfficeDesktopPath
    path = "F:\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\\FTMS\\DataProcessingScripts\\data\\"
    #HomeDesktopPath
    #path = "C:\\Users\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\\FTMS\\DataProcessingScripts\\data\\"

#define the colourmap of choice
glocmap = list(bp.Viridis10)
glocmap2 = list(bp.Inferno8)

glocmap.reverse()
glocmap2.reverse()

#path = "G:/DATA/FTICRMS/KEW-20160319/data/PORg CSV/"
inputpath = path +"OutputCSV/"
outputpath = path+"OutputHTML/"
FTPM.make_sure_path_exists(outputpath)
files = os.listdir(inputpath)

def intfileplot():
    filesA = os.listdir(inputpath)
    filesB = []
    for y in filesA:
        if y[-8:] =="hits.csv" and y[-10:] != "nohits.csv" and y[-11:] !="isohits.csv":
            filesB.append(y)
    for y in filesB:
        data,hetclassintdf = FTPM.mycsvreader(inputpath+y)
        isodata = FTPM.isocsvreader(inputpath+y[:-8]+"isohits.csv")
        nodata = FTPM.nohitsreader(inputpath+y[:-8]+"nohits.csv")
        intplotter(data,isodata,nodata,y,hetclassintdf)
        reset_output() #cleans up the cache which reduces file size
        
def intplotter(data,isodata,nodata,y,hetclassintdf):
    linewidth = 1.5
    source = ColumnDataSource(data)
    s2 = ColumnDataSource(data=dict(mz=data["mz"],Error=data["Error"],RA=data["RA"],Formula=data["Formula"],HeteroClass=data["HeteroClass"]))    
    isosource = ColumnDataSource(isodata)
    nosource = ColumnDataSource(nodata)
    url = "http://www.chemspider.com/Search.aspx?q=@Formula"
    TOOLS="crosshair,pan,wheel_zoom,box_zoom,reset,tap,previewsave,box_select,poly_select,lasso_select,hover"
    
    figdims = (900,500) #pixel dimensions for the normal figures
    msxlim = [200,700] #x limits in m/z for the mass spectra
    
    vkxlim = [0,1]
    vkylim = [0,2]
    p1 = figure(tools=TOOLS, title=y[:-9]+" - Van Krevelen",width=figdims[0], height=figdims[1], x_axis_label='O/C',y_axis_label='H/C',webgl=True,x_range=vkxlim,y_range=vkylim)
    color_mapper = LinearColorMapper(palette=glocmap, low=msxlim[0], high=msxlim[1])
    p1.scatter(x='OC', y='HC',source=source,size='VKsize', fill_color={'field': 'mz', 'transform': color_mapper}, fill_alpha=0.75,line_color=None) #use size not radius.
    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([('Formula', "@Formula"),('Mass',"@mz{1.11111}"),('Error (ppm)',"@Error{1.11}")])
    taptool = p1.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    
    color_bar = ColorBar(color_mapper=color_mapper, title="m/z", border_line_color=None, location=(0,0), scale_alpha=0.7)
                         #orientation='horizontal',location='top_left', scale_alpha=0.7)#,ticker=FixedTicker(ticks=[2,6,10,14,18]))
    p1.add_layout(color_bar,"right")
    

    dbexlim = [0,45]
    dbeylim = [0,40]
    cmax = max(data["Ono"])
    cmax = int(5 * round(float(cmax)/5))
    p2 = figure(tools=TOOLS, title=y[:-9]+" - DBE vs C# Plot",width=figdims[0], height=figdims[1], x_axis_label='C#',y_axis_label='DBE',webgl=True,x_range=dbexlim,y_range=dbeylim)
    color_mapper2 = LinearColorMapper(palette=glocmap2, low=0, high=cmax)
    p2.scatter(x='Cno', y='DBE',source=source,size='VKsize', fill_color={'field': 'Ono', 'transform': color_mapper2}, fill_alpha=0.75,line_color=None)
    hover = p2.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([('Formula', "@Formula"),('Mass',"@mz{1.11111}"),('Error (ppm)',"@Error{1.11}")])
    taptool = p2.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    
    color_bar2 = ColorBar(color_mapper=color_mapper2, title="O#", border_line_color=None, location=(0,0), scale_alpha=0.7,ticker=FixedTicker(ticks=[0,int(cmax/4),int(cmax/2),int(3*cmax/4),cmax]))
    p2.add_layout(color_bar2,"right")
    
    aixlim=[0,45]
    aiylim= [0,1]

    p3 = figure(tools=TOOLS, title=y[:-9]+" - AI(mod) vs C# Plot",width=figdims[0], height=figdims[1], x_axis_label='C#',y_axis_label='AI(mod)',webgl=True,x_range=aixlim,y_range=aiylim)
    color_mapper3 = LinearColorMapper(palette=glocmap2, low=0, high=cmax)    
    p3.scatter(x='Cno', y='AImod',source=source,size='VKsize', fill_color={'field': 'Ono', 'transform': color_mapper3}, fill_alpha=0.75,line_color=None)
    hover = p3.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([('Formula', "@Formula"),('Mass',"@mz{1.11111}"),('Error (ppm)',"@Error{1.11}")])
    taptool = p3.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    color_bar3 = ColorBar(color_mapper=color_mapper3,title="O#", border_line_color=None, location=(0,0), scale_alpha=0.7,ticker=FixedTicker(ticks=[0,int(cmax/4),int(cmax/2),int(3*cmax/4),cmax]))
                         #orientation='horizontal',location='top_left', scale_alpha=0.7)#,ticker=FixedTicker(ticks=[2,6,10,14,18]))
    p3.add_layout(color_bar3,"right")

    
    p4 = figure(tools=TOOLS, title=y[:-9]+" - Assigned Centroid MS",width=figdims[0], height=figdims[1], x_axis_label='m/z',y_axis_label='Abundance',y_range=[min(data["RA"]),max(data["RA"])],x_range=msxlim)
    p4.segment(x0=0,x1=800,y0=0,y1=0,line_width=1, line_color="black")
    p4.segment(x0='mz',y0=0,x1='mz',y1='RA',source=source,line_width=linewidth, line_color="black")
    p4.scatter(x='mz', y='RA',source=source,fill_color='black',line_color=None)
    hover = p4.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([('Formula', "@Formula"),('Mass',"@mz{1.11111}"),('Error (ppm)',"@Error{1.11}")])
    taptool = p4.select(type=TapTool)
    taptool.callback = OpenURL(url=url) 
    p4.yaxis[0].formatter = PrintfTickFormatter(format="%4.1e")
    
    """
    #this is me trying to plot a barplot of heteroatomic class distributions...
    
    p7 = figure(tools=TOOLS, title=y[:-9]+"",width=800, height=600, x_axis_label='HeteroClass',y_axis_label='Count',webgl=True)
    p7.quad(left="HetClassInts",y=hetclassdf[0],source=source,width=5,height=)
    
    t7 = layouts.Column(hist)
    tab7 = Panel(child=t7,title="test")
    """
    stretch = msxlim[0]*0.1
    p5 = figure(tools=TOOLS, title=y[:-9]+" - Assigned Centroid MS",width=1400, height=600, x_axis_label='m/z',y_axis_label='Abundance', y_range=[min(data["RA"]),max(data["RA"])],x_range=(msxlim[0]-stretch,msxlim[1]+stretch))
    p5.segment(x0=0,x1=800,y0=0,y1=0,line_width=1, line_color="black")
    no1 =p5.segment(x0='mz',y0=0,x1='mz',y1='RA',source=nosource,line_width=linewidth, line_color="red")
    no2 =p5.scatter(x='mz', y='RA',source=nosource,fill_color='red',line_color=None,legend="Unassigned Peaks")
    p5.scatter(x='mz', y='RA',source=source,fill_color='black',line_color=None,legend="Assigned Peaks")
    p5.segment(x0='mz',y0=0,x1='mz',y1='RA',source=source,line_width=linewidth, line_color="black")
    iso1 =p5.segment(x0='mz',y0=0,x1='mz',y1='RA',source=isosource,line_width=linewidth, line_color="green")
    iso2 =p5.scatter(x='mz', y='RA',source=isosource,fill_color='green',line_color=None, legend="Isotologue Peaks")

    hover = p5.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([('Formula', "@Formula"),('Mass',"@mz{1.11111}"),('Error (ppm)',"@Error{1.11}")])
    taptool = p5.select(type=TapTool)
    taptool.callback = OpenURL(url=url) 
    p5.yaxis[0].formatter = PrintfTickFormatter(format="%4.1e")


    js_code1 = "iso1.glyph.visible = false; iso2.glyph.visible = false; no1.glyph.visible = false; no2.glyph.visible = false;"
    cb1 = CustomJS(code=js_code1, args=dict(iso1=iso1,iso2=iso2,no1=no1,no2=no2))
    js_code2 = "iso1.glyph.visible = true; iso2.glyph.visible = true; no1.glyph.visible = true; no2.glyph.visible = true;"
    cb2 = CustomJS(code=js_code2, args=dict(iso1=iso1,iso2=iso2,no1=no1,no2=no2))
        
    toggleOn = Button(label="Hide", button_type="success",callback=cb1)
    toggleOff = Button(label="Show", button_type="success",callback=cb2)    
    
    top = layouts.Row(toggleOn,toggleOff)    
    t3 = layouts.Column(top,p5)
    tab3 = Panel(child=t3,title="Centroid MS with Isotopomers and No Hits")
    
 
    downloadbutton = Button(label="Download", button_type="success")
    downloadbutton.callback = CustomJS(args=dict(s2=s2), code="""
		var data = s2.get('data');
		var filetext = 'mz,Error,RA,Formula,HeteroClass\\n';
		for (i=0; i < data['mz'].length; i++) {
		var currRow = [data['mz'][i].toString(),
                             data['Error'][i].toString(),
        				data['RA'][i].toString(),
        				data['Formula'][i].toString(),
        				data['HeteroClass'][i].toString().concat('\\n')];
		var joined = currRow.join();
		filetext = filetext.concat(joined);
		}
		var filename = 'data_result.csv';
		var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

		//addresses IE
		if (navigator.msSaveBlob) {
			navigator.msSaveBlob(blob, filename);
		}

		else {
			var link = document.createElement("a");
			link = document.createElement('a');
			link.href = URL.createObjectURL(blob);
			link.download = filename;
			link.target = "_blank";
			link.style.visibility = 'hidden';
			link.dispatchEvent(new MouseEvent('click'))
		}
	""")       
 

    columns = [TableColumn(field="mz",title="m/z",formatter=NumberFormatter(format="0.00000")),
                TableColumn(field="Error", title="Error (ppm)",formatter=NumberFormatter(format="0.00")),
                TableColumn(field="RA",title="Abundance"),
                TableColumn(field="Formula",title="Formula"),
                TableColumn(field="HeteroClass",title="Heteroatomic Class")]    
    data_table = DataTable(source=s2, columns=columns, width=1400,row_headers=False,fit_columns=True)
    t4 = layouts.Column(data_table,downloadbutton)
    tab4=Panel(child=t4,title="Selected Data Table")
    
    source.callback = CustomJS(args=dict(s2=s2,dt=data_table), code="""
        var inds = cb_obj.get('selected')['1d'].indices;
        var d1 = cb_obj.get('data');
        var d2 = s2.get('data');
        if (inds.length == 0) {
            d2['mz'] = d1['mz']
            d2['Error'] = d1['Error']
            d2['RA'] = d1['RA']
            d2['Formula'] = d1['Formula']
            d2['HeteroClass'] = d1['HeteroClass']
        }
        else if (inds.length != 0) {
            d2['mz'] = []
            d2['Error'] = []
            d2['RA'] = []
            d2['Formula'] = []
            d2['HeteroClass'] = []
            for (i = 0; i < inds.length; i++) {
                d2['mz'].push(d1['mz'][inds[i]])
                d2['Error'].push(d1['Error'][inds[i]])
                d2['RA'].push(d1['RA'][inds[i]])
                d2['Formula'].push(d1['Formula'][inds[i]])
                d2['HeteroClass'].push(d1['HeteroClass'][inds[i]])
            }
        }
        s2.trigger('change');
        dt.trigger('change');
    """)
    
    """
    hetclasslist = hetclassintdf["HetClass"].tolist()
    hetclasslistnew = []
    for x in hetclasslist:
        hetclasslistnew.append([x,x])
    dropdown = Dropdown(label="Dropdown button", button_type="warning", menu=hetclasslistnew)
    """

    t1 = layouts.Row(p1,p4)
    t2 = layouts.Row(p2,p3)
    t12 = layouts.Column(t1,t2)
    tab1 = Panel(child=t12, title="Main")
    tabs = Tabs(tabs=[ tab1, tab3, tab4])
     
     
    for figs in [p1,p2,p3,p4,p5]:
        figs.xaxis.axis_label_text_font_size = "14pt"
        figs.xaxis.major_label_text_font_size = "14pt"  
        figs.yaxis.axis_label_text_font_size = "14pt"
        figs.yaxis.major_label_text_font_size = "14pt"  
        figs.title.text_font_size = "14pt"
        figs.toolbar_location = "above"
    
    
    with open(outputpath+'templates/index.html', 'r') as f:
        template = Template(f.read())
    
    #js_resources = JSResources(mode='inline')
    html = file_html(tabs,(CDN,CDN),"Interactive Van Krevelen Diagrams",template=template)
    output_file2 = outputpath+y[:-9]+'-plot.html'
    with open(output_file2, 'w') as f:
        f.write(html)
    view(output_file2)
    #show(tabs)  # open a browser
    
intfileplot()
