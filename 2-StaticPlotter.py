# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:20:10 2016

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

This script will read in a CSV of assigned formulae from a mass spectrum and plot a variety of scatter and hexbin plots.
This uses the matplotlib backend.
Examples include Van Krevelens and DBE vs C# plots, but the script is flexible to any 2D scatter plot.
The plots also encode colour and size, allowing up to four dimensions of information to be visualised on a single figure.

The script should be fairly self-evident how to customise to your tastes.
An example input file is included with this distribution - your input files must have the same format.

The script will ask you to confirm file locations, and then if you want titles, and if you want hexbin output as well as normal output. 

This tool was used in our recent paper on Scotch Whisky - https://link.springer.com/article/10.1007/s13361-016-1513-y 
"""
from __future__ import print_function    # Python 2 compatibility
from __future__ import absolute_import   # Python 2 compatibility


#Import a few modules
import os, sys
import matplotlib.pyplot as plt
from matplotlib import cm

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
   # path = "C:\\Users\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\\FTMS\\DataProcessingScripts\\data\\"


    
"""
Here we define the location of the input CSV(s) in the right format for us to use. 
These will with the headers = [,"Exp. m/z","Recal m/z","Theor. Mass","Error","Rel. Abundance","DBE","Cno","Hno","Nno","Ono","Sno","Formula","HeteroClass"]
Empty first column for the pandas index.
n rows for assigned formulae
Here we define the input folder location. This directory should contain CSV files output by petroorg. It should not matter if there are other files or directories present, as long as they dont end in .csv
#This section takes the user input to check the correct directory is being examined. In case it isnt, the user is prompted to define the correct data location.
"""
print(path)
isok = input("Is the path above correct for your data - Y or N? ")
if isok.upper() != "Y":
    newpath = input("What is the correct location of your data? ")
    while os.path.isdir(newpath) == False:
        print("The location does not exist at " +str(newpath))
        newpath = input("What is the correct location of your data? ")
    else:
        path = newpath


"""
This is the choice of colourmap. You should always select an appropriate colour scheme for figures.
see https://www.youtube.com/watch?v=xAoljeRJ3lU
"""
#glocmap = cm.inferno_r
glocmap = cm.viridis_r
#glocmap = cm.plasma_r
#glocmap = cm.magma_r


titlelogic = input("Do you want titles on the pictures - Y or N? ")
if titlelogic.upper() == "Y":
    titlelogic = True #Do you want titles? True or False
else:
    titlelogic = False
    
hexbinlogic = input("Do you want hexbin outputs - Y or N? ")
if hexbinlogic.upper() == "Y":
    hexbinlogic = True
else:
    hexbinlogic = False

#This function lists all the files within our directory, and parses out only the "hits.csv" files.
#It then reads in the data using the mycsvreader function.
#Then it passes the data for each sample in turn to the plotter.
def fileplotter():
    inputpath = path +"OutputCSV/"
    print("Looking for CSVs in " + inputpath)
    filesA = os.listdir(inputpath)
    filesB = []
    for y in filesA:
        if y[-8:] =="hits.csv" and y[-10:] != "nohits.csv" and y[-11:] !="isohits.csv":
            filesB.append(y)
    filenumbertotal = str(len(filesB))
    print("There are " + filenumbertotal +" CSVs to process")
    filenumber = 1
    for y in filesB:
        data,hetclassintdf = FTPM.mycsvreader(inputpath+y) # this reads the data using a csv reader function in the processing module.
        #this next bit passes the appropriate values to the plotting functions.
        produceplots(y[:-9],data["RA"], data["DBE"], data["AI"],
                     data["OC"], data["HC"], data["NC"], data["SC"], 
                     data["mz"],data["Error"], data["Cno"], data["Ono"])
        print("Processed file " + str(filenumber) +" of "+ filenumbertotal +".")
        filenumber = filenumber+1


#this function lets us quickly change the plot type - they're all constructed in the same way, with a scatter plot, but the choice of X, Y, labels, colour scale, etc can be made here.
def produceplots(sample,RA, DBE, AI, OC, HC, NC, SC, mz,error, Cno, Ono):
    RAnorm = [(float(i)/sum(RA))*30000 for i in RA] #this normalises the relative abundance to a better number for the sizes. Note that the 30000 is an arbitrary "fudge factor".
    #Plotter(X,Y,title,size,Cvar,CvarLab,Xlabel,Ylabel,xlim,ylim,clim,location,plottype,gridsize)
    #Below: X = O/C ratio, Y = H/C Ratio, X limit = 0-1.5, y limit = 0-2.5, colour limit = 0-800.
    Plotter(OC,HC,sample+" - Van Krevelen by mz",RAnorm,mz,"Mass","O/C","H/C",[0,1.5],[0,2.5],[0,800],"VanK","scatter",0)  #this makes a scatter van krevelen with mass for colour and abundance for size.
    Plotter(Cno,DBE,sample+" - DBE vs Carbon Number",RAnorm,Ono,"O #","C Number","DBE",[0,50],[0,35],[0,20],"DBE","scatter",0)  #this makes a DBE versus C# scatter with oxygen number for colour and abundance for size.
    if hexbinlogic:    
        Plotter(OC,HC,sample+" - Van Krevelen",0,0,"Count","O/C","H/C",[0,1.5],[0,2.5],35,"Hex-VanK","hexbin",(20,15)) #this makes a hexbin van krevelen
        Plotter(Cno,DBE,sample+" - DBE vs Carbon Number",0,0,"Count","C Number","DBE",[0,50],[0,35],40,"Hex-DBE","hexbin",(20,15)) #this makes a hexbin dbe plot.
                                        

#This takes areguments from "produceplots" and generates our images.     
#Arguments are: X (list/1D array for X axis), Y (same), title(string title),size(1D array variable to size the scatter points by),Cvar(1D array var for colour scale),
#CvarLab(String label for colourbar), XLabel(String),Ylabel(string),Xlim(list with two elements for x min and max),Ylim(same),clim(same for colour bar),
#location(string output path),plottype(string plot type, e.g. scatter, hexbin),gridsize(tuple for hexbin gridsize variable).
def Plotter(X,Y,title,size,Cvar,CvarLab,Xlabel,Ylabel,xlim,ylim,clim,location,plottype,gridsize):
    labelsize = 20 #this defines the text size of the labels for Xlabel and Ylabel
    widthmm, heightmm = 171, 233 #this defines the figure size in milimeters (useful for publication image size control)
    ratio = (widthmm/2.0)/(heightmm/4.0) #this calculates a ratio of the figure size. The 2x4 value is arbitrary but works
    widthinch = widthmm/25.40  #converts the value to inches as required by matplotlib
    fig = plt.figure(figsize=(widthinch,widthinch/ratio)) #this creates the figure object with the appropriate size
    ax1 = fig.add_subplot(111) #this adds a single subplot.
    #here we define if it is a scatter or hexbin plot we want. Possible to add further types should need be.
    if plottype == "scatter":
        ax1 = plt.scatter(X,Y,s=size,c=Cvar,cmap=glocmap,alpha=0.75,edgecolor='k',linewidth=0) #alpha is transparency and set at 75%. Edgecolour is black but we use a zero linewidth so there is no edge.
    elif plottype == "hexbin":
        #ax1 = plt.hexbin(X,Y,cmap=glocmap,bins="log",gridsize=gridsize,alpha=0.75,mincnt=1,vmax=math.log10(clim)) #this plot the hexbin on a log scale
        #below: mincnt=1 means zero values are white, not coloured. Vmax is clim, and linewidth is also zero.
        ax1 = plt.hexbin(X,Y,cmap=glocmap,bins=None,gridsize=gridsize,alpha=0.75,mincnt=1,vmax=clim,linewidth=0) #this plots the hexbin with mostly default variables
    plt.xlabel(Xlabel,size=labelsize) #add the xlabel
    plt.ylabel(Ylabel,size=labelsize) #add the ylabel
    imagefolder = "Images/"
    if titlelogic == True: #If true, the figures have titles and go into a title images directory.
        imagefolder = "TitleImages/"
        plt.title(title)
    axes = plt.gca() 
    axes.set_xlim(xlim) #this sets the xlim
    axes.set_ylim(ylim) #this sets the y limit
    plt.locator_params(axis='y',nbins=4) #this defines the number of ticks on y axis
    plt.locator_params(axis='x',nbins=5) #this defines the number of ticks on x axis
    plt.tick_params(axis="both",which="major",left="on",bottom="on",top="off",right="off",labelsize=labelsize) #this sets the major ticks and sizes.
    plt.tick_params(axis="both",which="minor",left="on",bottom="on",top="off",right="off") #this sets the minor ticks and sizes 
    axes.get_yaxis().set_tick_params(which='both', direction='out',width=1.25,length=8) #this sets tick size and direction y axis
    axes.get_xaxis().set_tick_params(which='both', direction='out',width=1.25,length=8) #this sets tick size and direction x axis
    
    for tick in axes.xaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize) #this would set a size, but it is unused.
        tick.label1.set_fontweight('bold') #this makes the tick labels bold on x axis
    for tick in axes.yaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold') #this makes the tick labels bold on y axis
    
    if plottype == "scatter": #sets a couple more variables for scatter plot colourbar
        cb = plt.colorbar(ax1)
        cb.set_label(CvarLab,size=labelsize) 
        cb.ax.tick_params(labelsize=labelsize) 
        cb.set_alpha(1)
        cb.set_clim(clim)
        cb.draw_all()
    elif plottype =="hexbin": #sets a couple more variables for scatter plot colourbar - currently identical to scatter, but kept separate in case of need to change.
        cb = plt.colorbar(ax1)
        cb.set_label(CvarLab,size=labelsize)
        cb.ax.tick_params(labelsize=labelsize) 
        cb.set_alpha(1)
        cb.draw_all()
    plt.tight_layout() #this makes sure the entire figure fits within the predefined image size.
    FTPM.make_sure_path_exists(path+imagefolder+location+"/") #Makes sure the output directory exists, and creates it if not.
    plt.savefig(path+imagefolder+location+"/"+title+".png",dpi=600) #saves the PNG (raster) at high DPI
    plt.savefig(path+imagefolder+location+"/"+title+".eps",dpi=600) #saves an EPS (vector) at high DPI.
    plt.show() #shows the image (useful if running in an IPython console.)
    
  
#This calls the entire script to start.
fileplotter()
