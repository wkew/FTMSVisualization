# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:42:36 2016

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

This script will read in an assigned peaklist (example input file included) and calculate the heteroatomic class distribution.
The output is a vbar plot of heteroamtic class versus count. You can also have the calculated numbers output in a format for replotting. 
This tool uses Seaborn - http://seaborn.pydata.org/

A number of (partially tested) other functions to plot output are included, though commented out. 

This tool was used in our recent paper on Scotch Whisky - https://link.springer.com/article/10.1007/s13361-016-1513-y 
The prompt for the user about whisky samples is thus borne from this - it also serves as an example of how to customise which classes to include. 
"""

from __future__ import print_function    # Python 2 compatibility
from __future__ import absolute_import   # Python 2 compatibility

import os, sys
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

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
    
    
whisky = input("Are these Whisky samples - Y or N?" )
if whisky.upper() == "Y":
    whisky = True
else:
    whisky = False
    
inputpath = path +"OutputCSV/"
outputpath = path + "Images/Classes/"
FTPM.make_sure_path_exists(outputpath) #this function checks the output directory exists; if it doesnt, it creates it.
print("Looking for CSVs in " + inputpath)
filesA = os.listdir(inputpath)
filesB = []
for y in filesA:
    if y[-8:] =="hits.csv" and y[-10:] != "nohits.csv" and y[-11:] !="isohits.csv":
        filesB.append(y)        
nfiles = len(filesB)

samplenames=[]
for x in filesB:
    samplenames.append(x[:-9])
    
heteroclasses=[]
for z in filesB:
    df1 = pd.read_csv(inputpath+z,index_col=0)
    hetclas = df1["HeteroClass"]
    hetclaslist = hetclas.tolist()
    heteroclasses.append(hetclaslist)
    
heteroclasses = [item for sublist in heteroclasses for item in sublist]
hetclasset = list(set(heteroclasses))

indexlist = []
for i in samplenames:
    for n in range(len(hetclasset)):
        indexlist.append(i)
    
###This section is relevant to my whisky samples
if whisky == True:        
    columnnames = ["Sample","Class","WoodType","Region","Age","Peated","HeteroClass","HeteroClassCount"]
    df4 = pd.read_csv(path+"SampleInfo-Dict.csv",index_col=0)
    df4 = df4.T
    dict4 = df4.to_dict()
    
    outputdata = pd.DataFrame(index = range(len(indexlist)), columns=columnnames)
    a = 0
    for y in filesB:
        df2 = pd.read_csv(inputpath+y,index_col=0)
        counter = Counter(df2["HeteroClass"])
        for x in counter:
            outputdata.iloc[a][0] = y[:-9]
            outputdata.iloc[a][1] = dict4[y[:-9]]["Class"]
            outputdata.iloc[a][2] = dict4[y[:-9]]["Total Wood"]
            outputdata.iloc[a][3] = dict4[y[:-9]]["Region"]
            outputdata.iloc[a][4] = dict4[y[:-9]]["Age"]
            outputdata.iloc[a][5] = dict4[y[:-9]]["Peated"]
            outputdata.iloc[a][6] = x  
            outputdata.iloc[a][7] = counter[x]
            a = a+1
    outputdata = outputdata.dropna(how="all",axis=0)
    
else:
    columnnames = ["Sample","Class","HeteroClass","HeteroClassCount"]
    outputdata = pd.DataFrame(index = range(len(indexlist)), columns=columnnames)
    a = 0
    for y in filesB:
        df2 = pd.read_csv(inputpath+y,index_col=0)
        counter = Counter(df2["HeteroClass"])
        for x in counter:
            outputdata.iloc[a][0] = y[:-9]
            outputdata.iloc[a][1] = y[:-9] #this is the Class variable, and should be defined as approrpriate for what you're plotting. In the case of single samples, it can be the sample name.
            outputdata.iloc[a][2] = x  
            outputdata.iloc[a][3] = counter[x]
            a = a+1
    outputdata = outputdata.dropna(how="all",axis=0)
    
pd.to_numeric(outputdata["HeteroClassCount"],errors="raise")

saveoutputdata = input("Do you want to save the output data in a text file for later re-processing - Y or N? ")
if saveoutputdata.upper() == "Y":
    outputdata.to_excel(inputpath+"HetClassByClass-longform.xlsx") #this saves the info out in a longform for plotting.
    
#outputdata = pd.read_excel(inputpath+"HetClassByClass-longform.xlsx") #this reads that data back in. Only necessary for manually re-running bits of script.

# This section creates a unique, naturally sorted list of heteroatom classes for plotting. Only really works for CHO formula. 
# If you have exotic heteroatoms, will need to refigure this yourself, or just hardcode the order you want. easy to do in Excel. 

order = outputdata["HeteroClass"].tolist()
order= list(set(order))                
order.sort(key=FTPM.natural_sort_key) # this natural sort function ensures a logical order to your barplot. 

if whisky == True:
    CHOorder = ["O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16","O17","O18","O19"]
    Fullorder = ["O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16","O17","O18",
    "O19","O1S1","O2S1","O3S1","O4S1","O5S1","O6S1","O7S1","O8S1","O9S1","O10S1","O11S1","O12S1"]
    CHOSorder =["O1S1","O2S1","O3S1","O4S1","O5S1","O6S1","O7S1","O8S1","O9S1","O10S1","O11S1","O12S1"]
    CHOSorderNew = ["O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16","O17","O18","O19","OnS"]
    labels = ["O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16","O17","O18","O19",r'O$\mathregular {_n}$S']

else:
    df = outputdata
    #colours = ["#a6cee3","#1f78b4","#b2df8a"] #colorblind and print friendly colours picked from http://colorbrewer2.org/
    colours = ["#1b9e77","#d95f02","#7570b3"] #as above, but brighter

def barplot():
    sns.set_style("white")
    sns.set_context("paper",font_scale=2)    
    ax = sns.barplot(x="HeteroClass",y="HeteroClassCount",hue="Class",
                     data=outputdata,order=order,palette=sns.color_palette(colours))
    ax.set(xlabel='Heteroatomic Class', ylabel='Count')
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) == 1:
        ax.legend_.remove()
    sns.despine()
    fig = ax.get_figure()
    plt.xticks(rotation=90)
    fig.set_size_inches(8, 6, forward=True)
    fig.savefig(outputpath+"Barplot.png",dpi=600,bbox_inches="tight")    
    fig.savefig(outputpath+"Barplot.eps",dpi=600,bbox_inches="tight")   

barplot() #plots a barplot.

""" 
# Here are some further examples of the Seaborn Plotting library applied to this problem.  
# Most of these rely on having many samples across a small number of classes you wish to compare 
def violinplot():
    sns.set_style("white")
    sns.set_context("paper",font_scale=2)    
    ax = sns.violinplot(x="HeteroClass",y="HeteroClassCount",hue="Class",data=outputdata,
                        order=order,
                        palette=sns.color_palette("bright"),
                        split=False,bw="silverman",scale_hue=True,scale="width",
                        cut=2,linewidth=1.5,inner="quartiles",saturation=1)
    ax.set(xlabel='Heteroatomic Class', ylabel='Count')
    sns.despine()
    fig = ax.get_figure()
    locs, labels = plt.xticks()
    plt.xticks(locs, labels, rotation=90)
    cur_ylim = ax.get_ylim()
    ax.set_ylim(0,cur_ylim[1])
    fig.set_size_inches((POPM.mm2inch(171,80)), forward=True)
    fig.savefig(outputpath+"violinplot-scalewidth.png",dpi=600,bbox_inches="tight")
    fig.savefig(outputpath+"violinplot-scalewidth.eps",dpi=600,bbox_inches="tight")
    
        
def boxplot():
    sns.set_style("white")
    sns.set_context("paper",font_scale=2)    
    ax = sns.boxplot(x="HeteroClass",y="HeteroClassCount",hue="Class",data=outputdata,order=order,palette=sns.color_palette("bright"))
    ax.set(xlabel='Heteroatomic Class', ylabel='Count')
    sns.despine()
    fig = ax.get_figure()
    plt.xticks(rotation=90)
    fig.set_size_inches(8, 6, forward=True)
    fig.savefig(outputpath+"Boxplot-comparison-CHO-only.png",dpi=300,bbox_inches="tight")

def swarmplot():
    sns.set_style("white")
    sns.set_context("paper",font_scale=2)    
    ax = sns.swarmplot(x="HeteroClass",y="HeteroClassCount",hue="Class",data=outputdata,order=order,palette=sns.color_palette("bright"))
    ax.set(xlabel='Heteroatomic Class', ylabel='Average Count')
    sns.despine()
    fig = ax.get_figure()
    plt.xticks(rotation=90)
    fig.set_size_inches(8, 6, forward=True)
    fig.savefig(outputpath+"swarmplot-comparison-CHO-only.png",dpi=300,bbox_inches="tight")
    
def stripplot():
    sns.set_style("white")
    sns.set_context("paper",font_scale=2)    
    ax = sns.stripplot(x="HeteroClass",y="HeteroClassCount",hue="Class",data=outputdata,order=order,palette=sns.color_palette("bright"),jitter=False,split=True)
    ax.set(xlabel='Heteroatomic Class', ylabel='Average Count')
    sns.despine()
    fig = ax.get_figure()
    plt.xticks(rotation=90)
    fig.set_size_inches(8, 6, forward=True)
    fig.savefig(outputpath+"striplot-comparison-CHO-only.png",dpi=300,bbox_inches="tight")
    
""" 
#EOF