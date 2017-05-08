# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:02:46 2016

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

This script reads in N CSV files (produced by PetroOrg) and outputs a reformatted CSV into a new directory. 
It also separates out monoisotopic assignments from isotopologues and "no hits".

This new CSV is of a more regular and interpretable format, making further data analysis or visualisation possible.

This script can be run interactively in an ipython console - i.e. in Spyder - or directly from the terminal.

Script now should work for Python 2 and 3. 
There is a processing script module which must be imported. Make sure to set the "scriptlocation" variable appropriately if running from ipython. 
If running as an executable, you shouldnt need to.

The PetroOrg assignments can include CHNOS elements. Addition of adducts work in progress - currently K or Na are questioned.
"""
from __future__ import print_function    # Python 2 compatibility
from __future__ import absolute_import   # Python 2 compatibility

#First we need to import a few modules
import csv, os, sys, re
import numpy as np
import pandas as pd


#We import also the POrgProcessingModule which contains a few useful functions.
#here we define where the scripts are stored. Make sure to change this to where you have saved these scripts.
try: #test if running in ipython
    __IPYTHON__
except NameError: #if not running in ipython....
    import FTMSVizProcessingModule as FTPM
    path  = os.getcwd()+"\\data\\" #example data location
else: #if running in ipython
    scriptlocation = "F:/Will/Dropbox/Documents/University/Edinburgh/Coding/Python3/FTMS/FTMSVisToolkit/Scripts/"
    sys.path.append(scriptlocation)
    import FTMSVizProcessingModule as FTPM
    
    path = "F:/Will/Dropbox/Documents/University/Edinburgh/Coding/Python3/FTMS/FTMSVisToolkit/"

    petrorgcsvloc = "POrgCSV/"
    
#Here we define the input folder location. This directory should contain CSV files output by petroorg. It should not matter if there are other files or directories present, as long as they dont end in .csv
#This section takes the user input to check the correct directory is being examined. In case it isnt, the user is prompted to define the correct data location.
print(path+petrorgcsvloc)
isok = input("Is the path above correct for your data - Y or N? ")
if isok.upper() != "Y":
    newpath = input("What is the correct location of your data? ")
    while os.path.isdir(newpath+petrorgcsvloc) == False:
        print("The location does not exist at " +str(newpath+petrorgcsvloc))
        newpath = input("What is the correct location of your data? ")
    else:
        path = newpath+petrorgcsvloc
        

#this function lists the files within the path we defined. It only includes ones ending in ".csv". It then runs the main function for each element in the list.
def fileconverter():
    print("Looking for CSVs in " + path+petrorgcsvloc)
    filesA = os.listdir(path+petrorgcsvloc)
    filesB = []
    for x in filesA:
        if x[-4:] == ".csv":
            filesB.append(x)
    filenumbertotal = str(len(filesB))
    print("There are " + filenumbertotal +" CSVs to process")
    filenumber = 1
    for x in filesB:
        sample = x[:-4]
        porgcsvreader(sample)
        print("Processed file " + str(filenumber) +" of "+ filenumbertotal +".")
        filenumber = filenumber+1

#This is the main body of the script, a function which parses the files and outputs the new files.
def porgcsvreader(sample):
    inputfile = path+petrorgcsvloc +"/" +sample+".csv"        #input file names
    data=[]                                     #create an empty list
    with open(inputfile,'r') as csvfile:       #read the input csv line by line
        linereader = csv.reader(csvfile,delimiter=',')
        for row in linereader:
            if row == []:
                continue                        #if a line is empty, ignore it.
            if FTPM.intchecker(row[0]):         #if a line begins with an integer, use this line. 
                data.append(row[1:])
    data2=[]                                    #PetroOrg CSVs have an empty column due to a trailing comma on every line. This removes those.
    for x in data:
        if x[-1] == '':
            data2.append(x[:-1])
        else:
            data2.append(x)   
            
    #Here we define the headers we need. Elemental headers are added on the fly if necessary, and so we DONT have to limit our elements anymore! 
    newheaders = ["Exp. m/z", "Recal m/z", "Theor. Mass", "Error", "Rel. Abundance", "Signal2Noise", "DBE"]
    nohitsheaders = ['Exp. m/z', 'Recal m/z', 'Theor. Mass', 'Rel. Abundance']
    
    nohits = pd.DataFrame(np.zeros((len(data2),len(nohitsheaders))),columns=nohitsheaders,dtype=object)
    hits = pd.DataFrame(np.zeros((len(data2),len(newheaders))),columns=newheaders,dtype=object)
    isolist = pd.DataFrame(np.zeros((len(data2),len(newheaders))),columns=newheaders,dtype=object) #now we create three output dataframes with the correct headers
    a, b, c = 0, 0, 0 #this is some counting
    for row in data2:
        if "C" not in row: #this will break if you have formula with no carbon. 
            nohits.loc[a,['Exp. m/z', 'Recal m/z', 'Rel. Abundance','Theor. Mass']] = row
            a+= 1
        elif "13C" in row: #this will break if you have non-13C isotopologues, i.e. 18O. PetroOrg doesn't seem to assign those, anyway.
            isolist.loc[b,["Exp. m/z", "Recal m/z", "Theor. Mass", "Error", "Rel. Abundance", "Signal2Noise", "DBE"]] = row[:7]
            for x in range(len(row[7:])):
                if FTPM.intchecker(row[7+x]):
                    isolist.loc[b,row[7+x-1]] = row[7+x]
            b+= 1
        elif "13C" not in row:
            hits.loc[c,["Exp. m/z", "Recal m/z", "Theor. Mass", "Error", "Rel. Abundance", "Signal2Noise", "DBE"]] = row[:7]
            for x in range(len(row[7:])):
                if FTPM.intchecker(row[7+x]):
                    hits.loc[c,row[7+x-1]] = row[7+x]
            c+= 1  
    nohits = FTPM.cleanupDF(nohits) #this function cleansup the dataframe (removing empty rows, re-sorting, etc)
    isolist = FTPM.cleanupDF(isolist)
    hits = FTPM.cleanupDF(hits)
    hits = hits.sort_index()
   
    formulaewithadducts, formulae, adducts = FTPM.porgformulator(hits)
   
    
    #this section generates a single string for the heteroatomic class
    heteroclasA = []
    elementclass = []
    for y in formulae.values:
        x = str(y)[2:-2]
        split = re.split('([H][\d]+)',x)
        heteroclasstemp = split[-1]
        if heteroclasstemp == '':
            heteroclasstemp = "CH"
        elementclasstemp = ''.join([i for i in heteroclasstemp if not i.isdigit()])
        elementclass.append(elementclasstemp)
        heteroclasA.append(heteroclasstemp)
    heteroclasA2 = pd.DataFrame(heteroclasA,columns=["HeteroClass"]) 
    elementclasses = pd.DataFrame(elementclass,columns=["ElementClass"])
    
    #this bit groups the resulting dataframes into one
    result = pd.concat([hits,formulae,formulaewithadducts,adducts,heteroclasA2,elementclasses],axis=1)
    result = result.sort_values("Exp. m/z")
    result.drop('Signal2Noise',axis=1,inplace=True) #this drops the SNR column as we dont use it because of how we input to PetroOrg.
    
    FTPM.make_sure_path_exists(path +"/OutputCSV/") #this function checks the output directory exists; if it doesnt, it creates it.
    #this then writes the output files 
    nohits.to_csv(path +"/OutputCSV/"+sample+"-nohits.csv")
    result.to_csv(path +"/OutputCSV/"+sample+"-hits.csv")
    isolist.to_csv(path +"/OutputCSV/"+sample+"-isohits.csv")
        
#this line runs the entire script.    
fileconverter()