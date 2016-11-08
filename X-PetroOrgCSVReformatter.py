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
    import POrgProcessingModule as POPM
    path  = os.getcwd()+"\\data\\" #example data location
else: #if running in ipython
    scriptlocation = "C:\\Users\\Will\\Dropbox\\Documents\\University\\Edinburgh\\Coding\\Python3\\FTMS\\DataProcessingScripts"
    sys.path.append(scriptlocation)
    import POrgProcessingModule as POPM
    path  = scriptlocation + "\\data\\" #This is the location of our example data. C

#Here we define the input folder location. This directory should contain CSV files output by petroorg. It should not matter if there are other files or directories present, as long as they dont end in .csv
#This section takes the user input to check the correct directory is being examined. In case it isnt, the user is prompted to define the correct data location.
print(path)
isok = input("Is the path above correct for your data - Y or N? ")
if isok.upper() != "Y":
    newpath = input("What is the correct location of your data? ")
    while os.path.isdir(newpath) == False:
        print("The location does not exist at " +str(newpath))
        newpath = input("What is the correct location of your data? ")
    else:
        path = newpath
        
isNa, isK = "n", "n"
ionisationmode = "Negative"#"Positive"

whationisationmode = input("Is this Negative or Positive? ")
while whationisationmode.upper()!= "NEGATIVE" and whationisationmode.upper()!= "POSITIVE":
    print("Please type negative or positive.")
    whationisationmode = input("Is this Negative or Positive? ")
    ionisationmode = "Negative"#"Positive"
else:
    if whationisationmode.upper()!= "NEGATIVE":
        ionisationmode = "Negative"#"Positive"
    elif whationisationmode.upper()!= "POSITIVE":
        ionisationmode = "Positive"
        isadducts =input("Have you included adducts such as Na, K - Y or N? ")
        if isadducts.upper() == "Y":
            isNa = input("Sodium - Y or N? ")
            isK = input("Potassium - Y or N? ")
    
#this are locations used historically

#this function lists the files within the path we defined. It only includes ones ending in ".csv". It then runs the main function for each element in the list.
def fileconverter():
    print("Looking for CSVs in " + path)
    filesA = os.listdir(path)
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
    inputfile = path +"/" +sample+".csv"        #input file names
    data=[]                                     #create an empty list
    with open(inputfile,'r') as csvfile:       #read the input csv line by line
        linereader = csv.reader(csvfile,delimiter=',')
        for row in linereader:
            if row == []:
                continue                        #if a line is empty, ignore it.
            if POPM.intchecker(row[0]):         #if a line begins with an integer, use this line. 
                data.append(row[1:])
    data2=[]                                    #PetroOrg CSVs have an empty column due to a trailing comma on every line. This removes those.
    for x in data:
        if x[-1] == '':
            data2.append(x[:-1])
        else:
            data2.append(x)   
    #here we define the list of headers we want in our output datafile. This list is exhaustive - we trim it down for the main output file. This is designed to work with CHNOS elemental assignments only. 
    #CHNOS and in negative mode. 
            
    newheaders = ["Exp. m/z", "Recal m/z", "Theor. Mass", "Error", "Rel. Abundance", "Signal2Noise", "DBE",
                  'Clet', 'Cno', 'Hlet', 'Hno', 'Nlet', 'Nno', 'Olet', 'Ono', 'Slet', 'Sno', '13C1let', '13C1no']
    #here we allow the addition of adducts
    if isNa.upper() == "Y":
        if isK.upper() != "Y":
            temp = newheaders[:-2]
            temp.extend(["Nalet","Nano",'13C1let', '13C1no'])
        elif isK.upper() == "Y":
            temp = newheaders[:-2]
            temp.extend(["Nalet","Nano","Klet","Kno",'13C1let', '13C1no']) #this assumes that potassium comes after sodium - check!
        newheaders = temp
    elif isK.upper() == "Y":
        temp = newheaders[:-2]
        temp.extend(["Klet","Kno",'13C1let', '13C1no'])
        newheaders = temp
        
    shortlist = newheaders[:3]+[newheaders[4]]
    shortlist = shortlist[:2]+[shortlist[3]]+[shortlist[2]] #shortlist is the no hits.
    shortlistlength = len(shortlist)
    midlist = newheaders[:-2] #mid list is for monoisotopic assignments, removes the 13C columns
    midlistlength = len(midlist)
    isolistlength = len(newheaders) #len(max(data,key=len))-1 #isolist is for the isotopologues - i.e. includes 13C assignments
    sizeB = len(data)
    nohits = np.zeros((sizeB,shortlistlength))
    nohits = pd.DataFrame(nohits,columns=shortlist,dtype=object)
    hits = np.zeros((sizeB,midlistlength))   
    hits = pd.DataFrame(hits,columns=midlist,dtype=object)
    isolist = np.zeros((sizeB,isolistlength))
    isolist = pd.DataFrame(isolist,columns=newheaders,dtype=object) #now we create three output dataframes with the correct headers
    a, b, c = 0, 0, 0 #this is some counting
    for row in data2:
        i,j,k = 0,0,0
        if len(row) ==shortlistlength:
            for x in range(len(row)):
                position=shortlist[i]
                nohits[position][a] = row[x]
                i += 1
            a+= 1
        elif len(row) == isolistlength:
            for x in range(len(row)):
                position=newheaders[j]
                isolist[position][b] = row[x]
                j += 1
            b+= 1
        elif len(row) == midlistlength:
            for x in range(len(row)):
                position=midlist[k]
                hits[position][c] = row[x]
                k += 1
            c+= 1  
    nohits = POPM.cleanupDF(nohits) #this function cleansup the dataframe (removing empty rows, re-sorting, etc)
    isolist = POPM.cleanupDF(isolist)
    hits = POPM.cleanupDF(hits)
    hits = hits.sort_index()
    
    #this section generates a single string for the formula.
    formulae = []
    elementcolumns = [x for x in list(hits.columns.values) if "no" in x]
    formulatorinput = hits.as_matrix(columns=elementcolumns)
    for x in formulatorinput:
        c = int(x[0])
        h = int(x[1])
        n = int(x[2])
        o = int(x[3])
        s = int(x[4])
        if isNa.upper() == "Y":
            na = int(x[5])
            if isK.upper() =="Y":
                k = int(x[6])
        elif isK.upper() =="Y":
            k =int(x[5])
        else:
            na = False
            k = False
        formula = POPM.formulator(c,h,n,o,s,na,k,ionisationmode)
        formulae.append(formula)
    formulae2 = pd.DataFrame(formulae,columns=["Formula"])
    
    #this section generates a single string for the heteroatomic class
    heteroclasA = []
    for x in formulae:
        split = re.split('([H][\d]+)',x)
        heteroclasA.append(split[-1])
    heteroclasA2 = pd.DataFrame(heteroclasA,columns=["HeteroClass"]) 
    
    #this bit groups the resulting dataframes into one
    result = pd.concat([hits,formulae2,heteroclasA2],axis=1)
    result = result.sort_values("Exp. m/z")
    result.drop('Signal2Noise',axis=1,inplace=True) #this drops the SNR column as we dont use it because of how we input to PetroOrg.
    
    POPM.make_sure_path_exists(path +"/OutputCSV/") #this function checks the output directory exists; if it doesnt, it creates it.
    #this then writes the output files 
    nohits.to_csv(path +"/OutputCSV/"+sample+"-nohits.csv")
    result.to_csv(path +"/OutputCSV/"+sample+"-hits.csv")
    isolist.to_csv(path +"/OutputCSV/"+sample+"-isolist.csv")
        
#this line runs the entire script.    
fileconverter()