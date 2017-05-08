# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 11:13:54 2016

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



This reads in data, and outputs a combined data matrix for statistical analysis. 

Input files will be the output "hits" from my POrgProcessor

Output data in a non-normalised form, and it a fully normalised way. 
Also outputs both with the na columns dropped, and others without.
"""
from __future__ import absolute_import   # Python 2 compatibility

import os, sys
import pandas as pd
import numpy as np

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



inputloc = path+"MergeInput/"
outputloc = path+"MergeOutput/"
FTPM.make_sure_path_exists(outputloc)

files = os.listdir(inputloc)
nfiles = len(files)

samplenames=[]
for x in files:
    samplenames.append(x[:-9])
    

formulae=[]
i = 1
for z in files:
    df1 = pd.read_csv(inputloc+z,index_col=0)
    formula = df1["Formula"]
    formulalist = formula.tolist()
    formulae.append(formulalist)
    formulalist2 = [item for sublist in formulae for item in sublist]
    i = i+1


formulalist = [item for sublist in formulae for item in sublist]

formulalist = list(set(formulalist))
totaluniqueformula = len(formulalist)

indexes = ["mz"]+list(samplenames)

outputdata = pd.DataFrame(index = indexes, columns=formulalist)

dfrand = pd.DataFrame(index = samplenames, columns=formulalist, data=np.random.randint(200000,high=1500000,size=[len(samplenames),totaluniqueformula]))
for y in files:
    df2 = pd.read_csv(inputloc+y,index_col=0)    
    for i in range(len(df2)):
        mass = df2.loc[i]["Theor. Mass"]
        form = df2.loc[i]["Formula"]
        try:
            RA = df2.loc[i]["Rel. Abundance"]
        except:
            RA = df2.loc[i]["Abundance"]
        outputdata.loc[y[:-9]][form] = RA
        if pd.isnull(outputdata.iloc[0][form]):
            outputdata.iloc[0][form] = mass

outputname = "output"
    
outputdata.to_excel(outputloc+outputname+".xlsx")

outputdatadropna = outputdata.dropna(axis=1,how="any",inplace=False)
outputdatadropna.to_excel(outputloc+outputname+"-dropna.xlsx",index_label="Sample")

outputdatadropna50 = outputdata.dropna(axis=1,how="any",inplace=False,thresh=int((len(outputdata)*0.50)+1))
outputdatadropna50.to_excel(outputloc+outputname+"-dropna50.xlsx",index_label="Sample")

outputdatadropna75 = outputdata.dropna(axis=1,how="any",inplace=False,thresh=int((len(outputdata)*0.75)+1))
outputdatadropna75.to_excel(outputloc+outputname+"-dropna75.xlsx",index_label="Sample")

outputdatadropna25 = outputdata.dropna(axis=1,how="any",inplace=False,thresh=int((len(outputdata)*0.25)+1))
outputdatadropna25.to_excel(outputloc+outputname+"-dropna25.xlsx",index_label="Sample")

totalformulacommontoall = len(outputdatadropna.T)
totalformulacommonto50 = len(outputdatadropna50.T)
totalformulacommonto75 = len(outputdatadropna75.T)
totalformulacommonto25 = len(outputdatadropna25.T)


statsdf = pd.DataFrame(index=["Count"],columns=["Unique","Common to 75%","Common to 50%", "Common to 25%","Common to All"])
statsdf["Unique"]["Count"] = totaluniqueformula
statsdf["Common to All"]["Count"] = totalformulacommontoall
statsdf["Common to 50%"]["Count"] = totalformulacommonto50
statsdf["Common to 75%"]["Count"] = totalformulacommonto75
statsdf["Common to 25%"]["Count"] = totalformulacommonto25


statsdf.to_excel(outputloc+outputname+"-statistics.xlsx")


outputnorm = outputdata.div(outputdata.sum(axis=1),axis=0)
outputnorm.to_excel(outputloc+"outputdata-normalised.xlsx",index_label="Sample")

outputnormdropna = outputnorm.dropna(axis=1,how="any",inplace=False)
outputnormdropna.to_excel(outputloc+"outputdata-normalised-dropna.xlsx",index_label="Sample")

outputfillna = outputdata
outputfillna[outputfillna.isnull()] = dfrand[outputfillna.isnull()]
outputfillna.to_excel(outputloc+"outputdata-NAFill.xlsx",index_label="Sample")

outputfillnanorm = outputfillna.div(outputfillna.sum(axis=1),axis=0)
outputfillnanorm.to_excel(outputloc+"outputdata-NAFill-normalised.xlsx",index_label="Sample")


  