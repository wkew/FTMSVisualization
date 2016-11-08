# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:06:31 2016

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

This file is a holder for a number of key or long functions. 
Many of these are used multiple times in different scripts, and keeping them here allows for easier maintenance.

"""
import os, errno, re, math
import numpy as np
import pandas as pd
from collections import Counter


#this function reads in a csv of the appropriate formatting
#calculates a few extra variables - i.e. H/C or O/C ratio. 
def mycsvreader(inputfile):
    data = pd.read_csv(inputfile,index_col=0)
    data=data.rename(columns = {'Exp. m/z':'mz', 'Recal m/z':'Recalmz','Theor. Mass':'TheorMass','Abundance':'RA'}) #Rename some columns for ease of use/avoid issues
    calcheaders = ["AI","AImod","OC","HC","NC","SC","VKsize"] #these are some values we need to calculate.
    vkfactor = (data["RA"].sum()/(5*data["RA"].mean()))/data["RA"].mean()
    calcdata = pd.DataFrame(columns=calcheaders,index=data.index)
    hetclassintdf = heteroclasstoints(data)
    for x in range(len(data)):
        C = data["Cno"][x]
        H = data["Hno"][x]
        N = data["Nno"][x]
        O = data["Ono"][x]
        S = data["Sno"][x]
        RA = data["RA"][x]
        AI = AIcalc(C,H,N,O,S)
        AImod = AImodcalc(C,H,N,O,S)
        calcdata["AI"][x]=AI
        calcdata["AImod"][x]=AImod
        calcdata["OC"][x] = float(O)/float(C)
        calcdata["HC"][x] = float(H)/float(C)
        calcdata["NC"][x] = float(N)/float(C)
        calcdata["SC"][x] = float(S)/float(C)
        calcdata["VKsize"][x] = areatoradii(float(RA)*vkfactor)
    calcdata = calcdata.apply(pd.to_numeric)
    data = data.join(calcdata)
    #data["HetClassInts"] = hetclassints
    return data, hetclassintdf

# This function reads in an isotopologue hit list    
def isocsvreader(inputfile):
    data = pd.read_csv(inputfile,index_col=0)
    data=data.rename(columns = {'Exp. m/z':'mz', 'Recal m/z':'Recalmz','Theor. Mass':'TheorMass','Abundance':'RA'}) #Rename some columns for ease of use/avoid issues
    calcheaders = ["AI","AImod","OC","HC","NC","SC","VKsize"] #these are some values we need to calculate.
    vkfactor = (data["RA"].sum()/(5*data["RA"].mean()))/data["RA"].mean()
    calcdata = pd.DataFrame(columns=calcheaders,index=data.index)
    for x in range(len(data)):
        C = data["Cno"][x]
        H = data["Hno"][x]
        N = data["Nno"][x]
        O = data["Ono"][x]
        S = data["Sno"][x]
        RA = data["RA"][x]
        AI = AIcalc(C,H,N,O,S)
        AImod = AImodcalc(C,H,N,O,S)
        calcdata["AI"][x]=AI
        calcdata["AImod"][x]=AImod
        calcdata["OC"][x] = float(O)/float(C)
        calcdata["HC"][x] = float(H)/float(C)
        calcdata["NC"][x] = float(N)/float(C)
        calcdata["SC"][x] = float(S)/float(C)
        calcdata["VKsize"][x] = areatoradii(float(RA)*vkfactor)
    calcdata = calcdata.apply(pd.to_numeric)
    data = data.join(calcdata)
    return data

#this function reads in an unassigned peak list
def nohitsreader(inputfile):
    data = pd.read_csv(inputfile,index_col=0)
    data=data.rename(columns = {'Exp. m/z':'mz', 'Recal m/z':'Recalmz','Theor. Mass':'TheorMass','Abundance':'RA'})
    return data
	
#Aromaticity Index calculation. Unlike DBE, this factors in O and S. Taken from DOI: 10.1002/rcm.2386   
def AIcalc(C,H,N,O,S):
    top = 1+C-O-S-(0.5*H)
    btm = C-O-S-N #should also be minus P, but P not an element we consider here
    if btm == 0:
        AI = 0
    else:
        AI = top/btm
    if AI < 0:
        AI = 0
    return AI
    
#This halves the number of oxygens, assuming only half are counting towards aromaticity. e.g. carboxylic acid. Same reference as above.
def AImodcalc(C,H,N,O,S):
    O = O/2
    top = 1+C-O-S-(0.5*H)
    btm = C-O-S-N #should also be minus P, but P not an element we consider here
    if btm == 0:
        AI = 0
    else:
        AI = top/btm
    if AI < 0:
        AI = 0
    return AI

   
def DBEcalc(C,H,N,ionisationmode):
    if ionisationmode == "negative":
        H = H+1
    DBE = (C+1-((H)/2)+(N/2))
    return DBE
    
def areatoradii(x):
    radii = math.sqrt(x/math.pi)
    return radii    

#This sorts heteroatomic class into logical orders.    
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', s)]

#This creates an integer number for each heteroatomic class, useful for specific plotting functions... # Unsure if still necessary    
def heteroclasstoints(data):
    hetclasses = data["HeteroClass"].tolist()
    cnt = Counter(hetclasses)
    cntdf = pd.DataFrame.from_dict(cnt,orient="index")
    cntdflist = cntdf.index.values.tolist()
    cntdflist.sort(key=natural_sort_key)
    cntdf = pd.DataFrame(columns=["HetClass","HetClassCount"],index=range(len(cnt)))
    cntdf["HetClass"] = cntdflist
    cntdf["HetClassInt"] = list(range(len(cnt)))
    hetclasscountlist = []
    for i in range(len(cntdf)):
        hetclasscountlist.append(cnt[cntdf["HetClass"][i]])
    cntdf["HetClassCount"] = hetclasscountlist
    return cntdf

#This determines the heteroclass    
def heteroclass(c,h,n,o,s):
    if o > 0:
        if n > 0:
            if s > 0:
                formula = "N"+str(n)+"O"+str(o)+"S"+str(s)  
            else:
                formula = "N"+str(n)+"O"+str(o)
        elif s > 0:
            formula = "O"+str(o)+"S"+str(s)
        else:            
            formula = "O"+str(o)
    elif n > 0:
        if s > 0:
            formula = "N"+str(n)+"S"+str(s)  
        else:
            formula = "N"+str(n)
    elif s > 0 :
        formula = "S"+str(s)
    else:
        formula = "C"+str(c)+"H"+str(h)
    return formula

#This important function makes sure a given path exists before a user attempts to save a file in a location.    
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
#This checks if a value is an integer - useful for converting petroorg output CSVs to our CSV format.
def intchecker(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

#This cleans up a dataframe
def cleanupDF(df):
    df = df.loc[~(df==0).all(axis=1)]
    droplist = [x for x in list(df.columns.values) if "let" in x]
    for dropA in droplist:
        df = df.drop(dropA,1)
    df = df.sort_values("Exp. m/z")
    return df
    
#This converts milimeters to inches for matplotlib figure sizes    
def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
        
#This function locates the positions for the start and end between certain values for mz.    
def mzfinder(fstart,fend,mz):
    tempmz = np.array(mz)
    idx = (tempmz>fstart)*(tempmz<fend)
    pos = np.where(idx)
    pstart = pos[0][0]
    pend = pos[0][-1]
    return pstart, pend

#This generates a string for a formula - it takes as input the atomic counts for the IONIC formula, and outputs a neutral molecular ion formula.     
def formulator(c,h,n,o,s,na,k,ionisationmode):
    if ionisationmode =="negative":
        h=h+1
    if o > 0:
        if n > 0:
            if s > 0:
                formula = "C"+str(c)+"H"+str(h)+"N"+str(n)+"O"+str(o)+"S"+str(s)  
            else:
                formula = "C"+str(c)+"H"+str(h)+"N"+str(n)+"O"+str(o)
        elif s > 0:
            formula = "C"+str(c)+"H"+str(h)+"O"+str(o)+"S"+str(s)
        else:            
            formula = "C"+str(c)+"H"+str(h)+"O"+str(o)
    elif n > 0:
        if s > 0:
            formula = "C"+str(c)+"H"+str(h)+"N"+str(n)+"S"+str(s)  
        else:
            formula = "C"+str(c)+"H"+str(h)+"N"+str(n)
    elif s > 0 :
        formula = "C"+str(c)+"H"+str(h)+"S"+str(s)
    else:
        formula = "C"+str(c)+"H"+str(h)
        """ #this section appends the Na and K count to our formula - however, we want molecular ion formula, not adduct formula. 
    if na != False and na > 0:
        if k != False and k > 0:
            formula = formula +"Na"+str(na)+"K"+str(k) #this is probably very unlikely unless you have doubly charged species, which we arent actually dealing with...
        else:
            formula = formula +"Na"+str(na)
    elif k != False and k > 0:
        formula = formula +"K"+str(k)
        """
    return formula

#This as above but for isotopes - generates a string for a formula - it takes as input the atomic counts for the IONIC formula, and outputs a neutral molecular ion formula.     
def isotopeformulator(c,h,n,o,s,na,k,ionisationmode,c13,o18):
    #if ionisationmode =="negative":
    #   h=h+1
    if c13 > 0:
        if o > 0:
            if n > 0:
                if s > 0:
                    formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" N"+str(n)+" O"+str(o)+" S"+str(s)  
                else:
                    formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" N"+str(n)+" O"+str(o)
            elif s > 0:
                formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" O"+str(o)+" S"+str(s)
            else:            
                formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" O"+str(o)
        elif n > 0:
            if s > 0:
                formula = "13C"+str(c13)+" 12C"+str(c)+"H"+str(h)+" N"+str(n)+" S"+str(s)  
            else:
                formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" N"+str(n)
        elif s > 0 :
                formula = "13C"+str(c13)+" 12C"+str(c)+" H"+str(h)+" S"+str(s)
        else:
            formula = "13C"+str(c13)+" 12C "+str(c)+" H"+str(h)
        if na != False and na > 0:
            if k != False and k > 0:
                formula = formula +" Na"+str(na)+" K"+str(k) #this is probably very unlikely unless you have doubly charged species, which we arent actually dealing with...
            else:
                formula = formula +" Na"+str(na)
        elif k != False and k > 0:
            formula = formula +" K"+str(k)
    if o18 > 0:
        if o > 0:
            if n > 0:
                if s > 0:
                    formula = "C"+str(c)+" H"+str(h)+" N"+str(n)+" 18O"+str(o18)+" 16O"+str(o)+" S"+str(s)  
                else:
                    formula = "C"+str(c)+" H"+str(h)+"N"+str(n)+" 18O"+str(o18)+" 16O"+str(o)
            elif s > 0:
                formula = "C"+str(c)+" H"+str(h)+" 18O"+str(o18)+" 16O"+str(o)+" S"+str(s)
            else:            
                formula = "C"+str(c)+" H"+str(h)+" 18O"+str(o18)+" 16O"+str(o)
        elif n > 0:
            if s > 0:
                formula = "C"+str(c)+" H"+str(h)+" N"+str(n)+" S"+str(s)  
            else:
                formula = "C"+str(c)+" H"+str(h)+" N"+str(n)
        elif s > 0 :
                formula = "C"+str(c)+" H"+str(h)+str(h)+"S"+str(s)
        else:
            formula = "C"+str(c)+" H"+str(h)
        if na != False and na > 0:
            if k != False and k > 0:
                formula = formula +" Na"+str(na)+" K"+str(k) #this is probably very unlikely unless you have doubly charged species, which we arent actually dealing with...
            else:
                formula = formula +" Na"+str(na)
        elif k != False and k > 0:
            formula = formula +" K"+str(k)
    return formula   
   
#This function splits a string formula into its constituent parts of elemtnal numbers, and then returns the heteroclass
def Form_To_Heteroclass(formula):        
    het = []
    for x in formula:
        c,h,n,o,s = 0,0,0,0,0
        if 'S' not in x:
            split = re.split('(\d+)',x)
            c = split[1]
            h = split[3]
            o = split[5]
        elif 'S' in x:
            split = re.split('(\d+)',x)
            c = split[1]
            h = split[3]
            o = split[5]
            s = split[7]
        heterocl = heteroclass(c,h,n,o,s)
        het.append(heterocl)
    return het
    

