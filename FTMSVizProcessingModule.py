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
    data=data.rename(columns = {'Cno':'C', 'Hno':'H','Sno':'S','Ono':'O',"Nno":"N","Pno":"P"}) #Rename some columns for ease of use/avoid issues
    data=data.rename(columns = {'Exp. m/z':'mz', 'Recal m/z':'Recalmz','Theor. Mass':'TheorMass','Abundance':'RA',"Rel. Abundance":"RA"}) #Rename some columns for ease of use/avoid issues
    calcheaders = ["AI","AImod","OC","HC","NC","SC","PC","VKsize"] #these are some values we need to calculate.
    RAsum = data.sum()["RA"]
    RAmean = data.mean()["RA"]
    vkfactor = (RAsum/(5*RAmean))/RAmean
    calcdata = pd.DataFrame(columns=calcheaders,index=data.index)
    hetclassintdf = heteroclasstoints(data)
    for x in range(len(data)):
        C = data["C"][x]
        H = data["H"][x]
        N = data["N"][x]
        O = data["O"][x]
        S = data["S"][x]
        try:
            P = data["P"][x]
        except:
            P = 0
        RA = data["RA"][x]
        AI = AIcalc(C,H,N,O,S,P)
        AImod = AImodcalc(C,H,N,O,S,P)
        calcdata["AI"][x]=AI
        calcdata["AImod"][x]=AImod
        calcdata["OC"][x] = float(O)/float(C)
        calcdata["HC"][x] = float(H)/float(C)
        calcdata["NC"][x] = float(N)/float(C)
        calcdata["SC"][x] = float(S)/float(C)
        calcdata["PC"][x] = float(P)/float(C)
        calcdata["VKsize"][x] = areatoradii(float(RA)*vkfactor)
    calcdata = calcdata.apply(pd.to_numeric)
    data = data.join(calcdata)
    #data["HetClassInts"] = hetclassints
    return data, hetclassintdf

# This function reads in an isotopologue hit list
def isocsvreader(inputfile):
    data = pd.read_csv(inputfile,index_col=0)
    data=data.rename(columns = {'Exp. m/z':'mz', 'Recal m/z':'Recalmz','Theor. Mass':'TheorMass','Abundance':'RA',"Rel. Abundance":"RA"}) #Rename some columns for ease of use/avoid issues
    calcheaders = ["AI","AImod","OC","HC","NC","SC","PC","VKsize"] #these are some values we need to calculate.
    vkfactor = (data["RA"].sum()/(5*data["RA"].mean()))/data["RA"].mean()
    calcdata = pd.DataFrame(columns=calcheaders,index=data.index)
    for x in range(len(data)):
        C = data["C"][x]
        H = data["H"][x]
        N = data["N"][x]
        O = data["O"][x]
        S = data["S"][x]
        try:
            P = data["P"][x]
        except:
            P = 0
        RA = data["RA"][x]
        AI = AIcalc(C,H,N,O,S,P)
        AImod = AImodcalc(C,H,N,O,S,P)
        calcdata["AI"][x]=AI
        calcdata["AImod"][x]=AImod
        calcdata["OC"][x] = float(O)/float(C)
        calcdata["HC"][x] = float(H)/float(C)
        calcdata["NC"][x] = float(N)/float(C)
        calcdata["SC"][x] = float(S)/float(C)
        calcdata["PC"][x] = float(P)/float(C)
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
def AIcalc(C,H,N,O,S,P):
    top = 1+C-O-S-(0.5*H)
    btm = C-O-S-N-P
    if btm == 0:
        AI = 0
    else:
        AI = top/btm
    if AI < 0:
        AI = 0
    return AI

#This halves the number of oxygens, assuming only half are counting towards aromaticity. e.g. carboxylic acid. Same reference as above.
def AImodcalc(C,H,N,O,S,P):
    O = O/2
    top = 1+C-O-S-(0.5*H)
    btm = C-O-S-N-P
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
    df = df.fillna(0)
    df = df.loc[~(df==0).all(axis=1)]
    columnlist = list(df.columns.values)
    for x in columnlist:
        df[x] = pd.to_numeric(df[x])
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
def oldformulator(c,h,n,o,s,p,na,k,ionisationmode):
    if ionisationmode =="negative": #this is if you have assigned formula using the provided tool, not petroorg.
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
    if p > 0:
        formula = formula +"P"+str(p)
    return formula

#This new formulator is much simpler 06/01/2017
def porgformulator(df):
    formulae = []
    adductfreeformulae = []
    adducts = []
    for y,x in df.iterrows():
        formula = []
        for i, j in x[7:].iteritems():
            if j > 0:
                formula.append(str(i)+str(int(j)))
        formulae.append(''.join(formula))
    formulae2 = pd.DataFrame(formulae,columns=["Formula with Adducts"])
    for y,x in formulae2.iterrows():
        if 'Na' in str(x.values):
            adduct = "Na"
            adductfreeformula = str(x.values)[2:-5]
        elif 'K' in str(x.values):
            adduct = "K"
            adductfreeformula = str(x.values)[2:-4]
        else:
            adduct = "H"
            adductfreeformula = str(x.values)[2:-2]
        adducts.append(adduct)
        adductfreeformulae.append(adductfreeformula)
    adducts2 = pd.DataFrame(adducts,columns=["Adduct"])
    adductfreeformulae2 = pd.DataFrame(adductfreeformulae,columns=["Formula"])
    return formulae2, adductfreeformulae2, adducts2

#This generates a string for a formula - it takes as input the atomic counts for the IONIC formula, and outputs a neutral molecular ion formula.
#For petrorg assigned formula, it is the neutral mass formula, not the ionic one.
def formulator(c,h,n,o,s,p,na,k,ionisationmode):
    if ionisationmode =="negative": #this is if you have assigned formula using the provided tool, not petroorg.
        h=h+1
    formula = "C"+str(c)+"H"+str(h)
    if n > 0:
        formula = formula + "N"+str(n)
    if o > 0:
        formula = formula + "O"+str(o)
    if s > 0:
        formula = formula + "S"+str(s)
    if p > 0:
        formula = formula + "P"+str(p)
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

#This generates a string for an elemental heteroclass, i.e. "O" covers all O counts, "NS" covers all N S counts.
#If petroorg assigned, it is neutral formula input.
def elementclass(c,h,n,o,s,p,na,k,ionisationmode):
    if ionisationmode =="negative": #this is if you have assigned formula using the provided tool, not petroorg.
        h=h+1
    formula = ''#"C"+"H"
    if n > 0:
        formula = formula + "N"
    if o > 0:
        formula = formula + "O"
    if s > 0:
        formula = formula + "S"
    if p > 0:
        formula = formula + "P"
    """ #this section appends the Na and K count to our formula - however, we want molecular ion formula, not adduct formula.
    if na != False and na > 0:
        if k != False and k > 0:
            formula = formula +"Na"+str(na)+"K"+str(k) #this is probably very unlikely unless you have doubly charged species, which we arent actually dealing with...
        else:
            formula = formula +"Na"+str(na)
    elif k != False and k > 0:
        formula = formula +"K"+str(k)
    """
    if formula == '':
        formula  = "CH"
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
