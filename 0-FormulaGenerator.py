# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 23:11:57 2015

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


Aspects based on 
Seven Golden Rules - doi:  10.1186/1471-2105-8-105

This tool will generate a list of all possible chemical ***IONIC*** formulae within specific constraints, notably:
    High and low mass limits
    Elemental (CHNOS Na K) limits
    Ionisation mode limits
    Seven golden rules, see above

The output are a series of csv files which contain lists of formulae and their exact masses. 
The formulae correspond to singly charged ionic species.
These may be protonated, sodium or potassium adducts (positive mode)
Or deprotonated (negative mode).

The generated lists can be used for:
    a) determining an appropriate error of assignment threshold at a given m/z
    b) automated assignment of peaklists

################
Changelog below
#####
The script defines two functions firstly - getmass and getabun. 
These are used to calculate the exact mass of a given formula & its percentage probability based on existing based on the natural abundance of the constituent isotopes.
So a compound of exactly 1 Carbon-12 will have a 98.9% probability, despite obviously not existing as a "real" molecule, and something C100 would be 33% (98.9% ^ 100)

Version 4 corrects to add the exact mass of an electron - taken from NIST.
The mass of this electron should *removed* from the final mass for a positive mode compound.
OR *added* to the final mass for a negative mode compound. 
Prior to this version, the calculator was producing neutral mass compounds. 

Version 4 is also more strict on elemental ratios (rules 4-6 from Seven Golden Rules, taking the "99.7% common range" rules.)

Version 5 
Adds ability to switch between negative and positive modes. 
Formula constraints tightened for what is known in Scotch Whisky (Negative ESI only shows CHOS, and max S1).
May need adjusting for other ionisation sources and modes. 
Negative mode formula must now have an odd number of protons - this means they have been ionised by loss of a proton. Radicals, etc, will not be determined. 
#####
"""

import numpy as np
import sys, os

"""
# We import also the FTMSVizProcessingModule which contains a few useful functions.
# here we define where the scripts are stored. 
# Make sure to change this to where you have saved these scripts.
"""
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
    
from datetime import datetime

#This is our list of masses - i.e. the center point of each output csv file. For ease/speed, we chunk our lists into 100 m/z blocks.
masses  = [150,250,350,450,550,650,750]
dictionarywindow = 50 # the window is plus or minus this value. So 50 = 100 m/z window. 


#This section checks what ionisation mode you wish to generate a dictionary for. ### Future versions, move this code to the ProcessingModule? -wk
mode = input("Do you want a dictionary of positive or negative mode ions? ")
while mode.lower() != "negative" and mode.lower() != "positive":
    print("Please enter either negative or positive")
    mode = input("Do you want a dictionary of positive or negative mode ions? ")
else:
    if mode.lower() == "negative":
        mode = "negative" 
    elif mode.lower() == "positive":
        mode = "positive"
        
startTime = datetime.now()

print("Elemental formulae limits are coded into the script. Please double check they are suitable for your application.")

#This calculates the mass of a given formulae - for an ION
def getmass(c,h,o,n,s,p,na,k):
    massC = chemdict['C'][0] * c
    massH = chemdict['H'][0] * h
    massO = chemdict['O'][0] * o
    massN = chemdict['N'][0] * n
    massS = chemdict['S'][0] * s
    massP = chemdict['P'][0] * p
    massNa = chemdict['Na'][0] * na
    massK = chemdict['K'][0] * k
    masse = chemdict['e'][0] * 1
    if mode == "negative":
        massTotal = massC + massH + massO + massN + massS +massP + masse
    elif mode == "positive":
        massTotal = massC + massH + massO + massN + massS + +massP + massNa + massK - masse
    return massTotal

#This calculates the natural abundance of a given ion, not yet used much. 
def getabun(c,h,o,n,s): # Dont need to update this section for Sodium as abundance = 100%
    if c > 0:
        abunC = chemdict['C'][1] ** c
    else: 
        abunC = 1

    if h > 0:
        abunH = chemdict['H'][1] ** h
    else: 
        abunH = 1

    if o > 0:
        abunO = chemdict['O'][1] ** o
    else: 
        abunO = 1

    if n > 0:
        abunN = chemdict['N'][1] ** n
    else: 
        abunN = 1

    if s > 0:
        abunS = chemdict['S'][1] ** s
    else: 
        abunS = 1
  
    abunTotal = abunC * abunH * abunO * abunN * abunS
    return abunTotal

#this determines which homologous series/general formula the formula provided fits into, returning a number for simplicity   
def homochecker(o,n,s,p): 
	homo = str(o) + str(n) + str(s) +str(p)
	homoval = o + n + s +p
	return homo, homoval

#checks if the formula assignment makes sense for a positive mode ion - based on max elemental limits of CHONSPKNa, where N rule is strictly followed. 
def pos_adduct_checker(h,n,na,k):
	logicstatement = False
	if n % 2 == 0:
		if h % 2 == 0:
			if na == 1 and k == 0:
				logicstatement = True
			elif na == 0 and k == 1:
				logicstatement = True
		elif h % 2 != 0:
			if na ==0 and k == 0:
				logicstatement = True
	elif n % 2 != 0:
		if h % 2 != 0:
			if na == 1 and k == 0:
				logicstatement = True
			elif na == 0 and k == 1:
				logicstatement = True
		elif h % 2 == 0:
			if na ==0 and k == 0:
				logicstatement = True
	return logicstatement

#Calculates formulae in positive mode for given limits. Includes possibility of Sodium and Potassium adducts. 
def pos_form_calc(maxC, maxH, maxO, maxN, maxS, maxNa, maxK, low, high):
	maxC = min((int(high) / 12), maxC) #here we say the max carbon count has to be the smaller of the total mass/12 or predefined maxC
	maxH = min((maxC * 4), maxH) #max hydrogen count is the smaller of 4 times the number of carbons or the predefined max hydrogen number.
	maxO = min((int(high) / 16), maxO) #here we say the max oxygen count has to be the smaller of the total mass/16 or predefined maxO
	maxN = maxN + 1
	maxS = maxS + 1
	maxNa = maxNa + 1
	maxK = maxK + 1
	allpossformula = []
	allposs = []
	for c in range(int(maxC))[1:]: #obviously our molecules contain at least 1 C and 1 H
		for h in range(int(maxH))[1:]: #Based on seven golden rules, a minimum/maximum H/C Ratio - should be 0.125 to 3.1 for 99.7% of molecules, but inc to 4 to represent ESI adduct possibilities
			hcrat = float(h)/float(c)
			if 0.2 < hcrat < 3.1:
				for o in range(int(maxO)):
					ocrat = float(o)/float(c)
					if ocrat < 1.2:
						for n in range(maxN):
							ncrat = float(n)/float(c)
							if ncrat < 1.3:
								for s in range(maxS):
									scrat = float(s)/float(c)
									if scrat < 0.8:
										for na in range(maxNa):
											for k in range(maxK):
												if pos_adduct_checker(h,n,na,k):
													mass = getmass(c,h,o,n,s,na,k) 
													homo,homoval = homochecker(o,n,s)
													if 0 < float(homoval) < (c*1.3) : #this checker ensures that there are a max 1.3*C heteroatoms. I.e. C10H20O13 is OK, but C10H20O14 is not OK. May need adjustment. 
														if low < mass < high:
															formula = "C%iH%iO%iN%iS%iNa%iK%i" % (c,h,o,n,s,na,k) 
															abundance = getabun(c,h,o,n,s)
															allposs.append([mass,abundance,c,h,o,n,s,na,k,homo,homoval])
															allpossformula.append([formula])
	return allposs   

#Checks N rule (and #H) to see if assigned formula is logical for a negative mode ion. 
def neg_nhchecker(h,n):
	if h % 2 != 0 and n % 2 == 0:
		logicstatement = True
	elif h % 2 == 0 and n % 2 != 0:
		logicstatement = True
	else:
		logicstatement = False
	return logicstatement

#Calculates the negative mode ion for given limits	
def neg_form_calc(maxC, maxH, maxO, maxN, maxS,maxP, low, high):
	maxC = min((int(high) / 12), maxC) #here we say the max carbon count has to be the smaller of the total mass/12 or predefined maxC
	maxH = min((maxC * 4), maxH) #max hydrogen count is the smaller of 4 times the number of carbons or the predefined max hydrogen number.
	maxO = min((int(high) / 16), maxO) #here we say the max oxygen count has to be the smaller of the total mass/16 or predefined maxO
	maxN = maxN + 1
	maxS = maxS + 1
	maxP = maxP + 1
	allpossformula = []
	allposs = []
	for c in range(int(maxC))[1:]: #obviously our molecules contain at least 1 C and 1 H
		for h in range(int(maxH))[1:]: #Based on seven golden rules, a minimum/maximum H/C Ratio - should be 0.125 to 3.1 for 99.7% of molecules, but inc to 4 to represent ESI adduct possibilities
			hcrat = float(h)/float(c)
			if 0.2 < hcrat < 3.1:
				for p in range(maxP):
					for o in range(int(maxO))[1:]:#we want species containing one oxygen. This is a Whisky/SRFA specific requirement, as in negative ESI we must see it.
						ocrat = float(o)/float(c)
						if ocrat < 1.2:
							for n in range(maxN):
								ncrat = float(n)/float(c)
								if ncrat < 1.3:
									if neg_nhchecker(h,n):
										for s in range(maxS):
											scrat = float(s)/float(c)
											if scrat < 0.8: 
												mass = getmass(c,h,o,n,s,p,0,0) 
												homo,homoval = homochecker(o,n,s,p)
												if 0 < float(homoval) < (c*1.3) : #this checker ensures that there are a max 1.3*C heteroatoms. I.e. C10H20O13 is OK, but C10H20O14 is not OK. May need adjustment. 
													if low < mass < high:
														formula = "C%iH%iO%iN%iS%iP%iNa%iK%i" % (c,h,o,n,s,p,0,0) 
														abundance = getabun(c,h,o,n,s)
														allposs.append([mass,abundance,c,h,o,n,s,p,0,0,homo,homoval])
														allpossformula.append([formula])
	return allposs    	


#atomic masses taken from Pure Appl. Chem. 2016; 88(3): 265–291, Atomic Weights of the Elements 2013, doi: 10.1515/pac-2015-0305
#Atomic Masses taken from AME2012 - Chinese Physics C 36 (2012)  1603-2014, Wang, Audi, Wapstra, Kondex, MacCormic, Xu, and Pfeiffer. doi: 10.1088/1674-1137/36/12/003
#isotoptic abundances from Pure Appl. Chem. 2016; 88(3): 293–306, Isotopic compositions of the elements 2013 (IUPAC Technical Report), doi: 10.1515/pac-2015-0503
#electron mass from NIST http://physics.nist.gov/cgi-bin/cuu/Value?meu|search_for=electron+mass
chemdict = {'H':(1.007825, 0.99984),
            'C':(12.000000, 0.98892),
            'N':(14.003074, 0.99634),
            'O':(15.994915, 0.99762),
            'Na':(22.989769, 1.0),
            'P':(30.973763,1.0),
            'S':(31.972071, 0.95041),
            'Cl':(34.968853, 0.75765),
            'K':(38.963706, 0.93258),
            'Br':(78.918338, 0.50686),
            'e':(0.0005485799, 1.0)} 

#################################################
# This section is important.
# Here you define your elemental limits.
# In short - if you have tight limits, you'll get faster results.
# Broad limits - longer calculations and risk of multiple possible assignments in the next script. Not tested fully at high limits.
# CHONS limits were derived from the literature, however the following negative mode limits have been  tightened for the samples our group looks at
# for example, fulvic acids and scotch whisky. 
# See, for example, Kew, W., Goodall, I., Clarke, D., Uhrin, D.
#                   "Chemical Diversity and Complexity of Scotch Whisky as Revealed by High-Resolution Mass Spectrometry" 
#                   J. Am. Soc. Mass Spectrom., 2016. doi:10.1007/s13361-016-1513-y
#################################################

def elementallimits(low,high):
    if mode == "negative":
        if high < 500:
            maxC = 29 #numbers based on paper
            maxH = 72
            maxO = 18
            maxN = 0
            maxS = 2#1 #max S found in Whisky was 1
            maxP = 0
            maxNa = 0
            maxK = 0
        elif 500 <= high <= 1000:
            maxC = 66
            maxH = 126
            maxO = 27
            maxN = 0
            maxS = 2#1 #max S found in Whisky was 1 
            maxP = 0
            maxNa = 0
            maxK = 0
            
    elif mode == "positive": #these numbers are quite broad, based on seven golden rules (etc). on  you may need to tailor to your appication - many possible formulae!
        if high < 500:
            maxC = 29
            maxH = 72
            maxO = 18
            maxN = 0#10
            maxS = 0#7
            maxP = 0
            maxNa = 1
            maxK = 0
        elif 500 <= high < 1000:
            maxC = 66
            maxH = 126
            maxO = 27
            maxN = 0#10
            maxS = 0#4
            maxP = 0
            maxNa = 1
            maxK = 0
    return maxC, maxH, maxO, maxN, maxS,maxP, maxNa, maxK
    
import pandas as pd
#Finally, this bit runs all of the code and saves the output dictionaries as csv files.     
for i in masses:
    low = i - dictionarywindow
    high = i + dictionarywindow
    print("Calculating formulae between " +str(low) + " and " +str(high) + " m/z")
    maxC, maxH, maxO, maxN, maxS,maxP, maxNa, maxK = elementallimits(low,high)
    maxlimstring = "C" +str(maxC) + " H"+str(maxH) + " N"+str(maxN) + " O"+str(maxO) +" S"+str(maxS)+" P"+str(maxP)
    if mode == "negative":
        allposs=neg_form_calc(maxC, maxH, maxO, maxN, maxS,maxP, low,high)
    elif mode == "positive":
        allposs=pos_form_calc(maxC, maxH, maxO, maxN, maxS,maxNa,maxK,low,high)
        maxlimstring = maxlimstring + " Na"+str(maxNa) +" K"+str(maxK)
    x = np.array(allposs)
    x = x[np.argsort(x[:,0])]
    filepathtosave = path+"FormulaDictionaries/"+str(mode[:3])
    df = pd.DataFrame(x,columns=["mass","abundance","C","H","O","N","S","P","Na","K","homo","homoval"])
    FTPM.make_sure_path_exists(filepathtosave) #Makes sure the output directory exists, and creates it if not.
    df.to_csv(filepathtosave+"\\"+"dict"+str(low)+".csv",index=False)
    #np.savetxt(filepathtosave+"\\"+"dict"+str(low)+".csv",x,delimiter=',',fmt="%s")
    if i == masses[-1]:
        print("Max Elemental Limits were:")
        print(maxlimstring)
        print("Time taken to calculate was " + str(datetime.now() - startTime))
