**FTMS Visualisation**
----
i-van Krevelen
------------------

**README**

*@author: Will Kew*
*will.kew@gmail.com*


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
	


This set of tools are python scripts for processing and visualising high-resolution FTMS data. 

The tools are easy to use, but will require the user to make minor changes for their system. 

All of the python scripts have been tested in Python 3, though should be largely compatible with Python 2.
The example notebook - 5-Interactive FTMS and Visualisation.ipynb - was originally written for Python 2 and may not work for Python 3.

For all files, you will need to correct the location of the scripts & the input data. An example of this location is included in the files.
File locations must be absolute (Python 3).

The files are as follows:

---------------------
0-FormulaGenerator.py
---------------------	
This script will generate lists of possible ion formulae according to strict rules & predefined elemental limits.
The user will need to check the elemental limits suit their application.
To do this, they will have to open the file in a text editor and change the limits as defined at the end of the file.
The output will be a set of files (known as dictXXX.csv), where XXX is the m/z region calculated. 
By default, the script will generate between 100 and 800 m/z, but this can be adjusted by the user.
Elements CHNOS and adducts K and Na are coded in. Additional elements will require more significant modification by the user.
	

----------------------
1-FormulaAssignment.py
----------------------
This script will assign peaklists molecular formula (by referencing the lists of possible formulae generated above).
This tool has been tested positive and negative mode data, and appears to work OK for a limited set. 
However, THIS TOOL HAS NOT BEEN VALIDATED THOROUGHLY AND WILL NEED RIGOROUS EXAMINATION BY A USER FOR CONFIDENCE ON THEIR OWN SAMPLES.
Specifically, CHO negative mode data can be confidently assigned due to limited possible hits.
CHNOS+Na or +K becomes significantly more difficult. This tool, currently, has no redundancy for assigning the same peak to multiple formulae.
Furthermore, the error limits will need to be tweaked depending on your data set & calibration.
It is provided without any guarantee and simply as a tool which users may wish to build upon for their own interest.

The output of this tool will be three files in "OutputCSV/":
	*XXX-hits.csv*
	*XXX-isohits.csv*
	*XXX-nohits.csv*
Where XXX = sample name. The hits list contains monoisotopic hits, the isohits contains confirmed isotopologues, and the nohits contains the remaining unassigned peaks. 
These are the input files for the remaining numbered scripts. 
No matter what formulae assignment tool you use, you will need to get your data into this format. 
As such, these example files are included for reference.
	
------------------
2-StaticPlotter.py
------------------	
This tool generates a number of van Krevelen, DBE, etc. plots. 
The output are publication quality PNG and EPS files. 
The backend is matplotlib, and it is very customisable for your own uses. 
Without customisation, the images should still be of interest.
The file will prompt you a few questions - do you want titles? do you want hexbin outputs?. These are self-explanatory.
	
-----------------------
3-HeteroClassPlotter.py
-----------------------	
This script leverages Seaborn plotting library to produce a heteroatomic class distribution plot. 
Without modification (aside from the location), the script should work. 
Note that it will build a plot for all the -hits.csv files in the OutputCSV, so for a single sample make sure to have only a single -hits.csv in the output folder, etc.
Examples of other types of plots are included in a commented out section of code.

-----------------------	
4-InteractivePlotter.py
-----------------------
This is the main tool as described in the research paper.
This tool will read in the -hits, -isohits, and -nohits files, and produce an interactive set of plots in a static HTML file.
The output is based on a template, and therefore you should keep the OutputHTML directory as is, unless you remove the template.
This whole script relies on the Bokeh Plotting library. This has been tested on version 0.12.3.
Find more - http://bokeh.pydata.org/en/latest/ 
	
--------------------------	
FTMSVizProcessingModule.py
--------------------------
This is a script full of functions which are used throughout the preceding scripts.
It must be loaded at the start of each script.
In python 3, in ipython, absolute locations are needed hence the start of each script.
In a future version, may be beneficial to recode the whole thing into a package?


---------------------------	
X-PetroOrgCSVReformatter.py
---------------------------
This script is will convert the output (Csv) from PetroOrg (The Florida State University) to formats usable by the plotting scripts.
It takes as input a .csv and outputs the hits, isohits, and nohits.
This script may need user review - it may well require modification depending on the elemental limits the user has used in PetroOrg.
It has been tested with PetroOrg S-10.2 for Negative mode data. It has not been tested in positive mode data. That should work, but it is not proven.
It is no longer actively maintained and is provided for reference only, with no guarantees

-----------------------
Required File Structure
-----------------------
Example input files are provided, along with example output files.
The folder hierarchy should be as follows:

/--- (this is the "scriptlocation" as defined in your script files)
	
	*scripts*
	
	data/--- (this is the data location as defined in your script files)
	
	InputPeaklist/---
	
	*.txt (if you want to do formula assignment)
	
	OutputCSV/--- (output from your formula assignment tool, or from the petroorg reformatter, or from the provided formula assignment tool)
	
		*-hits.csv
	
		*-isohits.csv
		
		*-nohits.csv
	
	OutputHTML/---
	
		templates/---
	
			index.html (template HTML file, example included)
	
		themes.yaml


All other folders and files should be created automatically as you run the scripts. 

If you have any problems, make sure all the paths are correct first.
	