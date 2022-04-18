## Epi Info Data Imports and Analyses
This package brings Epi Info statistical analysis routines to the Python environment. It has methods for converting datasets into Python objects and for passing those objects to analysis functions.
### Classes
#### ImportData
Contains functions for importing data sources to Python objects (lists of dictionaries) for consumption by analysis functions.<br><br>
```
from epiinfo.ImportData import *
```
##### CSV
```
dataset = eicsv('path/filename.csv')
```
Reads a CSV file and creates a list of dictionaries called 'dataset'.<br>
##### JSON
```
dataset = eijson('path/filenamejson')
```
Reads a JSON file and creates a list of dictionaries called 'dataset'.<br>
##### Sync File
```
initVector = 'D1041B94D49F66DF120AD40269E102EA'
passwordSalt = '15A4D4258940BD492EE2'
dataset = eisync('path/filename.epi7', initVector, passwordSalt, 'password')
```
Reads a sync file (encrypted XML exported from Epi Info mobile apps) and creates a list of dictionaries called 'dataset'. <i>initVector</i> is the Init Vector that was generated by the mobile device to encrypt the data on the mobile device. <i>passwordSalt</i> is the Salt that was generated by the mobile device to encrypt the data on the mobile device. <i>password</i> is the password that was used to encrypt the data on the mobile device.<br><br>
<i>This function can only be use with data that was encrypted using the Epi Info iOS app's Custom Keys option. See <a href="https://github.com/Epi-Info/Epi-Info-Python-Package/wiki/Mobile-Sync-Files">here</a> for details.</i><br>
#### LogisticRegression
Contains functions for Logistic Regression analysis. In addition to computing identical results to Epi Info Dashboard on identical datasets, this routine also computes interaction odds ratios and confidents limits for interaction effects that are not also main effects.<br>
```
from epiinfo.LogisticRegression import *
ivdict = {'Ill' : 'dependvar'} #, 'intercept' : True , 'includemissing' : False} # intercept and includemissing are optional keys
ivdict['exposureVariables'] = ['column0', 'column1', column2, 'column0*column1']
# ivdict['Group'] = 'matchgroupcolumn' # optional for matched case control studies
reg = LogisticRegression()
rslts = reg.doRegression(ivdict, dataset)
print(rslts.Variables)
print(rslts.Beta)
print(rslts.SE)
print(rslts.OR)
print(rslts.ORLCL)
print(rslts.ORUCL)
print(rslts.Z)
print(rslts.PZ)
print(rslts.Iterations)
print(rslts.MinusTwoLogLikelihood)
print(rslts.CasesIncluded)
print(rslts.Score, rslts.ScoreDF, rslts.ScoreP)
print(rslts.LikelihoodRatio, rslts.LikelihoodRatioDF, rslts.LikelihoodRatioP)
for ior in rslts.InteractionOR:
    print(ior)
```
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables and options.<br>
#### EICSTables
Contains the ComplexSampleTables class, which contains the functions ComplexSampleTables and ComplexSampleFrequencies.<br>
```
from epiinfo.EICSTables import *
csTables = ComplexSampleTables()
results = csTables.ComplexSampleTables(ivdict, dataset)
results = csTables.ComplexSampleFrequencies(ivdict, dataset)
```
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables.<br>
#### EICSMeans
Contains the ComplexSampleMeans class, which contains the ComplexSampleMeans function.<br>
```
from epiinfo.EICSMeans import *
csMeans = ComplexSampleMeans()
results = csMeans.ComplexSampleMeans(ivdict, dataset)
```
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables.<br>
