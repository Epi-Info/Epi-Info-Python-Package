## Epi Info Data Imports and Analyses
This package brings Epi Info statistical analysis routines, plus a log binomial regression analysis routine for computing adjusted relative risks of associations with binary outcome data, to the Python environment. It has methods for converting datasets into Python objects and for passing those objects to analysis functions.
#### Current version available in dist folder
```
pip install path/epiinfo-1.0.3-py3-none-any.whl
```
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
```
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables and options.<br>
```
rslts.show()
```
```
LOGISTIC REGRESSION RESULTS
Variable          Coefficient    Standard Error    Odds Ratio    Lower   Upper  
ChefSalad            1.145           0.3429          3.1424      1.6046  6.1539 
EggSaladSandwich     1.0418          0.3146          2.8343       1.53   5.2506 
CONSTANT            -0.7644          0.3602     

Number of Iterations             4
Minus Two Log-Likelihood    393.37
Number of Observations         309

Fit Test             Value      DF       P
Score              14.7777       2  0.0006
Likelihood Ratio   15.5999       2  0.0004
```
Individual result functions provide greater precision, as well as meaningful odds ratios for interaction terms.<br>
```
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
```
['ChefSalad', 'EggSaladSandwich', 'CONSTANT']
[1.1449910849486071, 1.041809442441627, -0.7644353230444337]
[0.34290905726619875, 0.31456430262135265, 0.3601766334566801]
[3.1424133419378486, 2.834340954032848]
[1.6046436032972575, 1.5300107079664744]
[6.153865937145245, 5.250609425070682]
[3.3390517418143193, 3.311912488988536, -2.1223901053991483]
[0.0008406490630090172, 0.0009266052827040479, 0.03380499507237897]
4
393.3736390464802
309
14.777703590233484 2 0.0006
15.599932822362803 2 0.0004
```
Results are from the fictional Salmonellosis dataset often used in Epi Info training sessions.<br>
#### LogBinomialRegression
Contains functions for Log-Binomial Regression analysis to provide adjusted risk ratios when the independent variable is binary.<br>
```
from epiinfo.LogBinomialRegression import *
ivdict = {'Ill' : 'dependvar'} #, 'intercept' : True , 'includemissing' : False} # intercept and includemissing are optional keys
ivdict['exposureVariables'] = ['column0', 'column1', column2, 'column0*column1']
#ivl['StartValues'] = [0.0, 0.0, 0.0, 0.0, -0.45] # Optional user-supplied starting beta values
reg = LogBinomialRegression()
rslts = reg.doRegression(ivdict, dataset)
```
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables and options.<br>
```
rslts.show()
```
```
LOG BINOMIAL REGRESSION RESULTS
Variable          Coefficient    Standard Error    Risk Ratio    Lower   Upper  
EggSaladSandwich     0.287           0.0863          1.3325      1.1251  1.578  
ChefSalad            0.3359          0.1156          1.3992      1.1155  1.755  
CONSTANT            -0.8579          0.1275     

Number of Iterations             7
Log-Likelihood             -197.81
Number of Observations         309
```
Individual result functions provide greater precision, as well as parameter iteration history and meaningful risk ratios for interaction terms.<br>
```
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
for ior in rslts.InteractionOR:
    print(ior)
# To see start values and subsequent beta parameter values in the iteration history.
print()
for ph in rslts.ParameterHistory:
    print(ph)
```
```
['EggSaladSandwich', 'ChefSalad', 'CONSTANT']
[0.28702656278093225, 0.33587378874771817, -0.8578882747952187]
[0.08629620218172233, 0.11560278700279132, 0.127540265427839]
[1.3324596068382382, 1.399162423625005]
[1.1251193436636948, 1.1154930317062357]
[1.5780091364122923, 1.7549688183079077]
[3.3260625094081475, 2.90541255497247, -6.726411238971455]
[0.0008808217459920762, 0.003667693292066409, 1.7389864495598986e-11]
7
-197.80804698840117
309

[0.015412272950627526, 0.018694620739870826, -0.1434745513345579]
[0.047175480510386114, 0.053313801454389464, -0.2716280030939807]
[0.12025288102095905, 0.13274826883584623, -0.4755035906059664]
[0.22117400493838188, 0.24761248723857526, -0.7059060377582919]
[0.2784656944889591, 0.3214894897321413, -0.8352770236314679]
[0.2868833735148343, 0.33549848896238826, -0.857374553660365]
[0.28702656278093225, 0.33587378874771817, -0.8578882747952187]
```
Results are from the fictional Salmonellosis dataset often used in Epi Info training sessions.<br>
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
