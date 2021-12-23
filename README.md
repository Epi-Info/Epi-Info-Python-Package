## Epi Info Analyses
This package brings Epi Info statistical analysis routines to the Python environment. It has methods for converting datasets into Python objects and for passing those objects to analysis functions.
### Classes
#### ImportData
Contains functions for importing data sources to Python objects (lists of dictionaries) for consumption by analysis functions.<br><br>
<code>from epiinfo.ImportData import *</code><br><br>
##### CSV
<code>dataset = eicsv('<i>path/filename</i>.csv')</code><br>
Reads a CSV file an creates a list of dictionaries called 'dataset'.<br>
#### EICSTables
Contains the ComplexSampleTables class, which contains the functions ComplexSampleTables and ComplexSampleFrequencies.<br><br>
<code>from epiinfo.EICSTables import *</code><br>
<code>csTables = ComplexSampleTables()</code><br>
<code>results = csTables.ComplexSampleTables(ivdict, dataset)</code><br>
<code>results = csTables.ComplexSampleFrequencies(ivdict, dataset)</code><br><br>
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables.<br>
#### EICSMeans
Contains the ComplexSampleMeans class, which contains the ComplexSampleMeans function.<br><br>
<code>from epiinfo.EICSMeans import *</code><br>
<code>csMeans = ComplexSampleMeans()</code><br>
<code>results = csMeans.ComplexSampleMeans(ivdict, dataset)</code><br><br>
'dataset' is the list of dictionaries containing the analysis data. 'ivdict' is a dictionary of analysis variables.<br>
