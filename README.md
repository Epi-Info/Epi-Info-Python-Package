## Epi Info Analyses
This package brings Epi Info statistical analysis routines to the Python environment. It has methods for converting datasets into Python objects and for passing those objects to analysis functions.
### Classes
#### ImportData
Contains functions for importing data sources to Python objects (lists of dictionaries) for consumption by analysis functions.<br><br>
<code>from epiinfo.ImportData import *</code><br>
##### CSV
<code>dataset = eicsv('<i>path/filename</i>.csv')</code><br>
Reads a CSV file an creates a list of dictionaries called 'dataset'.<br>
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
