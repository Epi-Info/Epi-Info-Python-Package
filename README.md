## Epi Info Data Imports and Analyses
This package brings Epi Info statistical analysis routines to the Python environment. It has methods for converting datasets into Python objects and for passing those objects to analysis functions.
### Classes
#### ImportData
Contains functions for importing data sources to Python objects (lists of dictionaries) for consumption by analysis functions.<br><br>
<code>from epiinfo.ImportData import *</code><br>
##### CSV
<code>dataset = eicsv('<i>path/filename</i>.csv')</code><br>
Reads a CSV file and creates a list of dictionaries called 'dataset'.<br>
##### JSON
<code>dataset = eijson('<i>path/filename</i>.json')</code><br>
Reads a JSON file and creates a list of dictionaries called 'dataset'.<br>
##### Sync File
<code>dataset = eisync('<i>path/filename</i>.epi7', '<i>password</i>')</code><br>
Reads a sync file (encrypted XML exported from Epi Info mobile apps) and creates a list of dictionaries called 'dataset'. <i>password</i> is the password that was used to encrypt the data on the mobile device.<br>
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
