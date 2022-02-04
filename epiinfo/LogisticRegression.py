from scipy.stats import t as tdist
import math
import time
from .randata import randata
from .RegressionUtilities import *

class LogisticRegression:
  def __init__(self):
    self.strataVar = None
    self.mainVar = None
    self.crosstabVar = None
    self.domainVar = None
    self.columnNames = None

    self.currentTable = None

    self.logisticResults = LogisticRegressionResults()

    self.mMatrixLikelihood = EIMatrix()

  def doLogistic(self, inputVariableList, dataTable):
    """ Executes the supporting functions to run the analysis
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
          dataTable (list(dict)): The analysis dataset
        Returns:
          self.logisticResuts (LogisticRegressionResults): This object contains a Rows property.
          It is a list of MeansRow objects, which have properties: Label, Count,
          Mean, StdErr, LCL, and UCL. These are the displayed output of the analysis.
          There is a MeansRow for TOTAL and one for each value of the crosstab
          variable, if present.
    """
    self.currentTable = dataTable
    return self.logisticResults
