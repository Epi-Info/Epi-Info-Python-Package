from scipy.stats import t as tdist
import math
import time
from .randata import randata
from .RegressionUtilities import *

class LogisticRegression:
  def __init__(self):
    self.mstrC = None
    self.mdblC = None
    self.mdblP = None
    self.mlngIter = None
    self.mdblConv = None
    self.mdblToler = None
    self.mboolIntercept = None
    self.mstraBoolean = None
    self.mstrMatchVar = None
    self.mstrWeightVar = None
    self.mstrDependVar = None
    self.mstrGroupVar = None
    self.mstraTerms = None
    self.mStrADiscrete = None
    self.terms = None
    self.discrete = None
    self.lintIntercept = None
    self.lintweight = None
    self.NumRows = None
    self.NumColumns = None

    self.currentTable = None
    self.lStrAVarNames = None

    self.logisticResults = LogisticRegressionResults()

    self.mMatrixLikelihood = EIMatrix()

  def zFromP(self, _p):
    """ Computs a Z statistic for a P value
        Parameters:
          _p (float): p-value
        Returns: float
    """
    ZFP = 0.0
    P0 = -0.322232431088
    P2 = -0.342242088547
    P3 = -0.0204231210245
    P4 = -4.53642210148E-05
    Q0 = 0.099348462606
    Q1 = 0.588581570495
    Q2 = 0.531103462366
    Q3 = 0.10353775285
    Q4 = 0.0038560700634
    F = _p
    if F >= 1.0:
      return ZFP
    if F >= 0.5:
      F = 1.0 - F
    if F == 0.5:
      return ZFP
    T = (math.log(1 / F ** 2.0)) ** 0.5
    T = T + ((((T * P4 + P3) * T + P2) * T - 1) * T + P0) / ((((T * Q4 + Q3) * T+ Q2) * T + Q1) * T + Q0)
    if _p >= 0.5:
      ZFP = -T
    else:
      ZFP = T
    return ZFP

  def createSettings(self, inputVariableList):
    """ Sets initial values of properties
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
        Returns: none
    """
    self.mstrC = "95"
    self.mdblC = 0.05
    self.mdblP = self.zFromP(self.mdblC * 0.5)
    self.mlngIter = 15
    self.mdblConv = 0.00001
    self.mdblToler = 0.000001
    self.mboolIntercept = True
    self.mstraBoolean = ['No', 'Yes', 'Missing']
    self.mstrMatchVar = ""
    self.mstrWeightVar = ""
    self.mstrDependVar = ""
    self.mstraTerms = []
    self.mStrADiscrete = []
    self.terms = 0
    self.discrete = 0

    for key in inputVariableList:
      if str(inputVariableList[key]).lower() == "term":
        self.mstraTerms.append(key)
        self.terms += 1
      if str(inputVariableList[key]).lower() == "discrete":
        self.mStrADiscrete.append(key)
        self.discrete += 1
        self.mstraTerms.append(key)
        self.terms += 1
      if str(inputVariableList[key]).lower() == "matchvar":
        self.mstrMatchVar = key
      if str(inputVariableList[key]).lower() == "weightvar":
        self.mstrWeightVar = key
      if str(inputVariableList[key]).lower() == "dependvar":
        self.mstrDependVar = key
      if str(inputVariableList[key]).lower() == "groupvar":
        self.mstrGroupVar = key
      if str(inputVariableList[key]).lower() == "intercept":
        self.mboolIntercept = inputVariableList[key]
      if str(inputVariableList[key]).lower() == "p":
        self.mdblP = float(inputVariableList[key])
        self.mdblC = 1.0 - self.mdblP
        self.mstrC = str(self.mdblP * 100.0)
        self.mdblP = self.zFromP(mdblC * 0.5)
      if str(inputVariableList[key]).lower() == "unsorted":
        self.mStrADiscrete.append(key)
        self.discrete += 1
        self.mstraTerms.append(key)
        self.terms += 1

  def removeRecordsWithNulls(self, currentTableMA):
    """ Removes records having null values in analysis variables
        Parameters:
          currentTableMA (list of dictionaries)
        Returns: none
    """
    i = len(currentTableMA) - 1
    while i >= 0:
      rowcopy = currentTableMA[i]
      for k in rowcopy:
        if str(rowcopy[k]).lower() in ['none','null','(null)','.',''] or str(rowcopy[k]).isspace():
          del currentTableMA[i]
          break
      i -= 1

  def outcomeOneZero(self, currentTableMA):
    """ Maps all outcome values to 1 or 0 if possible
        Parameters:
          currentTableMA (list of dictionaries)
        Returns: bool
    """
    isOneZero = True
    isYesNo = True
    isOneTwo = True
    isTrueFalse = True

    for lnsa in currentTableMA:
      for k, loutcome in lnsa.items():
        if isOneZero:
          if str(loutcome) not in ['1', '0']:
            isOneZero = False
        if isOneTwo:
          if str(loutcome) not in ['1', '2']:
            isOneTwo = False
        if isYesNo:
          if str(loutcome).lower() not in ['yes', 'no']:
            isYesNo = False
        if isTrueFalse:
          if str(loutcome).lower() not in ['true', 'false']:
            isTrueFalse = False
        break
        if isOneZero == False and isOneTwo == false and isYesNo == False and isTrueFalse == False:
          return False

    if isOneTwo:
      for lnsa in currentTableMA:
        for k, v in lnsa.items():
          if str(v) == '2':
            lnsa[k] = 0
          elif str(v) == '1':
            lnsa[k] = 1
          break
    elif isYesNo:
      for lnsa in currentTableMA:
        for k, v in lnsa.items():
          if str(v).lower == 'no':
            lnsa[k] = 0
          elif str(v).lower == 'yes':
            lnsa[k] = 1
          break
    elif isTrueFalse:
      for lnsa in currentTableMA:
        for k, v in lnsa.items():
          if str(v).lower == 'false':
            lnsa[k] = 0
          elif str(v).lower == 'true':
            lnsa[k] = 1
          break

    return True

  def checkIndependentVariables(self, currentTableMA, independentVariables):
    """ Creates dummy variables for independet variables when necessary
        Parameters:
          currentTableMA (list of dictionaries)
          independentVariables (list)
        Returns: bool
    """
    variablesNeedingDummies = []
    valuesForDummies = []
    rowOne = currentTableMA[0]

  def getCurrentTable(self, outcomeVariable, independentVariables):
    """ Creates an analysis dataset having a 0/1 dependent variable
        Parameters:
          outcomeVariable (str)
          IndependentVariables (list)
        Returns: bool
    """
    mutableCurrentTable = []
    lStrAVarNames = [iV for iV in independentVariables]
    columnsQueried = 1 + len(independentVariables) + 1
    for rowi in self.currentTable:
      row = {k : v for k, v in rowi.items() if k in [outcomeVariable] + independentVariables}
      row['RecStatus'] = 1
      mutableCurrentTable.append(row)

    if self.mstrGroupVar == None:
      pass
    else:
      # Matching later
      pass

    self.removeRecordsWithNulls(mutableCurrentTable)
    if self.outcomeOneZero(mutableCurrentTable) == False:
      return False
    self.checkIndependentVariables(mutableCurrentTable, independentVariables)
    self.currentTable = mutableCurrentTable

    return True

  def getRawData(self):
    """ Sets values for intercept and weight
        Parameters:
          none
        Returns: none
    """
    if self.mstrMatchVar is not None and len(str(self.mstrMatchVar)) > 0:
      self.mboolIntercept = False
    if self.mboolIntercept == True:
      self.lintIntercept = 1
    else:
      self.lintIntercept = 0
    if len(str(self.mstrWeightVar)) > 0:
      self.lintweight = 1
    else:
      self.lintweight = 0

    k = 0
    NumRows = len(self.currentTable)
    NumColumns = len(self.currentTable[0])

    if NumRows == 0:
      return

    lIntIsMatch = 0

    if len(mstrMatchVar) == 0:
      # Match Variable: write later
      i = 0
      lintnull = 0

    mVarArray = []
    if self.lStrAVarNames is not None:
      for i in range(len(self.lStrAVarNames)):
        mVarArray.append("")

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
    self.createSettings(inputVariableList)
    self.logisticResults = LogisticRegressionResults()
    if self.getCurrentTable(self.mstrDependVar, inputVariableList['exposureVariables']) == False:
      return
    self.getRawData

    return self.logisticResults
