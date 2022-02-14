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
    self.mboolFirst = None
    self.mstraBoolean = None
    self.mstrMatchVar = None
    self.mstrWeightVar = None
    self.mstrDependVar = None
    self.mstrGroupVar = None
    self.mstraTerms = None
    self.mStrADiscrete = None
    self.dummiesNSMA = []
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

  def pFromZ(self, _z):
    """ Computs a P value for a Z statistic
        Parameters:
          _z (float): Z statistic
        Returns: float
    """
    PFZ = 0.0
    UTZERO = 12
    CON = 1.28
    _x = _z
    if _z < 0.0:
      _x = _z * -1.0
    if _x > UTZERO:
      if _z < 0.0:
        PFZ = 1.0
      else:
        PFZ = 0.0
      return PFZ
    _y = _z ** 2.0 / 2.0
    if _x > CON:
      PFZ = _x - 0.151679116635 + 5.29330324926 / (_x + 4.8385912808 - 15.1508972451 / (_x + 0.742380924027 + 30.789933034 / (_x + 3.99019417011)))
      PFZ = _x + 0.000398064794 + 1.986158381364 / PFZ
      PFZ = _x - 0.000000038052 + 1.00000615302 / PFZ
      PFZ = 0.398942280385 * math.exp(-_y) / PFZ
    else:
      PFZ = _y / (_y + 5.75885480458 - 29.8213557808 / (_y + 2.624331121679 + 48.6959930692 / (_y + 5.92885724438)))
      PFZ = 0.398942280444 - 0.399903438504 * PFZ
      PFZ = 0.5 - _x * PFZ
    if _z < 0.0:
      PFZ = 1 - PFZ
    return PFZ

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

  def is_a_number(self, s):
    """ Determines whether a string is a number
        Parameters:
          s (str)
        Returns:
          bool
    """
    try:
      float(s)
      return True
    except ValueError:
      return False

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

  def makeDummies(self, currentTableMA, variablesNeedingDummies, valuesForDummies, independentVariables):
    """ Creates dummy variables for independent variables
        Parameters:
          currentTableMA (list of lists)
          variablesNeetingDummies (list)
          valuesForDummies (list)
          independentVariables (list)
        Returns: bool
    """
    if len(valuesForDummies) != len(variablesNeedingDummies):
      return
    newIndependentVariables = [iV for iV in independentVariables]
    i = 0
    for indexOfVariable in variablesNeedingDummies:
      valuesi = valuesForDummies[i]
      j = 0
      for nsmaij in currentTableMA:
        k = len(valuesi) - 1
        while k > 0:
          nsmaij.insert(indexOfVariable + 1, "")
          if j == 0:
            if indexOfVariable - 1 < len(newIndependentVariables):
              newIndependentVariables.insert(indexOfVariable, "")
            else:
              newIndependentVariables.append("")
          k -= 1
        k = len(valuesi) - 1
        while k >= 0:
          if nsmaij[indexOfVariable] == valuesi[k]:
            nsmaij[indexOfVariable + k] = '1'
          else:
            nsmaij[indexOfVariable + k] = '0'
          if j == 0:
            newIndependentVariables[indexOfVariable - 1 + k] = valuesi[k]
          k -= 1
        currentTableMA[j] = nsmaij
        j += 1
      i += 1
    independentVariables.clear()
    independentVariables += [iV for iV in newIndependentVariables]

  def checkIndependentVariables(self, currentTableMA, independentVariables):
    """ Checks whether independent variables need dummy variables
        Parameters:
          currentTableMA (list of lists)
          independentVariables (list)
        Returns: bool
    """
    variablesNeedingDummies = []
    valuesForDummies = []
    rowOne = currentTableMA[0]
    for j in range(1, len(rowOne)):
      isOneZero = True
      isYesNo = True
      isOneTwo = True
      isNumeric = True
      isTrueFalse = True
      for lsna in currentTableMA:
        loutcome = str(lsna[j])
        if isOneZero:
          if loutcome not in ['1','0']:
            isOneZero = False
        if isOneTwo:
          if loutcome not in ['1','2']:
            isOneTwo = False
        if isYesNo:
          if loutcome.lower() not in ['yes','no']:
            isYesNo = False
        if isTrueFalse:
          if loutcome.lower() not in ['true','false']:
            isTrueFalse = False
        if isNumeric:
          isNumeric = self.is_a_number(loutcome)
      if isOneTwo:
        for lsna in currentTableMA:
          if str(lsna[j]) == '1':
            lsna[j] = 1
          else:
            lsna[j] = 0
      elif isYesNo:
        for lsna in currentTableMA:
          if str(lsna[j]).lower() == 'yes':
            lsna[j] = 1
          else:
            lsna[j] = 0
      elif isTrueFalse:
        for lsna in currentTableMA:
          if str(lsna[j]).lower() == 'true':
            lsna[j] = 1
          else:
            lsna[j] = 0
      elif isNumeric == False or independentVariables[j - 1] in self.dummiesNSMA:
        variablesNeedingDummies.append(j)
        valuesForThisJ = []
        for lnsmai in currentTableMA:
          lnsmaij = str(lnsmai[j])
          if lnsmaij not in valuesForThisJ:
            valuesForThisJ.append(lnsmaij)
        if len(valuesForThisJ) > 1:
          valuesForThisJ.sort()
          valuesForDummies += valuesForThisJ[1:]
    if len(variablesNeedingDummies) > 0:
      self.makeDummies(currentTableMA, variablesNeedingDummies, valuesForDummies, independentVariables)

  def getCurrentTable(self, outcomeVariable, independentVariables):
    """ Creates an analysis dataset having a 0/1 dependent variable
        Converts analysis dataset from list of dictionaries to lise of lists
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
      # This section in Objective-C is for the value mapper,
      # which this routine does not have or need
      pass
    else:
      # This section in Objective-C is for the value mapper,
      # which this routine does not have or need
      pass

    self.removeRecordsWithNulls(mutableCurrentTable)
    if self.outcomeOneZero(mutableCurrentTable) == False:
      return False
    currentTableMutable = []
    for rd in mutableCurrentTable:
      rl = []
      for k, v in rd.items():
        rl.append(v)
      currentTableMutable.append(rl)
    self.checkIndependentVariables(currentTableMutable, independentVariables)
    self.currentTable = []
    for ctmr in currentTableMutable:
      self.currentTable.append(ctmr)

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
    self.NumRows = len(self.currentTable)
    self.NumColumns = len(self.currentTable[0])

    if self.NumRows == 0:
      return

    lIntIsMatch = 0

    if len(self.mstrMatchVar) == 0:
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
    self.getRawData()

    lintConditional = 0
    lintweight = 0
    ldblFirstLikelihood = 0.0
    ldblScore = 0.0

    self.mboolFirst = True
    self.mMatrixLikelihood = EIMatrix()
    self.mMatrixLikelihood.set_mboolFirst(self.mboolFirst)
    self.mMatrixLikelihood.set_mboolIntercept(self.mboolIntercept)
    if self.mstrGroupVar is not None:
      self.mMatrixLikelihood.set_mstrMatchVar(self.mstrGroupVar)
      lintConditional = 1
      self.NumColumns -= 1
      groupValues = []
      for i in range(len(self.currentTable)):
        groupValue = self.currentTable[i][1]
        if groupValue not in groupValues:
          groupValues.append(groupValue)
      self.mMatrixLikelihood.set_matchGroupValues(len(groupValues))
    else:
      self.mMatrixLikelihood.set_mstrMatchVar("")
    self.mMatrixLikelihood.MaximizeLikelihood(
                           self.NumRows,
                           self.NumColumns,
                           self.currentTable,
                           lintweight + lintConditional + 1,
                           self.NumColumns - (lintweight + lintConditional + 1),
                           self.mlngIter,
                           self.mdblToler,
                           self.mdblConv,
                           False)

    for ev in inputVariableList['exposureVariables']:
      self.logisticResults.Variables.append(ev)
    self.logisticResults.Variables.append('CONSTANT')
    mdblP = self.zFromP(0.025)
    i = 0
    for B in self.mMatrixLikelihood.get_mdblaB():
      self.logisticResults.Beta.append(B)
      self.logisticResults.SE.append(self.mMatrixLikelihood.get_mdblaInv()[i][i] ** 0.5)
      self.logisticResults.OR.append(math.exp(B))
      moe = mdblP * self.logisticResults.SE[i]
      self.logisticResults.ORLCL.append(math.exp(B - moe))
      self.logisticResults.ORUCL.append(math.exp(B + moe))
      self.logisticResults.Z.append(B / self.logisticResults.SE[i])
      self.logisticResults.PZ.append(2.0 * self.pFromZ(abs(self.logisticResults.Z[i])))
      i += 1

    return self.logisticResults
