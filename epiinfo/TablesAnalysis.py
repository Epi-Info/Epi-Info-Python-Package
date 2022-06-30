from scipy.stats import t as tdist
import math
import time
import itertools
from .randata import randata
from .RegressionUtilities import *

class TablesAnalysis:
  """ The TablesAnalysis class computes cell values for
      an MxN table that considers the relationship between
      two variables.

     Author: John Copeland
  """
  def __init__(self):
    self.includeMissing = False
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
    self.ColumnsAndValues = {}

    self.StartValues = []

    self.currentTable = None
    self.lStrAVarNames = None

    self.TableForInteractionTerms = []
    self.InteractionTerms = []

    self.logisticResults = LogisticRegressionResults()

    self.mMatrixLikelihood = EIMatrix()

  def PValFromChiSq(self, x, df):
    """ Computs a P value for a Chi Squared statistic
        with degrees of freedom.
        Parameters:
          x (float): Chi Squared statistic
          df (float): Degrees of Freedom
        Returns: float
    """
    j = 0.0
    k = 0.0
    l = 0.0
    m = 0.0
    pi = math.pi
    absx = x
    if x < 0.0:
      absx = -x
    if x < 0.000000001 or df < 1.0:
      return 1.0
    rr = 1.0
    ii = int(df * 1)
    while ii >= 2:
      rr *= float(ii * 1.0)
      ii -= 2
    k = math.exp(math.floor((df + 1.0) * 0.5) * math.log(absx) - x * 0.5) / rr
    if k < 0.00001:
      return 0.0
    if math.floor(df * 0.5) == df * 0.5:
      j = 1.0
    else:
      j = (2.0 / x / pi) ** 0.5
    l = 1.0
    m = 1.0
    if math.isnan(x) == False and math.isinf(x) == False:
      while m >= 0.00000001:
        df += 2.0
        m = m * x / df
        l = l + m
    return round(10000 * (1 - j * k * l)) / 10000

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

  def Run(self, inputVariableList, dataTable):
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
    if 'includeMissing' in inputVariableList:
      self.includeMissing = inputVariableList['includemissing']
    variablenames = []
    values = []
    tables = []
    for exposure in inputVariableList['exposureVariables']:
      outcomeValues = []
      exposureValues = []
      for tableRow in dataTable:
        if self.includeMissing == False and (tableRow[exposure] is None or tableRow[inputVariableList['outcomeVariable']] is None):
          continue
        if tableRow[inputVariableList['outcomeVariable']] not in outcomeValues:
          outcomeValues.append(tableRow[inputVariableList['outcomeVariable']])
        if tableRow[exposure] not in exposureValues:
          exposureValues.append(tableRow[exposure])
      if len(outcomeValues) == 2:
        if 1 in outcomeValues and 0 in outcomeValues and type(outcomeValues[0]) == int and type(outcomeValues[1]) == int:
          outcomeValues = [1, 0]
        elif True in outcomeValues and False in outcomeValues:
          outcomeValues = [True, False]
        elif 'Yes' in outcomeValues and 'No' in outcomeValues:
          outcomeValues = ['Yes', 'No']
        else:
          outcomeValues.sort()
      else:
        outcomeValues.sort()
      if len(exposureValues) == 2:
        if 1 in exposureValues and 0 in exposureValues and type(exposureValues[0]) == int and type(exposureValues[1]) == int:
          exposureValues = [1, 0]
        elif True in exposureValues and False in exposureValues:
          exposureValues = [True, False]
        elif 'Yes' in exposureValues and 'No' in exposureValues:
          exposureValues = ['Yes', 'No']
        else:
          exposureValues.sort()
      else:
        exposureValues.sort()
      onetable = [[0 for i in outcomeValues] for j in exposureValues]
      for tableRow in dataTable:
        if self.includeMissing == False and (tableRow[exposure] is None or tableRow[inputVariableList['outcomeVariable']] is None):
          continue
        rowindex = exposureValues.index(tableRow[exposure])
        colindex = outcomeValues.index(tableRow[inputVariableList['outcomeVariable']])
        onetable[rowindex][colindex] += 1
      tables.append(onetable)
      variablenames.append(inputVariableList['outcomeVariable'] + ' * ' + exposure)
      values.append([outcomeValues, exposureValues])
    print(variablenames)
    print(values)
    print(tables)
