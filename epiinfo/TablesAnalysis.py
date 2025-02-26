from scipy.stats import t as tdist
import math
import time
import itertools
from .randata import randata
from .BigDouble import BigDouble
from .RegressionUtilities import *

class TablesAnalysis:
  """ The TablesAnalysis class computes cell values for
      an MxN table that considers the relationship between
      two variables.

     Author: John Copeland
  """
  def __init__(self):
    self.includeMissing = False
    self.useCommonReference = False
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
    return round(100000000 * (1 - j * k * l)) / 100000000

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

  def RowTotals(self, table):
    """ Computes row totals for a table of counts
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list of ints
    """
    totals = [0] * len(table)
    for r in range(len(table)):
      for c in range(len(table[r])):
        totals[r] += table[r][c]
    return totals

  def ColumnTotals(self, table):
    """ Computes column totals for a table of counts
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list of ints
    """
    totals = [0] * len(table[0])
    for r in range(len(table)):
      for c in range(len(table[r])):
        totals[c] += table[r][c]
    return totals

  def EvalPoly(self, degC, r, bigC):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    y = 0.0
    bigY = BigDouble('doubleValue', 1.0)
    if r == 0.0:
      bigY.logValue = bigC[0].logValue
    elif r <= 1.0:
      bigY.logValue = bigC[int(degC)].logValue
      if r < 1.0:
        i = degC - 1
        while i >= 0:
          bigY.times(r)
          bigY.plusLog(bigC[i].logValue)
          i -= 1
      else:
        i = int(degC - 1)
        while i >= 0:
          bigY.plusLog(bigC[i].logValue)
          i -= 1
    elif r > 1.0:
      bigY.logValue = bigC[0].logValue
      r = 1.0 / r
      for i in range(1, int(degC) + 1):
        bigY.times(r)
        bigY.plusLog(bigC[i].logValue)
    return bigY.logValue

  def Func(self, r, sumA, value, degN, degD, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    func = 0.0
    logNumOverDenom = 0.0
    logNumer = self.EvalPoly(degN - 1, r, bigPolyN)
    logDenom = self.EvalPoly(degD - 1, r, bigPolyD)
    if r <= 1.0:
      logNumOverDenom = logNumer - logDenom
    else:
      logNumOverDenom = (logNumer - math.log10(r ** float((degD - 1) - (degN - 1)))) - logDenom
    func = 10 ** logNumOverDenom - value
    return func

  def BracketRoot(self, approx, x0, x1, f0, f1, sumA, value, degN, coordinates, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          none
    """
    iter = 0
    x1 = 0.5
    if approx > x1:
      x1 = approx;
    f0 = self.Func(x0, sumA, value, degN, degN, bigPolyD, bigPolyN)
    f1 = self.Func(x1, sumA, value, degN, degN, bigPolyD, bigPolyN)
    while f1 * f0 > 0.0 and iter < 10000:
      iter += 1
      x0 = x1
      f0 = f1
      x1 = x1 * 1.5 * iter
      f1 = self.Func(x1, sumA, value, degN, degN, bigPolyD, bigPolyN)
    coordinates[0] = x0
    coordinates[1] = x1
    coordinates[2] = f0
    coordinates[3] = f1

  def BracketRoots(self, approx, x0, x1, f0, f1, sumA, value, degN, coordinates, degD, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          none
    """
    iter = 0
    x1 = 0.5
    if approx > x1:
      x1 = approx;
    f0 = self.Func(x0, sumA, value, degN + 1, degD, bigPolyD, bigPolyN)
    f1 = self.Func(x1, sumA, value, degN + 1, degD, bigPolyD, bigPolyN)
    while f1 * f0 > 0.0 and iter < 10000:
      iter += 1
      x0 = x1
      f0 = f1
      x1 = x1 * 1.5 * iter
      f1 = self.Func(x1, sumA, value, degN + 1, degD, bigPolyD, bigPolyN)
    coordinates[0] = x0
    coordinates[1] = x1
    coordinates[2] = f0
    coordinates[3] = f1

  def Zero(self, x0, x1, f0, f1, sumA, value, degN, degD, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    found = False
    f2 = 0.0
    x2 = 0.0
    swap = 0.0
    iter = 0
    errorRenamed = 0
    absf0 = f0
    absf1 = f1
    if f0 < 0:
      absf0 = -f0
    if f1 < 0:
      absf1 = -f1
    if absf0 < absf1:
      swap = x0
      x0 = x1
      x1 = swap
      swap = f0
      f0 = f1
      f1 = swap
    found = (f1 == 0.0)
    if found == False and f0 * f1 > 0.0:
      errorRenamed = 1
    while found == False and iter < 10000 and errorRenamed == 0:
      iter += 1
      x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
      f2 = self.Func(x2, sumA, value, int(degN), int(degD), bigPolyD, bigPolyN)
      if f1 * f2 < 0.0:
        x0 = x1
        f0 = f1
      else:
        f0 = f0 * f1 / (f1 + f2)
      x1 = x2
      f1 = f2
      found = (abs(x1 / 0.0000001 - x0 / 0.0000001) * 0.0000001 < abs(x1 / 0.0000001) * 0.0000001 * 0.0000001 or f1 == 0.0)
    cmle = x1
    if found == False and iter > 10000:
      cmle = math.nan
    return cmle

  def Converge(self, approx, sumA, value, degN, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    coordinates = [0.0] * 4
    self.BracketRoot(approx, 0.0, 0.0, 0.0, 0.0, sumA, value, degN, coordinates, bigPolyD, bigPolyN)
    return self.Zero(coordinates[0], coordinates[1], coordinates[2], coordinates[3], sumA, value, degN, degN, bigPolyD, bigPolyN)

  def Converges(self, approx, sumA, value, degN, degD, bigPolyD, bigPolyN):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    coordinates = [0.0] * 4
    self.BracketRoots(approx, 0.0, 0.0, 0.0, 0.0, sumA, value, degN, coordinates, degD, bigPolyD, bigPolyN)
    return self.Zero(coordinates[0], coordinates[1], coordinates[2], coordinates[3], sumA, value, degN + 1, degD, bigPolyD, bigPolyN)

  def GetCmle(self, approx, minSumA, sumA, degN, bigPolyD):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    bigPolyN = []
    i = 0
    while i < degN:
      bigPolyN.append(BigDouble('logValue', bigPolyD[i].timesReturn(minSumA + i)))
      i += 1
    return self.Converge(approx, sumA, sumA, degN, bigPolyD, bigPolyN)

  def CalcCmle(self, approx, minSumA, sumA, maxSumA, degN, bigPolyD):
    """ Support function for exact odds ratios
        Parameters:
          numbers
        Returns:
          float
    """
    cmle = 0.0
    if minSumA < sumA and sumA < maxSumA:
      cmle = self.GetCmle(approx, minSumA, sumA, degN, bigPolyD)
    elif sumA == maxSumA:
      cmle = math.inf
    return cmle

  def addLogToLog(self, logValue, addedLogValue):
    if 10.0 ** logValue == 0:
        return addedLogValue
    power = round(logValue)
    return power + math.log10(10 ** (logValue - power) + 10.0 ** (addedLogValue - power))

  def ExactTests(self, minSumA, sumA, degD, ExactTests, bigPolyD, polyD):
    """ Support function for exact odds ratios
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list
    """
    diff = int(sumA - minSumA)
    bigUpTail = BigDouble('logValue', bigPolyD[int(degD)].logValue)
    bigUpTail = polyD[int(degD)]
    bigTwoTail = BigDouble('doubleValue', 1.0)
    bigTwoTail = float('-inf')
    if bigUpTail <= math.log10(1.000001) + polyD[diff]:
      bigTwoTail = bigUpTail
#    if bigUpTail.logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
#      bigTwoTail.plusLog(bigUpTail.logValue)
    i = int(degD) - 1
    while i >= diff:
      bigUpTail = self.addLogToLog(bigUpTail, polyD[i])
      if polyD[i] <= math.log10(1.000001) + polyD[diff]:
        bigTwoTail = self.addLogToLog(bigTwoTail, polyD[i])
#      bigUpTail.plusLog(bigPolyD[i].logValue)
#      if bigPolyD[i].logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
#        bigTwoTail.plusLog(bigPolyD[i].logValue)
      i -= 1
#    bigDenom = BigDouble('logValue', bigUpTail.logValue)
    bigDenom = bigUpTail
    i = diff - 1
    while i >= 0:
#      bigDenom.plusLog(bigPolyD[i].logValue)
      bigDenom = self.addLogToLog(bigDenom, polyD[i])
      if polyD[i] <= math.log10(1.000001) + polyD[diff]:
        bigTwoTail = self.addLogToLog(bigTwoTail, polyD[i])
#      if bigPolyD[i].logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
#        bigTwoTail.plusLog(bigPolyD[i].logValue)
      i -= 1
#    ExactTests[0] = 1.0 - (10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** (bigPolyD[diff].logValue - bigDenom.logValue))
#    ExactTests[1] = 10.0 ** (bigUpTail.logValue - bigDenom.logValue)
#    ExactTests[2] = 10.0 ** (bigTwoTail.logValue - bigDenom.logValue)
#    ExactTests[3] = 1.0 - (10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** ((math.log10(0.5) + bigPolyD[diff].logValue) - bigDenom.logValue))
#    ExactTests[4] = 10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** ((math.log10(0.5) + bigPolyD[diff].logValue) - bigDenom.logValue)
    ExactTests[0] = 1.0 - (10.0 ** (bigUpTail - bigDenom) - 10.0 ** (polyD[diff] - bigDenom))
    ExactTests[1] = 10.0 ** (bigUpTail - bigDenom)
    ExactTests[2] = 10.0 ** (bigTwoTail - bigDenom)
    ExactTests[3] = 1.0 - (10.0 ** (bigUpTail - bigDenom) - 10.0 ** ((math.log10(0.5) + polyD[diff]) - bigDenom))
    ExactTests[4] = 10.0 ** (bigUpTail - bigDenom) - 10.0 ** ((math.log10(0.5) + polyD[diff]) - bigDenom)

  def RRStats(self, table):
    """ Support function for exact odds ratios
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list
    """
    RRstats = [0.0] * 12
    a = table[0][0]
    b = table[0][1]
    c = table[1][0]
    d = table[1][1]
    n1 = table[0][0] + table[0][1]
    n0 = table[1][0] + table[1][1]
    m1 = table[0][0] + table[1][0]
    m0 = table[0][1] + table[1][1]
    n = m1 + m0
    re = table[0][0] / n1
    ru = table[1][0] / n0
    if ru < 0.00001:
      RRstats[0] = -1.0
      RRstats[1] = -1.0
      RRstats[2] = -1.0
    else:
      RRstats[0] = re / ru
      if re < 0.00001:
        RRstats[1] = -1.0
        RRstats[2] = -1.0
      else:
        RRstats[1] = math.exp(math.log((a / n1) / (c / n0)) - 1.96 * (d / (c * n0) + b / (n1 * a)) ** 0.5)
        RRstats[2] = math.exp(math.log((a / n1) / (c / n0)) + 1.96 * (d / (c * n0) + b / (n1 * a)) ** 0.5)
    RRstats[3] = (re - ru) * 100
    RRstats[4] = (re - ru - 1.96 * (re * (1 - re) / n1 + ru * (1 - ru) / n0) ** 0.5) * 100
    RRstats[5] = (re - ru + 1.96 * (re * (1 - re) / n1 + ru * (1 - ru) / n0) ** 0.5) * 100
    h3 = m1 * m0 * n1 * n0
    phi = (a * d - b * c) / h3 ** 0.5
    RRstats[6] = n * phi ** 2.0
    RRstats[7] = self.PValFromChiSq(RRstats[6], 1.0)
    RRstats[8] = (n - 1) / h3 * (a * d - b * c) ** 2.0
    RRstats[9] = self.PValFromChiSq(RRstats[8], 1.0)
    RRstats[10] = n / h3 * (max(0.0, float(abs(int(a * d - b * c))) - n * 0.5) ** 2.0)
    RRstats[11] = self.PValFromChiSq(RRstats[10], 1.0)
    return RRstats

  def CalcPoly(self, table):
    """ Support function for exact odds ratios
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list
    """
    ExactResults = [0.0, 0.0, 0.0, 0.0]
    n0 = float(table[1][0] + table[1][1])
    n1 = float(table[0][0] + table[0][1])
    m1 = float(table[0][0] + table[1][0])
    minA = 0.0
    if m1 - n0 > 0.0:
      minA = m1 - n0
    maxA = m1
    if n1 < m1:
      maxA = n1
    bigPolyD = [BigDouble('logValue', 0.0) for i in range(int(maxA - minA + 1))]
    polyDi = [None] * (round(maxA - minA) + 1)
    bigPolyD[0] = BigDouble('doubleValue', 1.0)
    polyDi[0] = 0.0
    aa = minA
    bb = m1 - minA + 1.0
    cc = n1 - minA + 1.0
    dd = n0 - m1 + minA
    for i in range(1, round(maxA - minA) + 1):
      bigPolyD[i] = BigDouble('logValue', bigPolyD[i - 1].timesReturn(((bb - i) / (aa + i)) * ((cc - i) / (dd + i))))
      polyDi[i] = float(polyDi[i - 1] * ((bb - i) / (aa + i)) * ((cc - i) / (dd + i)))
      polyDi[i] = polyDi[i - 1] + math.log10(((bb - i) / (aa + i)) * ((cc - i) / (dd + i)))
    ExactResults[0] = self.CalcCmle(1.0, minA, table[0][0], maxA, maxA - minA + 1, bigPolyD)
    ExactTestsList = [None] * 5
    self.ExactTests(minA, table[0][0], maxA - minA, ExactTestsList, bigPolyD, polyDi)
    ExactResults[1] = ExactTestsList[3]
    if ExactTestsList[4] < ExactTestsList[3]:
      ExactResults[1] = ExactTestsList[4]
    ExactResults[2] = ExactTestsList[0]
    if ExactTestsList[1] < ExactTestsList[0]:
      ExactResults[2] = ExactTestsList[1]
    ExactResults[3] = ExactTestsList[2]
    return ExactResults

  def GetExactLim(self, pbLower, pbFisher, approx, minSumA, sumA, degD, bigPolyD):
    """ Computes the single 2x2 table statistics
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          float
    """
    degN = int(sumA) - int(minSumA)
    pnConfLevel = 0.95
    value = 0.5 * (1.0 - pnConfLevel)
    if pbLower == True:
      value = 0.5 * (1 + pnConfLevel)
    if pbLower == True and pbFisher == True:
      degN = int(sumA) - int(minSumA) - 1
    bigPolyN = [BigDouble('logValue', 0.0) for i in range(int(degN + 1))]
    bigPolyNI = BigDouble('doubleValue', 1.0)
    for i in range(degN + 1):
      bigPolyN[i] = BigDouble('logValue', bigPolyD[i].logValue)
      bigPolyNI.logValue = bigPolyN[i].logValue
    if pbFisher == False:
      bigPolyN[degN].logValue = math.log10(0.5) + bigPolyD[degN].logValue
    limit = self.Converges(approx, sumA, value, degN, degD, bigPolyD, bigPolyN)
    return limit

  def CalcExactLim(self, pbLower, pbFisher, approx, minSumA, sumA, maxSumA, degD, bigPolyD):
    """ Computes the single 2x2 table statistics
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          float
    """
    limit = 0.0
    if minSumA < sumA and sumA < maxSumA:
      limit = self.GetExactLim(pbLower, pbFisher, approx, minSumA, sumA, degD, bigPolyD)
    elif sumA == minSumA:
      if pbLower == False:
        limit = self.GetExactLim(pbLower, pbFisher, approx, minSumA, sumA, degD, bigPolyD)
    elif sumA == maxSumA:
      if pbLower == True:
        limit = self.GetExactLim(pbLower, pbFisher, approx, minSumA, sumA, degD, bigPolyD)
      else:
        limit = math.inf
    return limit

  def FishOR(self, a, b, c, d, alpha, OR, plusA, one, minus, pbLower, pbFisher):
    """ Computes the single 2x2 table statistics
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          float
    """
    m1 = a + c
    n0 = c + d
    n1 = a + b
    maxK = m1
    if n1 < m1:
        maxK = n1
    minK = 0
    if m1 - n0 > 0:
      minK = m1 - n0
    bigPolyDD = [BigDouble('logValue', 0.0) for i in range(int(maxK - minK + 1))]
    bigPolyDD[0] = BigDouble('doubleValue', 1.0)
    bb = m1 - minK + 1
    cc = n1 - minK + 1
    dd = n0 - m1 + minK
    bigPolyDDI = BigDouble('logValue', bigPolyDD[0].logValue)
    degD = int(maxK - minK) + 1
    for i in range(1, degD):
      bigPolyDD[i] = BigDouble('logValue', bigPolyDD[i - 1].timesReturn((float(bb - i) / float(minK + i)) * (float(cc - i) / float(dd + i))))
      bigPolyDDI.logValue = bigPolyDD[i].logValue
    return self.CalcExactLim(pbLower, pbFisher, OR, minK, a, maxK, degD, bigPolyDD)

  def TwoX2Compute(self, table):
    """ Computes the single 2x2 table statistics
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          dict
    """
    stats = {}
    a = table[0][0]
    b = table[0][1]
    c = table[1][0]
    d = table[1][1]
    # Odds Ratio
    if b * c == 0:
      stats['OR'] = math.inf
      stats['ORLL'] = math.inf
      stats['ORUL'] = math.inf
    elif a * d == 0:
      stats['OR'] = 0.0
      stats['ORLL'] = 0.0
      stats['ORUL'] = 0.0
    else:
      stats['OR'] = ((float(a) * float(d)) / (float(b) * float(c)))
      stats['ORLL'] = math.exp(math.log((float(a) * float(d)) / (float(b) * float(c))) - 1.96 * (1 / float(a) + 1 / float(b) + 1 / float(c) + 1 / float(d)) ** 0.5)
      stats['ORUL'] = math.exp(math.log((float(a) * float(d)) / (float(b) * float(c))) + 1.96 * (1 / float(a) + 1 / float(b) + 1 / float(c) + 1 / float(d)) ** 0.5)
    ExactResults = self.CalcPoly(table)
    stats['MidPOR'] = ExactResults[0]
    stats['MidPORLL'] = self.FishOR(a, b, c, d, 0.05, ExactResults[0], 0, 0, 1, True, False)
    stats['MidPORUL'] = self.FishOR(a, b, c, d, 0.05, ExactResults[0], 1, 1, -1, False, False)
    stats['FisherORLL'] = self.FishOR(a, b, c, d, 0.05, ExactResults[0], 0, 0, 1, True, True)
    stats['FisherORUL'] = self.FishOR(a, b, c, d, 0.05, ExactResults[0], 1, 1, -1, False, True)
    stats['MidPExact1Tail'] = ExactResults[1]
    stats['FisherExact1Tail'] = ExactResults[2]
    stats['FisherExact2Tail'] = ExactResults[3]
    RRStatsResults = self.RRStats(table)
    stats['RiskRatio'] = RRStatsResults[0]
    stats['RiskRatioLL'] = RRStatsResults[1]
    stats['RiskRatioUL'] = RRStatsResults[2]
    stats['RiskDifference'] = RRStatsResults[3]
    stats['RiskDifferenceLL'] = RRStatsResults[4]
    stats['RiskDifferenceUL'] = RRStatsResults[5]
    stats['UncorrectedX2'] = RRStatsResults[6]
    stats['UncorrectedX2P'] = RRStatsResults[7]
    stats['MHX2'] = RRStatsResults[8]
    stats['MHX2P'] = RRStatsResults[9]
    stats['CorrectedX2'] = RRStatsResults[10]
    stats['CorrectedX2P'] = RRStatsResults[11]
    return stats

  def transposeMatrix(self, matrix):
    """ Computes the MxN table statistics
        Parameters:
          matrix (list(list)): A list of lists of ints
        Returns:
          list of lists
    """
    c = len(matrix[0])
    r = len(matrix)
    mtx = []
    for j in range(c):
      nsma0 = []
      for i in range(r):
        nsa = matrix[i]
        cellvalue = nsa[j]
        nsma0.append(cellvalue)
      mtx.append(nsma0)
    return mtx

  def iwork(self, iwkmax, iwkpt, number, itype):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          int
    """
    iwork = iwkpt[0]
    if itype == 2 or itype == 3:
      iwkpt[0] += number
    else:
      if iwork % 2 != 0:
        iwork += 1
      iwkpt[0] += 2 * number
      iwork /= 2
    return int(iwork)

  def f9xact(self, n, mm, ir, iroffset, fact):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          float
    """
    f9xact = fact[mm]
    for k in range(1, n + 1):
      f9xact -= fact[ir[k + iroffset - 1]]
    return f9xact

  def f6xact(self, nrow, irow, iflag, kyy, key, keyoffset, ldkey, last, ipn):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    kval = 0
    while True:
      last[0] += 1
      if last[0] <= ldkey:
        if key[last[0] + keyoffset - 1] < 0:
          continue
        kval = key[last[0] + keyoffset - 1]
        key[last[0] + keyoffset - 1] = -9999
        j = nrow 
        while j >= 2:
          irow[j] = int(kval / kyy[j])
          kval -= irow[j] * kyy[j]
          j -= 1
        irow[1] = kval
        ipn[0] = last[0]
        break
      else:
        last[0] = 0
        iflag[0] = 3
        break

  def f8xact(self, irow, irowoffset, iz, i1, izero, noo, noooffset):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    i = 0
    for i in range(1, i1):
      noo[i + noooffset - 1] = irow[i + irowoffset - 1]
    resetI = True
    for i in range(i1, izero):
      if iz >= irow[i + irowoffset - 1 + 1]:
        resetI = False
        break
      noo[i + noooffset - 1] = irow[i + irowoffset - 1 + 1]
    if resetI == True:
      i = izero
    noo[i + noooffset - 1] = iz
    while True:
      i += 1
      if i > izero:
        return
      noo[i + noooffset - 1] = irow[i + irowoffset - 1]

  def alogam(self, x, ifault):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          float
    """
    alogam = 0.0
    a1 = 0.918938533204673
    a2 = 0.00595238095238
    a3 = 0.00793650793651
    a4 = 0.002777777777778
    a5 = 0.08333333333333
    half = 0.5
    zero = 0.0
    one = 1.0
    seven = 7.0
    alogam = zero
    ifault = 1
    if x < zero:
      return alogam
    ifault = 0
    y = x
    f = zero
    if y > seven:
      z = one / (y * y)
      alogam = f + (y - half) * log(y) - y + a1 + (((-a2 * z + a3) * z - a4) * z + a5) / y
      return alogam
    f = y
    while True: # Fortran line 10
      y = y + one
      if y >= seven:
        pass
      else:
        f = f * y
        continue
      f = -1.0 * log(f)
      z = one / (y * y)
      alogam = f + (y - half) * log(y) - y + a1 + (((-a2 * z + a3) * z - a4) * z + a5) / y
      return alogam

  def gamds(self, y, p, ifault):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          float
    """
    gammds = 0.0
    ifail = 1
    e = 1.0e-6
    zero = 0.0
    one = 1.0
    ifault = 1
    gammds = zero
    if y <= 0 or p <= 0:
      return gammds
    ifault = 2
    f = math.exp(p * math.log(y) - self.alogam(p + one, ifail) - y)
    if f == zero:
      return gammds
    ifault = 0
    c = one
    gammds = one
    a = p
    while True: # Fortran line 10
      a = a + one
      c = c * y / a
      gammds = gammds + c
      if c / gammds > e:
        continue
      break
    gammds = gammds * f
    return gammds

  def f11act(self, irow, irowoffset, i1, i2, noo, noooffset):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    for i in range(1, i1):
      noo[i + noooffset - 1] = irow[i + irowoffset - 1]
    for i in range(i1, i2 + 1):
      noo[i + noooffset - 1] = irow[i + irowoffset - 1 + 1]

  def f10act(self, nrow, irow, irowoffset, ncol, icol, icoloffset, val, xmin, fact, nd, ne, m):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    for i in range(1, nrow):
      nd[i] = 0
    iz = int(icol[1 + icoloffset - 1] / nrow)
    ne[1] = iz
    ix = int(icol[1 + icoloffset - 1] - nrow * iz)
    m[1] = ix
    if ix != 0:
      nd[ix] += 1
    for i in range(2, ncol + 1):
      ix = int(icol[i + icoloffset - 1] / nrow)
      ne[i] = ix
      iz += ix
      ix = int(icol[i + icoloffset - 1] - nrow * ix)
      m[i] = ix
      if ix != 0:
        nd[ix] += 1
    i = nrow - 2
    while i >= 1:
      nd[i] += nd[i + 1]
      i -= 1
    ix = 0
    nrw1 = int(nrow + 1)
    i = nrow
    while i >= 2:
      ix += int(iz + nd[nrw1 - i] - irow[i + irowoffset - 1])
      if ix < 0:
        return
      i -= 1
    for i in range(1, ncol + 1):
      ix = int(ne[i])
      iz = int(m[i])
      val[0] += iz * fact[ix + 1] + (nrow - iz) * fact[ix]
    xmin[0] = True

  def f4xact(self, nrow, irow, irowoffset, ncol, icol, icoloffset, dsp, fact, icstkk, ncstk, lstk, mstk, nstk, nrstk, irstkk, ystk, tol):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    i = 1
    j = 1
    irstk = [[0 for r in range(0, nrow + ncol + 1)] for c in range(0, nrow + ncol + 1)]
    icstk = [[0 for r in range(0, nrow + ncol + 1)] for c in range(0, nrow + ncol + 1)]
    if nrow == 1:
      for i in range(1, ncol + 1):
        dsp[0] -= fact[icol[i + icoloffset - 1]]
      return
    if ncol == 1:
      for i in range(1, nrow + 1):
        dsp[0] -= fact[irow[i + irowoffset - 1]]
      return
    if nrow * ncol == 4:
      if irow[2 + irowoffset - 1] <= icol[2 + icoloffset - 1]:
        dsp[0] = dsp[0] - fact[irow[2 + irowoffset - 1]] - fact[icol[1 + icoloffset - 1]] - fact[icol[2 + icoloffset - 1] - irow[2 + irowoffset - 1]]
      else:
        dsp[0] = dsp[0] - fact[icol[2 + icoloffset - 1]] - fact[irow[1 + irowoffset - 1]] - fact[irow[2 + irowoffset - 1] - icol[2 + icoloffset - 1]]
      return
    for i in range(1, nrow + 1):
      irstk[1][i] = irow[nrow - i + 1 + irowoffset - 1]
    for j in range(1, ncol + 1):
      icstk[1][j] = icol[ncol - j + 1 + icoloffset - 1]
    nro = nrow;
    nco = ncol
    nrstk[1] = nro
    ncstk[1] = nco
    ystk[1] = 0.0
    y = 0.0
    istk = 1
    l = 1
    amx = 0.0
    while True: # Fortran line 50
      if True:
        goto50bool = False
        ir1 = irstk[istk][1]
        ic1 = icstk[istk][1]
        m = 0
        n = 0
        k = 0
        mn = 0
        if ir1 > ic1:
          if nro >= nco:
            m = nco - 1
            n = 2
          else:
            m = nro
            n = 1
        elif ir1 < ic1:
          if nro <= nco:
            m = nro - 1
            n = 1
          else:
            m = nco
            n = 2
        else:
          if nro <= nco:
            m = nro - 1
            n = 1
          else:
            m = nco - 1
            n = 2
        while True: # Fortran line 60
          if True:
            goto60bool = False
            if n == 1:
              i = l
              j = 1
            else:
              i = 1
              j = l
            irt = irstk[istk][i]
            ict = icstk[istk][j]
            mn = irt
            if mn > ict:
              mn = ict
            y += fact[mn]
            if irt == ict:
              nro -= 1
              nco -= 1
              irstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA[zfj] = int(irstk[istk][zfj])
              irstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA1[zfj] = int(irstk[istk + 1][zfj])
              self.f11act(irstkIA, 1, i, nro, irstkIA1, 1)
              icstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA[zfj] = int(icstk[istk][zfj])
              icstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA1[zfj] = int(icstk[istk + 1][zfj])
              self.f11act(icstkIA, 1, j, nco, icstkIA1, 1)
              for tk in range(0, len(icstkIA)):
                icstk[istk][tk] = icstkIA[tk]
                icstk[istk + 1][tk] = icstkIA1[tk]
              for tk in range(0, len(icstkIA)):
                irstk[istk][tk] = irstkIA[tk]
                irstk[istk + 1][tk] = irstkIA1[tk]
            elif irt >= ict:
              nco -= 1
              irstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA[zfj] = int(irstk[istk][zfj])
              irstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA1[zfj] = int(irstk[istk + 1][zfj])
              icstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA[zfj] = int(icstk[istk][zfj])
              icstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA1[zfj] = int(icstk[istk + 1][zfj])
              self.f11act(icstkIA, 1, j, nco, icstkIA1, 1)
              self.f8xact(irstkIA, 1, irt - ict, i, nro, irstkIA1, 1)
              for tk in range(0, len(icstkIA)):
                icstk[istk][tk] = icstkIA[tk]
                icstk[istk + 1][tk] = icstkIA1[tk]
              for tk in range(0, len(icstkIA)):
                irstk[istk][tk] = irstkIA[tk]
                irstk[istk + 1][tk] = irstkIA1[tk]
            else:
              nro -= 1
              irstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA[zfj] = int(irstk[istk][zfj])
              irstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                irstkIA1[zfj] = int(irstk[istk + 1][zfj])
              icstkIA = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA[zfj] = int(icstk[istk][zfj])
              icstkIA1 = [0] * (nrow + ncol + 1)
              for zfj in range(0, nrow + ncol + 1):
                icstkIA1[zfj] = int(icstk[istk + 1][zfj])
              self.f11act(irstkIA, 1, i, nro, irstkIA1, 1)
              self.f8xact(icstkIA, 1, ict - irt, j, nco, icstkIA1, 1)
              for tk in range(0, len(icstkIA)):
                icstk[istk][tk] = icstkIA[tk]
                icstk[istk + 1][tk] = icstkIA1[tk]
              for tk in range(0, len(icstkIA)):
                irstk[istk][tk] = irstkIA[tk]
                irstk[istk + 1][tk] = irstkIA1[tk]
            goto90bool = False
            if nro == 1:
              for k in range(1, nco + 1):
                y += fact[icstk[istk + 1][k]]
              goto90bool = True
            if goto90bool == False and nco == 1:
              for k in range(1, nro + 1):
                y += fact[irstk[istk + 1][k]]
              goto90bool = True
            if goto90bool == False:
              lstk[istk] = l
              mstk[istk] = m
              nstk[istk] = n
              istk += 1
              nrstk[istk] = nro
              ncstk[istk] = nco
              ystk[istk] = y
              l = 1
              goto50bool = True
              break
            if y >  amx:
              amx = y
              if dsp[0] - amx <= tol:
                dsp[0] = 0.0
                return
            while True: # Fortran line 100
              goto100bool = False
              istk -= 1
              if istk == 0:
                dsp[0] -= amx
                if dsp[0] - amx <= tol:
                  dsp[0] = 0.0
                return
              l = int(lstk[istk] + 1)
              while True: # Fortran line 110
                if l > mstk[istk]:
                  goto100bool = True
                  break
                n = nstk[istk]
                nro = nrstk[istk]
                nco = ncstk[istk]
                y = ystk[istk]
                if n == 1:
                  if irstk[istk][l] < irstk[istk][l - 1]:
                    goto60bool = True
                    break
                elif n == 2:
                  if icstk[istk][l] < icstk[istk][l - 1]:
                    goto60bool = True
                    break
                l += 1
                continue
              if goto100bool:
                continue
              break
            if goto60bool:
              continue
            break
        if goto50bool:
          continue
        break

  def f3xact(self, nrow, irow, irowoffset, ncol, icol, icoloffset, dlp, mm, fact, ico, iro, it, lb, nr, nt, nu, itc, ist, stv, alen, tol):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    n11 = 0
    n12 = 0
    nro = 0
    nco = 0
    val = [0.0]
    nn = 0
    xmin = False
    nitc = 0
    nst = 0
    nn1 = 0
    nc1 = 0
    ic1 = 0
    ic2 = 0
    ii = 0
    key = 0
    ipn = 0
    itp = 0
    for i in range(0, ncol + 1):
      alen[i] = 0.0
    for i in range(1, 401):
      ist[i] = -1
    if nrow <= 1:
      if nrow > 0:
        dlp[0] -= fact[icol[1]]
        for i in range(2, ncol + 1):
          dlp[0] -= fact[icol[i + icoloffset - 1]]
      return
    if ncol <= 1:
      if ncol > 0:
        dlp[0] -= fact[irow[1]] - fact[irow[2]]
        for i in range(3, nrow + 1):
          dlp[0] -= fact[irow[i + irowoffset - 1]]
      return
    if nrow * ncol == 4:
      n11 = int((irow[1 + irowoffset - 1] + 1) * (icol[1 + icoloffset - 1] + 1) / (mm[0] + 2))
      n12 = int(irow[1 + irowoffset - 1] - n11)
      dlp[0] = dlp[0] - fact[n11] - fact[n12] - fact[icol[1 + icoloffset - 1] - n11] - fact[icol[2 + icoloffset - 1] - n12]
      return
    val = [0.0]
    xmin = [False]
    if irow[nrow + irowoffset - 1] <= irow[1 + irowoffset - 1] + ncol:
      self.f10act(nrow, irow, irowoffset, ncol, icol, icoloffset, val, xmin, fact, lb, nu, nr)
    if xmin[0] == False:
      if icol[ncol + icoloffset - 1] <= icol[1 + icoloffset - 1] + nrow:
        self.f10act(ncol, icol, icoloffset, nrow, irow, irowoffset, val, xmin, fact, lb, nu, nr)
    if xmin[0] == True:
      dlp[0] -= val[0]
      return
    nn = mm[0]
    if nrow >= ncol:
      nro = nrow
      nco = ncol
      for i in range(1, nrow + 1):
        iro[i] = irow[i + irowoffset - 1]
      ico[1] = icol[1 + icoloffset - 1]
      nt[1] = nn - ico[1]
      for i in range(2, ncol + 1):
        ico[i] = icol[i + icoloffset - 1]
        nt[i] = nt[i - 1] - ico[i]
    else:
      nro = ncol
      nco = nrow
      ico[1] = irow[1 + irowoffset - 1]
      nt[1] = nn - ico[1]
      for i in range(2, nrow + 1):
        ico[i] = irow[i + irowoffset - 1]
        nt[i] = nt[i - 1] - ico[i]
      for i in range(1, ncol + 1):
        iro[i] = icol[i + icoloffset - 1]
    vmn = 1.0e10
    nc1s = nco - 1
    irl = 1
    ks = 0
    ldst = 200
    k = ldst
    kyy = ico[nco] + 1
    bool90 = False # to goto 100
    while True: # Fortran line 90
      if True:
        goto90bool = False
        if bool90:
          xmin = [False]
          if iro[nro - 1] <= iro[irl - 1] + nco:
            self.f10act(nro, iro, irl, nco, ico, 1, val, xmin, fact, lb, nu, nr)
          if xmin[0] == False:
            if ico[nco - 1] <= ico[0] + nro:
              self.f10act(nco, ico, 1, nro, iro, irl, val, xmin, fact, lb, nu, nr)
          if xmin[0] == True:
            if val[0] < vmn:
              vmn = val[0]
            while True: # Fortran line 200
              if nitc > 0:
                itp = itc[nitc + k] + k
                nitc -= 1
                val[0] = stv[itp]
                key = int(ist[itp])
                ist[itp] = -1
                i = nco
                while i >= 2:
                  ico[i] = key % kyy
                  key = int(key / kyy)
                  i -= 1
                ico[i] = key
                nt[i] = nn - ico[1]
                for i in range(2, nco + 1):
                  nt[i] = nt[i - 1] - ico[i]
                goto90bool = True
                break
              elif nro > 2 and nst > 0:
                nitc = nst
                nst = 0
                k = ks
                ks = ldst - ks
                nn -= iro[irl]
                irl += 1
                nro -= 1
                continue
              dlp[0] -= vmn
              return
        if goto90bool:
          continue
        while True: # Fortran line 100
          if True:
            lev = 1
            nr1 = nro - 1
            nrt = iro[irl]
            nct = ico[1]
            lb[1] = int(float((nrt + 1) * (nct + 1)) / float(nn + nr1 * nc1s + 1) - tol) - 1
            nu[1] = int(float((nrt + nc1s) * (nct + nr1)) / float(nn + nr1 * nc1s)) - lb[1] + 1
            nr[1] = nrt - lb[1]
            while True: # Fortran line 110
              if True:
                goto110bool = False
                nu[lev] -= 1
                if nu[lev] == 0:
                  if lev == 1:
                    while True: # Fortran line 200
                      if nitc > 0:
                        itp = itc[nitc + k] + k
                        nitc -= 1
                        val[0] = stv[itp]
                        key = int(ist[itp])
                        ist[itp] = -1
                        i = nco
                        while i >= 2:
                          ico[i] = key % kyy
                          key = int(key / kyy)
                          i -= 1
                        ico[i] = key
                        nt[i] = nn - ico[1]
                        for i in range(2, nco + 1):
                          nt[i] = nt[i - 1] - ico[i]
                        goto90bool = True
                        break
                      elif nro > 2 and nst > 0:
                        nitc = nst
                        nst = 0
                        k = ks
                        ks = ldst - ks
                        nn -= iro[irl]
                        irl += 1
                        nro -= 1
                        continue
                      dlp[0] -= vmn
                      return
                    if goto90bool:
                      break
                  lev -= 1
                  continue
                lb[lev] += 1
                nr[lev] -= 1
                while True: # Fortran line 120
                  if True:
                    alen[lev] = alen[lev - 1] + fact[lb[lev]]
                    if lev < nc1s:
                      nn1 = int(nt[lev])
                      nrt = nr[lev]
                      lev += 1
                      nc1 = nco - lev
                      nct = ico[lev]
                      lb[lev] = int(float((nrt + 1) * (nct + 1)) / float(nn1 + nr1 * nc1 + 1) - tol)
                      nu[lev] = int(float((nrt + nc1) * (nct + nr1)) / float(nn1 + nr1 + nc1) - lb[lev] + 1)
                      nr[lev] = nrt - lb[lev]
                      continue
                    alen[nco] = alen[lev] + fact[nr[lev]]
                    lb[nco] = nr[lev]
                    v = val[0] + alen[nco]
                    if nro == 2:
                      v += fact[ico[1] - lb[1]] + fact[ico[2] - lb[2]]
                      for i in range(3, nco + 1):
                        v += fact[ico[i] - lb[i]]
                      if v < vmn:
                        vmn = v
                    elif nro == 3 and nco == 2:
                      nn1 = int(nn - iro[irl] + 2)
                      ic1 = int(ico[1] - lb[1])
                      ic2 = int(ico[2] - lb[2])
                      n11 = int((iro[irl + 1] + 1) * (ic1 + 1) / nn1)
                      n12 = int(iro[irl + 1] - n11)
                      v += fact[n11] + fact[n12] + fact[ic1 - n11] + fact[ic2 - n12]
                      if v < vmn:
                        vmn = v
                    else:
                      for i in range(1, nco + 1):
                        it[i] = ico[i] - lb[i]
                      if nco == 2:
                        if it[1] > it[2]:
                          ii = it[1]
                          it[1] = it[2]
                          it[2] = ii
                      elif nco == 3:
                        ii = it[1]
                        if ii > it[3]:
                          if ii > it[2]:
                            if it[2] > it[3]:
                              it[1] = it[3]
                              it[3] = ii
                            else:
                              it[1] = it[2]
                              it[2] = it[3]
                              it[3] = ii
                          else:
                            it[1] = it[3]
                            it[3] = it[2]
                            it[2] = ii
                        elif ii > it[2]:
                          it[1] = it[2]
                          it[2] = ii
                        elif it[2] > it[3]:
                          ii = it[2]
                          it[2] = it[3]
                          it[3] = ii
                      else:
                        it.sort()
                      key = int(it[1] * kyy + it[2])
                      for i in range(3, nco + 1):
                        key = int(it[i] + key * kyy)
                      ipn = int(key % ldst + 1)
                      skipto180 = False
                      goto190 = False
                      ii = ks + ipn
                      for itp in range(ipn, ldst + 1):
                        if ist[ii] < 0:
                          skipto180 = True
                          goto190 = False
                          break
                        elif ist[ii] == key:
                          skipto180 = True
                          goto190 = True
                          break
                        ii += 1
                      if skipto180 == False:
                        ii = ks + 1
                        for itp in range(1, ipn ):
                          if ist[ii] < 0:
                            skipto180 = True
                            goto190 = False
                            break
                          elif ist[ii] == key:
                            skipto180 = True
                            goto190 = True
                            break
                          ii += 1
                      if goto190 == False:
                        ist[ii] = key # Fortran line 180
                        stv[ii] = v
                        nst += 1
                        ii = nst + ks
                        itc[ii] = itp
                        goto110bool = True
                        break
                      stv[ii] = min(v, stv[ii]) # Fortran line 190
                    goto110bool = True
                    break
                if goto110bool:
                  continue
                break
            break
        if goto90bool:
          bool90 = True
          continue
        break

  def f7xact(self, nrow, imax, idif, k, ks, iflag):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    m = 0
    k1 = 0
    mm = 0
    iflag[0] = 0
    if ks[0] == 0:
      while True: # Fortran line 10
        ks[0] += 1
        if idif[ks[0]] == imax[ks[0]]:
          continue
        break
    if idif[k[0]] > 0 and k[0] > ks[0]:
      idif[k[0]] = idif[k[0]] - 1
      while True: # Fortran line 30
        k[0] -= 1
        if imax[k[0]] == 0:
          continue
        break
      m = k[0]
      while True: # Fortran line 40
        if idif[m] >= imax[m]:
          m -= 1
          continue
        break
      idif[m] = idif[m] + 1
      if m == ks[0]:
        if idif[m] == imax[m]:
          ks[0] = k[0]
    else:
      while True: # Fortran line 50
        goto50bool = False
        while True: # Also for Fortran line 50
          goto70bool = False
          for k1 in range(k[0] + 1, nrow + 1):
            if idif[k1] > 0:
              goto70bool = True
              break
          if goto70bool:
            break
          iflag[0] = 1
          return
        mm = 1
        for i in range(1, k[0] + 1):
          mm += idif[i]
          idif[i] = 0
        k[0] = k1
        while True: # Fortran line 90
          k[0] -= 1
          m = min(mm, imax[k[0]])
          idif[k[0]] = m
          mm -= m
          if mm > 0 and k[0] != 1:
            continue
          if mm > 0:
            if k1 != nrow:
              k[0] = k1
              goto50bool = True
              break
            iflag[0] = 1
            return
          break
        if goto50bool:
          continue
        idif[k1] = idif[k1] - 1
        ks[0] = 0
        while True: # Fortran line 100
          ks[0] += 1
          if ks[0] > k[0]:
            return
          if idif[ks[0]] >= imax[ks[0]]:
            continue
          return
        break

  def f5xact(self, pastp, tol, kval, key, keyoffset, ldkey, ipoin, ipoinoffset, stp, stpoffset, ldstp, ifrq, ifrqoffset, npoin, npoinoffset, nr, nroffset, nl, nloffset, ifreq, itop, ipsh, itp):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    ird = 1
    ipn = 1
    goto40bool = False
    if ipsh:
      goto30bool = False
      ird = int(kval % ldkey) + 1
      for itp[0] in range(ird, ldkey + 1):
        if key[itp[0] + keyoffset - 1] == kval:
          goto40bool = True
          break
        if key[itp[0] + keyoffset - 1] < 0:
          goto30bool = True
          break
      if goto40bool == False and goto30bool == False:
        for itp[0] in range(1, ird):
          if key[itp[0] + keyoffset - 1] == kval:
            goto40bool = True
            break
          if key[itp[0] + keyoffset - 1] < 0:
            break
      if goto40bool == False:
        key[itp[0] + keyoffset - 1] = kval
        itop[0] += 1
        ipoin[itp[0] + ipoinoffset - 1] = itop[0]
        npoin[itop[0] + npoinoffset - 1] = -1
        nr[itop[0] + nroffset - 1] = -1
        nl[itop[0] + nloffset - 1] = -1
        stp[itop[0] + stpoffset - 1] = pastp
        ifrq[itop[0] + ifrqoffset - 1] = ifreq
        return
    ipn = ipoin[itp[0] + ipoinoffset - 1]
    test1 = pastp - tol
    test2 = pastp + tol
    while True: # Fortran line 50
      if stp[ipn + stpoffset - 1] < test1:
        ipn = nl[ipn + nloffset - 1]
        if ipn > 0:
          continue
      elif stp[ipn + stpoffset - 1] > test2:
        ipn = nr[ipn + nroffset - 1]
        if ipn > 0:
          continue
      else:
        ifrq[ipn + ifrqoffset - 1] = ifrq[ipn + ifrqoffset - 1] + ifreq
        return
      itop[0] += 1
      ipn = ipoin[itp[0] + ipoinoffset - 1]
      itmp = ipn
      break
    while True: # Fortran line 60
      if stp[ipn + stpoffset - 1] < test1:
        itmp = ipn
        ipn = nl[ipn + nloffset - 1]
        if ipn > 0:
          continue
        else:
          nl[itmp + nloffset - 1] = itop[0]
      elif stp[ipn + stpoffset - 1] > test2:
        itmp = ipn
        ipn = nr[ipn + nroffset - 1]
        if ipn > 0:
          continue
        else:
          nr[itmp + nroffset - 1] = itop[0]
      break
    npoin[itop[0] + npoinoffset - 1] = npoin[itmp + npoinoffset - 1]
    npoin[itmp + npoinoffset - 1] = itop[0]
    stp[itop[0] + stpoffset - 1] = pastp
    ifrq[itop[0] + ifrqoffset - 1] = ifreq
    nl[itop[0] + nloffset - 1] = -1
    nr[itop[0] + nroffset - 1] = -1

  def f2xact(self, nrow, ncol, table, ldtabl, expect, percnt, emin, prt, pre, fact, ico, iro, kyy, idif, irn, key, ldkey, ipoin, ldstp, stp, ifrq, dlp, dsp, tm, key2, iwk, rwk):
    """ Supports the MxN table statistics
        Parameters:
          several
        Returns:
          none
    """
    f5itp = [0]
    for i in range(1, 2 * ldkey + 1):
      key[i] = -9999
      key2[i] = -9999
    preops = 0
    ncell = 0
    ifault = 1
    iflag = [1]
    tmp = 0.0
    pv = 0.0
    df = 0.0
    obs2 = 0.0
    obs3 = 0.0
    chisq = False
    pre[0] = 0.0
    itop = [0]
    emn = 0.0
    emx = 1.0e30
    tol = 3.45254e-7
    amiss = -12345.0
    imax = 2147483647
    if expect > 0.0:
      emn = emin
    else:
      emn = emx
    k = ncol
    # Variables for f3xact
    i31 = 1
    i32 = i31 + k
    i33 = i32 + k
    i34 = i33 + k
    i35 = i34 + k
    i36 = i35 + k
    i37 = i36 + k
    i38 = i37 + k
    # Variables for f4xact
    k = nrow + ncol + 1
    i41 = 1
    i42 = i41 + k
    i43 = i42 + k
    i44 = i43 + k
    i45 = i44 + k
    i46 = i45 + k
    if nrow > ldtabl:
      return
    if ncol <= 1:
      return
    ntot = [0]
    for r in range(1, nrow + 1):
      iro[r] = 0
      for c in range(1, ncol + 1):
        if table[r][c] < -0.0001:
          return
        iro[r] += table[r][c]
        ntot[0] += table[r][c]
    riro = [0] * (nrow + 1)
    for r in range(1, nrow + 1):
      riro[r] = iro[r]
    iro = riro
    if ntot[0] == 0:
      prt[0] = amiss;
      pre[0] = amiss;
      return
    for c in range(1, ncol + 1):
      ico[c] = 0
      for r in range(1, nrow + 1):
        ico[c] += table[r][c]
    rico = [0] * (ncol + 1)
    for c in range(1, ncol + 1):
      rico[c] = ico[c]
    ico = rico
    iro.sort()
    ico.sort()
    nro = nrow
    nco = ncol
    rkyy = [0] * (nro + 1)
    rkyy[1] = 1
    mj = ncol
    for r in range(2, nro + 1):
      rkyy[r] = rkyy[r - 1] * (iro[r - 1] + 1)
      mj /= rkyy[r - 1]
    kyy = rkyy
    kmax = (iro[nro] + 1) * kyy[nro - 1]
    fact[0] = 0.0
    fact[1] = 0.0
    fact[2] = math.log(2.0)
    i = 3
    while i <= ntot[0]:
      fact[i] = fact[i - 1] + math.log(float(i))
      mj = i + 1
      if mj <= ntot[0]:
        fact[mj] = fact[i] + fact[2] + fact[int(mj / 2)] - fact[int(mj / 2) - 1]
      i += 2
    obs = tol
    ntot[0] = 0
    dd = 0.0
    for mj in range(1, nco + 1):
      dd = 0.0
      for r in range(1, nro + 1):
        dd += fact[table[r][mj]]
        ntot[0] += table[r][mj]
      obs += fact[ico[mj]] - dd
    dro = self.f9xact(nro, ntot[0], iro, 1, fact)
    prt[0] = math.exp(obs - dro)
    # Initialize pointers
    k = nco
    last = [ldkey + 1]
    jkey = ldkey + 1
    jstp = ldstp + 1
    jstp2 = 3 * ldstp + 1
    jstp3 = 4 * ldstp + 1
    jstp4 = 5 * ldstp + 1
    ikkey = 0
    ikstp = 0
    ikstp2 = 2 * ldstp
    ipo = [1]
    ipoin[1] = 1
    stp[1] = 0.0
    ifrq[1] = 1
    ifrq[ikstp2 + 1] = -1
    while True: # Fortran line 110
      if True:
        goto110bool = False
        kb = nco - k + 1
        ks = [0]
        n = ico[kb] #Ends up being the lowest column total
        kd = [nro + 1]
        kmax = nro
        for i in range(1, nro + 1):
          idif[i] = 0
        while True: # Fortran line 130
          if True:
            kd[0] -= 1 # So kd is now highest index of row totals vector
            ntot[0] = min(n, iro[kd[0]]) # The lowest column total or the highest row total??
            idif[kd[0]] = ntot[0]
            if idif[kmax] == 0:
              kmax -= 1
            n -= ntot[0]
            if n > 0 and kd[0] != 1:
              continue
            k1 = 0
            if n != 0:
              while True: # Fortran line 310
                iflag = [1]
                self.f6xact(nro, iro, iflag, kyy, key, ikkey + 1, ldkey, last, ipo)
                if iflag[0] == 3:
                  k = k - 1
                  itop[0] = 0
                  ikkey = jkey - 1
                  ikstp = jstp - 1
                  ikstp2 = jstp2 - 1
                  jkey = ldkey - jkey + 2
                  jstp = ldstp - jstp + 2
                  jstp2 = 2 * ldstp + jstp
                  for f in range(1,  2 * ldkey + 1):
                    key2[f] = -9999
                  if k >= 2:
                    continue
                  return
                else:
                  goto110bool = True
                  break
              if goto110bool:
                break
            k1 = k - 1
            n = ico[kb]
            ntot[0] = 0
            # kb began as 1 less than the FORTRAN value so this is the same as in FORTRAN
            for i in range(kb + 1, nco + 1):
              ntot[0] += ico[i]
            while True: # Fortran line 150
              if True:
                goto150bool = False
                for i in range(1, nro + 1):
                  irn[i] = iro[i] - idif[i]
                nrb = 0;
                nro2 = 0
                ii = 0
                if k1 > 1:
                  i = nro
                  if nro == 2:
                    if irn[1] > irn[2]:
                      ii = irn[1]
                      irn[1] = irn[2]
                      irn[2] = ii
                  elif nro == 3:
                    ii = irn[1]
                    if ii > irn[3]:
                      if ii > irn[2]:
                        if irn[2] > irn[3]:
                          irn[1] = irn[3]
                          irn[3] = ii
                        else:
                          irn[1] = irn[1]
                          irn[2] = irn[2]
                          irn[3] = ii
                      else:
                        irn[1] = irn[3]
                        irn[3] =  irn[2]
                        irn[2] = ii
                    elif ii > irn[2]:
                      irn[1] = irn[2]
                      irn[2] = ii
                    elif irn[2] > irn[3]:
                      ii = irn[2]
                      irn[2] = irn[3]
                      irn[3] = ii
                  else:
                    for j in range(2, nro + 1):
                      i = j - 1
                      ii = irn[j]
                      while True:
                        if ii < irn[i]:
                          irn[i + 1] = irn[i]
                          i -= 1
                          if i > 0:
                            continue
                        irn[i + 1] = ii
                        break
                  for i in range(1, nro + 1):
                    if irn[i] != 0:
                      break
                  nrb = i
                  nro2 = nro - i + 1
                else:
                  nrb = 1
                  nro2 = nro
                ddf = self.f9xact(nro, n, idif, 1, fact)
                drn = self.f9xact(nro2, ntot[0], irn, nrb, fact) - dro + ddf
                itp = 0
                kval = 1
                if k1 > 1:
                  kval = irn[1] + irn[2] * kyy[2]
                  i = 2
                  for i in range(3, nro + 1):
                    kval += irn[i] * kyy[i]
                  i = kval % (2 * ldkey) + 1
                  bool240 = False
                  for itp in range(i, 2 * ldkey + 1):
                    ii = key2[itp]
                    if ii == kval:
                      bool240 = True
                      break
                    elif ii < 0:
                      key2[itp] = kval
                      dlp[itp] = 1.0
                      dsp[itp] = 1.0
                      bool240 = True
                      break
                  if bool240 is False:
                    for itp in range(1, i):
                      ii = key2[itp]
                      if ii == kval:
                        bool240 = True
                        break
                      elif ii < 0:
                        key2[itp] = kval
                        dlp[itp] = 1.0
                        bool240 = True
                        break
                while True: # Fortran line 240
                  if True:
                    ipsh = True
                    ipn = ipoin[ipo[0] + ikkey]
                    pastp = stp[ipn + ikstp]
                    ifreq = ifrq[ipn + ikstp]
                    if k1 > 1:
                      obs2 = obs - fact[ico[kb + 1]] - fact[ico[kb + 2]] - ddf
                      for i in range(3, k1 + 1):
                        obs2 -= fact[ico[kb + i]]
                      if dlp[itp] > 0.0:
                        dspt = obs - obs2 - ddf
                        dlp[itp] = 0.0
                        iwk0 = []
                        iwk1 = []
                        iwk2 = []
                        iwk3 = []
                        iwk4 = []
                        iwk5 = []
                        iwk6 = []
                        iwk7 = []
                        iwk8 = []
                        for iwki in iwk:
                          iwk0.append(iwki)
                          iwk1.append(iwki)
                          iwk2.append(iwki)
                          iwk3.append(iwki)
                          iwk4.append(iwki)
                          iwk5.append(iwki)
                          iwk6.append(iwki)
                          iwk7.append(iwki)
                          iwk8.append(iwki)
                        rwk0 = []
                        rwk1 = []
                        for rwki in rwk:
                          rwk0.append(rwki)
                          rwk1.append(rwki)
                        dlpitp = [dlp[itp]]
                        self.f3xact(nro2, irn, nrb, k1, ico, kb + 1, dlpitp, ntot, fact, iwk0, iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, iwk7, iwk8, rwk0, rwk1, tol)
                        dlp[itp] = dlpitp[0]
                        dlp[itp] = min(0.0, dlp[itp])
                        dsp[itp] = dspt
                        dspitp = [dsp[itp]]
                        self.f4xact(nro2, irn, nrb, k1, ico, kb + 1, dspitp, fact, iwk0, iwk1, iwk2, iwk3, iwk4, iwk5, iwk6, rwk0, tol)
                        dsp[itp] = dspitp[0]
                        dsp[itp] = min(0.0, dsp[itp] - dspt)
                        # Use chi-squared approximation
                        if float((irn[nrb] * ico[kb + 1]) / float(ntot[0])) > emn:
                          ncell = 0
                          for j in range(1, nro2 + 1):
                            for l in range(1, k1 + 1):
                              if irn[nrb + j - 1] * ico[kb + l] >= ntot[0] * expect:
                                ncell += 1
                          if ncell * 100 >= k1 * nro2 * percnt:
                            tmp = 0.0
                            for j in range(1, nro2 + 1):
                              tmp += fact[irn[nrb + j - 1]] - fact[irn[nrb + j - 1] - 1]
                            tmp *= (k1 - 1)
                            for j in range(1, k1 + 1):
                              tmp += (nro2 - 1) * (fact[ico[kb + j]] - fact[ico[kb + j] - 1])
                            df = (nro2 - 1) * (k1 - 1)
                            tmp = tmp + df * 1.83787706640934548356065947281
                            tmp = tmp - (nro2 * k1 - 1) * (fact[ntot[0]] - fact[ntot[0] - 1])
                            tm[itp] = -2.0 * (obs - dro) - tmp
                          else:
                            tm[itp] = -9876.0
                        else:
                          tm[itp] = -9876.0
                      obs3 = obs2 - dlp[itp]
                      obs2 -= dsp[itp]
                      if tm[itp] == -9876.0:
                        chisq = False
                      else:
                        chisq = True
                        tmp = tm[itp]
                    else:
                      obs2 = obs - drn - dro
                      obs3 = obs2
                    while True: # Fortran line 300
                      if pastp <= obs3:
                        pre[0] += ifreq * math.exp(pastp + drn)
                        preops += 1
                        if preops == 106 or preops == 13:
                          checkpoint = True
                      elif pastp < obs2:
                        if chisq:
                          df = (nro2 - 1) * (k1 - 1)
                          pv = self.gamds(max(0.0, tmp + 2.0 * (pastp + drn)) / 2.0, df / 2.0, ifault)
                          pre[0] += ifreq * math.exp(pastp + drn) * pv
                        else:
                          self.f5xact(pastp + ddf, tol, kval, key, jkey, ldkey, ipoin, jkey, stp, jstp, ldstp, ifrq, jstp, ifrq, jstp2, ifrq, jstp3, ifrq, jstp4, ifreq, itop, ipsh, f5itp)
                          ipsh = False
                      ipn = ifrq[ipn + ikstp2]
                      if ipn > 0:
                        pastp = stp[ipn + ikstp]
                        ifreq = ifrq[ipn + ikstp]
                        continue
                      self.f7xact(kmax, iro, idif, kd, ks, iflag)
                      if iflag[0] != 1:
                        goto150bool = True
                        break
                      while True: # Fortran line 310
                        iflag[0] = 1
                        self.f6xact(nro, iro, iflag, kyy, key, ikkey + 1, ldkey, last, ipo)
                        if iflag[0] == 3:
                          k -= 1
                          itop[0] = 0
                          ikkey = jkey - 1
                          ikstp = jstp - 1
                          ikstp2 = jstp2 - 1
                          jkey = ldkey - jkey + 2
                          jstp = ldstp - jstp + 2
                          jstp2 = 2 * ldstp + jstp
                          for f in range(1, 2 * ldkey + 1):
                            key2[f] = -9999
                          if k >= 2:
                            continue
                        else:
                          goto110bool = True
                          break
                        break
                      break
                      if goto110bool:
                        break
                    if goto150bool:
                      break
                    if goto110bool:
                      break
                    break
                if goto150bool:
                  continue
                if goto110bool:
                  break
                break
            if goto110bool:
              break
            break
        if goto110bool:
          continue
        break

  def FEXACT(self, SortedRows):
    """ Computes the MxN table statistics
        Parameters:
          SortedRows (list(list)): A list of lists of ints
        Returns:
          float
    """
    table0 = []
    for i in range(len(SortedRows) + 1):
      nsma = [0]
      if i == 0:
        for j in range(len(SortedRows[0])):
          nsma.append(0)
        table0.append(nsma)
        continue
      dr0 = SortedRows[i - 1]
      for d in dr0:
        nsma.append(d)
      table0.append(nsma)
    table = table0
    if len(table) > len(table[0]):
      table = self.transposeMatrix(table0)
    ncol = len(table[0]) - 1
    nrow = len(table) - 1
    ldtabl = nrow
    expect = 0.0
    percent = 90.0
    emin = 1.0
    prt = [0.0]
    pre = [0.0]
    iwkmax = 200000
    mult = 30
    ireal = 4
    iwkpt = [1]
    ntot = 0
    for r in range(1, nrow + 1):
      for c in range(1, ncol + 1):
        ntot += table[r][c]
    nco = ncol
    nro = nrow
    k = nrow + ncol + 1
    kk = k * ncol
    ldkey = 0
    ldstp = 0
    numb = 0
    i1 = self.iwork(iwkmax, iwkpt, ntot + 1, ireal) - 1
    i2 = self.iwork(iwkmax, iwkpt, nco, 2) - 1
    i3 = self.iwork(iwkmax, iwkpt, nco, 2) - 1
    i3a = self.iwork(iwkmax, iwkpt, nco, 2) - 1
    i3b = self.iwork(iwkmax, iwkpt, nro, 2) - 1
    i3c = self.iwork(iwkmax, iwkpt, nro, 2) - 1
    iiwk = self.iwork(iwkmax, iwkpt, max(5 * k + 2 * kk, 800 + 7 * ncol), 2) - 1
    irwk = self.iwork(iwkmax, iwkpt, max(400 + ncol + 1, k), ireal) - 1
    if ireal == 4:
      numb = 18 + 10 * mult
      ldkey = int((iwkmax - iwkpt[0] + 1) / numb)
    else:
      numb = 12 * 8 * mult
      ldkey = int((iwkmax - iwkpt[0] + 1) / numb)
    ldstp = mult * ldkey
    i4 = self.iwork(iwkmax, iwkpt, 2 * ldkey, 2) - 1
    i5 = self.iwork(iwkmax, iwkpt, 2 * ldkey, 2) - 1
    i6 = self.iwork(iwkmax, iwkpt, 2 * ldstp, ireal) - 1
    i7 = self.iwork(iwkmax, iwkpt, 6 * ldstp, 2) - 1
    i8 = self.iwork(iwkmax, iwkpt, 2 * ldkey, ireal) - 1
    i9 = self.iwork(iwkmax, iwkpt, 2 * ldkey, ireal) - 1
    i9a = self.iwork(iwkmax, iwkpt, 2 * ldkey, ireal) - 1
    i10 = self.iwork(iwkmax, iwkpt, 2 * ldkey, 2) - 1
    i1array = [0.0] * irwk
    i2array = [0] * (i3 - i2 + 1)
    i3array = [0] * (i3a - i3 + 1)
    i3aarray = [0] * (i3b - i3a + 1)
    i3barray = [0] * (i3c - i3b + 1)
    i3carray = [0] * (iiwk - i3c + 1)
    i4array = [0] * (i5 - i4 + 1)
    i5array = [0] * (i7 - i5 + 1)
    i6array = [0.0] * (i8 - i6 + 1)
    i7array = [0] * (i10 - i7 + 1)
    i8array = [0.0] * (i9 - i8 + 1)
    i9array = [0.0] * (i9a - i9 + 1)
    i9aarray = [0.0] * (iwkmax - i9a + 1)
    i10array = [0] * int(2 * ldkey + 1)
    iiwkarray = [0] * (i4 - iiwk + 1)
    irwkarray = [0.0] * (i6 - irwk + 1)
    self.f2xact(nrow, ncol, table, ldtabl, expect, percent, emin, prt, pre, i1array, i2array, i3array, i3aarray, i3barray, i3carray,
                i4array, int(ldkey), i5array, ldstp, i6array, i7array, i8array, i9array, i9aarray, i10array, iiwkarray, irwkarray)
    return pre[0]

  def MXNCompute(self, table):
    """ Computes the MxN table statistics
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          dict
    """
    stats = {}
    if len(table) < 2 or len(table[0]) < 2:
      return stats
    chiSq = 0.0
    lowExpectation = False
    totals = [None] * (len(table[0]) + 1)
    totals[len(table[0])] = 0.0
    rowtotals = [None] * len(table)
    for i in range(len(table)):
      rowtotals[i] = 0.0
      for j in range(len(table[0])):
        if i == 0:
          totals[j] = 0.0
        totals[j] += table[i][j]
        totals[len(table[0])] += table[i][j]
        rowtotals[i] += table[i][j]
    ps = [0.0] * len(table[0])
    for i in range(len(table[0])):
      ps[i] = totals[i] / totals[len(table[0])]
    for i in range(len(table)):
      for j in range(len(table[0])):
        observed = table[i][j]
        expected = rowtotals[i] * ps[j]
        OminusESqOverE = (observed - expected) ** 2.0 / expected;
        chiSq += OminusESqOverE
        if expected < 5.0:
          lowExpectation = True
    chiSqP = self.PValFromChiSq(chiSq, (len(table) - 1) * (len(table[0]) - 1))
    stats['ChiSq'] = chiSq
    stats['ChiSqDF'] = (len(table) - 1) * (len(table[0]) - 1)
    stats['ChiSqP'] = chiSqP
    stats['FishersExact'] = self.FEXACT(table)
    return stats

  def Subset(self, dataTable, column, values):
    """ Executes the supporting functions to run the analysis
        Parameters:
          dataTable (list(dict)): The analysis dataset
          column (string): The value on which to subset
          values (list): The values to keep
        Returns:
          None. Stores results as class variables.
    """
    newdt = []
    for row in dataTable:
      if row[column] in values:
        newdt.append(row)
    return newdt

  def ZSElnOR(self, a, b, c, d):
    p1, p2, p3, p4, p5, ZValue = 0.0, 0.0, 0.0, 0.0, 0.0, 1.96
    for i in range(len(a)):
        p1 += ((a[i] + d[i]) / (a[i] + b[i] + c[i] + d[i])) * (a[i] * d[i] / (a[i] + b[i] + c[i] + d[i]))
        p2 += ((a[i] + d[i]) / (a[i] + b[i] + c[i] + d[i])) * (b[i] * c[i] / (a[i] + b[i] + c[i] + d[i]))
        p2 += ((b[i] + c[i]) / (a[i] + b[i] + c[i] + d[i])) * (a[i] * d[i] / (a[i] + b[i] + c[i] + d[i]))
        p3 += ((b[i] + c[i]) / (a[i] + b[i] + c[i] + d[i])) * (b[i] * c[i] / (a[i] + b[i] + c[i] + d[i]))
        p4 += (a[i] * d[i] / (a[i] + b[i] + c[i] + d[i]))
        p5 += (b[i] * c[i] / (a[i] + b[i] + c[i] + d[i]))
    SElnOR = (p1 / (2 * p4 * p4) + p2 / (2 * p4 * p5) + p3 / (2 * p5 * p5)) ** 0.5
    return ZValue * SElnOR

  def ZSElnRR(self, a, b, c, d):
    numerator, denom1, denom2, ZValue = 0.0, 0.0, 0.0, 1.96
    for i in range(len(a)):
        numplus = ((a[i] + c[i]) * (a[i] + b[i]) * (c[i] + d[i]) - a[i] * c[i] * (a[i] + b[i] + c[i] + d[i]))
        numplus /= ((a[i] + b[i] + c[i] + d[i]) * (a[i] + b[i] + c[i] + d[i]))
        numerator += numplus
        denom1 += (a[i] * (c[i] + d[i])) / (a[i] + b[i] + c[i] + d[i])
        denom2 += (c[i] * (a[i] + b[i])) / (a[i] + b[i] + c[i] + d[i])
    SElnOR = (numerator / (denom1 * denom2)) ** 0.5
    return ZValue * SElnOR

  def strat2x2(self, a, b, c, d):
    numerator, denominator = 0.0, 0.0
    for i in range(len(a)):
        numerator += (a[i] * d[i]) / (a[i] + b[i] + c[i] + d[i])
        denominator += (b[i] * c[i]) / (a[i] + b[i] + c[i] + d[i])
    return numerator / denominator

  def ComputedRR(self, a, b, c, d):
    numerator, denominator = 0.0, 0.0
    for i in range(len(a)):
        numerator += (a[i] * (c[i] + d[i])) / (a[i] + b[i] + c[i] + d[i])
        denominator += (c[i] * (a[i] + b[i])) / (a[i] + b[i] + c[i] + d[i])
    return numerator / denominator

  def exactorln(self, aas, bbs, ccs, dds):
    M1, M0, N1, N0 = [0.0] * len(aas), [0.0] * len(aas), [0.0] * len(aas), [0.0] * len(aas)
    ls, us, xs = [0.0] * len(aas), [0.0] * len(aas), []
    x, l, u, maxuldiff = 0, 0, 0, 0.0
    for i in range(len(aas)):
        xs.append(aas[i])
        if aas[i] > 999999 or bbs[i] > 999999 or ccs[i] > 999999 or dds[i] > 999999:
            return float('nan')
    for i in range(len(aas)):
        M1[i] = aas[i] + ccs[i]
        M0[i] = bbs[i] + dds[i]
        N1[i] = aas[i] + bbs[i]
        N0[i] = ccs[i] + dds[i]
        ls[i] = max(0, N1[i] - M0[i])
        us[i] = min(M1[i], N1[i])
        x += xs[i]
        l += ls[i]
        u += us[i]
    Cs = []
    dimC2 = 0
    C2 = []
    for i in range(len(aas)):
        dimC2 += us[i] - ls[i]
        Cs.append([])
    for i in range(dimC2 + 1):
        C2.append(0.0)
    maxuldiff = 0
    maxusi = us[0]
    for i in range(len(aas)):
        if us[i] > maxusi:
            maxusi = us[i]
    for i in range(len(aas)):
        if us[i] - ls[i] > maxuldiff:
            maxuldiff = us[i] - ls[i]
            Cs[i] = [0] * (maxusi + 1)
        for s in range(ls[i], us[i] + 1):
            Cs[i][s - ls[i]] = self.choosey(M1[i], s) * self.choosey(M0[i], N1[i] - s)
    Y = [0.0] * (int(u) - int(l) + 1)
    for j in range(us[0] - ls[0] + 1):
        for k in range(us[1] - ls[1] + 1):
            C2[j + k] += Cs[0][j] * Cs[1][k]
    bound = 0
    for i in range(2, len(aas)):
        for j in range(0, u - l + 1):
            Y[j] = C2[j]
            C2[j] = 0.0
        bound = 0
        for j in range(i - 1 + 1):
            bound += us[i] - ls[i]
        for j in range(u - l - 1):
            for k in range(us[i] - ls[i] + 1):
                if j + k <= u - l:
                    C2[j + k] += Y[j] * Cs[i][k]
    R = 0.0
    Ds = [0.0] * len(aas)
    d2 = 1.0
    FR = 1.0
    adder = 0.0
    while FR > 0.975:
        for i in range(len(aas)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R += 1
        for j in range(len(aas)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        for i in range((x - 1) - l + 1):
            adder = C2[i]
            for j in range(len(Ds) - 1 + 1):
                adder /= Ds[j]
            adder *= R ** (i + l)
            FR += adder
    aa = R - 1.0
    bb = R + 0.5
    precision = 0.00001
    while bb - aa > precision:
        for i in range(len(aas)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R = (bb + aa) / 2.0
        for j in range(len(aas)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        for i in range((x - 1) - l + 1):
            adder = C2[i]
            for j in range(len(Ds) - 1 + 1):
                adder /= Ds[j]
            adder *= R ** (i + l)
            FR += adder
        if FR > 0.975:
            aa = R
        else:
            bb = R
    return R

  def exactorun(self, aas, bbs, ccs, dds, minimum):
    M1, M0, N1, N0 = [0.0] * len(aas), [0.0] * len(aas), [0.0] * len(aas), [0.0] * len(aas)
    ls, us, xs = [0.0] * len(aas), [0.0] * len(aas), []
    x, l, u, maxuldiff = 0, 0, 0, 0.0
    for i in range(len(aas)):
        xs.append(aas[i])
        if aas[i] > 999999 or bbs[i] > 999999 or ccs[i] > 999999 or dds[i] > 999999:
            return float('nan')
    for i in range(len(aas)):
        M1[i] = aas[i] + ccs[i]
        M0[i] = bbs[i] + dds[i]
        N1[i] = aas[i] + bbs[i]
        N0[i] = ccs[i] + dds[i]
        ls[i] = max(0, N1[i] - M0[i])
        us[i] = min(M1[i], N1[i])
        x += xs[i]
        l += ls[i]
        u += us[i]
    Cs = []
    dimC2 = 0
    C2 = []
    for i in range(len(aas)):
        dimC2 += us[i] - ls[i]
        Cs.append([])
    for i in range(dimC2 + 1):
        C2.append(0.0)
    maxuldiff = 0
    maxusi = us[0]
    for i in range(len(aas)):
        if us[i] > maxusi:
            maxusi = us[i]
    for i in range(len(aas)):
        if us[i] - ls[i] > maxuldiff:
            maxuldiff = us[i] - ls[i]
            Cs[i] = [0] * (maxusi + 1)
        for s in range(ls[i], us[i] + 1):
            Cs[i][s - ls[i]] = self.choosey(M1[i], s) * self.choosey(M0[i], N1[i] - s)
    Y = [0.0] * (int(u) - int(l) + 1)
    for j in range(us[0] - ls[0] + 1):
        for k in range(us[1] - ls[1] + 1):
            C2[j + k] += Cs[0][j] * Cs[1][k]
            if C2[j + k] == float('inf'):
                return float('nan')
    bound = 0
    for i in range(2, len(aas)):
        for j in range(0, u - l + 1):
            Y[j] = C2[j]
            C2[j] = 0.0
        bound = 0
        for j in range(i - 1 + 1):
            bound += us[i] - ls[i]
        for j in range(u - l - 1):
            for k in range(us[i] - ls[i] + 1):
                if j + k <= u - l:
                    C2[j + k] += Y[j] * Cs[i][k]
                if j + k <= u - l:
                    if C2[j + k] == float('inf'):
                        return float('nan')
    R = minimum - 0.1
    Ds = [0.0] * len(aas)
    d2 = 1.0
    FR = 0.0
    MiddlePart = 1.0
    while True:
        for i in range(len(aas)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R += 0.1
        for j in range(len(aas)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        MiddlePart = 1 / Ds[0]
        for i in range(x - l + 1):
            if i + l > 0:
                MiddlePart = R / Ds[0]
            for j in range(2, i + l + 1):
                MiddlePart *= R
            for j in range(1, len(aas)):
                MiddlePart /= Ds[j]
            FR += C2[i] * MiddlePart
        if FR <= 0.025:
            break
    aa = R - 1.0
    bb = R + 0.5
    precision = 0.00001
    while bb - aa > precision:
        for i in range(len(aas)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R = (bb + aa) / 2.0
        for j in range(len(aas)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        MiddlePart = 1 / Ds[0]
        for i in range(x - l + 1):
            if i + l > 0:
                MiddlePart = R / Ds[0]
            for j in range(2, i + l + 1):
                MiddlePart *= R
            for j in range(1, len(aas)):
                MiddlePart /= Ds[j]
            FR += C2[i] * MiddlePart
        if FR > 0.025:
            aa = R
        else:
            bb = R
    return R

  def choosey(self, chooa, choob):
    ccccc = chooa - choob
    if choob < chooa / 2:
        choob = ccccc
    choosey = 1.0
    for i in range(choob + 1, chooa + 1):
                   choosey = (choosey * i) / (chooa - (i - 1))
    return choosey

  def ucestimaten(self, a, b, c, d):
    M1, M0, N1, N0, ls, us, xs, x, l, u, maxuldiff = [], [], [], [], [], [], [], 0, 0, 0, 0.0
    for i in range(len(a)):
        if a[i] > 999999 or b[i] > 999999 or c[i] > 999999 or d[i] > 999999:
            return float('nan')
    for i in range(len(a)):
        xs.append(a[i])
    for i in range(len(a)):
        M1.append(a[i] + c[i])
        M0.append(b[i] + d[i])
        N1.append(a[i] + b[i])
        N0.append(c[i] + d[i])
        ls.append(max(0, N1[i] - M0[i]))
        us.append(min(M1[i], N1[i]))
        x += xs[i]
        l += ls[i]
        u += us[i]
    Cs = []
    dimC2 = 0
    C2 = []
    for i in range(len(a)):
        dimC2 += us[i] - ls[i]
        Cs.append([])
    for i in range(dimC2 + 1):
        C2.append(0.0)
    maxuldiff = 0
    maxusi = us[0]
    for i in range(len(a)):
        if us[i] > maxusi:
            maxusi = us[i]
    for i in range(len(a)):
        if us[i] - ls[i] > maxuldiff:
            maxuldiff = us[i] - ls[i]
            Cs[i] = [0] * (maxusi + 1)
        for s in range(ls[i], us[i] + 1):
            Cs[i][s - ls[i]] = self.choosey(M1[i], s) * self.choosey(M0[i], N1[i] - s)
    Y = [0.0] * (u - l + 1)
    for j in range(us[0] - ls[0] + 1):
        for k in range(us[1] - ls[1] + 1):
            C2[j + k] += Cs[0][j] * Cs[1][k]
    bound = 0
    for i in range(2, len(a)):
        for j in range(0, u - l + 1):
            Y[j] = C2[j]
            C2[j] = 0.0
        bound = 0
        for j in range(i - 1 + 1):
            bound += us[i] - ls[i]
        for j in range(u - l - 1):
            for k in range(us[i] - ls[i] + 1):
                if j + k <= u - l:
                    C2[j + k] += Y[j] * Cs[i][k]
    R = 0.0
    Ds = [0.0] * len(a)
    d2 = 1.0
    FR = 0.0
    adder = 0.0
    while FR < x:
        for i in range(len(a)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R += 1
        for j in range(len(a)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        for i in range(u - l + 1):
            adder = (i + l) * C2[i]
            for j in range(len(Ds) - 1 + 1):
                adder /= Ds[j]
            adder *= R ** (i + l)
            FR += adder
    aa = R - 1.0
    bb = R + 0.5
    precision = 0.00001
    while bb - aa > precision:
        for i in range(len(a)):
            Ds[i] = 0.0
        d2 = 1.0
        FR = 0.0
        R = (bb + aa) / 2.0
        for j in range(len(a)):
            for i in range(us[j] - ls[j] + 1):
                Ds[j] += Cs[j][i] * R ** (ls[j] + i)
            d2 *= Ds[j]
        for i in range(u - l + 1):
            adder = (i + l) * C2[i]
            for j in range(len(Ds) - 1 + 1):
                adder /= Ds[j]
            adder *= R ** (i + l)
            FR += adder
        if FR < x:
            aa = R
        else:
            bb = R
    return R

  def ComputeUnChisq(self, a, b, c, d):
    numerator, denominator = 0.0, 0.0
    for i in range(len(a)):
        numerator += (a[i] * d[i] - b[i] * c[i]) / (a[i] + b[i] + c[i] + d[i])
        dplus = ((a[i] + b[i]) * (c[i] + d[i]) * (a[i] + c[i]) * (b[i] + d[i]))
        dplus /= (((a[i] + b[i] + c[i] + d[i]) - 1) * (a[i] + b[i] + c[i] + d[i]) * (a[i] + b[i] + c[i] + d[i]))
        denominator += dplus
    return (numerator * numerator) / denominator

  def ComputeCorrChisq(self, a, b, c, d):
    numerator, denominator = 0.0, 0.0
    for i in range(len(a)):
        numerator += (a[i] * d[i] - b[i] * c[i]) / (a[i] + b[i] + c[i] + d[i])
        dplus = ((a[i] + b[i]) * (c[i] + d[i]) * (a[i] + c[i]) * (b[i] + d[i]))
        dplus /= (((a[i] + b[i] + c[i] + d[i]) - 1) * (a[i] + b[i] + c[i] + d[i]) * (a[i] + b[i] + c[i] + d[i]))
        denominator += dplus
    return ((abs(numerator) - 0.5) * (abs(numerator) - 0.5)) / denominator

  def bdOR(self, yy, yn, ny, nn):
    bd = 0.0
    bor = [0.0] * len(yy)
    w = [0.0] * len(yy)
    for i in range(len(yy)):
        bor[i] = (yy[i] / yn[i]) / (ny[i] / nn[i])
        w[i] = 1 / (1 / yy[i] + 1 / yn[i] + 1 / ny[i] + 1 / nn[i])
    numerator, denominator = 0.0, 0.0
    for i in range(len(yy)):
        numerator += w[i] * math.log(bor[i])
        denominator += w[i]
    orDirect = math.exp(numerator / denominator)
    for i in range(len(yy)):
        bd += (math.log(bor[i]) - math.log(orDirect)) ** 2.0 * w[i]
    return bd

  def bdtOR(self, yy, yn, ny, nn, mleOR):
    bd, sumYY, sumAk, sumVar = 0.0, 0.0, 0.0, 0.0
    for i in range(0, len(yy)):
        N1k = yy[i] + ny[i]
        N0k = yn[i] + nn[i]
        tk = yy[i] + yn[i]
        Nk = yy[i] + yn[i] + ny[i] + nn[i]
        a = 0.0
        b = min(tk, N1k)
        Ak = (a + b) / 2.0
        precision = 0.00001
        psi = 0.0
        while True:
            Ak = (a + b) / 2.0
            psi = (Ak * (N0k - tk + Ak)) / ((N1k - Ak) * (tk - Ak))
            if psi < 0.0 or psi < mleOR:
                a = Ak
            else:
                b = Ak
            if abs(mleOR - psi) < precision:
                break
        var = 1 / (1 / Ak + 1 / (N1k - Ak) + 1 / (tk - Ak) + 1 / (N0k - tk + Ak))
        bd += ((yy[i] - Ak) * (yy[i] - Ak)) / var
        sumYY += yy[i]
        sumAk += Ak
        sumVar += var
    correction = ((sumYY - sumAk) * (sumYY - sumAk)) / sumVar
    bd -= correction
    return bd

  def bdRR(self, yy, yn, ny, nn):
    bd = 0.0
    rr = [0.0] * len(yy)
    w = [0.0] * len(yy)
    for i in range(len(yy)):
        rr[i] = (yy[i] / (yy[i] + yn[i])) / (ny[i] / (ny[i] + nn[i]))
        w[i] = 1 / (yn[i] / (yy[i] * (yy[i] + yn[i])) + nn[i] / (ny[i] * (ny[i] + nn[i])))
    numerator, denominator = 0.0, 0.0
    for i in range(len(yy)):
        numerator += w[i] * math.log(rr[i])
        denominator += w[i]
    rrDirect = math.exp(numerator / denominator)
    for i in range(len(yy)):
        bd += (math.log(rr[i]) - math.log(rrDirect)) ** 2.0 * w[i]
    return bd

  def Summarize(self, yyArray, ynArray, nyArray, nnArray):
    ret = {}
    cumulativeYY = sum(yyArray)
    cumulativeYN = sum(ynArray)
    cumulativeNY = sum(nyArray)
    cumulativeNN = sum(nnArray)
    ret['computedOddsRatio'] = self.strat2x2(yyArray, ynArray, nyArray, nnArray)
    ret['computedOddsRatioMHLL'] = ret['computedOddsRatio'] *  math.exp(-self.ZSElnOR(yyArray, ynArray, nyArray, nnArray))
    ret['computedOddsRatioMHUL'] = ret['computedOddsRatio'] *  math.exp(self.ZSElnOR(yyArray, ynArray, nyArray, nnArray))
    ret['computedRR'] = self.ComputedRR(yyArray, ynArray, nyArray, nnArray)
    ret['computedRRMHLL'] = ret['computedRR'] *  math.exp(-self.ZSElnRR(yyArray, ynArray, nyArray, nnArray))
    ret['computedRRMHUL'] = ret['computedRR'] *  math.exp(self.ZSElnRR(yyArray, ynArray, nyArray, nnArray))
    ret['mleOR'] = float('nan')
    ret['ExactORLL'] = float('nan')
    ret['ExactORUL'] = float('nan')
    if cumulativeYN == 0.0 or cumulativeNY == 0.0:
        ret['mleOR'] = float('inf')
        ret['ExactORLL'] = self.exactorln(yyArray, ynArray, nyArray, nnArray)
        ret['ExactORUL'] = float('inf')
    elif cumulativeYY == 0.0 or cumulativeNN == 0.0:
        ret['mleOR'] = 0.0
        ret['ExactORLL'] = 0.0
        ret['ExactORUL'] = self.exactorun(yyArray, ynArray, nyArray, nnArray, ret['mleOR'])
    else:
        ret['mleOR'] = self.ucestimaten(yyArray, ynArray, nyArray, nnArray)
        ret['ExactORLL'] = self.exactorln(yyArray, ynArray, nyArray, nnArray)
        ret['ExactORUL'] = self.exactorun(yyArray, ynArray, nyArray, nnArray, ret['mleOR'])
    ret['uncorrecedChiSquare'] = self.ComputeUnChisq(yyArray, ynArray, nyArray, nnArray)
    ret['uncorrectedChiSquareP'] = self.PValFromChiSq(float(ret['uncorrecedChiSquare']), 1)
    ret['corrChisq'] = self.ComputeCorrChisq(yyArray, ynArray, nyArray, nnArray)
    ret['corrChisqP'] = self.PValFromChiSq(float(ret['corrChisq']), 1)
    ret['bdtOR'] = self.bdtOR(yyArray, ynArray, nyArray, nnArray, float(ret['computedOddsRatio']))
    ret['bdtORP'] = self.PValFromChiSq(float(ret['bdtOR']), len(yyArray) - 1)
    ret['bdOR'] = self.bdOR(yyArray, ynArray, nyArray, nnArray)
    ret['bdORP'] = self.PValFromChiSq(float(ret['bdOR']), len(yyArray) - 1)
    ret['bdRR'] = self.bdRR(yyArray, ynArray, nyArray, nnArray)
    ret['bdRRP'] = self.PValFromChiSq(float(ret['bdRR']), len(yyArray) - 1)
    return ret

  def Run(self, inputVariableList, dataTable):
    """ Executes the supporting functions to run the analysis
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
          dataTable (list(dict)): The analysis dataset
        Returns:
          None. Stores results as class variables.
    """
    if 'includeMissing' in inputVariableList:
      self.includeMissing = inputVariableList['includemissing']
    try:
      self.useCommonReference = inputVariableList['usecommonreference']
    except KeyError:
      pass
    variablenames = []
    values = []
    tables = []
    rowtotals = []
    coltotals = []
    statistics = []
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
      if len(outcomeValues) == 2 and len(exposureValues) == 2:
        tables.append(onetable)
        rowtotals.append(self.RowTotals(onetable))
        coltotals.append(self.ColumnTotals(onetable))
        statistics.append(self.TwoX2Compute(onetable))
        variablenames.append(inputVariableList['outcomeVariable'] + ' * ' + exposure)
        values.append([outcomeValues, exposureValues])
      elif len(outcomeValues) == 2 and len(exposureValues) > 2 and self.useCommonReference:
        ivl0 = {}
        for ivlk in inputVariableList:
          ivl0[ivlk] = inputVariableList[ivlk]
        ivl0['exposureVariables'] = [exposure]
        for ev in exposureValues[1:]:
          crresult = self.Run(ivl0, self.Subset(dataTable, exposure, [exposureValues[0], ev]))
          variablenames += crresult['Variables']
          values += crresult['VariableValues']
          tables += crresult['Tables']
          rowtotals += crresult['RowTotals']
          coltotals += crresult['ColumnTotals']
          statistics += crresult['Statistics']
      else:
        tables.append(onetable)
        rowtotals.append(self.RowTotals(onetable))
        coltotals.append(self.ColumnTotals(onetable))
        statistics.append(self.MXNCompute(onetable))
        variablenames.append(inputVariableList['outcomeVariable'] + ' * ' + exposure)
        values.append([outcomeValues, exposureValues])
    return {'Variables' : variablenames, 'VariableValues' : values, 'Tables' : tables, 'RowTotals' : rowtotals, 'ColumnTotals' : coltotals, 'Statistics' : statistics}
