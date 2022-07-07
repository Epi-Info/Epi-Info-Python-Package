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
      f2 = self.Func(x2, sumA, value, degN, degD, bigPolyD, bigPolyN)
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

  def ExactTests(self, minSumA, sumA, degD, ExactTests, bigPolyD):
    """ Support function for exact odds ratios
        Parameters:
          table (list(list)): A list of lists of ints
        Returns:
          list
    """
    diff = int(sumA - minSumA)
    bigUpTail = BigDouble('logValue', bigPolyD[int(degD)].logValue)
    bigTwoTail = BigDouble('doubleValue', 1.0)
    if bigUpTail.logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
      bigTwoTail.plusLog(bigUpTail.logValue)
    i = int(degD) - 1
    while i >= diff:
      bigUpTail.plusLog(bigPolyD[i].logValue)
      if bigPolyD[i].logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
        bigTwoTail.plusLog(bigPolyD[i].logValue)
      i -= 1
    bigDenom = BigDouble('logValue', bigUpTail.logValue)
    i = diff - 1
    while i >= 0:
      bigDenom.plusLog(bigPolyD[i].logValue)
      if bigPolyD[i].logValue <= math.log10(1.000001) + bigPolyD[diff].logValue:
        bigTwoTail.plusLog(bigPolyD[i].logValue)
      i -= 1
    ExactTests[0] = 1.0 - (10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** (bigPolyD[diff].logValue - bigDenom.logValue))
    ExactTests[1] = 10.0 ** (bigUpTail.logValue - bigDenom.logValue)
    ExactTests[2] = 10.0 ** (bigTwoTail.logValue - bigDenom.logValue)
    ExactTests[3] = 1.0 - (10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** ((math.log10(0.5) + bigPolyD[diff].logValue) - bigDenom.logValue))
    ExactTests[4] = 10.0 ** (bigUpTail.logValue - bigDenom.logValue) - 10.0 ** ((math.log10(0.5) + bigPolyD[diff].logValue) - bigDenom.logValue)

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
    bigPolyD[0] = BigDouble('doubleValue', 1.0)
    aa = minA
    bb = m1 - minA + 1.0
    cc = n1 - minA + 1.0
    dd = n0 - m1 + minA
    for i in range(1, int(maxA - minA) + 1):
      bigPolyD[i] = BigDouble('logValue', bigPolyD[i - 1].timesReturn(((bb - i) / (aa + i)) * ((cc - i) / (dd + i))))
    ExactResults[0] = self.CalcCmle(1.0, minA, table[0][0], maxA, maxA - minA + 1, bigPolyD)
    ExactTestsList = [None] * 5
    self.ExactTests(minA, table[0][0], maxA - minA, ExactTestsList, bigPolyD)
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
    return stats

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
      tables.append(onetable)
      rowtotals.append(self.RowTotals(onetable))
      coltotals.append(self.ColumnTotals(onetable))
      if len(outcomeValues) == 2 and len(exposureValues) == 2:
        statistics.append(self.TwoX2Compute(onetable))
      else:
        statistics.append(self.MXNCompute(onetable))
      variablenames.append(inputVariableList['outcomeVariable'] + ' * ' + exposure)
      values.append([outcomeValues, exposureValues])
    return {'Variables' : variablenames, 'VariableValues' : values, 'Tables' : tables, 'RowTotals' : rowtotals, 'ColumnTotals' : coltotals, 'Statistics' : statistics}
