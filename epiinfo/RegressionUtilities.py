""" This file defines several supporting classes for the
    Epi Info Regression Routines
"""
class EIMatrix:
  def __init__(self):
    self._mboolFirst = True
    self._mboolIntercept = True
    self._mstrMatchVar = None # str
    self._matchGroupValues = 0
    self._mdblaJacobian = [] # list
    self._mdblaInv = [] # list
    self._mdblaB = [] # list
    self._mmdblaF = [] # list
    self._mboolConverge = True
    self._mboolErrorStatus = True
    self._mdblllfst = 0.0
    self._mdbllllast = 0.0
    self._mdblScore = 0.0
    self._mintIterations = 0
    self._lstrError = None # str

  def get_mboolFirst(self):
    return self._mboolFirst
  def set_mboolFirst(self, v):
    self._mboolFirst = v

  def get_mboolIntercept(self):
    return self._mboolIntercept
  def set_mboolIntercept(self, v):
    self._mboolIntercept = v

  def get_mstrMatchVar(self):
    return self._mstrMatchVar
  def set_mstrMatchVar(self, v):
    self._mstrMatchVar = v

  def get_matchGroupValues(self):
    return self._matchGroupValues
  def set_matchGroupValues(self, v):
    self._matchGroupValues = v

  def get_mdblaJacobian(self):
    return self._mdblaJacobian
  def set_mdblaJacobian(self, v):
    self._mdblaJacobian = v

  def get_mdblaInv(self):
    return self._mdblaInv
  def set_mdblaInv(self, v):
    self._mdblaInv = v

  def get_mdblaB(self):
    return self._mdblaB
  def set_mdblaB(self, v):
    self._mdblaB = v

  def get_mdblaF(self):
    return self._mdblaF
  def set_mdblaF(self, v):
    self._mdblaF = v

  def get_mboolConverge(self):
    return self._mboolConverge
  def set_mboolConverge(self, v):
    self._mboolConverge = v

  def get_mboolErrorStatus(self):
    return self._mboolErrorStatus
  def set_mboolErrorStatus(self, v):
    self._mboolErrorStatus = v

  def get_mdblllfst(self):
    return self._mdblllfst
  def set_mdblllfst(self, v):
    self._mdblllfst = v

  def get_mdbllllast(self):
    return self._mdbllllast
  def set_mdbllllast(self, v):
    self._mdbllllast = v

  def get_mdblScore(self):
    return self._mdblScore
  def set_mdblScore(self, v):
    self._mdblScore = v

  def get_mintIterations(self):
    return self._mintIterations
  def set_mintIterations(self, v):
    self._mintIterations = v

  def get_lstrError(self):
    return self._lstrError
  def set_lstrError(self, v):
    self._lstrError = v

  def lubksb(self, a, n, indx, B):
    """ Performs calculations and stores the results in a list.
        Parameters:
          a (list)
          n (int)
          indx (list)
          B (list)
    """
    ii = 0
    ip = 0
    sum = 0.0
    for i in range(0, n):
      ip = indx[i + 1]
      sum = B[ip]
      B[ip] = B[i + 1]
      if ii > 0:
        for j in range(ii, i + 1):
          sum -= float(a[i][j-1]) * B[j]
      elif sum != 0.0:
        ii = i + 1
    i = n - 1
    while i >= 0:
      sum = B[i + 1]
      for j in range(i + 1, n):
        sum -= float(a[i][j]) * B[j + 1]
      B[i + 1] = sum / float(a[i][i])
      i -= 1

  def fabs(self, a):
    """ Computes and returns absolute value
        Parameters:
          a (int of float))
        Returns: Absolute value of a
    """
    if a < 0:
      return -a
    return a

  def ludcmp(self, a, n, indx, d):
    """ Performs calculations and stores the results in a list.
        Parameters:
          a (list)
          n (int)
          indx (list)
          d (list): A mutable double address in ObjC; a list in Python
        Returns: 1 or -1
    """
    d[0] = 44.5
    TINY = 1.0e-20

    imax = 0
    dum = 0.0
    big = 0.0
    sum = 0.0
    vv = []

    d[0] = 1.0
    for i in range(0, n):
      big = 0.0
      for j in range(0, n):
        absaij = self.fabs(float(a[i][j]))
        if absaij > big:
          big = absaij
      if big == 0.0:
        return -1
      vv.append(1.0 / big)
    for j in range(0, n):
      for i in range(0, j):
        sum = float(a[i][j])
        for k in range(0, i):
          sum -= float(a[i][k]) * float(a[k][j])
        a[i][j] = sum
      big = 0.0
      for i in range(j, n):
        sum = float(a[i][j])
        for k in range(0, j):
          sum -= float(a[i][k]) * float(a[k][j])
        a[i][j] = sum
        dum = vv[i] * self.fabs(sum)
        if dum >= big:
          big = dum
          imax = i
      if j != imax:
        for k in range(0, n):
          dum = float(a[imax][k])
          a[imax][k] = float(a[j][k])
          a[j][k] = dum
        d[0] = -1 * d[0]
        vv[imax] = vv[j]
      indx[j] = imax + 1
      if float(a[i][j]) == 0.0:
        a[j][j] = TINY
      if j != n:
        dum = 1.0 / float(a[j][j])
        for i in range(j + 1, n):
          oldValue = float(a[i][j])
          newValue = oldValue * dum
          a[i][j] = newValue
    return 1

  def inv(self, a, invA):
    """ Inverts a matrix and stores the results in a list.
        Parameters:
          a (list of lists)
          invA (list of lists)
        Returns: none
    """
    n = len(a[0])
    indx = []
    col = []
    for i in range(0, n + 2):
      indx.append(0)
      col.append(0.0)
    d = [0.0]
    self.ludcmp(a, n, indx, d)
    for j in range(0, n):
      for i in range(0, n):
        col[i] = 0.0
      col[j] = 1.0
      indxShifted = []
      colShifted = []
      for i in range(0, n + 3):
        indxShifted.append(0)
        colShifted.append(0.0)
      for k in range(0, n + 2):
        indxShifted[k + 1] = indx[k]
        colShifted[k + 1] = col[k]
      self.lubksb(a, n, indxShifted, colShifted)
      for i in range(0, n):
        inbA[i][j] = float(colShifted[i + 1])

  def Conditional(self, lintOffset, ldblaDataArray, ldblaJacobian, ldblB, ldblaF, nRows):
    """ Conditional Logistic Regression stores results in a list
        Parameters:
          lintOffset (int)
          ldblaDataArray (list)
          ldblaJacobian (list)
          ldblB (list)
          ldblF (list)
          nRows (int)
        Returns: float
    """
    conditional = 0.0
    x = []
    lDblAParamSum = []
    t = []
    c = []
    for i in range(0, len(ldblB)):
      x.append(0.0)
      lDblAParamSum.append(0.0)
      t.append(0.0)
      c.append([])
      for j in range(0, len(ldblB)):
        c[i].append(0.0)
    IthLikelihood = 0.0
    LogLikelihood = 0.0
    likelihood = 1.0
    ldblweight = 1.0
    lIntRow = 0
    lLevels = 0
    lLeveldata = 0.0
    lDblT0 = 0.0
    cases = 0
    count = 0

    for s in range(0, self.get_matchGroupValues()):
      lIntRow += lLevels
      lLevels = 1
      lLeveldata = float(ldblaDataArray[lIntRow][1])
      if s + 1 == self.get_matchGroupValues():
        lLevels = nRows - lIntRow + 0
      else:
        while lLeveldata == float(ldblaDataArray[lIntRow + lLevels - 0][1]):
          lLevels += 1
      lDblT0 = 0.0
      cases = 0
      count = 0
      for i in range(lIntRow - 0, lLevels + lIntRow - 0):
        if lintOffset == 3:
          ldblweight = float(ldblaDataArray[i][lintOffset])
        count += ldblweight
        for j in range(0, len(ldblB)):
          x[j] = float(ldblaDataArray[i][j + lIntOffset])
        if float(ldblaDataArray[i][0]) > 0:
          for j in range(0, len(ldblB)):
            lDblAParamSum[j] += x[j]
          cases += ldblweight
        IthLikelihood = 0.0
        for j in range(0, len(ldblB)):
          IthLikelihood += x[j] * float(ldblB[j])
        IthLikelihood = math.exp(IthLikelihood)
        IthLikelihood *= ldblweight
        lDblT0 += IthLikelihood
        for k in range(0, len(ldblB)):
          t[k] += IthLikelihood * x[k]
        for k in range(0, len(ldblB)):
          for j in range(0, len(ldblB)):
            c[j][k] += x[j] * x[k] * IthLikelihood

      IthLikelihood = 0.0
      for i in range(0, len(ldblB)):
        IthLikelihood += lDblAParamSum[i] * float(ldblB[i])
      contrast = Yes
      if cases == count or cases == 0:
        contrast = No
      if contrast == Yes:
        conditional += IthLikelihood - math.log(lDblT0)
        for i in range(0, len(ldblB)):
          ldblaFi = float(ldblaF[i])
          ldblaFi += lDblAParamSum[i] - t[i] / lDblT0
          ldblaF[i] = ldblaFi
        for i in range(0, len(ldblB)):
          for k in range(0, len(ldblB)):
            ldblaJacobianik = float(ldblaJacobian[i][k])
            ldblaJacobianik = ldblaJacobianik + c[i][k] / lDblT0 - t[i] * t[k] / (lDblT0 * lDblT0)
            ldblaJacobian[i][k] = ldblaJacobianik
      for i in range(0, len(ldblB)):
        lDblAParamSum[i] = 0.0
        t[i] = 0.0
        for k in range(0, len(ldblB)):
          c[i][k] = 0.0

    return conditional

  def UnConditional(self, lintOffset, ldblaDataArray, ldblaJacobian, ldblB, ldblaF, nRows):
    """ UnConditional Logistic Regression stores results in a list
        Parameters:
          lintOffset (int)
          ldblaDataArray (list)
          ldblaJacobian (list)
          ldblB (list)
          ldblF (list)
          nRows (int)
        Returns: float
    """
    unconditional = 0.0
    x = []
    for i in range(0, len(ldblB)):
      x.append(0.0)
    ldblIthLikelihood = 0.0
    ldblIthContribution = 0.0
    ldblweight = 1.0
    for i in range(0, nRows):
      for j in range(0, len(ldblB)):
        x[j] = float(ldblaDataArray[i][j + lintOffset])
      if lintOffset == 2:
        ldblweight = float(ldblaDataArray[i][1])
      ldblIthLikelihood = 0.0
      for j in range(0, len(ldblB)):
        ldblIthLikelihood += x[j] * float(ldblB[j])
      ldblIthLikelihood = 1 / (1 + math.exp(-ldblIthLikelihood))
      if float(ldblaDataArray[i][0]) == 0:
        ldblIthContribution = 1.0 - ldblIthLikelihood
      else:
        ldblIthContribution = ldblIthLikelihood
      for k in range(0, len(ldblB)):
        oldldblaF = 0.0
        if len(ldblaF) > k:
          oldldblaF = float(ldblaF[k])
        if float(ldblaDataArray[i][0]) > 0.0:
          newldblaF = oldldblaF + (1 - ldblIthLikelihood) * x[k] * ldblweight
          ldblaF[k] = ldblaF
        else:
          newldblaF = oldldblaF + (0 - ldblIthLikelihood) * x[k] * ldblweight
          ldblaF[k] = ldblaF
        for j in range(0, len(ldblB)):
          oldldblaJacobianjk = 0.0
          if len(ldblaJacobian) > j:
            oldldblaJacobianj = ldblaJacobian[j]
            if len(oldldblaJacobianj) > k:
              oldldblaJacobianjk = float(oldldblaJacobianj[k])
            else:
              ldblaJacobian.append([None] * len(ldblB))
            newldblaJacobianjk = oldldblaJacobianjk + (x[k] * x[j] * (1 - ldblIthLikelihood) * ldblIthLikelihood) * ldblweight
            ldblaJacobian[j][k] = newldblaJacobianjk
      unconditional = unconditional + log(ldblIthContribution) * ldblweight
    return unconditional

  def CalcLikelihood(self, lintOffset, ldblA, ldblB, ldblaJacobian, ldblaF, nRows, likelihood, strError, booStartAtZero):
    """ Computes likelihood and stores the result as a float in a list
        Parameters:
          lintOffset (int)
          ldblA (list)
          ldblB (list)
          ldblaJacobian (list)
          ldblF (list)
          nRows (int)
          likelihood (list of float): A float in a list so it will be mutable
          strError (list of str): A string in a list so it will be mutable
          booStartAtZero (bool)
        Returns: none
    """
    i = 0
    k = False
    ncases = 0.0
    nrecs = 0.0
    if len(strError) > 0:
      strError[0] = ""
    else:
      strError.append("")

    if self.get_mboolFirst() == True and self.get_mboolIntercept() == True:
      self.set_mboolFirst(False)
      ncases = 0.0
      for i in range(0, len(ldblA)):
        if int(ldblA[i][0]) == 1:
          ncases += 1.0
        nrecs += 1.0
      if ncases > 0.0 and nrecs - ncases > 0.0:
        if booStartAtZero == True:
          ldblB[len(ldblB) - 1] = 0.0
        else:
          ldblB[len(ldblB) - 1] = math.log(ncases / (nrecs - ncases))
      elif ncases == 0.0:
        strError[0] = "Dependent variable contains no cases."
        return
      elif nrecs - ncases == 0.0:
        strError[0] = "Dependent variable contains no controls."
        return
    if self.get_mstrMatchVar() is not None and len(self.get_mstrMatchVar()) > 0: 
      for i in range(0, len(ldblB)):
        ldblaF.append(0.0)
        arrayfori = []
        for j in range(0, len(ldblB)):
          arrayfori.append(0.0)
        ldblaJacobian.append(arrayfori)
      likelihood[0] = self.Conditional(lintOffset, ldblA, ldblaJacobian, ldblB, ldblaF, nRows)
    else:
      likelihood[0] = self.UnConditional(lintOffset, ldblA, ldblaJacobian, ldblB, ldblaF, nRows)

class LogisticRegressionResults:
  """ Class for Logistic Regression results
      which includs a list of LogisticRow objects
  """
  def __init__(self):
    self._Rows = []
    self._ErrorMessage = ''
  def get_Rows(self):
    return self._Rows
  def set_Rows(self, v):
    self._Rows = v
  def get_ErrorMessage(self):
    return self._ErrorMessage
  def set_ErrorMessage(self, v):
    self._ErrorMessage = v
