import math
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
    self._lstrError = [""]
    self._startValues = []
    self._parameterHistory = []

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

  def get_StartValues(self):
    return self._startValues
  def set_StartValues(self, v):
    self._startValues = v

  def get_ParameterHistory(self):
    return self._parameterHistory
  def set_ParameterHistory(self, v):
    self._parameterHistory = v

  def setListItem(self, thelist, theitem, theindex):
    """ Sets and item at a specific list index; grows list if necessary
        Parameters:
          thelist (list)
          theitem (anything)
          theindex (int)
    """
    if len(thelist) > theindex:
      thelist[theindex] = theitem
      return
    for i in range(1 + theindex - len(thelist)):
      thelist.append(None)
    thelist[theindex] = theitem

  def setMatrixItem(self, thematrix, theitem, theiindex, thejindex):
    """ Sets and item at a specific list of lists index; grows list(s) if necessary
        Parameters:
          thematrix (list of lists)
          theitem (anything)
          theiindex (int)
          thejindex (int)
    """
    if len(thematrix) > theiindex:
      self.setListItem(thematrix[theiindex], theitem, thejindex)
      return
    self.setListItem(thematrix, [], theiindex)
    self.setListItem(thematrix[theiindex], theitem, thejindex)

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
      B[i + 1] = sum
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
          a (int or float)
        Returns: Absolute value of a
    """
    if a < 0:
      return -a
    return a

  def trans(self, aRows, aCols, a, b):
    """ Transforms a matrix
        Parameters:
          aRows (int)
          aCols (int)
          a (list)
          b (list)
        Returns: none
    """
    for i in range(aRows):
      for j in range(aCols):
        b[j][i] = a[i][j]

  def mul(self, rowA, colA, rowB, colB, a, b, c):
    """ Populates one matrix from two others
        Parameters:
          aRows (int)
          aCols (int)
          a (list)
          b (list)
        Returns: none
    """
    if colA != rowB:
      return
    for i in range(rowA):
      for k in range(colB):
        c[i][k] = 0
        for j in range(colA):
          c[i][k] += a[i][j] * b[j][k]

  def ludcmpforlinear(self, n, d, a, indx):
    """ Performs calculations and stores the results in a list.
        Parameters:
          n (int)
          d (float)
          a (list of lists)
          indx (list)
        Returns: none
    """
    TINY = 1.0e-20
    imax = 0
    dum = 0.0
    big = 0.0
    sum = 0.0
    vv = []

    d = 1.0
    for i in range(0, n):
      big = 0.0
      for j in range(0, n):
        absaij = self.fabs(float(a[i][j]))
        if absaij > big:
          big = absaij
      if big == 0.0:
        return
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
        d = -1 * d
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
    return

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
        invA[i][j] = float(colShifted[i + 1])

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
          x[j] = float(ldblaDataArray[i][j + lintOffset])
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
      contrast = True
      if cases == count or cases == 0:
        contrast = False
      if contrast == True:
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
          self.setListItem(ldblaF, newldblaF, k)
        else:
          newldblaF = oldldblaF + (0 - ldblIthLikelihood) * x[k] * ldblweight
          self.setListItem(ldblaF, newldblaF, k)
        for j in range(0, len(ldblB)):
          oldldblaJacobianjk = 0.0
          if len(ldblaJacobian) > j:
            oldldblaJacobianj = ldblaJacobian[j]
            if len(oldldblaJacobianj) > k:
              oldldblaJacobianjk = 0.0 if oldldblaJacobianj[k] is None else float(oldldblaJacobianj[k])
          else:
            ldblaJacobian.append([None] * len(ldblB))
          newldblaJacobianjk = oldldblaJacobianjk + (x[k] * x[j] * (1 - ldblIthLikelihood) * ldblIthLikelihood) * ldblweight
          ldblaJacobian[j][k] = newldblaJacobianjk
      unconditional = unconditional + math.log(ldblIthContribution) * ldblweight
    return unconditional

  def UnConditionalLB(self, lintOffset, ldblaDataArray, ldblaJacobian, ldblB, ldblaF, nRows):
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
    ldblweight = 1.0
    self._parameterHistory.append([])
    for beta in ldblB:
      self._parameterHistory[len(self.get_ParameterHistory()) - 1].append(beta)
    for i in range(0, nRows):
      for j in range(0, len(ldblB)):
        x[j] = float(ldblaDataArray[i][j + lintOffset])
      if lintOffset == 2:
        ldblweight = float(ldblaDataArray[i][1])
      ldblIthLikelihood = 0.0
      for j in range(0, len(ldblB)):
        ldblIthLikelihood += x[j] * float(ldblB[j])
      ldblIthLikelihood = math.exp(ldblIthLikelihood)
      if float(ldblaDataArray[i][0]) == 0:
        ldblIthContribution = 1.0 - ldblIthLikelihood
      else:
        ldblIthContribution = ldblIthLikelihood
      for k in range(0, len(ldblB)):
        oldldblaF = 0.0
        if len(ldblaF) > k:
          oldldblaF = float(ldblaF[k])
        ysubi = 0.0
        if float(ldblaDataArray[i][0]) > 0.0:
          newldblaF = oldldblaF + (1 - 0) * x[k] * ldblweight
          self.setListItem(ldblaF, newldblaF, k)
          ysubi = 1.0
        else:
          newldblaF = oldldblaF + ((0 - ldblIthLikelihood) / (1 - ldblIthLikelihood)) * x[k] * ldblweight
          self.setListItem(ldblaF, newldblaF, k)
        for j in range(0, len(ldblB)):
          oldldblaJacobianjk = 0.0
          if len(ldblaJacobian) > j:
            oldldblaJacobianj = ldblaJacobian[j]
            if len(oldldblaJacobianj) > k:
              oldldblaJacobianjk = 0.0 if oldldblaJacobianj[k] is None else float(oldldblaJacobianj[k])
          else:
            ldblaJacobian.append([None] * len(ldblB))
          newldblaJacobianjk = oldldblaJacobianjk + ((1 - ysubi) * x[k] * x[j] * (1 / (1 - ldblIthLikelihood) ** 2.0) * ldblIthLikelihood) * ldblweight
          ldblaJacobian[j][k] = newldblaJacobianjk
      if ldblIthContribution <= 0.0:
        self.set_lstrError(['Could not maximize likelihood. Try different starting values.'])
        return 0.0
      unconditional = unconditional + math.log(ldblIthContribution) * ldblweight
    return unconditional

  def CalcLikelihood(self, lintOffset, ldblA, ldblB, ldblaJacobian, ldblaF, nRows, likelihood, strError, booStartAtZero):
    """ Computes likelihood and stores the result as a float in a list
        Parameters:
          lintOffset (int)
          ldblA (list)
          ldblB (list)
          ldblaJacobian (list)
          ldblF (list)
          nRows (list of int)
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

  def doLinear(self, currentTable):
    """ Computes the starting betas for log-binomial regression
        Uses weighted log linear regression.
        Parameters:
          ldblA (list)
        Returns: list
    """
    results = []
    lintWRows = 0
    lintweight = 1
    lintrowCount = 0
    ldblMagic = 0.0
    NumRows = len(currentTable)
    NumColumns = len(currentTable[0])
    for i in range(NumRows):
      lintWRows += 1
      if lintweight == 1:
        lintrowCount += currentTable[i][1]
      else:
        lintrowCount = lintWRows
    y = []
    x = []
    xx = []
    invxx = []
    xy = []
    tx = []
    B = []
    yhat = []
    resid = []
    sse = 0.0
    mse = 0.0
    rmse = 0
    indx = []
    d = 0.0
    fvalue = []
    covb = []
    probf = []
    stdb = []
    for i in range(NumColumns - 1 - lintweight):
      fvalue.append(None)
      probf.append(None)
      stdb.append(None)
      xy.append([None])
      B.append([None])
      indx.append(None)
      covbrow = []
      xxrow = []
      invxxrow = []
      for j in range(NumColumns - 1 - lintweight):
        covbrow.append(None)
        xxrow.append(None)
        invxxrow.append(None)
      covb.append(covbrow)
      xx.append(xxrow)
      invxx.append(invxxrow)
      txrow = []
      for j in range(NumRows):
        txrow.append(None)
      tx.append(txrow)
    for i in range(NumRows):
      y.append([None])
      yhat.append([None])
      resid.append(None)
      xrow = []
      for j in range(NumColumns - 1 - lintweight):
        xrow.append(None)
      x.append(xrow)
    coeff = []
    for i in range(4):
      coeffrow = []
      for j in range(NumColumns - 1 - lintweight):
        coeffrow.append(None)
      coeff.append(coeffrow)
    coeff[0][0] = 0.0
    meanY = 0.0
    ra2 = 0.0
    r2 = 0.0
    ssy = 0.0
    ftest = 0.0
    p = 0
    k = 0
    for i in range(NumRows):
      if lintweight == 1:
        y[k][0] = float(currentTable[i][0]) * float(currentTable[i][1]) ** 0.5
        for j in range(1 + lintweight, NumColumns):
          x[k][j - 1 - lintweight] = float(currentTable[i][j]) * float(currentTable[i][1]) ** 0.5
        k += 1
      else:
        # Unnecessary because only weighted analysis is called here
        k += 1
    mboolFirst = True
    mboolIntercept = True
    self.trans(NumRows, NumColumns - 1 - lintweight, x, tx)
    self.mul(NumColumns - 1 - lintweight, NumRows, NumRows, NumColumns - 1 - lintweight, tx, x, xx)
    self.mul(NumColumns - 1 - lintweight, NumRows, NumRows, 1, tx, y, xy)
    for e in range(NumColumns - 1 - lintweight):
      for c in range(NumColumns - 1 - lintweight):
        invxx[e][c] = xx[e][c]
    self.ludcmpforlinear(NumColumns - 1 - lintweight, d, invxx, indx)
    d = 1
    for i in range(NumColumns - 1 - lintweight):
      d *= invxx[i][i]
      if d == 0:
        print('Colinear Data; using log(p) as start value')
        return B
    if self.fabs(d) < 0.000001:
      print('Linear matrix tolerance exceeded; using log(p) as start value')
      return B
    self.inv(xx, invxx)
    self.mul(NumColumns - 1 - lintweight, NumColumns - 1 - lintweight, NumColumns - 1 - lintweight, 1, invxx, xy, B)
    return B

  def GetStartValues(self, ldblA):
    """ Computes the starting betas for log-binomial regression
        Uses weighted log linear regression.
        Parameters:
          ldblA (list)
        Returns: list
    """
    startvalues = []
    datamatrix = []
    for r in ldblA:
      dmr = [0]
      for dmrc in r:
        dmr.append(dmrc)
      y = 0.9 * dmr[1] + 0.1 * (1 - dmr[1])
      dmr[0] = math.log(y)
      dmr[1] = (y / (1 - y))
      datamatrix.append(dmr)
    regresults = self.doLinear(datamatrix)
    for rr in regresults:
      if type(rr) is list:
        startvalues.append(rr[0])
      elif type(rr) is float:
        startvalues.append(rr)
      else:
        startvalues.append(0.0)
    return startvalues

  def CalcLikelihoodLB(self, lintOffset, ldblA, ldblB, ldblaJacobian, ldblaF, nRows, likelihood, strError, booStartAtZero):
    """ Computes likelihood and stores the result as a float in a list
        Parameters:
          lintOffset (int)
          ldblA (list)
          ldblB (list)
          ldblaJacobian (list)
          ldblF (list)
          nRows (list of int)
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
          ldblB[len(ldblB) - 1] = math.log(ncases / nrecs)
          startvalues = self.GetStartValues(ldblA)
          if len(self.get_StartValues()) > 0:
            startvalues = self.get_StartValues()
            while len(startvalues) < len(ldblB):
              startvalues.insert(0, 0.0)
          if startvalues[len(startvalues) - 1] is not None:
            for stv in range(len(startvalues)):
              ldblB[stv] = startvalues[stv]
      elif ncases == 0.0:
        strError[0] = "Dependent variable contains no cases."
        return
      elif nrecs - ncases == 0.0:
        strError[0] = "Dependent variable contains no controls."
        return
    if False: # No Conditional regression with Log-Binomial
      for i in range(0, len(ldblB)):
        ldblaF.append(0.0)
        arrayfori = []
        for j in range(0, len(ldblB)):
          arrayfori.append(0.0)
        ldblaJacobian.append(arrayfori)
      likelihood[0] = self.Conditional(lintOffset, ldblA, ldblaJacobian, ldblB, ldblaF, nRows)
    else:
      likelihood[0] = self.UnConditionalLB(lintOffset, ldblA, ldblaJacobian, ldblB, ldblaF, nRows)

  def MaximizeLikelihood(self, nRows, nCols, dataArray, lintOffset, lintMatrixSize, llngIters, ldblToler, ldblConv, booStartAtZero):
    """ Maximizes Likelihood
        Parameters:
          nRows (int)
          nCols (int)
          dataArray (list of lists)
          lintOffset (int)
          lintMatrixSize (list)
          llngIters (int)
          ldblToler (float)
          ldblConv (float)
          booStartAtZero (bool)
        Returns: none
    """
    strCalcLikelihoodError = [""]
    self.set_lstrError(strCalcLikelihoodError)
    ldbllfst = [0.0]
    self.set_mdblScore(0.0)
    oldmdblScore = 0.0
    self.set_mboolConverge(True)
    self.set_mboolErrorStatus(False)
    self.set_mdblaJacobian([])
    oldmdblaJacobian = [[None] * lintMatrixSize for j in range(lintMatrixSize)]
    self.set_mdblaInv([])
    self.set_mdblaF([])
    oldmdblaF = [None] * lintMatrixSize
    self.set_mdblaB([0.0] * lintMatrixSize)

    self.CalcLikelihood(lintOffset, dataArray, self.get_mdblaB(), self.get_mdblaJacobian(), self.get_mdblaF(), nRows, ldbllfst, strCalcLikelihoodError, booStartAtZero)

    for i in range(len(self.get_mdblaB())):
      forInv = []
      for j in range(len(self.get_mdblaB())):
        forInv.append(0.0)
        self.setMatrixItem(oldmdblaJacobian, float(self.get_mdblaJacobian()[i][j]), i, j)
        oldmdblaJacobian[i][j] = float(self.get_mdblaJacobian()[i][j])
      self.get_mdblaInv().append(forInv)
      oldmdblaF[i] = float(self.get_mdblaF()[i])

    if len(strCalcLikelihoodError[0]) > 0:
      print(strCalcLikelihoodError[0])
      return

    self.set_mintIterations(1)
    ldbloldll = ldbllfst[0]
    ldbll = []
    ldbll.append(ldbllfst[0])
    if ldbllfst[0] > 0:
      self.set_mboolConverge(False)
      strCalcLikelihoodError[0] = "Positive Log-Likelihood, regression is diverging"
      print("Positive Log-Likelihood, regression is diverging")
      return

    self.inv(self.get_mdblaJacobian(), self.get_mdblaInv())

    oldmdblaB = []
    oldmdblaInv = []
    for pop in range(len(self.get_mdblaB())):
      oldmdblaB.append(self.get_mdblaB()[pop])
    for pop in range(len(self.get_mdblaInv())):
      oldmdblaInv.append(self.get_mdblaInv()[pop])
    ldblDet = 1.0
    for i in range(len(self.get_mdblaB())):
      ldblDet *= float(self.get_mdblaJacobian()[i][i])
    if self.fabs(ldblDet) < ldblToler:
      self.set_mboolConverge(False)
      strCalcLikelihoodError[0] = "Matrix Tolerance Exceeded"
      self.set_lstrError(strCalcLikelihoodError)
      print("Matrix Tolerance Exceeded")
      return

    ldblaScore = [0.0] * lintMatrixSize
    for i in range(len(self.get_mdblaB())):
      for k in range(len(self.get_mdblaB())):
        mdblainfik = float(self.get_mdblaInv()[i][k])
        mdblafk = float(self.get_mdblaF()[k])
        onetimestheother = mdblafk * mdblainfik
        ldblaScore[i] = ldblaScore[i] + onetimestheother
        self.get_mdblaJacobian()[i][k] = 0.0
      self.get_mdblaB()[i] = float(self.get_mdblaB()[i]) + ldblaScore[i]
      self.set_mdblScore(self.get_mdblScore() + ldblaScore[i] * float(self.get_mdblaF()[i]))

    ridge = 0.0

    self.set_mintIterations(2)
    while self.get_mintIterations() < llngIters:
      for i in range(len(self.get_mdblaF())):
        self.get_mdblaF()[i] = 0.0
      self.CalcLikelihood(lintOffset, dataArray, self.get_mdblaB(), self.get_mdblaJacobian(), self.get_mdblaF(), nRows, ldbll, strCalcLikelihoodError, booStartAtZero)
      doThisStuff = True
      if ldbloldll - ldbll[0] > ldblConv:
        if ridge > 0.0 and ridge < 1000.0:
          self.set_mintIterations(self.get_mintIterations() - 1)
          ridge *= 4.0
          for i in range(len(self.get_mdblaB())):
            for j in range(len(self.get_mdblaB())):
              iequalsj = 0.0
              if i == j:
                iequalsj = 1.0
              self.get_mdblaJacobian()[i][j] = oldmdblaJacobian[i][j] * (1 + iequalsj * ridge)
            self.get_mdblaB()[i] = float(oldmdblaB[i])
            self.get_mdblaF()[i] = float(oldmdblaF[i])
          doThisStuff = False
        if doThisStuff == True:
          self.set_mboolConverge(False)
          self.set_mdblllfst(ldbllfst)
          strCalcLikelihoodError[0] = "Regression not converging"
      elif ldbll[0] - ldbloldll < ldblConv:
        self.set_mdblaB([])
        for i in range(len(oldmdblaB)):
          self.get_mdblaB().append(oldmdblaB[i])
        self.set_mintIterations(self.get_mintIterations() - 1)
        self.set_mdblllfst(ldbllfst)
        self.set_mdbllllast(ldbll)
        return

      if doThisStuff == True:
        oldmdblaInv = []
        for i in range(len(self.get_mdblaInv())):
          oldmdblaInv.append(self.get_mdblaInv()[i])
        oldmdblaB = []
        for i in range(len(self.get_mdblaB())):
          oldmdblaB.append(self.get_mdblaB()[i])
        for i in range(len(self.get_mdblaB())):
          for j in range(len(self.get_mdblaB())):
            oldmdblaJacobian[i][j] = float(self.get_mdblaJacobian()[i][j])
          oldmdblaF[i] = float(self.get_mdblaF()[i])
        ridge = 0.0
        ldbloldll = ldbll[0]

      self.inv(self.get_mdblaJacobian(), self.get_mdblaInv())
      ldblDet = 1.0
      for i in range(len(self.get_mdblaB())):
        ldblDet *= self.get_mdblaJacobian()[i][i]
      if self.fabs(ldblDet) < ldblToler:
        self.set_mboolConverge(False)
        strCalcLikelihoodError[0] = "Matrix Tolerance Exceeded"
        self.set_lstrError(strCalcLikelihoodError)
        print("Matrix Tolerance Exceeded")
        return
      for i in range(len(self.get_mdblaB())):
        for k in range(len(self.get_mdblaB())):
          self.get_mdblaB()[i] = float(self.get_mdblaB()[i]) + float(self.get_mdblaF()[k]) * float(self.get_mdblaInv()[i][k])
          self.get_mdblaJacobian()[i][k] = 0.0

      self.set_mintIterations(self.get_mintIterations() + 1)

    self.set_mdblllfst(ldbllfst)
    self.set_mdbllllast(ldbll)

  def MaximizeLikelihoodLB(self, nRows, nCols, dataArray, lintOffset, lintMatrixSize, llngIters, ldblToler, ldblConv, booStartAtZero):
    """ Maximizes Likelihood
        Parameters:
          nRows (int)
          nCols (int)
          dataArray (list of lists)
          lintOffset (int)
          lintMatrixSize (list)
          llngIters (int)
          ldblToler (float)
          ldblConv (float)
          booStartAtZero (bool)
        Returns: none
    """
    strCalcLikelihoodError = [""]
    self.set_lstrError(strCalcLikelihoodError)
    ldbllfst = [0.0]
    self.set_mdblScore(0.0)
    oldmdblScore = 0.0
    self.set_mboolConverge(True)
    self.set_mboolErrorStatus(False)
    self.set_mdblaJacobian([])
    oldmdblaJacobian = [[None] * lintMatrixSize for j in range(lintMatrixSize)]
    self.set_mdblaInv([])
    self.set_mdblaF([])
    oldmdblaF = [None] * lintMatrixSize
    self.set_mdblaB([0.0] * lintMatrixSize)

    self.CalcLikelihoodLB(lintOffset, dataArray, self.get_mdblaB(), self.get_mdblaJacobian(), self.get_mdblaF(), nRows, ldbllfst, strCalcLikelihoodError, booStartAtZero)

    if self.get_lstrError()[0] == 'Could not maximize likelihood. Try different starting values.':
      print(self.get_lstrError()[0])
      return

    for i in range(len(self.get_mdblaB())):
      forInv = []
      for j in range(len(self.get_mdblaB())):
        forInv.append(0.0)
        self.setMatrixItem(oldmdblaJacobian, float(self.get_mdblaJacobian()[i][j]), i, j)
        oldmdblaJacobian[i][j] = float(self.get_mdblaJacobian()[i][j])
      self.get_mdblaInv().append(forInv)
      oldmdblaF[i] = float(self.get_mdblaF()[i])

    if len(strCalcLikelihoodError[0]) > 0:
      print(strCalcLikelihoodError[0])
      return

    self.set_mintIterations(1)
    ldbloldll = ldbllfst[0]
    ldbll = []
    ldbll.append(ldbllfst[0])
    if ldbllfst[0] > 0:
      self.set_mboolConverge(False)
      strCalcLikelihoodError[0] = "Positive Log-Likelihood, regression is diverging"
      print("Positive Log-Likelihood, regression is diverging")
      return

    self.inv(self.get_mdblaJacobian(), self.get_mdblaInv())

    oldmdblaB = []
    oldmdblaInv = []
    for pop in range(len(self.get_mdblaB())):
      oldmdblaB.append(self.get_mdblaB()[pop])
    for pop in range(len(self.get_mdblaInv())):
      oldmdblaInv.append(self.get_mdblaInv()[pop])
    ldblDet = 1.0
    for i in range(len(self.get_mdblaB())):
      ldblDet *= float(self.get_mdblaJacobian()[i][i])
    if self.fabs(ldblDet) < ldblToler:
      self.set_mboolConverge(False)
      strCalcLikelihoodError[0] = "Matrix Tolerance Exceeded"
      self.set_lstrError(strCalcLikelihoodError)
      print("Matrix Tolerance Exceeded")
      return

    ldblaScore = [0.0] * lintMatrixSize
    for i in range(len(self.get_mdblaB())):
      for k in range(len(self.get_mdblaB())):
        mdblainfik = float(self.get_mdblaInv()[i][k])
        mdblafk = float(self.get_mdblaF()[k])
        onetimestheother = mdblafk * mdblainfik
        ldblaScore[i] = ldblaScore[i] + onetimestheother
        self.get_mdblaJacobian()[i][k] = 0.0
      self.get_mdblaB()[i] = float(self.get_mdblaB()[i]) + ldblaScore[i]
      self.set_mdblScore(self.get_mdblScore() + ldblaScore[i] * float(self.get_mdblaF()[i]))

    ridge = 0.0

    self.set_mintIterations(2)
    while self.get_mintIterations() < llngIters:
      for i in range(len(self.get_mdblaF())):
        self.get_mdblaF()[i] = 0.0
      self.CalcLikelihoodLB(lintOffset, dataArray, self.get_mdblaB(), self.get_mdblaJacobian(), self.get_mdblaF(), nRows, ldbll, strCalcLikelihoodError, booStartAtZero)

      if self.get_lstrError()[0] == 'Could not maximize likelihood. Try different starting values. See the epiinfo package README for further instructions.':
        print(self.get_lstrError()[0])
        return

      doThisStuff = True
      if ldbloldll - ldbll[0] > ldblConv:
        if ridge > 0.0 and ridge < 1000.0:
          self.set_mintIterations(self.get_mintIterations() - 1)
          ridge *= 4.0
          for i in range(len(self.get_mdblaB())):
            for j in range(len(self.get_mdblaB())):
              iequalsj = 0.0
              if i == j:
                iequalsj = 1.0
              self.get_mdblaJacobian()[i][j] = oldmdblaJacobian[i][j] * (1 + iequalsj * ridge)
            self.get_mdblaB()[i] = float(oldmdblaB[i])
            self.get_mdblaF()[i] = float(oldmdblaF[i])
          doThisStuff = False
        if doThisStuff == True:
          self.set_mboolConverge(False)
          self.set_mdblllfst(ldbllfst)
          strCalcLikelihoodError[0] = "Regression not converging"
      elif ldbll[0] - ldbloldll < ldblConv:
        #self.set_mdblaB([])
        #for i in range(len(oldmdblaB)):
        #  self.get_mdblaB().append(oldmdblaB[i])
        #self.set_mintIterations(self.get_mintIterations() - 1)
        self.set_mdblllfst(ldbllfst)
        self.set_mdbllllast(ldbll)
        self.inv(self.get_mdblaJacobian(), self.get_mdblaInv())
        return

      if doThisStuff == True:
        oldmdblaInv = []
        for i in range(len(self.get_mdblaInv())):
          oldmdblaInv.append(self.get_mdblaInv()[i])
        oldmdblaB = []
        for i in range(len(self.get_mdblaB())):
          oldmdblaB.append(self.get_mdblaB()[i])
        for i in range(len(self.get_mdblaB())):
          for j in range(len(self.get_mdblaB())):
            oldmdblaJacobian[i][j] = float(self.get_mdblaJacobian()[i][j])
          oldmdblaF[i] = float(self.get_mdblaF()[i])
        ridge = 0.0
        ldbloldll = ldbll[0]

      self.inv(self.get_mdblaJacobian(), self.get_mdblaInv())
      ldblDet = 1.0
      for i in range(len(self.get_mdblaB())):
        ldblDet *= self.get_mdblaJacobian()[i][i]
      if self.fabs(ldblDet) < ldblToler:
        self.set_mboolConverge(False)
        strCalcLikelihoodError[0] = "Matrix Tolerance Exceeded"
        self.set_lstrError(strCalcLikelihoodError)
        print("Matrix Tolerance Exceeded")
        return
      for i in range(len(self.get_mdblaB())):
        for k in range(len(self.get_mdblaB())):
          self.get_mdblaB()[i] = float(self.get_mdblaB()[i]) + float(self.get_mdblaF()[k]) * float(self.get_mdblaInv()[i][k])
          self.get_mdblaJacobian()[i][k] = 0.0

      self.set_mintIterations(self.get_mintIterations() + 1)

    self.set_mdblllfst(ldbllfst)
    self.set_mdbllllast(ldbll)

class LogisticRegressionResults:
  """ Class for Logistic Regression results
      which includs a list of LogisticRow objects
  """
  def __init__(self):
    self._Rows = []
    self._ErrorMessage = ''
    self.Variables = []
    self.Beta = []
    self.SE = []
    self.OR = []
    self.ORLCL = []
    self.ORUCL = []
    self.RR = []
    self.RRLCL = []
    self.RRUCL = []
    self.Z = []
    self.PZ = []
    self.Score = None
    self.ScoreDF = None
    self.ScoreP = None
    self.LikelihoodRatio = None
    self.LikelihoodRatioDF = None
    self.LikelihoodRatioP = None
    self.LogLikelihood = None
    self.MinusTwoLogLikelihood = None
    self.Iterations = None
    self.CasesIncluded = None
    self.InteractionOR = []
    self.InteractionRR = []
    self.ParameterHistory = None

  def get_Rows(self):
    return self._Rows
  def set_Rows(self, v):
    self._Rows = v
  def get_ErrorMessage(self):
    return self._ErrorMessage
  def set_ErrorMessage(self, v):
    self._ErrorMessage = v

  def show(self):
    if len(self.Beta) == 0:
      print('Nothing to show')
      return
    widthofcoeff = 16
    showFitTests = True
    for varname in self.Variables:
      if len(varname) > widthofcoeff:
        widthofcoeff = len(varname)
    regType = 'LOGISTIC'
    testType = 'Minus Two Log-Likelihood'
    ratioType = 'Odds'
    ratios = self.OR
    ratiolls = self.ORLCL
    ratiouls = self.ORUCL
    likelihood = self.MinusTwoLogLikelihood
    interactions = self.InteractionOR
    if len(self.OR) == 0:
      showFitTests = False
      regType = 'LOG BINOMIAL'
      testType = 'Log-Likelihood'
      ratioType = 'Risk'
      ratios = self.RR
      ratiolls = self.RRLCL
      ratiouls = self.RRUCL
      likelihood = self.LogLikelihood
      interactions = self.InteractionRR
    print(regType, 'REGRESSION RESULTS')
    # Saving this next line to record this other way of accomplishing the formatting
    # print('{0:<{6}s}{1:^16s}{2:^16s}{3:^16s}{4:^16s}{5:^16s}'.format('Variable','Coefficient','Standard Error',ratioType+' Ratio','Lower','Upper',widthofcoeff))
    print('{:<{width}s}{:^16s}{:^16s}{:^16s}{:^8s}{:^8s}'.format('Variable','Coefficient','Standard Error',ratioType+' Ratio','Lower','Upper',width=widthofcoeff))
    for i in range(len(self.Variables)):
      if i < len(self.Variables) - 1:
        print('{:<{width}s}{:^16s}{:^16s}{:^16s}{:^8s}{:^8s}'.format(self.Variables[i],str(round(self.Beta[i], 4)),str(round(self.SE[i], 4)),str(round(ratios[i], 4)),str(round(ratiolls[i], 4)),str(round(ratiouls[i], 4)),width=widthofcoeff))
    print('{:<{width}s}{:^16s}{:^16s}\n'.format(self.Variables[i],str(round(self.Beta[i], 4)),str(round(self.SE[i], 4)),width=widthofcoeff))
    print('{:<{width}s}{:>8s}'.format('Number of Iterations',str(round(self.Iterations, 0)),width=26))
    print('{:<{width}s}{:>8s}'.format(testType,str(round(likelihood, 2)),width=26))
    print('{:<{width}s}{:>8s}'.format('Number of Observations',str(round(self.CasesIncluded, 0)),width=26))
    if showFitTests:
      print()
      print('{:<{width}s}{:>8s}{:>8s}{:>8s}'.format('Fit Test','Value','DF','P',width=18))
      print('{:<{width}s}{:>8s}{:>8s}{:>8s}'.format('Score',str(round(self.Score, 4)),str(round(self.ScoreDF, 4)),str(round(self.ScoreP, 4)),width=18))
      print('{:<{width}s}{:>8s}{:>8s}{:>8s}'.format('Likelihood Ratio',str(round(self.LikelihoodRatio, 4)),str(round(self.LikelihoodRatioDF, 4)),str(round(self.LikelihoodRatioP, 4)),width=18))

  def showHTML(self):
    if len(self.Beta) == 0:
      print('Nothing to show<br>')
      return
    widthofcoeff = 16
    showFitTests = True
    for varname in self.Variables:
      if len(varname) > widthofcoeff:
        widthofcoeff = len(varname)
    regType = 'LOGISTIC'
    testType = 'Minus Two Log-Likelihood'
    ratioType = 'Odds'
    ratios = self.OR
    ratiolls = self.ORLCL
    ratiouls = self.ORUCL
    likelihood = self.MinusTwoLogLikelihood
    interactions = self.InteractionOR
    if len(self.OR) == 0:
      showFitTests = False
      regType = 'LOG BINOMIAL'
      testType = 'Log-Likelihood'
      ratioType = 'Risk'
      ratios = self.RR
      ratiolls = self.RRLCL
      ratiouls = self.RRUCL
      likelihood = self.LogLikelihood
      interactions = self.InteractionRR
    print('<style type="text/css">\ntd.PyStats\n{\nborder-style:none;\nborder-width:0px;\npadding-left:5px;\npadding-right:5px;\npadding-bottom:0px;\n}\n</style>\n\n')
    print(regType, 'REGRESSION RESULTS', '<br>')
    # Saving this next line to record this other way of accomplishing the formatting
    # print('{0:<{6}s}{1:^16s}{2:^16s}{3:^16s}{4:^16s}{5:^16s}'.format('Variable','Coefficient','Standard Error',ratioType+' Ratio','Lower','Upper',widthofcoeff), '<br>')
    # print('{:<{width}s}{:^16s}{:^16s}{:^16s}{:^8s}{:^8s}'.format('Variable','Coefficient','Standard Error',ratioType+' Ratio','Lower','Upper',width=widthofcoeff), '<br>')
    print('<br>\n<table>\n<tr><td class="PyStats"><strong>Variable</strong></td><td class="PyStats"><strong>Coefficient</strong></td><td class="PyStats"><strong>Standard Error</strong></td><td class="PyStats"><strong>',ratioType+' Ratio</strong></td><td class="PyStats"><strong>Lower</strong></td><td class="PyStats"><strong>Upper</strong></td></tr>')
    for i in range(len(self.Variables)):
      if i < len(self.Variables) - 1:
        print('<tr><td class="PyStats"><strong>', self.Variables[i], '</strong></td><td class="PyStats" align="right">', '{:.4f}'.format(round(self.Beta[i], 4)), '</td><td class="PyStats" align="right">', '{:.4f}'.format(round(self.SE[i], 4)), '</td><td class="PyStats" align="right">', '{:.4f}'.format(round(ratios[i], 4)), '</td><td class="PyStats">', '{:.4f}'.format(round(ratiolls[i], 4)), '</td><td class="PyStats">', '{:.4f}'.format(round(ratiouls[i], 4)), '</td></tr>\n')
    print('<tr><td class="PyStats"><strong>', self.Variables[i], '</strong></td><td class="PyStats" align="right">', '{:.4f}'.format(round(self.Beta[i], 4)), '</td><td class="PyStats" align="right">', '{:.4f}'.format(round(self.SE[i], 4)), '</td><td class="PyStats">&nbsp;</td><td class="PyStats">&nbsp;</td><td class="PyStats">&nbsp;</td></tr>\n</table>\n<br>\n')
    print('<table>\n<tr><td class="PyStats"><strong>', 'Number of Iterations', '</strong></td><td class="PyStats" align="right">', str(round(self.Iterations, 0)), '</td></tr>\n')
    print('<tr><td class="PyStats"><strong>', testType, '</strong></td><td class="PyStats" align="right">', str(round(likelihood, 2)), '</td></tr>\n')
    print('<tr><td class="PyStats"><strong>', 'Number of Observations', '</strong></td><td class="PyStats" align="right">', str(round(self.CasesIncluded, 0)), '</td></tr>\n</table>\n')
    if showFitTests:
      print('<br>\n')
      print('<table>\n<tr><td class="PyStats"><strong>', 'Fit Test', '</strong></td><td class="PyStats" align="center"><strong>', 'Value', '</strong></td><td class="PyStats"><strong>', 'DF', '</strong></td><td class="PyStats" align="center"><strong>', 'P', '</strong></td></tr>\n')
      print('<tr><td class="PyStats"><strong>', 'Score', '</strong></td><td class="PyStats">', str(round(self.Score, 4)), '</td><td class="PyStats" align="center">', str(round(self.ScoreDF, 4)), '</td><td class="PyStats">', str(round(self.ScoreP, 4)), '</td></tr>\n')
      print('<tr><td class="PyStats"><strong>', 'Likelihood Ratio', '</strong></td><td class="PyStats">', str(round(self.LikelihoodRatio, 4)), '</td><td class="PyStats" align="center">', str(round(self.LikelihoodRatioDF, 4)), '</td><td class="PyStats">', str(round(self.LikelihoodRatioP, 4)), '</td></tr>\n</table>\n<br>\n')
