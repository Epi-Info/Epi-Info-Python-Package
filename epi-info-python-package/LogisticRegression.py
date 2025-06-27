from scipy.stats import t as tdist
import math
import time
import itertools
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
    self.ColumnsAndValues = {}

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
    self.mstrGroupVar = None
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
        self.mstrGroupVar = key
      if str(inputVariableList[key]).lower() == "weightvar":
        self.mstrWeightVar = key
      if str(inputVariableList[key]).lower() == "dependvar":
        self.mstrDependVar = key
      if str(inputVariableList[key]).lower() == "groupvar":
        self.mstrGroupVar = key
      if str(key).lower() == "intercept":
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
      if key == "dummies":
        self.dummiesNSMA = []
        for dum in inputVariableList[key]:
          self.dummiesNSMA.append(dum)

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

  def makeDummies(self, currentTableMA, variablesNeedingDummies, valuesForDummies, referencesForDummies, independentVariables):
    """ Creates dummy variables for independent variables
        Parameters:
          currentTableMA (list of lists)
          variablesNeetingDummies (list)
          valuesForDummies (list)
          referencesForDummies (list)
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
            newIndependentVariables[indexOfVariable - 1 + k] = valuesi[k] + ':' + referencesForDummies[i]
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
    dummiesMap = {}
    variablesNeedingDummies = []
    valuesForDummies = []
    referencesForDummies = []
    rowOne = currentTableMA[0]
    startCol = 1
    if self.mstrMatchVar is not None and len(str(self.mstrMatchVar)) > 0:
      startCol = 2
    for j in range(startCol, len(rowOne) - 1):
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
      dummiesMap[independentVariables[j - startCol]] = [independentVariables[j - startCol]]
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
      elif isNumeric == False or independentVariables[j - startCol] in self.dummiesNSMA:
        variablesNeedingDummies.append(j)
        valuesForThisJ = []
        for lnsmai in currentTableMA:
          lnsmaij = str(lnsmai[j])
          if lnsmaij not in valuesForThisJ:
            valuesForThisJ.append(lnsmaij)
        if len(valuesForThisJ) > 1:
          valuesForThisJ.sort()
          valuesForDummies.append(valuesForThisJ[1:])
          dummiesMap[independentVariables[j - 1]] = valuesForThisJ[1:]
          referencesForDummies.append(valuesForThisJ[0])
    if len(variablesNeedingDummies) > 0:
      self.makeDummies(currentTableMA, variablesNeedingDummies, valuesForDummies, referencesForDummies, independentVariables)

    return dummiesMap

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
    # Ensure interaction terms are included
    allCovariates = []
    soloCovariates = []
    interactionCovariates = []
    interactionsList = []
    for iV in independentVariables:
      if '*' not in iV:
        if iV not in allCovariates:
          allCovariates.append(iV)
        if iV not in soloCovariates:
          soloCovariates.append(iV)
      else:
        interactionsList.append(iV.replace(" ", ""))
        iVList = interactionsList[-1].split("*")
        for ivli in iVList:
          if ivli not in allCovariates:
            allCovariates.append(ivli)
          if ivli not in interactionCovariates:
            interactionCovariates.append(ivli)
    for rowi in self.currentTable:
      row = {k : v for k, v in rowi.items() if k in [outcomeVariable] + allCovariates}
      if self.mstrGroupVar is not None and len(self.mstrGroupVar) > 0:
        row[self.mstrGroupVar] = rowi[self.mstrGroupVar]
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
    interactionTableMutable = []
    for rd in mutableCurrentTable:
      rl = []
      rlindep = []
      rl.append(rd[outcomeVariable])
      rlindep.append(rd[outcomeVariable])
      if self.mstrGroupVar is not None and len(self.mstrGroupVar) > 0:
        rl.append(rd[self.mstrGroupVar])
        rlindep.append(rd[self.mstrGroupVar])
      #for iV in independentVariables:
      for iV in allCovariates:
        #if '*' in iV:
        #  iVList = iV.replace(" ", "").split("*")
        #  iVProd = float(rd[iVList[0]])
        #  for ivli in iVList[1:]:
        #    iVProd *= float(rd[ivli])
        #  rl.append(iVProd)
        #else:
        #  rl.append(rd[iV])
        if iV in soloCovariates:
          rl.append(rd[iV])
        if iV in interactionCovariates:
          rlindep.append(rd[iV])
      rl.append(rd['RecStatus'])
      rlindep.append(rd['RecStatus'])
      currentTableMutable.append(rl)
      interactionTableMutable.append(rlindep)
    soloDummiesMap = self.checkIndependentVariables(currentTableMutable, soloCovariates)
    interactionDummiesMap = self.checkIndependentVariables(interactionTableMutable, interactionCovariates)
    interactionCovariatesWithoutRefs = []
    for iC in interactionCovariates:
      if ':' in iC:
        interactionCovariatesWithoutRefs.append(iC.split(':')[0])
      else:
        interactionCovariatesWithoutRefs.append(iC)
    tableWithInteractions = []
    expandedInteractionsList = []
    valueslists = []
    for iLi in interactionsList:
      valueslist = []
      for iLiItem in iLi.split('*'):
        valueslist.append(interactionDummiesMap[iLiItem])
      valueslists.append(valueslist)
    combos = []
    for vl in valueslists:
      combos.append(list(itertools.product(*vl)))
    for c in combos:
      eI = []
      for t in c:
        expandedInteractionsList.append('*'.join(list(t)))
    interactedTable = []
    addtoindexofcolumn1 = 1
    if self.mstrMatchVar is not None and len(str(self.mstrMatchVar)) > 0:
      addtoindexofcolumn1 = 2
    allOfTheInteractingTerms = [] # Added PYDATE 20220309
    allOfTheInteractingTermsDataTable = [] # Added PYDATE 20220309
    newInteractionCovariatesWithoutRefs = []
    for covariateFromAllCovariates in allCovariates:
      if covariateFromAllCovariates in interactionCovariatesWithoutRefs:
        newInteractionCovariatesWithoutRefs.append(covariateFromAllCovariates)
    for itmi in interactionTableMutable:
      itr = []
      aotitdtr = []
      for eili in expandedInteractionsList:
        datum = 1.0
        eilil = eili.split('*')
        for eilili in eilil:
          indexofcolumn = interactionCovariatesWithoutRefs.index(eilili) + addtoindexofcolumn1
          datum *= float(itmi[indexofcolumn])
          if len(allOfTheInteractingTermsDataTable) == 0:
            allOfTheInteractingTerms.append(eilili) # Added PYDATE 20220309
          aotitdtr.append(itmi[indexofcolumn]) # Added PYDATE 20220309
        itr.append(datum)
      interactedTable.append(itr + [1])
      allOfTheInteractingTermsDataTable.append(aotitdtr) # Added PYDATE 20220309
    interactedColumns = []
    for eili in expandedInteractionsList:
      eilil = eili.split('*')
      icol = ""
      for eilili in eilil:
        icol += interactionCovariates[interactionCovariatesWithoutRefs.index(eilili)] + "*"
      interactedColumns.append(icol[:-1])
    independentVariables.clear()
    for sC in soloCovariates:
      independentVariables.append(sC)
    for sC in interactedColumns:
      independentVariables.append(sC)
    self.currentTable = []
    ctmrownumber = 0
    for ctmr in currentTableMutable:
      if self.mboolIntercept:
        self.currentTable.append(ctmr[:-1] + interactedTable[ctmrownumber])
      else:
        self.currentTable.append(ctmr[:-1] + interactedTable[ctmrownumber][:-1])
      ctmrownumber += 1

    # Set the values in ColumnsAndValues
    colindex = 0
    addtocolindex = 1
    if self.mstrMatchVar is not None and len(str(self.mstrMatchVar)) > 0:
      addtocolindex = 2
    for col in independentVariables:
      cd = {'number' : colindex}
      if '*' not in col:
        if ':' in col:
          cd['ref'] = col.split(':')[1]
          cd['compare'] = col.split(':')[0]
        else:
          uniquevalues = []
          for datarow in self.currentTable:
            if len(uniquevalues) > 2:
              break
            if datarow[colindex + addtocolindex] not in uniquevalues:
              uniquevalues.append(datarow[colindex + addtocolindex])
          if len(uniquevalues) == 2:
            uniquevalues.sort()
            cd['ref'] = uniquevalues[0]
            cd['compare'] = uniquevalues[1]
          else:
            cd['ref'] = None
            cd['compare'] = None
      else:
        cd['ref'] = None
        cd['compare'] = None
      self.ColumnsAndValues[col] = cd
      colindex += 1
    for col in independentVariables:
      if '*' in col:
        colcols = col.split('*')
        for colcol in colcols:
          if colcol not in self.ColumnsAndValues:
            # self.InteractionTerms.append(colcol) # Added PYDATE 20220309
            cd = {'number' : -1}
            uniquevalues = []
            mctindex = 0 # Added PYDATE 20220310
            for datarow in mutableCurrentTable:
              # if len(self.InteractionTerms) == 1: # Added PYDATE 20220309
              #   self.TableForInteractionTerms.append([]) # Added PYDATE 20220309
              # self.TableForInteractionTerms[mctindex].append(datarow[colcol]) # Added PYDATE 20220309
              # mctindex += 1 # Added PYDATE 20220309
              if len(uniquevalues) > 2:
                break
              #   continue # PYDATE 20220309: changed break to continue
              # if datarow[colcol] not in uniquevalues: # Commented out PYDATE 20220310
              if allOfTheInteractingTermsDataTable[mctindex][newInteractionCovariatesWithoutRefs.index(colcol)] not in uniquevalues: # Added PYDATE 20220310
                # uniquevalues.append(datarow[colcol]) # Commented out PYDATE 20220310
                uniquevalues.append(allOfTheInteractingTermsDataTable[mctindex][newInteractionCovariatesWithoutRefs.index(colcol)]) # Added PYDATE 20220310
              mctindex += 1 # Added PYDATE 20220310
            if len(uniquevalues) == 2:
              uniquevalues.sort()
              cd['ref'] = uniquevalues[0]
              cd['compare'] = uniquevalues[1]
            else:
              cd['ref'] = None
              cd['compare'] = None
            self.ColumnsAndValues[colcol] = cd

    self.InteractionTerms = allOfTheInteractingTerms # Added PYDATE 20220309
    self.TableForInteractionTerms = allOfTheInteractingTermsDataTable # Added PYDATE 20220309

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

  def getInteractionIndexes(self, interactions, iaTerms, bLabels):
    """ Returns the index of a variable
        Parameters:
          interactions (int)
          iaTerms (int)
          bLabels (list)
        Returns: list
    """
    interactionIndexes = []
    if iaTerms == 2 and interactions == 1:
      i = len(bLabels) - 1
      while i >= 0:
        if '*' in bLabels[i]:
          if len(interactionIndexes) == 0:
            interactionIndexes.append(i)
          else:
            interactionIndexes.insert(0, i)
        i -= 1
    return interactionIndexes

  def getContinuousIndex(self, var, bLabels):
    """ Returns the index of a variable
        Parameters:
          var (str)
          bLabels (list)
        Returns: int
    """
    contIdx = 0
    i = len(bLabels) - 1
    while i >= 0:
      if bLabels[i] == var:
        contIdx = i
      i -= 1
    return contIdx

  def getColumnMean(self, columnNumber, DataArray):
    """ Computes the mean of a column
        Parameters:
          columnNumber (int)
          DataArray (list of lists)
        Returns: float
    """
    columnSum = 0.0
    for row in DataArray:
      columnSum += float(row[columnNumber])
    return columnSum / float(len(DataArray))

  def getRef(self, var, bLabels):
    """ Returns a variable's reference value
        Parameters:
          var (str)
          BLabels (list)
        Returns: float
    """
    refValue = None
    return refValue

  def DummyFirst(self, cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          cm (list of lists)
          bLabels (list)
          B (list)
          lastVar1 (str)
          lastVar2 (str)
          interactions (int)
          iaTerms (int)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    Z = self.zFromP(0.025)
    ref1 = self.ColumnsAndValues[lastVar1]['ref']
    column2 = 1
    i = 0
    for bLabel in bLabels:
      if bLabel == lastVar2:
        column2 += i
      i += 1
    ref2 = self.getColumnMean(column2, DataArray)
    otherValues1 = [self.ColumnsAndValues[lastVar1]['number'], self.ColumnsAndValues[lastVar1]['compare']]
    interactionIndexes = []
    for bLabel in bLabels:
      if '*' in bLabel and lastVar1 in bLabel and lastVar2 in bLabel:
        interactionIndexes.append(self.ColumnsAndValues[bLabel]['number'])
    for i in range(int(len(otherValues1) / 2)):
      est = B[int(otherValues1[2 * i])] + ref2 * B[interactionIndexes[i]]
      variance = cm[int(otherValues1[2 * i])][int(otherValues1[2 * i])] +\
                 ref2 ** 2.0 * cm[interactionIndexes[i]][interactionIndexes[i]] +\
                 2 * ref2 * cm[int(otherValues1[2 * i])][interactionIndexes[i]]
      lcl = est - Z * variance ** 0.5
      ucl = est + Z * variance ** 0.5
      iorOut.append([lastVar1, str(otherValues1[2 * i + 1]) + ' vs ' + str(ref1) + ' at ' + lastVar2 + '=' + str(ref2), math.exp(est), math.exp(lcl), math.exp(ucl)])
    return iorOut

  def DummyLast(self, cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          cm (list of lists)
          bLabels (list)
          B (list)
          lastVar1 (str)
          lastVar2 (str)
          interactions (int)
          iaTerms (int)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    Z = self.zFromP(0.025)
    ref1 = self.ColumnsAndValues[lastVar1]['ref']
    singleIndex = self.getContinuousIndex(lastVar2, bLabels)
    interactionIndexes = []
    for bLabel in bLabels:
      if '*' in bLabel and lastVar1 in bLabel and lastVar2 in bLabel:
        interactionIndexes.append(self.ColumnsAndValues[bLabel]['number'])
    otherValues1 = [self.ColumnsAndValues[lastVar1]['number'], self.ColumnsAndValues[lastVar1]['compare']]
    est = B[singleIndex]
    variance0 = cm[singleIndex][singleIndex]
    lcl0 = est - Z * variance0 ** 0.5
    ucl0 = est + Z * variance0 ** 0.5
    iorOut.append([lastVar2 + '*' + lastVar1, lastVar2 + ' at ' + lastVar1 + '=' + str(ref1), math.exp(est), math.exp(lcl0), math.exp(ucl0)])
    for i in range(int(len(otherValues1) / 2)):
      est = B[singleIndex] + B[interactionIndexes[i]]
      variance = cm[singleIndex][singleIndex] + cm[interactionIndexes[i]][interactionIndexes[i]] + 2 * cm[singleIndex][interactionIndexes[i]]
      lcl = est - Z * variance ** 0.5
      ucl = est + Z * variance ** 0.5
      iorOut.append([lastVar2 + '*' + lastVar1, lastVar2 + ' at ' + lastVar1 + '=' + str(otherValues1[2 * i + 1]), math.exp(est), math.exp(lcl), math.exp(ucl)])
    return iorOut

  def TwoDummyVariables(self, cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          cm (list of lists)
          bLabels (list)
          B (list)
          lastVar1 (str)
          lastVar2 (str)
          interactions (int)
          iaTerms (int)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    Z = self.zFromP(0.025)
    ref1 = self.ColumnsAndValues[lastVar1]['ref']
    ref2 = self.ColumnsAndValues[lastVar2]['ref']
    otherValues1 = [self.ColumnsAndValues[lastVar1]['number'], self.ColumnsAndValues[lastVar1]['compare']]
    otherValues2 = [self.ColumnsAndValues[lastVar2]['number'], self.ColumnsAndValues[lastVar2]['compare']]
    est = B[int(otherValues1[0])]
    lcl = est - Z * cm[int(otherValues1[0])][int(otherValues1[0])] ** 0.5
    ucl = est + Z * cm[int(otherValues1[0])][int(otherValues1[0])] ** 0.5
    iorOut.append([lastVar1, str(otherValues1[1]) + ' vs ' + str(ref1) + ' at ' + lastVar2 + '=' + str(ref2), math.exp(est), math.exp(lcl), math.exp(ucl)])
    # REVISIT: section for more than one dummy variable for an initial variable
    ####
    interactionIndexes = []
    multiple = int(len(otherValues2) / 2)
    for bLabel in bLabels:
      if '*' in bLabel and lastVar1 in bLabel and lastVar2 in bLabel:
        interactionIndexes.append(self.ColumnsAndValues[bLabel]['number'])
    for k in range(int(len(otherValues2) / 2)):
      for i in range(int(len(otherValues1) / 2)):
        est = (B[int(otherValues1[2 * i])] + B[int(interactionIndexes[multiple * i + k])])
        variance = cm[int(otherValues1[2 * i])][int(otherValues1[2 * i])] + \
                   cm[interactionIndexes[multiple * i + k]][interactionIndexes[multiple * i + k]] + \
                   2 * cm[int(otherValues1[2 * i])][interactionIndexes[multiple * i + k]]
        lcl = est - Z * variance ** 0.5
        ucl = est + Z * variance ** 0.5
        iorOut.append([lastVar1, str(otherValues1[2 * i + 1]) + ' vs ' + str(ref1) + ' at ' + lastVar2 + '=' + str(otherValues2[2 * k + 1]), math.exp(est), math.exp(lcl), math.exp(ucl)])
    return iorOut

  def NoDummyVariables(self, cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          cm (list of lists)
          bLabels (list)
          B (list)
          lastVar1 (str)
          lastVar2 (str)
          interactions (int)
          iaTerms (int)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    Z = self.zFromP(0.025)
    column2 = 1
    i = 0
    for bLabel in bLabels:
      if bLabel is not None and bLabel == lastVar2:
        column2 += i
      i += 1
    ref2 = self.getColumnMean(column2, DataArray)
    singleIndex = self.getContinuousIndex(lastVar1, bLabels)
    interactionIndexes = []
    for bLabel in bLabels:
      if '*' in bLabel and lastVar1 in bLabel and lastVar2 in bLabel:
        interactionIndexes.append(self.ColumnsAndValues[bLabel]['number'])
    est = B[singleIndex] + ref2 * B[interactionIndexes[0]]
    variance = cm[singleIndex][singleIndex] + ref2 ** 2.0 * cm[interactionIndexes[0]][interactionIndexes[0]] + 2.0 * ref2 * cm[singleIndex][interactionIndexes[0]]
    lcl = est - Z * variance ** 0.5
    ucl = est + Z * variance ** 0.5
    iorOut.append([lastVar1 + '*' + lastVar2, lastVar1 + ' at ' + lastVar2 + ' = ' + str(ref2), math.exp(est), math.exp(lcl), math.exp(ucl)])
    return iorOut

  def IORNotMainEffects(self, lastVar1, lastVar2, cm, bLabels, B, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          lastVar1 (str): the first interaction term
          lastVar2 (str): the second interaction term
          cm (list of lists)
          bLabels (list)
          B (list)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    Z = self.zFromP(0.025)
    oneIsDummy = self.ColumnsAndValues[lastVar1]['ref'] is not None
    twoIsDummy = self.ColumnsAndValues[lastVar2]['ref'] is not None
    if lastVar1 not in bLabels and lastVar2 not in bLabels:
      if oneIsDummy and twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        iOR = math.exp(beta * ref1)
        iLCL = math.exp((beta - Z * iSE) * ref1)
        iUCL = math.exp((beta + Z * iSE) * ref1)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref1), \
                       iOR, \
                       iLCL, \
                       iUCL])
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif oneIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar2), self.TableForInteractionTerms)
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        iOR = math.exp(beta * ref1)
        iLCL = math.exp((beta - Z * iSE) * ref1)
        iUCL = math.exp((beta + Z * iSE) * ref1)
        iorOut.append([lastVar1, 'At ' + lastVar2 + ' = ' + str(ref1), iOR, iLCL, iUCL])
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, 'At ' + lastVar2 + ' = ' + str(ref2), iOR, iLCL, iUCL])
      else:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar2), self.TableForInteractionTerms)
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, 'At ' + lastVar2 + ' = ' + str(ref2), iOR, iLCL, iUCL])
    elif lastVar1 not in bLabels:
      if oneIsDummy and twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        iOR = math.exp(beta * ref1)
        iLCL = math.exp((beta - Z * iSE) * ref1)
        iUCL = math.exp((beta + Z * iSE) * ref1)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref1), \
                       iOR, \
                       iLCL, \
                       iUCL])
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif oneIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar1), self.TableForInteractionTerms)
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        iOR = math.exp(beta * ref1)
        iLCL = math.exp((beta - Z * iSE) * ref1)
        iUCL = math.exp((beta + Z * iSE) * ref1)
        iorOut.append([lastVar1, 'At ' + lastVar2 + ' = ' + str(ref1), iOR, iLCL, iUCL])
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, 'At ' + lastVar2 + ' = ' + str(ref2), iOR, iLCL, iUCL])
      else:
        iNumber2 = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber2]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber2][iNumber2] ** 0.5
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar1), self.TableForInteractionTerms)
        iOR = math.exp(beta * ref2)
        iLCL = math.exp((beta - Z * iSE) * ref2)
        iUCL = math.exp((beta + Z * iSE) * ref2)
        iorOut.append([lastVar1, \
                       'At ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
    elif lastVar2 not in bLabels:
      if oneIsDummy and twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1]['number']
        iNumber2 = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        beta2 = beta + B[iNumber2]
        iSE = self.mMatrixLikelihood.get_mdblaInv()[iNumber][iNumber] ** 0.5
        variance = cm[iNumber][iNumber] + cm[iNumber2][iNumber2] + 2 * cm[iNumber][iNumber2]
        iSEt = variance ** 0.5
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        iOR = math.exp(beta)
        iLCL = math.exp((beta - Z * iSE))
        iUCL = math.exp((beta + Z * iSE))
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref1), \
                       iOR, \
                       iLCL, \
                       iUCL])
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        iOR = math.exp(beta2 * ref2)
        iLCL = math.exp((beta2 - Z * iSEt) * ref2)
        iUCL = math.exp((beta2 + Z * iSEt) * ref2)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif oneIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1]['number']
        iNumber2 = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        beta2 = B[iNumber2]
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar2), self.TableForInteractionTerms)
        beta3 = beta + B[iNumber2] * ref2
        variance = cm[iNumber][iNumber] +\
                   ref2 ** 2.0 * cm[iNumber2][iNumber2] +\
                   2 * ref2 * cm[iNumber][iNumber2]
        iSEt = variance ** 0.5
        iOR = math.exp(beta3)
        iLCL = math.exp(beta3 - Z * iSEt)
        iUCL = math.exp(beta3 + Z * iSEt)
        iorOut.append([lastVar1, \
                       str(self.ColumnsAndValues[lastVar1]['compare']) + ' vs ' + str(self.ColumnsAndValues[lastVar1]['ref']) + ' at ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      elif twoIsDummy:
        iNumber = self.ColumnsAndValues[lastVar1]['number']
        iNumber2 = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        beta2 = B[iNumber2]
        ref1 = self.ColumnsAndValues[lastVar2]['ref']
        ref2 = self.ColumnsAndValues[lastVar2]['compare']
        beta3 = beta + B[iNumber2] * ref2
        variance = cm[iNumber][iNumber]
        iSEt = variance ** 0.5
        iOR = math.exp(beta)
        iLCL = math.exp(beta - Z * iSEt)
        iUCL = math.exp(beta + Z * iSEt)
        iorOut.append([lastVar1, \
                       'At ' + lastVar2 + ' = ' + str(ref1), \
                       iOR, \
                       iLCL, \
                       iUCL])
        variance = cm[iNumber][iNumber] +\
                   ref2 ** 2.0 * cm[iNumber2][iNumber2] +\
                   2 * ref2 * cm[iNumber][iNumber2]
        iSEt = variance ** 0.5
        iOR = math.exp(beta3)
        iLCL = math.exp(beta3 - Z * iSEt)
        iUCL = math.exp(beta3 + Z * iSEt)
        iorOut.append([lastVar1, \
                       'At ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
      else:
        iNumber = self.ColumnsAndValues[lastVar1]['number']
        iNumber2 = self.ColumnsAndValues[lastVar1 + '*' + lastVar2]['number']
        beta = B[iNumber]
        beta2 = B[iNumber2]
        ref2 = self.getColumnMean(self.InteractionTerms.index(lastVar2), self.TableForInteractionTerms)
        beta3 = beta + B[iNumber2] * ref2
        variance = cm[iNumber][iNumber] +\
                   ref2 ** 2.0 * cm[iNumber2][iNumber2] +\
                   2 * ref2 * cm[iNumber][iNumber2]
        iSEt = variance ** 0.5
        iOR = math.exp(beta3)
        iLCL = math.exp(beta3 - Z * iSEt)
        iUCL = math.exp(beta3 + Z * iSEt)
        iorOut.append([lastVar1, \
                       'At ' + lastVar2 + ' = ' + str(ref2), \
                       iOR, \
                       iLCL, \
                       iUCL])
    return iorOut

  def IOR(self, cm, bLabels, B, DataArray):
    """ Computes odds ratios and CIs for interaction
        terms holding one value fixed.
        Parameters:
          cm (list of lists)
          bLabels (list)
          B (list)
          DataArray (list of lists)
        Returns: list
    """
    iorOut = []
    noInteractions = True
    for bLabel in bLabels:
      if noInteractions and bLabel is not None and '*' in bLabel:
        noInteractions = False
    if noInteractions:
      return iorOut

    iaTerms = 1
    for bLabel in bLabels:
      if bLabel is not None and '*' in bLabel:
        iat = bLabel.count('*')
        if iat > iaTerms:
          iaTerms = iat
    if iaTerms > 1:
      return iorOut
    iaTerms += 1

    lastVar1 = ""
    lastVar2 = ""
    interactions = 0
    for bLabel in bLabels:
      if bLabel is not None and '*' in bLabel:
        lastVar1 = bLabel.split('*')[0]
        lastVar2 = bLabel.split('*')[1]
        interactions += 1
    if lastVar1 not in bLabels or lastVar2 not in bLabels:
      return self.IORNotMainEffects(lastVar1, lastVar2, cm, bLabels, B, DataArray)
    oneIsDummy = self.ColumnsAndValues[lastVar1]['ref'] is not None
    twoIsDummy = self.ColumnsAndValues[lastVar2]['ref'] is not None
    if oneIsDummy and twoIsDummy:
      iorOut += self.TwoDummyVariables(cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray)
    elif oneIsDummy:
      iorOut += self.DummyFirst(cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray)
    elif twoIsDummy:
      iorOut += self.DummyLast(cm, bLabels, B, lastVar2, lastVar1, interactions, iaTerms, DataArray)
    else:
      iorOut.append(self.NoDummyVariables(cm, bLabels, B, lastVar1, lastVar2, interactions, iaTerms, DataArray))

    return iorOut

  def depVarIsBinary(self, dataTable):
    """ Checks that the dependent variable is binary
        Parameters:
          dataTable (list(dict)): The analysis dataset
        Returns:
          bool
    """
    depvalues = set([])
    for row in dataTable:
      depvalues.add(row[0])
      if len(depvalues) > 2:
        return False
    if len(depvalues) == 2:
      return True
    else:
      return False

  def doRegression(self, inputVariableList, dataTable):
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
    if self.depVarIsBinary(self.currentTable) == False:
      print('Dependent variable must have exactly two values.')
      return self.logisticResults
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

    if len(self.mMatrixLikelihood.get_lstrError()) > 0 and self.mMatrixLikelihood.get_lstrError()[0] == "Matrix Tolerance Exceeded":
      return self.logisticResults

    for ev in inputVariableList['exposureVariables']:
      self.logisticResults.Variables.append(ev)
    self.logisticResults.Variables.append('CONSTANT')
    self.logisticResults.InteractionOR = self.IOR(self.mMatrixLikelihood.get_mdblaInv(), self.logisticResults.Variables, self.mMatrixLikelihood.get_mdblaB(), self.currentTable)
    if self.mboolIntercept == False or (self.mstrGroupVar is not None and len(self.mstrGroupVar) > 0):
      del self.logisticResults.Variables[-1]
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
    if self.logisticResults.Variables[-1] == 'CONSTANT':
      del self.logisticResults.OR[-1]
      del self.logisticResults.ORLCL[-1]
      del self.logisticResults.ORUCL[-1]

    self.logisticResults.Score = self.mMatrixLikelihood.get_mdblScore()
    self.logisticResults.ScoreDF = len(self.mMatrixLikelihood.get_mdblaB())
    if self.mboolIntercept:
      self.logisticResults.ScoreDF -= 1
    self.logisticResults.ScoreP = self.PValFromChiSq(self.logisticResults.Score, self.logisticResults.ScoreDF)
    self.logisticResults.LikelihoodRatio = 2.0 * (self.mMatrixLikelihood.get_mdbllllast()[0] - self.mMatrixLikelihood.get_mdblllfst()[0])
    self.logisticResults.LikelihoodRatioDF = self.logisticResults.ScoreDF
    self.logisticResults.LikelihoodRatioP = self.PValFromChiSq(self.logisticResults.LikelihoodRatio, self.logisticResults.LikelihoodRatioDF)
    self.logisticResults.MinusTwoLogLikelihood = -2.0 * self.mMatrixLikelihood.get_mdbllllast()[0]
    self.logisticResults.Iterations = self.mMatrixLikelihood.get_mintIterations()
    self.logisticResults.CasesIncluded = self.NumRows

    return self.logisticResults
