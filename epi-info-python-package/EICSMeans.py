from scipy.stats import t as tdist
import math
import time
from .randata import randata
from .CSUtilities import *

class ComplexSampleMeans:
  def __init__(self):
    self.strataVar = None
    self.mainVar = None
    self.crosstabVar = None
    self.domainVar = None
    self.psuVar = None
    self.weightVar = None
    self.columnNames = None
    self.domain1 = None
    self.domain2 = None

    self.validCases = 0

    self.tableName = None
    self.booleanLabels = None
    self.outputLevel = None
    self.percents = None
    self.booleanValues = None

    self.sortedTable = [{}]
    self.distinctTable = [{}]

    self.mis = None
    self.first = CSMeansTotal()
    self.last = CSMeansTotal()
    self.com = None

    self.row = None

    self.outcome = CSField()
    self.domain = None
    self.strata = CSField()
    self.psu = CSField()
    self.weight = CSField()
    self.crossTab = None

    self.varT = None
    self.csOutputBuffer = None
    self.cnOutputLevel = None
    self.cbIncludePercents = None
    self.cbStandalone = None

    self.isDeleted = None
    self.isVerified = None

    self.varianceMultiplier = None
    self.errorMessage = None
    self.numErrors = None

    self.meansResults = CSMeansResults()

    self.currentTable = [{}]
    
    self.confidenceLevel = None

  def CreateSettings(self, inputVariableList):
    """ Initializes objectes necessary for the analysis.
        Checks for existance of Stratify, Weight, and Crosstab
        variables. Initializes class variables with the column
        names of the analysis variables in the dataset.
        Parameters:
          inputVariableList (list): A list of dictionaries sorted by strata and PSU
        Returns:
          bool
    """
    self.com = False
    self.outputLevel = 3
    self.booleanLabels = "Yes;No;Missing"
    self.percents = True
    self.booleanValues = False
    self.domain1 = ''
    self.domain2 = ''
   
    for kvp in inputVariableList:
      if kvp.lower() == "percents":
        self.percents = inputVariableList[kvp]

      if kvp.lower() == "stratavar" or kvp.lower() == "stratvarlist":
        self.strataVar = inputVariableList[kvp]

      if kvp.lower() == "numeric_variable" or kvp.lower() == "mainvar" or kvp.lower() == "identifier":
        self.mainVar = inputVariableList[kvp]

      if kvp.lower() == "cross_tabulation_variable" or kvp.lower() == "crosstabvar"or kvp.lower() == "identifier2":
        self.crosstabVar = inputVariableList[kvp]

      if kvp.lower() == "psuvar":
        self.psuVar = inputVariableList[kvp]

      if kvp.lower() == "weightvar":
        self.weightVar = inputVariableList[kvp]

      if kvp.lower() == "tablename":
        self.tableName = inputVariableList[kvp]

    self.cnOutputLevel = 3 #self.outputLevel
    self.cbIncludePercents = self.percents
    
    if self.psuVar is None or len(self.psuVar) == 0:
        self.errorMessage = 'PSU variable is missing'
        self.numErrors += 1
        return False
    
    if self.mainVar is None or len(self.mainVar) == 0:
        self.errorMessage = 'Main variable is missing'
        self.numErrors += 1
        return False

  def Init(self):
    """ Creates the analysis dataset (a list of dicts) by
        subsetting the input dataset to dicts having nonmissing
        values for the analysis variables; then creating a
        list of dicts sorted by strata and PSU.
        Initializes the class variables that hold analysis
        values with the values from the first dict in the
        sorted dataset.
        Parameters:
          None: uses class variables.
        Returns:
          bool
    """
    numRows = -1
    numCats = -1
    numStrata = -1
    self.columnNames = []
    self.isDeleted = False
    self.isVerified = False
    columnNamesArray = []
    
    self.meansResults = CSMeansResults()
    self.meansResults.set_Rows([])
    validCases = 0
    
    self.hasStrataVar = False
    self.hasPsuVar = False
    self.hasMainVar = False
    self.hasCrosstabVar = False
    self.hasWeightVar = False
    if self.strataVar is not None and len(self.strataVar) > 0:
      self.columnNames.append(self.strataVar)
      self.strata = CSField()
      self.strata.set_FieldLabel(self.strataVar)
      self.hasStrataVar = True
    else:
      self.strata = CSField()
      self.strata.set_FieldLabel("None")
      self.strata.set_FieldEntry(1)
    
    if self.weightVar is not None and len(self.weightVar) > 0:
      self.columnNames.append(self.weightVar)
      self.weight = CSField()
      self.weight.set_FieldLabel(self.weightVar)
      self.hasWeightVar = True
    
    if self.mainVar is not None and len(self.mainVar) > 0:
      self.columnNames.append(self.mainVar)
      self.outcome = CSField()
      self.outcome.set_FieldLabel(self.mainVar)
      self.hasMainVar = True
    
    if self.crosstabVar is not None and len(self.crosstabVar) > 0:
      self.columnNames.append(self.crosstabVar)
      self.domain = CSField()
      self.domain.set_FieldLabel(self.crosstabVar)
      self.hasCrosstabVar = True
    
    if self.psuVar is not None and len(self.psuVar) > 0:
      self.columnNames.append(self.psuVar)
      self.psu = CSField()
      self.psu.set_FieldLabel(self.psuVar)
      self.hasPsuVar = True
    if self.confidenceLevel is None or self.confidenceLevel == 0:
      self.confidenceLevel = 0.975 # REVISIT: check if t stat matches epi info's; OK, it does

    for cn in self.columnNames:
      columnNamesArray.append(cn)

    sortVariables = []
    keepVariables = []

    if self.hasStrataVar:
      sortVariables.append(self.strataVar)
      keepVariables.append(self.strataVar)
    if self.hasPsuVar:
      sortVariables.append(self.psuVar)
      keepVariables.append(self.psuVar)
    if self.hasMainVar:
      keepVariables.append(self.mainVar)
    if self.hasCrosstabVar:
      keepVariables.append(self.crosstabVar)
    if self.hasWeightVar:
      keepVariables.append(self.weightVar)

    self.row = 0

    unsortedTable = []
    for d in self.currentTable:
      appendd = True
      dsub = {}
      for v in keepVariables:
        if v not in d:
          appendd = False
          continue
        dsub[v] = d[v]
      if appendd:
        if 'RECSTATUS' in d:
          dsub['RECSTATUS'] = d['RECSTATUS']
        unsortedTable.append(dsub)
   
    self.sortedTable = sorted(unsortedTable, key = lambda ust: ([ust[sv] for sv in sortVariables]))
   
    self.listOfDicts = self.sortedTable
    numRows = len(self.sortedTable)
    
    rowRow = self.sortedTable[0]
    if self.hasMainVar:
      self.outcome.set_FieldEntry(rowRow[self.mainVar])
    if self.hasStrataVar:
      self.strata.set_FieldEntry(rowRow[self.strataVar])
    if self.hasWeightVar:
      self.weight.set_FieldEntry(rowRow[self.weightVar])
    if self.hasPsuVar:
      self.psu.set_FieldEntry(rowRow[self.psuVar])
    if self.hasCrosstabVar:
      self.domain.set_FieldEntry(rowRow[self.crosstabVar])

    if numRows <= 0:
      self.errorMessage = 'No Data available to load'
      self.numErrors += 1
      return False

    sortedTableRD = randata(self.sortedTable)
    if self.hasStrataVar:
      numCats = sortedTableRD.countdistinct([self.psuVar, self.strataVar])
      numStrata = sortedTableRD.countdistinct(self.strataVar)
    else:
      psuVarArray = [self.psuVar] # psuVarArray occurs only here in EICSTables.vb; what is the point? REVISIT?
      numCats = sortedTableRD.countdistinct(self.psuVar)
      numStrata = 1

    if numCats <= 1:
      self.varianceMultiplier = 1.96
    else:
      self.varianceMultiplier = tdist.ppf(self.confidenceLevel, numCats - numStrata)

    return True

  def GetNextRow(self):
    """ Iterates through the analysis dataset and writes the
        dict values to the class analysis variables.
        Parameters:
          None
        Returns:
          bool
    """
    if len(self.listOfDicts) == self.row:
      return False   
    lodRow = self.listOfDicts[self.row]
    if self.hasStrataVar:
      self.strata.set_FieldEntry(lodRow[self.strataVar])
      if len(str(self.strata.get_FieldEntry())) <= 0:
        self.strata.set_cbMissing(True)
      else:
        self.strata.set_cbMissing(False)
    else:
      self.strata.set_FieldEntry = 1
      self.strata.set_cbMissing(False)
    
    if self.hasWeightVar:
      self.weight.set_FieldEntry(lodRow[self.weightVar])
      self.weight.set_FieldReal(float(lodRow[self.weightVar]))
      if len(str(self.weight.get_FieldEntry())) <= 0:
        self.weight.set_cbMissing(True)
      else:
        self.weight.set_cbMissing(False)
   
    if self.hasMainVar:
      self.outcome.set_FieldEntry(lodRow[self.mainVar])
      self.outcome.set_FieldReal(float(lodRow[self.mainVar]))
      if len(str(self.outcome.get_FieldEntry())) <= 0:
        self.outcome.set_cbMissing(True)
      else:
        self.outcome.set_cbMissing(False)
   
    if self.hasCrosstabVar:
      self.domain.set_FieldEntry(lodRow[self.crosstabVar])
      if len(str(self.domain.get_FieldEntry())) <= 0:
        self.domain.set_cbMissing(True)
      else:
        self.domain.set_cbMissing(False)

    if self.hasPsuVar:
      self.psu.set_FieldEntry(lodRow[self.psuVar])
      if len(str(self.psu.get_FieldEntry())) <= 0:
        self.psu.set_cbMissing(True)
      else:
        self.psu.set_cbMissing(False)
    
    if 'RECSTATUS' in self.sortedTable[self.row]:
      recstatus = 1
      recstatus = self.sortedTable[self.row]['RECSTATUS']
      if recstatus < 1:
        self.isDeleted = True

    self.isDeleted = False # What is the point of seting this to True only to immediately set it back to False?
    self.isVerified = True

    self.row += 1
    
    return True

  def NewTot(self, dom):
    """ Creates a new CSMeansTotal object
        Parameters:
          dom (str): the object's crosstab level and/or 'TOTAL'
        Returns:
          CSMeansTotal
    """
    PROC_Name = "clsCMeans::NewTot"

    Ptr = CSMeansTotal()
    Ptr.set_Domain(dom)
    Ptr.set_YE(0)
    Ptr.set_SumW(0)
    Ptr.set_N(0)
    Ptr.set_Min(self.outcome.get_FieldReal())
    Ptr.set_Max(Ptr.get_Min())
    Ptr.set_NextTotal(None)
    return Ptr

  def ResetReader(self):
    """ Resets the dataset row iterator to zero
        Parameters:
          none
        Returns:
          none
    """
    self.row = 0

  def ValidCase(self):
    """ Checks a row for missing values in analysis variables
        Parameters:
          none
        Returns:
          bool: False if any analysis variable value is missing;
                True otherwise
    """
    PROC_Name = "clsCMeans::ValidCase"
   
    ValidCase = True
   
    if self.outcome is not None and ValidCase:
      if self.outcome.get_cbMissing():
        ValidCase = False
   
    if self.strata is not None and ValidCase:
      if self.strata.get_cbMissing():
        ValidCase = False
   
    if self.psu is not None and ValidCase:
      if self.psu.get_cbMissing():
        ValidCase = False
   
    if self.weight is not None and ValidCase:
      if self.weight.get_cbMissing():
        ValidCase = False
   
    if self.domain is not None and ValidCase:
      if self.domain.get_cbMissing():
        ValidCase = False
    
    if ValidCase:
      self.validCases += 1
   
    return ValidCase

  def GetWeight(self):
    """ Returns the value of the weight variable, if present,
        for the current row
        Parameters:
          none
        Returns:
          float
    """
    PROC_Name = "clsCMeans::GetWeight"

    if self.weight is not None:
      return self.weight.get_FieldReal()

    return 1.0

  def AccumYE(self, P):
    """ Uses the current data row's main and weight values
        to adjust the property values of a CSMeansTotal object
        Parameters:
          P (CSMeansTotal)
        Returns:
          none
    """
    PROC_Name = "clsCMeans::AccumYE"
    # P is CSMeansTotal type
    
    Value = self.outcome.get_FieldReal()
    P.set_YE(P.get_YE() + Value * self.GetWeight())
    P.set_SumW(P.get_SumW() + self.GetWeight())
    P.set_N(P.get_N() + 1)
    if P.get_N() == 1:
      P.set_Min(Value)
      P.set_Max(Value)
    else:
      if P.get_Min() > Value:
        P.set_Min(Value)
      elif P.get_Max() < Value:
        P.set_Max(Value)

  def AddTot(self, dom):
    """ Adds additional CSMeansTotal objects for crosstab variable
        values
        Parameters:
          dom (str)
        Returns:
          none
    """
    P = CSMeansTotal()
    inserted = False
    Ptr = self.NewTot(dom)
    self.AccumYE(Ptr)
    if self.first.get_NextTotal() is None:
      self.first.set_NextTotal(Ptr)
      self.last = Ptr
    else:
      P = self.first.get_NextTotal()
      if P.get_Domain() > dom:
        Ptr.set_NextTotal(P)
        self.first.set_NextTotal(Ptr)
      else:
        while P.get_NextTotal() is not None and not inserted:
          if P.get_NextTotal().get_Domain() > dom:
            Ptr.set_NextTotal(P.get_NextTotal())
            P.set_NextTotal(Ptr)
            inserted = True
          else:
            P = P.get_NextTotal()
        if not inserted:
          self.last.set_NextTotal(Ptr)
          self.last = Ptr

  def FindTotal(self, dom):
    """ Returns the CSMeansTotal for the value of dom
        Parameters:
          dom (str)
        Returns:
          CSMeansTotal
    """
    Ptr = self.first
    found = False
    while not found and Ptr is not None:
      if str(Ptr.get_Domain()) == str(dom):
        found = True
      else:
        Ptr = Ptr.get_NextTotal()
    return Ptr

  def FirstPass(self):
    """ The first loop through the dataset and adds the weighted outcome values
        Parameters:
          none
        Returns:
          bool
    """
    PROC_Name = "clsCMeans::FirstPass"
   
    P = CSMeansTotal()
   
    FirstPass = False
    while self.GetNextRow():
      if self.ValidCase() and not self.isDeleted:
        if self.domain is not None:
          P = self.FindTotal(self.domain.get_FieldEntry())
          if self.com:
            if P is not None:
              self.AccumYE(P)
              self.AccumYE(self.first)
          else:
            if P is None:
              self.AddTot(self.domain.get_FieldEntry())
            else:
              self.AccumYE(P)
            self.AccumYE(self.first)
        else:
          self.AccumYE(self.first)
      else:
        self.mis += 1

    return True

  def AccumVar(self, ah):
    """ Computes the variance for each crosstab value
        Parameters:
          ah (int)
        Returns:
          none
    """
    Ptr = self.first
    while Ptr is not None:
      if ah > 1:
        Ptr.set_VarT(Ptr.get_VarT() + (ah * Ptr.get_Sumqha2() - (Ptr.get_Sumqha() ** 2)) / (ah - 1))
      else:
        Ptr.set_VarT(-9999999.0)
      Ptr = Ptr.get_NextTotal()

  def AccumSumq(self):
    """ Computes components of variance
        Parameters:
          none
        Returns:
          none
    """
    Ptr = self.first
    while Ptr is not None:
      Ptr.set_Sumqha(Ptr.get_Sumqha() + Ptr.get_qha())
      Ptr.set_Sumqha2(Ptr.get_Sumqha2() + Ptr.get_qha() ** 2)
      Ptr = Ptr.get_NextTotal()

  def Accumqha(self, P):
    """ Computes components of variance
        Parameters:
          none
        Returns:
          none
    """
    Qhab = None
    if P.get_SumW() > 0:
      Qhab = (self.outcome.get_FieldReal() - (P.get_YE() / P.get_SumW())) * (self.GetWeight() / P.get_SumW())
    else:
      Qhab = 0.0
    P.set_qha(P.get_qha() + Qhab)

  def Qhab(self, P):
    """ Computes components of variance
        Parameters:
          none
        Returns:
          bool
    """
    Qhab = None
    if P.get_SumW() > 0:
      Qhab = (self.outcome.get_FieldReal() - (P.get_YE() / P.get_SumW())) * (self.GetWeight() / P.get_SumW())
    else:
      Qhab = 0.0
    return Qhab

  def SumqInit(self):
    Ptr = self.first
    while Ptr is not None:
      Ptr.set_Sumqha(0.0)
      Ptr.set_Sumqha2(0.0)
      Ptr = Ptr.get_NextTotal()

  def QhaInit(self):
    """ Initializes qha at zero for all CSMeansTotal objects
        Parameters:
          none
        Returns:
          none
    """
    Ptr = self.first
    while Ptr is not None:
      Ptr.set_qha(0.0)
      Ptr.set_qha2(0.0)
      Ptr = Ptr.get_NextTotal()

  def VarTInit(self):
    """ Initializes variance at zero for all CSMeansTotal objects
        Parameters:
          none
        Returns:
          none
    """
    Ptr = self.first
    while Ptr is not None:
      Ptr.set_VarT(0.0)
      Ptr = Ptr.get_NextTotal()

  def FieldColl(self, p1, s):
    """ Compares the strata or PSU values in two items of data
        Parameters:
          p1 (str, float, or int): a value from the strata or PUS variable
          s  (str, float, or int): a value from the strata or PUS variable
        Returns:
          int indicating greater than, less than, or equal
    """
    ft = None # integer
    i = None #  integer
    R = None #  double
    R2 = None # double
    FieldColl = 0
    
    if str(p1.get_FieldEntry()).isnumeric():
      if float(p1.get_FieldEntry()) % 1 == 0:
        i = int(s)
        FieldColl = p1.get_FieldInt() - i
      else:
        R = float(s)
        R2 = p1.get_FieldReal()
        if R2 > R:
          FieldColl = 1
        elif R2 < R:
          FieldColl = -1
        else:
          FieldColl = 0
    else:
      if float(p1.get_FieldEntry()) > float(s):
        FieldColl = 1
      elif float(p1.get_FieldEntry()) < float(s):
        FieldColl = -1
      else:
        FieldColl = 0
    return FieldColl

  def SecondPass(self):
    """ Loops over the analysis dataset computing results
        Parameters:
          none
        Returns:
          bool
    """
    PROC_Name = "clsCMeans::SecondPass"
    P = CSMeansTotal()
    Valid = False
    ah = 0
    Rec = True
    NowStrat = ""
    NowPSU = ""
    qha = None
    qha2 = None
    Sumqha = None
    Sumqha2 = None
    bContinue = True
    bHadValidPSU = True
    
    SecondPass = False
    self.VarTInit()
    self.varT = 0
    
    while Rec and not Valid:
      Rec = self.GetNextRow()
      if Rec:
        Valid = self.ValidCase()
        if Valid and not self.isDeleted:
          if self.strata is not None:
            NowStrat = self.strata.get_FieldEntry()
          if self.psu is not None:
            NowPSU = self.psu.get_FieldEntry()
    
    while True:
      self.SumqInit()
      Sumqha = 0.0
      Sumqha2 = 0.0
      ah = 0
      while True:
        qha = 0.0
        qha2 = 0.0
        self.QhaInit()
        bHadValidPSU = False
        while True:
          if self.ValidCase() and not self.isDeleted:
            bHadValidPSU = True
            if not self.com:
              if self.domain is not None:
                P = self.FindTotal(self.domain.get_FieldEntry())
                self.Accumqha(P)
              self.Accumqha(self.first)
            else:
              P = self.FindTotal(self.domain.get_FieldEntry())
              if P == self.first.get_NextTotal():
                qha += self.Qhab(P)
              else:
                if P == self.last:
                  qha -= self.Qhab(P)
              if P is not None:
                self.Accumqha(self.first)
                self.Accumqha(P)
          Rec = self.GetNextRow()
          if self.psu is not None:
            if self.psu.get_FieldEntry() != NowPSU:
              bContinue = True
            elif self.strata is not None and self.strata.get_FieldEntry() != NowStrat:
              bContinue = True
            else:
              bContinue = False
          else:
            bContinue = True
          if bContinue or Rec == False:
            break
        if self.psu is not None:
          if Rec == False or self.FieldColl(self.psu, NowPSU) > 0:
            NowPSU = self.psu.get_FieldEntry()
          elif self.strata is not None:
            if self.strata.get_FieldEntry() != NowStrat:
              NowPSU = self.psu.get_FieldEntry()
            else:
              pass
              #self.errorMessage = "File is not sorted!"
          else:
            self.errorMessage = "File is not sorted!"
            return False
        if bHadValidPSU:
          ah += 1
          self.AccumSumq()
          if self.com:
            Sumqha += qha
            Sumqha2 += qha ** 2
        if self.strata is not None:
          if self.strata.get_FieldEntry() != NowStrat:
            bContinue = True
          else:
            bContinue = False
        else:
          bContinue = True
        if bContinue or Rec == False:
          break
      if self.strata is not None:
        if Rec == False or self.FieldColl(self.strata, NowStrat) > 0:
          NowStrat = self.strata.get_FieldEntry()
        else:
          self.errorMessage = "File is not sorted!"
          SecondPass = False
          self.numErrors += 1
          return SecondPass
      self.AccumVar(ah)
      if ah > 1 and self.com:
        self.varT = self.varT + (ah * Sumqha2 - (Sumqha ** 2)) / (ah - 1)
      if Rec == False:
        break
    
    SecondPass = True
    return SecondPass

  def PrintValues(self, errorMessage):
    """ Computes the final statistical output of the analysis
        and stores the results for TOTAL and each crosstab value
        in a list of MeansRow objects
        Parameters:
          errorMessage (str)
        Returns:
          none: The resulting list is a property of the meansResults
                class variable.
    """
    Ptr = CSMeansTotal()
    Lo = None
    Up = None
    Diff = None
    i = 0
    sOutline = ''
    nOutfile = 0
    
    if self.cnOutputLevel > 0:
      if self.cbStandalone:
        # This is just building text output in VB
        pass
      # Lots of text output building in VB
      if self.domain is not None:
        Ptr = self.first.get_NextTotal()
      else:
        Ptr = self.first
      while Ptr is not None:
        mRow = MeansRow()
        mRow.set_Label(Ptr.get_Domain())
        if self.cnOutputLevel > 1:
          mRow.set_Count(float(Ptr.get_N()))
        if self.cnOutputLevel > 0:
          if Ptr.get_SumW() > 0:
            mRow.set_Mean(Ptr.get_YE() / Ptr.get_SumW())
          else:
            mRow.set_Mean(None)
        if self.cnOutputLevel > 2:
          if Ptr.get_VarT() > 0:
            mRow.set_StdErr(Ptr.get_VarT() ** 0.5)
          else:
            mRow.set_SteErr(None)
        if self.cnOutputLevel > 1 and self.cbIncludePercents:
          if Ptr.get_SumW() > 0 and Ptr.get_VarT() > 0:
            Lo = (Ptr.get_YE() / Ptr.get_SumW()) - (self.varianceMultiplier * Ptr.get_VarT() ** 0.5)
            mRow.set_LCL(Lo)
            Up = (Ptr.get_YE() / Ptr.get_SumW()) + (self.varianceMultiplier * Ptr.get_VarT() ** 0.5)
            mRow.set_UCL(Up)
        if self.cnOutputLevel:
          mRow.set_Min(Ptr.get_Min())
          mRow.set_Max(Ptr.get_Max())
        if Ptr == self.first:
          Ptr = None
        else:
          Ptr = Ptr.get_NextTotal()
          if Ptr is None:
            Ptr = self.first
        self.meansResults.get_Rows().append(mRow)
      if self.com and self.cnOutputLevel > 2:
        dRow = MeansRow()
        dRow.set_Label("Difference")
        if self.first.get_NextTotal().get_SumW() > 0 and self.last.get_SumW() > 0:
          Diff = (self.first.get_NextTotal().get_YE() / self.first.get_NextTotal().get_SumW()) - \
                 (self.last.get_YE() / self.last.get_SumW())
          dRow.set_Mean(Diff)
          if self.varT > 0:
            dRow.set_StdErr(self.varT ** 0.5)
            Lo = Diff - (self.varianceMultiplier * self.varT ** 0.5)
            Up = Diff + (self.varianceMultiplier * self.varT ** 0.5)
            dRow.set_LCL(Lo)
            dRow.set_UCL(Up)
          else:
            dRow.set_StdErr(None)
        else:
          dRow.set_Mean(None)
          dRow.set_StdErr(None)
        self.meansResults.get_Rows().append(dRow)

  def ComplexSampleMeans(self, inputVariableList, dataTable):
    """ Executes the supporting functions to run the analysis
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
          dataTable (list(dict)): The analysis dataset
        Returns:
          self.meansResuts (CSMeansResults): This object contains a Rows property.
          It is a list of MeansRow objects, which have properties: Label, Count,
          Mean, StdErr, LCL, and UCL. These are the displayed output of the analysis.
          There is a MeansRow for TOTAL and one for each value of the crosstab
          variable, if present.
    """
    csfstarttime = time.time()
    self.currentTable = dataTable

    self.CreateSettings(inputVariableList)

    self.errorMessage = ''

    self.numErrors = 0
    output = []

    self.meansResults.set_ErrorMessage('')

    if self.Init() == False:
      self.meansResults.set_ErrorMessage('There was a problem initializing the statistics.')
      return self.meansResults

    self.mis = 0
    self.first = self.NewTot("TOTAL")
    self.last = self.first

    self.first.set_NextTotal(None)
    result = None
   
    if self.com:
      self.GetNextRow()
      if self.Domain1 < self.Domain2:
        self.first.set_NextTotal(self.NewTot(self.Domain1))
        self.last = self.NewTot(self.Domain2)
      else:
        self.first.set_NextTotal(self.NewTot(self.Domain2))
        self.last = self.NewTot(self.Domain1)
      self.first.get_NextTotal().set_NextTotal(self.last)
      self.ResetReader()
    
    result = self.FirstPass()
    
    if self.errorMessage is not None and len(self.errorMessage) > 0:
      self.meansResults.set_ErrorMessage(self.errorMessage)
      return self.meansResults
    
    if self.first.get_NextTotal() is not None:
      if self.first.get_NextTotal().get_NextTotal() is not None and self.first.get_NextTotal().get_NextTotal() == self.last:
        com = True
    
    self.ResetReader()
    
    if result:
      self.errorMessage = ''
      result = self.SecondPass()
      if (self.errorMessage is not None and len(self.errorMessage) > 0) or self.numErrors > 0:
        self.meansResults.set_ErrorMessage(self.errorMessage)
        return self.meansResults
    
    if result:
      if self.cnOutputLevel > 0:
        self.errorMessage = ''
        self.PrintValues(self.errorMessage)
        if (self.errorMessage is not None and len(self.errorMessage) > 0) or self.numErrors > 0:
          self.meansResults.set_ErrorMessage(self.errorMessage)
          return self.meansResults
    
    return self.meansResults
