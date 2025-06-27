from scipy.stats import t as tdist
import math
import time
from .randata import randata, csvToRandata
from .CSUtilities import *
class ComplexSampleTables:
  def __init__(self):
    self.columnNames = None
    self.com = None
    self.outputLevel = None
    self.booleanLabels = None
    self.percents = None
    self.booleanValues = None
    self.domain1 = None
    self.domain2 = None
    self.identifiers = None
    self.mainVar = None
    self.strataVar = None
    self.domainVar = None
    self.psuVar = None
    self.weightVar = None
    self.tableName = None
    self.vntLabels = None
    self.cnOutputLevel = None
    self.cbIncludePercents = None
    self.isDeleted = None
    self.isVerified = None
    self.outcome = None
    self.domain = None
    self.strata = None
    self.psu = None
    self.weight = None
    self.confidenceLevel = None
    self.Row = None
    self.currentTable = None
    self.varianceMultiplier = None
    self.FirstDom = None
    self.LastDom = None
    self.TotalDom = None
    self.FirstCat = None
    self.LastCat = None
    self.Mis = None
    self.Pdm = None
    self.sortedTable = None
    self.collectedSortedTable = None
    self.distinctTable = None
    self.T22 = None
    self.vntResultsArray = []
    self.columnPercents = []
  
  def ResetReader(self):
    """ Resets the dataset row iterator to zero
        Parameters:
          none
        Returns:
          none
    """
    self.Row = 0
  
  def Init(self, errorMessage):
    """ Creates the analysis dataset (a list of dictionaries) by
        subsetting the input dataset to dictionaries having nonmissing
        values for the analysis variables; then creating a
        list of dicts sorted by strata and PSU.
        Initializes the class variables that hold analysis
        values with the values from the first dictionary in the
        sorted dataset.
        Parameters:
          errorMessage (str)
        Returns:
          bool
    """
    initstarttime = time.time()
    validCases = 0 # REMOVE
    numRows = -1
    numCats = -1
    numStrata = -1
    self.columnNames = []
    self.isDeleted = False
    self.isVerified = False
    columnNamesArray = []
    
    if self.strataVar is not None and len(self.strataVar) > 0:
      self.columnNames.append(self.strataVar)
      self.strata = CSField()
      self.strata.set_FieldLabel(self.strataVar)
    else:
      self.strata = CSField()
      self.strata.set_FieldLabel("None")
      self.strata.set_FieldEntry(1)
    
    if self.weightVar is not None and len(self.weightVar) > 0:
      self.columnNames.append(self.weightVar)
      self.weight = CSField()
      self.weight.set_FieldLabel(self.weightVar)
    
    if self.mainVar is not None and len(self.mainVar) > 0:
      self.columnNames.append(self.mainVar)
      self.outcome = CSField()
      self.outcome.set_FieldLabel(self.mainVar)
    
    if self.domainVar is not None and len(self.domainVar) > 0:
      self.columnNames.append(self.domainVar)
      self.domain = CSField()
      self.domain.set_FieldLabel(self.domainVar)
    
    if self.psuVar is not None and len(self.psuVar) > 0:
      self.columnNames.append(self.psuVar)
      self.psu = CSField()
      self.psu.set_FieldLabel(self.psuVar)
      if self.confidenceLevel is None or self.confidenceLevel == 0:
        self.confidenceLevel = 0.975 # REVISIT: check if t stat matches epi info's; OK, it does
    
    for cn in self.columnNames:
      columnNamesArray.append(cn)
    
    sortVariables = []
    
    if self.strataVar is not None and len(self.strataVar) > 0:
      sortVariables.append(self.strataVar)
    if self.psuVar is not None and len(self.psuVar) > 0:
      sortVariables.append(self.psuVar)
    if self.mainVar is not None and len(self.mainVar) > 0:
      sortVariables.append(self.mainVar)
    if self.domainVar is not None and len(self.domainVar) > 0:
      sortVariables.append(self.domainVar)
    if self.weightVar is not None and len(self.weightVar) > 0:
      sortVariables.append(self.weightVar)
    
    self.Row = 0
    
    whereClause = self.mainVar + ' is not null and ' + self.psuVar + ' is not null'
    if self.weightVar is not None and len(self.weightVar) > 0:
      whereClause = whereClause + ' and ' + self.weightVar + ' is not null'
    if self.strataVar is not None and len(self.strataVar) > 0:
      whereClause = whereClause + ' and ' + self.strataVar + ' is not null'
    
    unsortedTable = self.currentTable
    #unsortedTable.registerTempTable('unsortedTable')
    self.sortedTable = sorted(self.currentTable, key = lambda eics: ([eics[sv] for sv in sortVariables]))
    sql_Statement = 'select (row_number() over (order by "row_index")) - 1 as _row_number_'
    if self.mainVar is not None and len(self.mainVar) > 0:
      sql_Statement += ', ' + self.mainVar
    if self.strataVar is not None and len(self.strataVar) > 0:
      sql_Statement += ', ' + self.strataVar
    if self.weightVar is not None and len(self.weightVar) > 0:
      sql_Statement += ', ' + self.weightVar
    if self.psuVar is not None and len(self.psuVar) > 0:
      sql_Statement += ', ' + self.psuVar
    if self.domainVar is not None and len(self.domainVar) > 0:
      sql_Statement += ', ' + self.domainVar
    sql_Statement += ' from psdf'
    #self.sortedTable = spark.sql(sql_Statement)
    self.collectedSortedTable = self.sortedTable
    self.listOfDicts = self.sortedTable
    self.mainVarSet = []
    self.strataVarSet = []
    self.weightVarSet = []
    self.psuVarSet = []
    self.domainVarSet = []
    for dct in self.listOfDicts:
        if self.mainVar in sortVariables:
            self.mainVarSet.append(dct[self.mainVar])
        if self.strataVar in sortVariables:
            self.strataVarSet.append(dct[self.strataVar])
        if self.weightVar in sortVariables:
            self.weightVarSet.append(dct[self.weightVar])
        if self.psuVar in sortVariables:
            self.psuVarSet.append(dct[self.psuVar])
        if self.domainVar in sortVariables:
            self.domainVarSet.append(dct[self.domainVar])
    numRows = len(self.sortedTable)
    
    rowRow = self.sortedTable[0]
    if self.mainVar is not None and len(self.mainVar) > 0:
      self.outcome.set_FieldEntry(rowRow[self.mainVar])
    if self.strataVar is not None and len(self.strataVar) > 0:
      self.strata.set_FieldEntry(rowRow[self.strataVar])
    if self.weightVar is not None and len(self.weightVar) > 0:
      self.weight.set_FieldEntry(rowRow[self.weightVar])
    if self.psuVar is not None and len(self.psuVar) > 0:
      self.psu.set_FieldEntry(rowRow[self.psuVar])
    if self.domainVar is not None and len(self.domainVar) > 0:
      self.domain.set_FieldEntry(rowRow[self.domainVar])
    
    if numRows <= 0:
      errorMessage = 'No Data available to load'
      return False
    
    sortedTableRD = randata(self.sortedTable)
    if self.strataVar is not None and len(self.strataVar) > 0:
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
    
  def CreateSettings(self, inputVariableList):
    """ Initializes objectes necessary for the analysis.
        Checks for existance of Stratify, Weight, Exposure, and
        Outcome variables. Initializes class variables with the
        column names of the analysis variables in the dataset.
        Parameters:
          inputVariableList (list): A list of dictionaries sorted by strata and PSU
        Returns:
          none
    """
    self.com = False
    self.outputLevel = 3
    self.booleanLabels = "Yes;No;Missing"
    self.percents = True
    self.booleanValues = False
    self.domain1 = ''
    self.domain2 = ''
    
    if "IdentifierList" in inputVariableList:
      self.identifiers = inputVariableList["IdentifierList"].split(",")
      self.mainVar = self.identifiers[0]

    for kvp in inputVariableList:
      if kvp.lower() == "percents":
        self.percents = inputVariableList[kvp]

      if kvp.lower() == "stratavar" or kvp.lower() == "stratvarlist":
        self.strataVar = inputVariableList[kvp]

      if kvp.lower() == "exposure_variable" or kvp.lower() == "mainvar" or kvp.lower() == "identifier" or \
         kvp.lower() == "identifier1" or kvp.lower() == "exposurevar":
        self.domainVar = inputVariableList[kvp]

      if kvp.lower() == "outcome_variable" or kvp.lower() == "crosstabvar"or kvp.lower() == "identifier2" or kvp.lower() == "outcomevar":
        self.mainVar = inputVariableList[kvp]

      if kvp.lower() == "psuvar":
        self.psuVar = inputVariableList[kvp]

      if kvp.lower() == "weightvar":
        self.weightVar = inputVariableList[kvp]

      if kvp.lower() == "tablename":
        self.tableName = inputVariableList[kvp]

    self.vntLabels = []
    self.vntLabels.append('')
    self.vntLabels.append('')
    self.vntLabels.append('')

    # TODO: Get from configuration
    self.vntLabels[0] = "Yes" #context.SetProperties("RepresentationOfYes")
    self.vntLabels[1] = "No" #context.SetProperties("RepresentationOfNo")
    self.vntLabels[2] = "Missing" #context.SetProperties("RepresentationOfMissing")

    self.cnOutputLevel = self.outputLevel
    self.cbIncludePercents = self.percents
  
  def GetNextRow(self):
    """ Iterates through the analysis dataset and writes the
        dict values to the class analysis variables.
        Parameters:
          None
        Returns:
          bool
    """
    if len(self.sortedTable) == self.Row:
      return False
    
    if self.strataVar is not None and len(self.strataVar) > 0:
      self.strata.set_FieldEntry(self.listOfDicts[self.Row][self.strataVar])
      if len(str(self.strata.get_FieldEntry())) <= 0:
        self.strata.set_cbMissing(True)
      else:
        self.strata.set_cbMissing(False)
    else:
      self.strata.set_FieldEntry = 1
      self.strata.set_cbMissing(False)
    
    if self.weightVar is not None and len(self.weightVar) > 0:
      self.weight.set_FieldEntry(self.listOfDicts[self.Row][self.weightVar])
      if len(str(self.weight.get_FieldEntry())) <= 0:
        self.weight.set_cbMissing(True)
      else:
        self.weight.set_cbMissing(False)
    
    if self.mainVar is not None and len(self.mainVar) > 0:
      self.outcome.set_FieldEntry(self.mainVarSet[self.Row])
      if len(str(self.outcome.get_FieldEntry())) <= 0:
        self.outcome.set_cbMissing(True)
      else:
        self.outcome.set_cbMissing(False)
    
    if self.domainVar is not None and len(self.domainVar) > 0:
      self.domain.set_FieldEntry(self.listOfDicts[self.Row][self.domainVar])
      if len(str(self.domain.get_FieldEntry())) <= 0:
        self.domain.set_cbMissing(True)
      else:
        self.domain.set_cbMissing(False)
    
    if self.psuVar is not None and len(self.psuVar) > 0:
      self.psu.set_FieldEntry(self.psuVarSet[self.Row])
      if len(str(self.psu.get_FieldEntry())) <= 0:
        self.psu.set_cbMissing(True)
      else:
        self.psu.set_cbMissing(False)
    
    if 'RECSTATUS' in self.sortedTable[self.Row]:
      recstatus = 1
      recstatus = self.listOfDicts[self.Row]['RECSTATUS']
      if recstatus < 1:
        self.isDeleted = True
    
    self.isDeleted = False # What is the point of conditionally seting this to True only to immediately set it back to False?
    self.isVerified = True
    
    self.Row += 1
    
    return True
  
  def ValidCase(self):
    """ Checks a row for missing values in analysis variables
        Parameters:
          none
        Returns:
          bool: False if any analysis variable value is missing;
                True otherwise
    """
    PROC_Name = "clsCTables::ValidCase"
    
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
    
    return ValidCase
  
  def NewTot(self, dom, cat):
    """ Creates a new CSTotal object
        Parameters:
          dom (str): the object's exposure variable level and/or 'TOTAL'
          cat (str): the object's outcome variable level
        Returns:
          CSTotal
    """
    PROC_Name = "clsCTables::NewTot"
    
    Ptr = CSTotal()
    Ptr.set_Domain(dom)
    Ptr.set_Category(cat)
    Ptr.set_YE(0)
    Ptr.set_N(0)
    Ptr.set_NextDom(None)
    Ptr.set_NextCat(None)
    return Ptr
  
  def NewCat(self, cat):
    """ Creates a new CSCategory object
        Parameters:
          cat (str): the object's outcome variable level
        Returns:
          CSCategory
    """
    Pd = CSCategory()
    Pd.set_Category(cat)
    Pd.set_NextCat(None)
    return Pd
  
  def WhereIns(self, dom):
    """ Attempts to return a CSDomain object for the
        value of dom
        Parameters:
          dom (str): the object's exposure variable level and/or 'TOTAL'
        Returns:
          CSDomain
    """
    PROC_Name = "clsCTables::WhereIns"
    
    Pd = CSDomain()
    found = False
    
    WhereIns = Pd
    
    if str(dom) == 'TOTAL' or str(dom) == 'Difference':
      WhereIns = self.LastDom
    else:
      if self.FirstDom.get_Domain() > dom: # REVISIT: this assumes dom is a number but above it is a string
        WhereIns = None
      else:
        Pd = self.FirstDom
        while Pd.get_NextDom() is not None and not found:
          if Pd.get_NextDom().get_Domain() > dom:
            WhereIns = Pd
            found = True
          else:
            Pd = Pd.get_NextDom()
        if not found:
          WhereIns = self.LastDom
    
    return WhereIns
  
  def WhereInsC(self, cat):
    """ Attempts to return a CSCategory object for the
        value of cat
        Parameters:
          cat (str): the object's outcome variable level
        Returns:
          CSCategory
    """
    PROC_Name = "clsCTables::WhereInsC"
    
    Pc = CSCategory()
    found = False
    WhereInsC = None
    
    if self.FirstCat.get_Category() > cat:
      WhereInsC = None
    else:
      Pc = self.FirstCat
      while Pc.get_NextCat() is not None and not found:
        if Pc.get_NextCat().get_Category() > cat:
          WhereInsC = Pc
          found = True
        else:
          Pc = Pc.get_NextCat()
      if not found:
        WhereInsC = self.LastCat
    
    return WhereInsC
  
  def FindDom(self, dom):
    """ Attempts to find the CSDomain object for the
        value of dom
        Parameters:
          dom (str): the object's exposure variable level and/or 'TOTAL'
        Returns:
          CSDomain
    """
    PROC_Name = "clsCTables::FindDom"
    
    Pd = CSDomain()
    found = False
    
    Pd = self.FirstDom
    
    while Pd is not None and not found:
      if Pd.get_Domain() == dom:
        found = True
      else:
        Pd = Pd.get_NextDom()
    
    return Pd
  
  def AddDom(self, dom):
    """ Creates a new CSDomain object for the
        value of dom
        Parameters:
          dom (str): the object's exposure variable level and/or 'TOTAL'
        Returns:
          CSDomain
    """
    PROC_Name = "clsCTables::AddDom"
    
    Pd = CSDomain()
    Pc = CSCategory()
    R2 = CSCategory()
    Pt = CSDomain()
    P = CSTotal()
    Q = CSTotal()
    R = CSTotal()
    
    Pd.set_Domain(dom)
    Pd.set_SumW(self.GetWeight())
    Pd.set_N(1)
    Pd.set_NextDom(None)
    
    if self.LastDom is None:
      Pc = self.NewCat(self.outcome.get_FieldEntry())
      self.FirstCat = Pc
      self.LastCat = Pc
      P = self.NewTot(dom, self.outcome.get_FieldEntry())
      Pd.set_FirstCat(P)
      Pc.set_FirstDom(P)
      self.FirstDom = Pd
      self.LastDom = Pd
    else:
      Pt = self.WhereIns(dom)
      if Pt is None:
        Pd.set_NextDom(self.FirstDom)
        R2 = self.FirstCat
        P = self.NewTot(dom, R2.get_Category())
        P.set_NextDom(R2.get_FirstDom())
        Pd.set_FirstCat(P)
        R2.set_FirstDom(P)
        R2 = R2.get_NextCat()
        Q = P
        while R2 is not None:
          P = self.NewTot(dom, R2.get_Category())
          P.set_NextDom(R2.get_FirstDom())
          Q.set_NextCat(P)
          R2.set_FirstDom(P)
          R2 = R2.get_NextCat()
          Q = P
        self.FirstDom = Pd
      else:
        Pd.set_NextDom(Pt.get_NextDom())
        Pt.set_NextDom(Pd)
        R = Pt.get_FirstCat()
        P = self.NewTot(dom, R.get_Category())
        P.set_NextDom(R.get_NextDom())
        Pd.set_FirstCat(P)
        R.set_NextDom(P)
        R = R.get_NextCat()
        Q = P
        while R is not None:
          P = self.NewTot(dom, R.get_Category())
          P.set_NextDom(R.get_NextDom())
          Q.set_NextCat(P)
          R.set_NextDom(P)
          R = R.get_NextCat()
          Q = P
        if Pt == self.LastDom:
          self.LastDom = Pd
    AddDom = Pd
    return AddDom
  
  def GetWeight(self):
    """ Returns the value of the weight variable, if present,
        for the current row
        Parameters:
          none
        Returns:
          float
    """
    PROC_Name = "clsCTables::GetWeight"
    
    if self.weight is not None:
      return self.weight.get_FieldEntry()
    
    return 1.0
  
  def FindCat(self, P, cat):
    """ Attempts to find the CSCategory object for the
        value of cat
        Parameters:
          P (CSTotal)
          cat (str): the object's outcome variable level
        Returns:
          CSCategory
    """
    PROC_Name = "clsCTables::FindCat"
    
    found = False
    while P is not None and not found:
      if P.get_Category() == cat:
        found = True
      else:
        P = P.get_NextCat()
    FindCat = P
    return FindCat
  
  def AddCat(self, cat, dom):
    """ Creates a new CSCategory object for the
        value of cat
        Parameters:
          dom (str): the object's exposure variable level and/or 'TOTAL'
          cat (str): the object's outcome variable level
        Returns:
          CSCategory
    """
    PROC_Name = "clsCTables::AddCat"
    
    Pd = CSDomain()
    R2 = CSDomain()
    Pc = CSCategory()
    Pt = CSCategory()
    P = CSTotal()
    Q = CSTotal()
    R = CSTotal()
    
    Pc.set_Category(cat)
    Pc.set_NextCat(None)
    Pc.set_FirstDom(None)
    
    Pt = self.WhereInsC(cat)
    
    if Pt is None:
      Pc.set_NextCat(self.FirstCat)
      R2 = self.FirstDom
      P = self.NewTot(R2.get_Domain(), cat)
      if R2.get_Domain() == dom:
        self.AccumYE(P)
      P.set_NextCat(R2.get_FirstCat())
      Pc.set_FirstDom(P)
      R2.set_FirstCat(P)
      R2 = R2.get_NextDom()
      Q = P
      while R2 is not None:
        P = self.NewTot(R2.get_Domain(), cat)
        if R2.get_Domain() == dom:
          self.AccumYE(P)
        P.set_NextCat(R2.get_FirstCat())
        Q.set(NextDom(P))
        R2.set_FirstCat(P)
        R2 = R2.get_NextDom()
        Q = P
      self.FirstCat = Pc
    else:
      Pc.set_NextCat(Pt.get_NextCat())
      Pt.set_NextCat(Pc)
      R = Pt.get_FirstDom()
      P = self.NewTot(R.get_Domain(), cat)
      if R.get_Domain() == dom:
        self.AccumYE(P)
      P.set_NextCat(R.get_NextCat())
      Pc.set_FirstDom(P)
      R.set_NextCat(P)
      R = R.get_NextDom()
      Q = P
      while R is not None:
        P = self.NewTot(R.get_Domain(), cat)
        if R.get_Domain() == dom:
          self.AccumYE(P)
        P.set_NextCat(R.get_NextCat())
        Q.set_NextDom(P)
        R.set_NextCat(P)
        R = R.get_NextDom()
        Q = P
      if Pt == self.LastCat:
        self.LastCat = Pc
  
  def AccumYE(self, P):
    """ Uses the current data row's main and weight values
        to adjust the property values of a CSTotal object
        Parameters:
          P (CSTotal)
        Returns:
          none
    """
    PROC_Name = "clsCTables::AccumYE"
    P.set_YE(float(P.get_YE()) + float(self.GetWeight()))
    P.set_N(P.get_N() + 1)
  
  def FirstPassCom(self, errorMessage):
    """ The first loop through the dataset and adds the weighted outcome counts
        Parameters:
          errorMessage (str)
        Returns:
          bool
    """
    PROC_Name = "clsCTables::FirstPassCom"
    
    Pd = CSDomain()
    P = CSTotal()
    N = None
    Valid = None
    doma = None
    cate = None
    
    FirstPassCom = False
    
    while self.GetNextRow():
      Valid = self.ValidCase()
      doma = self.Domain.get_FieldEntry()
      if not Valid:
        self.Mis += 1
      if Valid and ((doma == self.domain1) or (doma == self.domain2)):
        cate = self.outcome.get_FieldEntry()
        Pd = self.FindDom(doma)
        if Pd is None:
          Pd = self.AddDom(doma)
        else:
          Pd.set_SumW(Pd.get_SumW() + self.GetWeight())
          Pd.set_N(pd.get_N() + 1)
        P = self.FindCat(Pd.get_FirstCat(), cate)
        if P is None:
          self.AddCat(cate, doma)
        else:
          self.AccumYE(P)
    
    FirstPassCom = True
    return FirstPassCom
  
  def FirstPass(self, errorMessage):
    """ The first loop through the dataset and adds the weighted outcome counts
        Parameters:
          errorMessage (str)
        Returns:
          bool
    """
    PROC_Name = "clsCTables::FirstPass"
    
    Pd = CSDomain()
    P = CSTotal()
    cate = None
    doma = None
    
    FirstPass = False
    
    while self.GetNextRow():
      if self.ValidCase():
        if self.domain is not None:
          doma = self.domain.get_FieldEntry()
          Pd = self.FindDom(doma)
          if Pd is None:
            Pd = self.AddDom(doma)
          else:
            Pd.set_SumW(Pd.get_SumW() + self.GetWeight())
            Pd.set_N(Pd.get_N() + 1)
        else:
          doma = "TOTAL"
          Pd = self.FirstDom
          if Pd is None:
            Pd = self.AddDom(doma)
          else:
            Pd.set_SumW(Pd.get_SumW() + self.GetWeight())
            Pd.set_N(Pd.get_N() + 1)
            self.TotalDom = Pd
        cate = self.outcome.get_FieldEntry()
        P = self.FindCat(Pd.get_FirstCat(), cate)
        if P is None:
          self.AddCat(cate, doma)
        else:
          self.AccumYE(P)
      else:
        self.Mis += 1
    
    return True
  
  def LnOf(self, X):
    """ Computes the natural logarithm of a number
        Parameters:
          X (float)
        Returns:
          float
    """
    PROC_Name = "clsCTables::LnOf"
    Y = (float(X) - 1.0) / (float(X) + 1.0)
    Y2 = Y * Y
    LnOf = float(0.0)
    num = Y
    Den = 1.0
    while True:
      LastLn = LnOf
      LnOf += num / Den
      num *= Y2
      Den += 2
      if LnOf == LastLn:
        break
    LnOf *= 2.0
    return LnOf
  
  def ExpOf(self, X):
    """ Computes e to a power
        Parameters:
          X (float)
        Returns:
          float
    """
    PROC_Name = "clsCTables::ExpOf"
    Ex = float(X) + 1.0
    Term = float(X)
    N = 1.0
    while True:
      N += 1.0
      Term *= X
      Term /= N
      LastEx = Ex
      Ex += Term
      if Ex == LastEx:
        break
    ExpOf = Ex
    return ExpOf
  
  def PrintRisk(self, RR, ODD, RD, VarRR, VarOD, VarRD, VarOR, VarLnRR):
    """ Stores Odds and Risk statistics in a class list of lists
        Parameters:
          RR      (float)
          ODD     (float)
          RD      (float)
          VarRR   (float)
          VarOD   (float)
          VarRD   (float)
          VarOR   (float)
          VarLnRR (float)
        Returns:
          none
    """
    PROC_Name = "clsCTables::PrintRisk"
    s1, s2, s3, s4 = ('html', 'html', 'html', 'html')
    
    self.vntResultsArray[0][8] = "Odds Ratio (OR)"
    self.vntResultsArray[0][9] = "Standard Error for OR"
    self.vntResultsArray[0][10] = "Lower 95% Confidence Limit for OR"
    self.vntResultsArray[0][11] = "Upper 95% Confidence Limit for OR"
    self.vntResultsArray[0][12] = "Risk Ratio (RR)"
    self.vntResultsArray[0][13] = "Standard Error for RR"
    self.vntResultsArray[0][14] = "Lower 95% Confidence Limit for RR"
    self.vntResultsArray[0][15] = "Upper 95% Confidence Limit for RR"
    self.vntResultsArray[0][16] = "Risk Difference (RD%)"
    self.vntResultsArray[0][17] = "Standard Error for RD(%)"
    self.vntResultsArray[0][18] = "Lower 95% Confidence Limit for RD(%)"
    self.vntResultsArray[0][19] = "Upper 95% Confidence Limit for RD(%)"
    
    if ODD > 0.0:
      self.vntResultsArray[1][8] = ODD
      self.vntResultsArray[1][9] = VarOD ** 0.5
      self.vntResultsArray[1][10] = self.ExpOf(self.LnOf(ODD) - (self.varianceMultiplier * VarOR ** 0.5))
      self.vntResultsArray[1][11] = self.ExpOf(self.LnOf(ODD) + (self.varianceMultiplier * VarOR ** 0.5))
    else:
      if self.cnOutputLevel > 1:
        self.vntResultsArray[1][10] = None
        self.vntResultsArray[1][11] = None
    
    if RR > 0.0:
      self.vntResultsArray[1][12] = RR
      self.vntResultsArray[1][13] = VarRR ** 0.5
      self.vntResultsArray[1][14] = self.ExpOf(self.LnOf(RR) - (self.varianceMultiplier * VarLnRR ** 0.5))
      self.vntResultsArray[1][15] = self.ExpOf(self.LnOf(RR) + (self.varianceMultiplier * VarLnRR ** 0.5))
      self.vntResultsArray[1][16] = RD * 100.0
      self.vntResultsArray[1][17] = VarRD ** 0.5 * 100.0
      self.vntResultsArray[1][18] = (RD - (self.varianceMultiplier * VarRD ** 0.5)) * 100.0
      self.vntResultsArray[1][19] = (RD + (self.varianceMultiplier * VarRD ** 0.5)) * 100.0
    else:
      self.vntResultsArray[1][12] = None
      self.vntResultsArray[1][13] = None
      self.vntResultsArray[1][14] = None
      self.vntResultsArray[1][15] = None
      self.vntResultsArray[1][16] = None
      self.vntResultsArray[1][17] = None
      self.vntResultsArray[1][18] = None
      self.vntResultsArray[1][19] = None
  
  def AccumVar(self, ah):
    """ Computes the variance for each exposure value
        Parameters:
          ah (int)
        Returns:
          none
    """
    PROC_Name = "clsCTables::AccumVar"
    Pd = self.FirstDom
    Pc = CSCategory()
    P = CSTotal()
    while Pd is not None:
      P = Pd.get_FirstCat()
      while P is not None:
        if ah > 1:
          P.set_VarT(P.get_VarT() + ((ah * P.get_Sumqha2()) - (P.get_Sumqha() ** 2)) / (ah - 1))
        P = P.get_NextCat()
      Pd = Pd.get_NextDom()
  
  def AccumSumq(self):
    """ Computes components of variance
        Parameters:
          none
        Returns:
          none
    """
    PROC_Name = "clsCTables::AccumSumq"
    Pd = self.FirstDom
    Pc = CSCategory()
    P = CSTotal()
    while Pd is not None:
      P = Pd.get_FirstCat()
      while P is not None:
        P.set_Sumqha(P.get_Sumqha() + P.get_qha())
        P.set_Sumqha2(P.get_Sumqha2() + P.get_qha()**2)
        P = P.get_NextCat()
      Pd = Pd.get_NextDom()
  
  def FieldColl(self, p1, s):
    """ Compares the strata or PSU values in two items of data
        Parameters:
          p1 (str, float, or int): a value from the strata or PUS variable
          s  (str, float, or int): a value from the strata or PUS variable
        Returns:
          int indicating greater than, less than, or equal
    """
    PROC_Name = "clsCTables::FieldColl"
    
    ft = None
    i = None
    R = None
    R2 = None
    
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
      if p1.get_FieldEntry() > s:
        FieldColl = 1
      elif p1.get_FieldEntry() < s:
        FieldColl = -1
      else:
        FieldColl = 0
    return FieldColl
  
  def Accumqha(self, P, Pd):
    """ Computes components of variance
        Parameters:
          none
        Returns:
          none
    """
    PROC_Name = "clsCTables::Accumqha"
    
    Qhab = 0.0
    Ptr = CSTotal()
    
    if float(Pd.get_SumW()) > 0:
      Qhab = (1 - (P.get_YE() / Pd.get_SumW())) * (self.GetWeight() / Pd.get_SumW())
      P.set_qha(P.get_qha() + Qhab)
      Ptr = Pd.get_FirstCat()
      while Ptr is not None:
        if Ptr != P:
          Ptr.set_qha(Ptr.get_qha() - (Ptr.get_YE() / Pd.get_SumW()) * (self.GetWeight() / Pd.get_SumW()))
        Ptr = Ptr.get_NextCat()
  
  def QhaInit(self):
    """ Initializes qha at zero for all CSTotal objects
        Parameters:
          none
        Returns:
          none
    """
    PROC_Name = "clsCTables::QhaInit"
    
    Pd = CSDomain()
    Pc = CSCategory()
    P = CSTotal()
    
    Pd = self.FirstDom
    
    while Pd is not None:
      P = Pd.get_FirstCat()
      while P is not None:
        P.set_qha(0)
        P.set_qha2(0)
        P = P.get_NextCat()
      Pd = Pd.get_NextDom()
  
  def SumQInit(self):
    """ Initializes Sumqha at zero for all CSTotal objects
        Parameters:
          none
        Returns:
          none
    """
    PROC_Name = "clsCTables::SumQInit"
    
    Pd = CSDomain()
    Pc = CSCategory()
    P = CSTotal()
    
    Pd = self.FirstDom
    
    while Pd is not None:
      P = Pd.get_FirstCat()
      while P is not None:
        P.set_Sumqha(0)
        P.set_Sumqha2(0)
        P = P.get_NextCat()
      Pd = Pd.get_NextDom()
  
  def VarTInit(self):
    """ Initializes variance at zero for all CSTotal objects
        Parameters:
          none
        Returns:
          none
    """
    PROC_Name = "clsCTables::VarTInit"
    
    Pd = CSDomain()
    Pc = CSCategory()
    P = CSTotal()
    
    Pd = self.FirstDom
    
    while Pd is not None:
      P = Pd.get_FirstCat()
      while P is not None:
        P.set_VarT(0)
        P = P.get_NextCat()
      Pd = Pd.get_NextDom()
  
  def ComputeTot(self, errorMessage):
    """ Sums the total counts
        Parameters:
          errorMessage (str)
        Returns:
          none
    """
    PROC_Name = "clsCTables::ComputeTot"
    
    Pd = CSDomain()
    Pt = CSDomain()
    Pc = CSCategory()
    P = CSTotal()
    Ptr = CSTotal()
    
    Pd = self.AddDom("TOTAL")
    Pd.set_SumW(0)
    Pd.set_N(0)
    self.TotalDom = Pd
    
    Pt = self.FirstDom
    while Pt != Pd:
      Pd.set_SumW(float(Pd.get_SumW()) + float(Pt.get_SumW()))
      Pd.set_N(Pd.get_N() + Pt.get_N())
      Pt = Pt.get_NextDom()
    Ptr = Pd.get_FirstCat()
    Pc = self.FirstCat
    while Pc is not None:
      P = Pc.get_FirstDom()
      while P != Ptr:
        Ptr.set_YE(Ptr.get_YE() + P.get_YE())
        Ptr.set_N(Ptr.get_N() + P.get_N())
        P = P.get_NextDom()
      Pc = Pc.get_NextCat()
      Ptr = Ptr.get_NextCat()
  
  def DesignEffect(self, Pd):
    """ Computes the design effect
        Parameters:
          Pd (CSDomain)
        Returns:
          float
    """
    PROC_Name = "clsCTables::DesignEffect"
    
    denominator = 0
    P = Pd.get_FirstCat()
    if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
      denominator = P.get_YE() / Pd.get_SumW()
    denominator = denominator * (1 - denominator)
    if denominator != 0:
      denominator = denominator / (Pd.get_N() - 1)
    if denominator != 0:
      return P.get_VarT() / denominator
    else:
      return -1.0
  
  def SecondPass(self, errorMessage):
    """ Loops over the analysis dataset computing results
        Parameters:
          errorMessage (str)
        Returns:
          bool
    """
    PROC_Name = "clsCTables::SecondPass"
    P = CSTotal
    p2 = CSTotal
    Pd = CSDomain
    Pc = CSCategory
    Valid = False
    risk = False
    RDRRok = False
    ORok = False
    ah = None
    Rec = True
    N = 0
    NowStrat = ''
    NowPSU = ''
    doma = None
    cate = None
    RR = None
    qhaRR = None
    qha2RR = None
    SumqhaRR = None
    Sumqha2RR = None
    VarRR = 0.0
    ODD = None
    qhaOD = None
    qha2OD = None
    SumqhaOD = None
    Sumqha2OD = None
    VarOD = 0.0
    R1 = None
    R2 = None
    RD = None
    qhaRD = None
    qha2RD = None
    SumqhaRD = None
    Sumqha2RD = None
    VarRD = 0.0
    qhaOR = 0.0
    SumqhaOR = 0.0
    Sumqha2OR = 0.0
    VarOR = 0.0
    qhaLnRR = 0.0
    SumqhaLnRR = 0.0
    Sumqha2LnRR = 0.0
    VarLnRR = 0.0
    a = None
    b = None
    c = None
    d = None
    bContinue = None
    bHadValidPSU = None
    
    SecondPass = False
    self.T22 = False
    self.VarTInit()
    
    if self.LastCat == self.FirstCat.get_NextCat():
      if self.FirstDom.get_NextDom() is not None:
        if self.TotalDom == self.FirstDom.get_NextDom().get_NextDom():
          self.T22 = True
          a = self.FirstCat.get_FirstDom().get_YE()
          d = self.LastCat.get_FirstDom().get_NextDom().get_YE()
          b = self.LastCat.get_FirstDom().get_YE()
          c = self.FirstCat.get_FirstDom().get_NextDom().get_YE()
          ODD = -1
          if a > 0 and b > 0 and c > 0 and d > 0:
            risk = True
            ORok = True
            ODD = (a * d) / (b * c)
          RR = -1
          if a > 0 and c > 0:
            risk = True
            RDRRok = True
            R1 = a / (a + b)
            R2 = c / (c + d)
            RD = R1 - R2
            RR = (a * (c + d)) / (c * (a + b))
    
    while not Valid and Rec:
      Rec = self.GetNextRow()
      Valid = self.ValidCase()
      if Valid:
        if self.strata is not None:
          NowStrat = self.strata.get_FieldEntry()
        if self.psu is not None:
          NowPSU = self.psu.get_FieldEntry()
    while True:
      self.SumQInit()
      SumqhaRR = 0
      Sumqha2RR = 0
      SumqhaOD = 0
      Sumqha2OD = 0
      SumqhaRD = 0
      Sumqha2RD = 0
      SumqhaOR = 0.0
      Sumqha2OR = 0.0
      SumqhaLnRR = 0.0
      Sumqha2LnRR = 0.0
      ah = 0
      while True:
        self.QhaInit()
        qhaRR = 0.0
        qha2RR = 0.0
        qhaOD = 0.0
        qha2OD = 0.0
        qhaRD = 0.0
        qha2RD = 0.0
        qhaOR = 0.0
        qhaLnRR = 0.0
        bHadValidPSU = False
        while True:
          if self.domain is not None:
            doma = self.domain.get_FieldEntry()
            Valid = False
            if self.com:
              if str(doma) != '':
                if self.ValidCase() and (doma == self.domain1 or doma == self.domain2):
                  Valid = True
              else:
                Valid = self.ValidCase()
            else:
              Valid = self.ValidCase()
          else:
            Valid = self.ValidCase()
          
          if Valid:
            bHadValidPSU = True
            cate = self.outcome.get_FieldEntry()
            if self.domain is not None:
              Pd = self.FindDom(doma)
              P = self.FindCat(Pd.get_FirstCat(), cate)
              self.Accumqha(P, Pd)
            
            if self.T22 and risk:
              if Pd == self.FirstDom:
                if P == self.FirstCat.get_FirstDom():
                  if RDRRok:
                    qhaRR = qhaRR + RR * ((1 - a / (a + b)) / a) * self.GetWeight()
                    qhaLnRR = qhaLnRR + (1 / a - 1 / (a + b)) * self.GetWeight()
                    qhaRD = qhaRD + R1 * ((1 - R1) / a) * self.GetWeight()
                  if ORok:
                    qhaOR = qhaOR + ((1 / a) * self.GetWeight())
                    qhaOD = qhaOD + (1 / a) * self.GetWeight()
                else:
                  if RDRRok:
                    qhaRR = qhaRR - RR * (1 / (a + b)) * self.GetWeight()
                    qhaLnRR = qhaLnRR - (1 / (a + b)) * self.GetWeight()
                    qhaRD = qhaRD - (R1 * R1 / a) * self.GetWeight()
                  if ORok:
                    qhaOR = qhaOR - (1 / b) * self.GetWeight()
                    qhaOD = qhaOD - ODD * (1 / b) * self.GetWeight()
              else:
                if P == self.FirstDom.get_NextDom().get_FirstCat():
                  if RDRRok:
                    qhaRR = qhaRR - RR * ((1 - c / (c + d)) / c) * self.GetWeight()
                    qhaLnRR = qhaLnRR + (1 / (c + d) - 1 / c) * self.GetWeight()
                    qhaRD = qhaRD - R2 * ((1 - R2) / c) * self.GetWeight()
                  if ORok:
                    qhaOR = qhaOR - (1 / c) * self.GetWeight()
                    qhaOD = qhaOD - ODD * (1 / c) * self.GetWeight()
                else:
                  if RDRRok:
                    qhaRR = qhaRR + RR * (1 / (c + d)) * self.GetWeight()
                    qhaLnRR = qhaLnRR + (1 / (c + d)) * self.GetWeight()
                    qhaRD = qhaRD + (R2 * R2 / c) * self.GetWeight()
                  if ORok:
                    qhaOR = qhaOR + (1 / d) * self.GetWeight()
                    qhaOD = qhaOD + ODD * (1 / d) * self.GetWeight()
            
            P = self.FindCat(self.TotalDom.get_FirstCat(), cate)
            self.Accumqha(P, self.TotalDom)
          
          Rec = self.GetNextRow()
          
          if self.psu is None:
            bContinue = True
          else:
            if self.psu.get_FieldEntry() != NowPSU:
              bContinue = True
            elif self.strata.get_FieldEntry() != NowStrat:
              bContinue = True
            else:
              bContinue = False
          
          if not Rec or bContinue:
            #print('exiting inner loop')
            break
        if self.psu is not None:
          if not Rec or self.FieldColl(self.psu, NowPSU) > 0:
            NowPSU = self.psu.get_FieldEntry()
          elif self.strata is not None:
            if self.strata.get_FieldEntry() != NowStrat:
              NowPSU = self.psu.get_FieldEntry()
            else:
              return False
          else:
            return False
        
        if bHadValidPSU:
          ah += 1
          self.AccumSumq()
          SumqhaRR = SumqhaRR + qhaRR
          Sumqha2RR = Sumqha2RR + (qhaRR ** 2)
          SumqhaOD = SumqhaOD + qhaOD
          Sumqha2OD = Sumqha2OD + (qhaOD ** 2)
          SumqhaOR = SumqhaOR + qhaOR
          Sumqha2OR = Sumqha2OR + (qhaOR ** 2)
          SumqhaLnRR = SumqhaLnRR + qhaLnRR
          Sumqha2LnRR = Sumqha2LnRR + (qhaLnRR ** 2)
          SumqhaRD = SumqhaRD + qhaRD
          Sumqha2RD = Sumqha2RD + (qhaRD ** 2)
        
        if self.strata is not None:
          if self.strata.get_FieldEntry() != NowStrat:
            bContinue = True
          else:
            bContinue = False
        else:
          bContinue = False
        
        if not Rec or bContinue:
          #print('exiting middle loop')
          break
      
      if self.strata is not None:
        if not Rec or self.FieldColl(self.strata, NowStrat) > 0:
          NowStrat = self.strata.get_FieldEntry()
        else:
          return False
      
      self.AccumVar(ah)
      if ah > 1:
        VarRR = VarRR + (ah * Sumqha2RR - (SumqhaRR ** 2)) / (ah - 1)
        VarOD = VarOD + (ah * Sumqha2OD - (SumqhaOD ** 2)) / (ah - 1)
        VarRD = VarRD + (ah * Sumqha2RD - (SumqhaRD ** 2)) / (ah - 1)
        VarOR = VarOR + (ah * Sumqha2OR - (SumqhaOR ** 2)) / (ah - 1)
        VarLnRR = VarLnRR + (ah * Sumqha2LnRR - (SumqhaLnRR ** 2)) / (ah - 1)
      
      if not Rec:
        #print('exiting outer loop')
        break
    
    self.vntResultsArray = []
    if self.T22 and risk and self.cnOutputLevel > 1:
      for i in range(0, 21):
        self.vntResultsArray.append([None] * 20)
      self.PrintRisk(RR, ODD, RD, VarRR, VarOD, VarRD, VarOR, VarLnRR)
    else:
      for i in range(0, 2):
        self.vntResultsArray.append([None] * 8)
    
    return SecondPass
  
  def ResultsArrayWithTotals(self):
    """ Builds and returns a list of analysis results
        Parameters:
          none
        Returns:
          list
    """
    PROC_Name = "clsCTables::ResultsArrayWithTotals"
    
    vntOutTable = []
    Pc = CSCategory()
    Pd = CSDomain()
    P = CSTotal()
    Total = CSTotal()
    Pct = None
    varDE = None
    Cl = None
    nCat = int()
    nDom = int()
    nRows = int()
    nCurrentRow = int()
    
    if self.outcome is None:
      vntOutTable = [[None]]
    elif self.weight is not None and self.strata is None and self.psu is None:
      vntOutTable = [[None]]
    else:
      Pd = self.FirstDom
      nCat = 0
      P = Pd.get_FirstCat()
      while P is not None:
        nCat += 1
        P = P.get_NextCat()
    nDom = 0
    while Pd is not None:
      nDom += 1
      Pd = Pd.get_NextDom()
    nRows = nCat * nDom
    vntOutTable = []
    for i in range(0, nRows):
      vntOutTable.append([None] * 9)
    Pd = self.FirstDom
    nCurrentRow = 0
    if Pd == self.TotalDom:
      P = Pd.get_FirstCat()
      while P is not None:
        vntOutTable[nCurrentRow][0] = P.get_Category()
        vntOutTable[nCurrentRow][1] = None
        vntOutTable[nCurrentRow][2] = P.get_N()
        if nDom == 1:
          vntOutTable[nCurrentRow][3] = 100.0
          if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
            vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Pd.get_SumW())
          else:
            vntOutTable[nCurrentRow][4] = None
        else:
          if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
            vntOutTable[nCurrentRow][3] = 100 * (P.get_YE() / Pd.get_SumW())
          else:
            vntOutTable[nCurrentRow][3] = None
          Total = P
          while Total.get_NextDom() is not None:
            Total = Total.get_NextDom()
          if P.get_YE() > 0.0 and Total.get_YE() > 0.0:
            vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Total.get_YE())
          else:
            vntOutTable[nCurrentRow][4] = None
        if P.get_VarT() >= 0.0:
          vntOutTable[nCurrentRow][5] = 100 * P.get_VarT() ** 0.5
        else:
          vntOutTable[nCurrentRow][5] = None
        Pct = 0.0
        if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
          Pct = 100 * (P.get_YE() / Pd.get_SumW())
        if Pct > 0.0:
          vntOutTable[nCurrentRow][6] = 100 * ((P.get_YE() / Pd.get_SumW()) + (-self.varianceMultiplier * P.get_VarT() ** 0.5))
          vntOutTable[nCurrentRow][7] = 100 * ((P.get_YE() / Pd.get_SumW()) + (self.varianceMultiplier * P.get_VarT() ** 0.5))
        else:
          vntOutTable[nCurrentRow][6] = None
          vntOutTable[nCurrentRow][7] = None
        varDE = self.DesignEffect(Pd)
        if varDE != -1.0:
          vntOutTable[nCurrentRow][8] = varDE
        else:
          vntOutTable[nCurrentRow][8] = None
        P = P.get_NextCat()
        nCurrentRow += 1
    else:
      while Pd is not None:
        P = Pd.get_FirstCat()
        while P is not None:
          vntOutTable[nCurrentRow][0] = P.get_Category()
          vntOutTable[nCurrentRow][1] = Pd.get_Domain()
          vntOutTable[nCurrentRow][2] = P.get_N()
          if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
            vntOutTable[nCurrentRow][3] = 100 * (P.get_YE() / Pd.get_SumW())
          else:
            vntOutTable[nCurrentRow][3] = None
          Total = P
          while Total.get_NextDom() is not None:
            Total = Total.get_NextDom()
          if P.get_YE() > 0.0 and Total.get_YE() > 0.0:
            vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Total.get_YE())
          else:
            vntOutTable[nCurrentRow][4] = None
          if P.get_VarT() >= 0.0:
            vntOutTable[nCurrentRow][5] = 100 * P.get_VarT() ** .5
          else:
            vntOutTable[nCurrentRow][5] = None
          Pct = 0.0
          if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
            Pct = 100 * (P.get_YE() / Pd.get_SumW())
          if Pct > 0.0:
            vntOutTable[nCurrentRow][6] = 100 * ((P.get_YE() / Pd.get_SumW()) + (-self.varianceMultiplier * P.get_VarT() ** 0.5))
            vntOutTable[nCurrentRow][7] = 100 * ((P.get_YE() / Pd.get_SumW()) + (self.varianceMultiplier * P.get_VarT() ** 0.5))
          else:
            vntOutTable[nCurrentRow][6] = None
            vntOutTable[nCurrentRow][7] = None
          varDE = self.DesignEffect(Pd)
          if varDE != -1.0:
            vntOutTable[nCurrentRow][8] = varDE
          else:
            vntOutTable[nCurrentRow][8] = None
          P = P.get_NextCat()
          nCurrentRow += 1
        Pd = Pd.get_NextDom()
    return vntOutTable
  
  def ResultsArray(self):
    """ Builds and returns a list of analysis results
        Parameters:
          none
        Returns:
          list
    """
    PROC_Name = "clsCTables::ResultsArray"
    
    vntOutTable = []
    vntVarNames = [None] * 10
    vntVarPrompts = [None] * 10
    Pc = CSCategory()
    Pd = CSDomain()
    P = CSTotal()
    Total = CSTotal()
    Pct = None
    varDE = None
    Cl = None
    nCat = int()
    nDom = int()
    nRows = int()
    nCurrentRow = int()
    
    if self.outcome is None:
      vntOutTable = [[None]]
    elif self.weight is not None and self.strata is None and self.psu is None:
      vntOutTable = [[None]]
    else:
      Pd = self.FirstDom
      nCat = 0
      P = Pd.get_FirstCat()
      while P is not None:
        nCat += 1
        P = P.get_NextCat()
      if Pd == self.TotalDom:
        nDom = 1
      else:
        nDom = 0
        while Pd is not None:
          nDom += 1
          Pd = Pd.get_NextDom()
      nRows = nCat * nDom
      vntOutTable = []
      for i in range(0, nRows):
        vntOutTable.append([None] * 12)
      Pd = self.FirstDom
      nCurrentRow = 0
      if Pd == self.TotalDom:
        P = Pd.get_FirstCat()
        while P is not None:
          vntOutTable[nCurrentRow][0] = P.get_Category()
          vntOutTable[nCurrentRow][1] = None
          vntOutTable[nCurrentRow][2] = P.get_N()
          vntOutTable[nCurrentRow][9] = P.get_YE()
          if nDom == 1:
            vntOutTable[nCurrentRow][3] = 100.0
            if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
              vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Pd.get_SumW())
            else:
              vntOutTable[nCurrentRow][4] = None
          else:
            if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
              vntOutTable[nCurrentRow][3] = 100 * (P.get_YE() / Pd.get_SumW())
            else:
              vntOutTable[nCurrentRow][3] = None
            Total = P
            while Total.get_NextDom() is not None:
              Total = Total.get_NextDom()
            if P.get_YE() > 0.0 and Total.get_YE() > 0.0:
              vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Total.get_YE())
            else:
              vntOutTable[nCurrentRow][4] = None
          if P.get_VarT() >= 0.0:
            vntOutTable[nCurrentRow][5] = 100 * P.get_VarT() ** 0.5
          else:
            vntOutTable[nCurrentRow][5] = None
          Pct = 0.0
          if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
            Pct = 100 * (P.get_YE() / Pd.get_SumW())
          if Pct > 0.0:
            vntOutTable[nCurrentRow][6] = 100 * ((P.get_YE() / Pd.get_SumW()) + (-self.varianceMultiplier * P.get_VarT() ** 0.5))
            vntOutTable[nCurrentRow][7] = 100 * ((P.get_YE() / Pd.get_SumW()) + (self.varianceMultiplier * P.get_VarT() ** 0.5))
            phat = P.get_YE() / Pd.get_SumW()
            logitLCL = math.log(phat / (1.0 - phat)) - self.varianceMultiplier * (P.get_VarT() ** 0.5 / (phat * (1.0 - phat)))
            logitUCL = math.log(phat / (1.0 - phat)) + self.varianceMultiplier * (P.get_VarT() ** 0.5 / (phat * (1.0 - phat)))
            vntOutTable[nCurrentRow][10] = 100 * (math.exp(logitLCL) / (1 + math.exp(logitLCL)))
            vntOutTable[nCurrentRow][11] = 100 * (math.exp(logitUCL) / (1 + math.exp(logitUCL)))
          else:
            vntOutTable[nCurrentRow][6] = None
            vntOutTable[nCurrentRow][7] = None
          varDE = self.DesignEffect(Pd)
          if varDE != -1.0:
            vntOutTable[nCurrentRow][8] = varDE
          else:
            vntOutTable[nCurrentRow][8] = None
          P = P.get_NextCat()
          nCurrentRow += 1
      else:
        while Pd is not None:
          P = Pd.get_FirstCat()
          while P is not None:
            vntOutTable[nCurrentRow][0] = P.get_Category()
            vntOutTable[nCurrentRow][1] = Pd.get_Domain()
            vntOutTable[nCurrentRow][2] = P.get_N()
            if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
              vntOutTable[nCurrentRow][3] = 100 * (P.get_YE() / Pd.get_SumW())
            else:
              vntOutTable[nCurrentRow][3] = None
            Total = P
            while Total.get_NextDom() is not None:
              Total = Total.get_NextDom()
            if P.get_YE() > 0.0 and Total.get_YE() > 0.0:
              vntOutTable[nCurrentRow][4] = 100 * (P.get_YE() / Total.get_YE())
            else:
              vntOutTable[nCurrentRow][4] = None
            if P.get_VarT() >= 0.0:
              vntOutTable[nCurrentRow][5] = 100 * P.get_VarT() ** .5
            else:
              vntOutTable[nCurrentRow][5] = None
            Pct = 0.0
            if P.get_YE() > 0.0 and Pd.get_SumW() > 0.0:
              Pct = 100 * (P.get_YE() / Pd.get_SumW())
            if Pct > 0.0:
              vntOutTable[nCurrentRow][6] = 100 * ((P.get_YE() / Pd.get_SumW()) + (-self.varianceMultiplier * P.get_VarT() ** 0.5))
              vntOutTable[nCurrentRow][7] = 100 * ((P.get_YE() / Pd.get_SumW()) + (self.varianceMultiplier * P.get_VarT() ** 0.5))
            else:
              vntOutTable[nCurrentRow][6] = None
              vntOutTable[nCurrentRow][7] = None
            varDE = self.DesignEffect(Pd)
            if varDE != -1.0:
              vntOutTable[nCurrentRow][8] = varDE
            else:
              vntOutTable[nCurrentRow][8] = None
            P = P.get_NextCat()
            nCurrentRow += 1
          Pd = Pd.get_NextDom()
    if self.outcome is not None:
      vntVarNames[1] = self.outcome.get_FieldLabel()
    else:
      vntVarNames[1] = None
    if self.domain is not None:
      vntVarNames[2] = self.domain.get_FieldLabel()
    else:
      vntVarNames[2] = None
    vntVarNames[3] = 'Count'
    vntVarNames[4] = 'RowPct'
    vntVarNames[5] = "ColPct"
    vntVarNames[6] = 'StdErr'
    vntVarNames[7] = "LCL"
    vntVarNames[8] = 'UCL'
    vntVarNames[9] = 'DesignEff'
    if self.outcome is not None:
      vntVarPrompts[1] = self.outcome.get_FieldLabel()
    else:
      vntVarPrompts[1] = None
    if self.domain is not None:
      vntVarPrompts[2] = self.domain.get_FieldLabel()
    else:
      vntVarPrompts[2] = None
    vntVarPrompts[3] = 'Count'
    vntVarPrompts[4] = 'Row Percent'
    vntVarPrompts[5] = "Column Percent"
    vntVarPrompts[6] = 'Standard Error Percent'
    vntVarPrompts[7] = "Lower Confidence Limit"
    vntVarPrompts[8] = 'Upper Confidence Limit'
    vntVarPrompts[9] = 'Design Effect'
    self.vntResultsArray[0][0] = 'Errors'
    self.vntResultsArray[1][0] = ''
    self.vntResultsArray[0][1] = 'VarNames'
    self.vntResultsArray[1][1] = vntVarNames
    self.vntResultsArray[0][2] = 'VarPrompts'
    self.vntResultsArray[1][2] = vntVarPrompts
    self.vntResultsArray[0][3] = 'OutTable'
    self.vntResultsArray[1][3] = vntOutTable
    self.vntResultsArray[0][4] = 'Excluded'
    self.vntResultsArray[1][4] = self.Mis
    self.vntResultsArray[0][5] = 'Weight Variable'
    if self.weight is not None:
      self.vntResultsArray[1][5] = self.weight.get_FieldLabel()
    else:
      self.vntResultsArray[1][5] = "None"
    self.vntResultsArray[0][6] = 'PSU Variable'
    if self.psu is not None:
      self.vntResultsArray[1][6] = self.psu.get_FieldLabel()
    else:
      self.vntResultsArray[1][6] = "None"
    self.vntResultsArray[0][7] = 'Stratification Variable'
    if self.strata is not None:
      self.vntResultsArray[1][7] = self.strata.get_FieldLabel()
    else:
      self.vntResultsArray[1][7] = "None"
    
    ResultsArray = self.vntResultsArray
    return ResultsArray

  def ComplexSampleFrequencies(self, inputVariableList, dataTable):
    """ Executes the supporting functions to run the CS Frequencies analysis
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
          dataTable (list(dict)): The analysis dataset
        Returns:
          csFrequencyResults (CSFrequencyResults): This object contains a Rows property.
          It is a list of CSRow objects, which have properties: Value, Domain, Count,
          WeightedCount, RowPercent, ColPercent, SE, LCL, UCL, DesignEffect, LogitLCL,
          and LogitUCL. These are the displayed output of the analysis. There is a CSRow
          for TOTAL and one for each value of the frequency variable, if present.
    """
    self.currentTable = dataTable
    self.CreateSettings(inputVariableList)
    
    errorMessage = ''
    
    output = []
    
    csFrequencyResults = CSFrequencyResults()
    csFrequencyResults.set_ErrorMessage('')
    csFrequencyResults.set_Rows([])
    
    if self.Init(errorMessage) == False:
      csFrequencyResults.set_ErrorMessage(errorMessage)
      return csFrequencyResults
    
    self.Mis = 0
    result = None
    
    if self.com:
      errorMessage = ''
      result = self.FirstPassCom(errorMessage)
      if errorMessage is not None and len(errorMessage) > 0:
        csFrequencyResults.set_ErrorMessage(errorMessage)
        return csFrequencyResults
    else:
      errorMessage = ''
      result = self.FirstPass(errorMessage)
      if errorMessage is not None and len(errorMessage) > 0:
        csFrequencyResults.set_ErrorMessage(errorMessage)
        return csFrequencyResults
    
    if result:
      if self.FirstDom.get_NextDom() == self.LastDom:
        self.com = True
        self.domain1 = self.FirstDom.get_Domain()
        self.domain2 = self.LastDom.get_Domain()
      
      if self.domain is not None:
        self.ComputeTot(errorMessage)
        if errorMessage is not None and len(errorMessage) > 0:
          csFrequencyResults.set_ErrorMessage(errorMessage)
          return csFrequencyResults
      
      self.ResetReader()
      
      result = self.SecondPass(errorMessage)
      
      if errorMessage is not None and len(errorMessage) > 0:
        csFrequencyResults.set_ErrorMessage(errorMessage)
        return csFrequencyResults
      
    vntResultsArray = self.ResultsArray()
    if vntResultsArray is None:
      return csFrequencyResults
    vntOutTable = vntResultsArray[1][3]
    for i in range(0, len(vntOutTable)):
      fRow = CSRow()
      fRow.set_Value(vntOutTable[i][0])
      fRow.set_Domain(vntOutTable[i][1])
      fRow.set_Count(vntOutTable[i][2])
      fRow.set_WeightedCount(vntOutTable[i][9])
      fRow.set_RowPercent(vntOutTable[i][3])
      fRow.set_ColPercent(vntOutTable[i][4])
      fRow.set_SE(vntOutTable[i][5])
      fRow.set_LCL(vntOutTable[i][6])
      fRow.set_UCL(vntOutTable[i][7])
      fRow.set_DesignEffect(vntOutTable[i][8])
      fRow.set_LogitLCL(vntOutTable[i][10])
      fRow.set_LogitUCL(vntOutTable[i][11])
      csFrequencyResults.get_Rows().append(fRow)
    
    return csFrequencyResults

  def ComplexSampleTables(self, inputVariableList, dataTable):
    """ Executes the supporting functions to run the CS Tables analysis
        Parameters:
          inputVariableList (dict): Indicates the names of the analysis variables
          dataTable (list(dict)): The analysis dataset
        Returns:
          tablesResults (CSTablesResults): This object contains a Rows property.
          It is a list of CSRow objects, which have properties: Value, Domain, Count,
          WeightedCount, RowPercent, ColPercent, SE, LCL, UCL, DesignEffect, LogitLCL, 
          and LogitUCL. These are the displayed output of the analysis. There is a CSRow
          for TOTAL and one for each value of the frequency variable, if present.
    """
    self.currentTable = dataTable
    
    self.CreateSettings(inputVariableList)
    
    errorMessage = ''
    
    output = []
    
    self.columnPercents = []
    tablesResults = CSTablesResults()
    tablesResults.set_ErrorMessage('')
    tablesResults.set_Rows([])
    
    if self.Init(errorMessage) == False:
      tablesResults.set_ErrorMessage(errorMessage)
      return tablesResults
    
    self.Mis = 0
    self.FirstDom = None
    self.LastDom = None
    self.TotalDom = None
    self.FirstCat = None
    self.LastCat = None
    self.Pdm = None
    result = None
    
    if self.com:
      errorMessage = ''
      result = self.FirstPassCom(errorMessage)
      if errorMessage is not None and len(errorMessage) > 0:
        tablesResults.set_ErrorMessage(errorMessage)
        return tablesResults
    else:
      errorMessage = ''
      result = self.FirstPass(errorMessage)
      if errorMessage is not None and len(errorMessage) > 0:
        tablesResults.set_ErrorMessage(errorMessage)
        return tablesResults
    
    if result:
      if self.FirstDom.get_NextDom() == self.LastDom:
        self.com = True
        self.domain1 = self.FirstDom.get_Domain()
        self.domain2 = self.LastDom.get_Domain()
      
      if self.domain is not None:
        self.ComputeTot(errorMessage)
        if errorMessage is not None and len(errorMessage) > 0:
          tablesResults.set_ErrorMessage(errorMessage)
          return tablesResults
      
      self.ResetReader()
      
      result = self.SecondPass(errorMessage) # TODO: code SecondPass (line 3450)
      if errorMessage is not None and len(errorMessage) > 0:
        tablesResults.set_ErrorMessage(errorMessage)
        return tablesResults
      
    vntResultsArray = self.ResultsArray()
    sw1 = self.FirstDom.get_SumW()
    sw2 = self.FirstDom.get_NextDom().get_SumW()
    vntOutTable = self.ResultsArrayWithTotals()
    
    if len(vntResultsArray[1]) > 9 and len(vntResultsArray[1]) < 22:
      tablesResults.set_StandardErrorOR(vntResultsArray[1][9])
      tablesResults.set_StandardErrorRR(vntResultsArray[1][13])
      tablesResults.set_StandardErrorRD(vntResultsArray[1][17])
      tablesResults.set_OddsRatio(vntResultsArray[1][8])
      tablesResults.set_RiskRatio(vntResultsArray[1][12])
      tablesResults.set_RiskDifference(vntResultsArray[1][16])
      tablesResults.set_LCLOR(vntResultsArray[1][10])
      tablesResults.set_LCLRR(vntResultsArray[1][14])
      tablesResults.set_LCLRD(vntResultsArray[1][18])
      tablesResults.set_UCLOR(vntResultsArray[1][11])
      tablesResults.set_UCLRR(vntResultsArray[1][15])
      tablesResults.set_UCLRD(vntResultsArray[1][19])
      for v in vntResultsArray[1][3]:
        tRow = CSTablesRow()
        tRow.set_Outcome(v[0])
        tRow.set_Exposure(v[1])
        tRow.set_Count(v[2])
        tRow.set_RowPercent(v[3])
        tRow.set_ColPercent(v[4])
        tRow.set_SE(v[5])
        tRow.set_LCL(v[6])
        tRow.set_UCL(v[7])
        tRow.set_DesignEffect(v[8])
        tablesResults.get_Rows().append(tRow)
    else:
      tablesResults.set_StandardErrorOR(None)
      tablesResults.set_StandardErrorRR(None)
      tablesResults.set_StandardErrorRD(None)
      tablesResults.set_OddsRatio(None)
      tablesResults.set_RiskRatio(None)
      tablesResults.set_RiskDifference(None)
      tablesResults.set_LCLOR(None)
      tablesResults.set_LCLRR(None)
      tablesResults.set_LCLRD(None)
      tablesResults.set_UCLOR(None)
      tablesResults.set_UCLRR(None)
      tablesResults.set_UCLRD(None)

      currRow = str(vntOutTable[i][1])
    
    return tablesResults
