""" This file defines several supporting classes for the
    Epi Info Complex Sample Analysis Routines
"""
class CSUtilities:
  def __init__(self):
    self._Exists = True
  def get_Exists(self):
    return self._Exists
  def set_Exists(self, v):
    self._Exists = v

class CSMeansTotal:
  """ Class to hold intermediate statistics during the
      CS Means Analysis
  """
  def __init__(self, dom=''):
    self._Domain = dom 
    self._YE = 0.0
    self._NextTotal = None # CSMeansTotal type
    self._SumW = 0.0
    self._N = 0
    self._Min = 0.0
    self._Max = 0.0
    self._qha = None
    self._qha2 = None
    self._Sumqha = None
    self._Sumqha2 = None
    self._VarT = None
  def get_Domain(self):
    return self._Domain
  def set_Domain(self, v):
    self._Domain = v
  def get_YE(self):
    return self._YE
  def set_YE(self, v):
    self._YE = v
  def get_NextTotal(self):
    return self._NextTotal
  def set_NextTotal(self, v):
    self._NextTotal = v
  def get_SumW(self):
    return self._SumW
  def set_SumW(self, v):
    self._SumW = v
  def get_N(self):
    return self._N
  def set_N(self, v):
    self._N = v
  def get_Min(self):
    return self._Min
  def set_Min(self, v):
    self._Min = v
  def get_Max(self):
    return self._Max
  def set_Max(self, v):
    self._Max = v
  def get_qha(self):
    return self._qha
  def set_qha(self, v):
    self._qha = v
  def get_qha2(self):
    return self._qha2
  def set_qha2(self, v):
    self._qha2 = v
  def get_Sumqha(self):
    return self._Sumqha
  def set_Sumqha(self, v):
    self._Sumqha = v
  def get_Sumqha2(self):
    return self._Sumqha2
  def set_Sumqha2(self, v):
    self._Sumqha2 = v
  def get_VarT(self):
    return self._VarT
  def set_VarT(self, v):
    self._VarT = v

class CSRow:
  """ Class for CS Frequency results
  """
  def __init__(self):
    self._Value = ''
    self._Domain = ''
    self._Count = None
    self._WeightedCount = None
    self._RowPercent = None
    self._ColPercent = None
    self._SE = None
    self._LCL = None
    self._UCL = None
    self._DesignEffect = None
    self._LogitLCL = None
    self._LogitUCL = None
  def get_Value(self):
    return self._Value
  def set_Value(self, v):
    self._Value = v
  def get_Domain(self):
    return self._Domain
  def set_Domain(self, v):
    self._Domain = v
  def get_Count(self):
    return self._Count
  def set_Count(self, v):
    self._Count = v
  def get_WeightedCount(self):
    return self._WeightedCount
  def set_WeightedCount(self, v):
    self._WeightedCount = v
  def get_RowPercent(self):
    return self._RowPercent
  def set_RowPercent(self, v):
    self._RowPercent = v
  def get_ColPercent(self):
    return self._ColPercent
  def set_ColPercent(self, v):
    self._ColPercent = v
  def get_SE(self):
    return self._SE
  def set_SE(self, v):
    self._SE = v
  def get_LCL(self):
    return self._LCL
  def set_LCL(self, v):
    self._LCL = v
  def get_UCL(self):
    return self._UCL
  def set_UCL(self, v):
    self._UCL = v
  def get_DesignEffect(self):
    return self._DesignEffect

class CSTablesRow:
  """ Class for CS Frequency results
  """
  def __init__(self):
    self._Outcome = ''
    self._Exposure = ''
    self._Count = None
    self._RowPercent = None
    self._ColPercent = None
    self._SE = None
    self._LCL = None
    self._UCL = None
    self._DesignEffect = None
    self._LogitLCL = None
    self._LogitUCL = None
  def get_Outcome(self):
    return self._Outcome
  def set_Outcome(self, v):
    self._Outcome = v
  def get_Exposure(self):
    return self._Exposure
  def set_Exposure(self, v):
    self._Exposure = v
  def get_Count(self):
    return self._Count
  def set_Count(self, v):
    self._Count = v
  def get_RowPercent(self):
    return self._RowPercent
  def set_RowPercent(self, v):
    self._RowPercent = v
  def get_ColPercent(self):
    return self._ColPercent
  def set_ColPercent(self, v):
    self._ColPercent = v
  def get_SE(self):
    return self._SE
  def set_SE(self, v):
    self._SE = v
  def get_LCL(self):
    return self._LCL
  def set_LCL(self, v):
    self._LCL = v
  def get_UCL(self):
    return self._UCL
  def set_UCL(self, v):
    self._UCL = v
  def get_DesignEffect(self):
    return self._DesignEffect
  def set_DesignEffect(self, v):
    self._DesignEffect = v
  def get_LogitLCL(self):
    return self._LogitLCL
  def set_LogitLCL(self, v):
    self._LogitLCL = v
  def get_LogitUCL(self):
    return self._LogitUCL
  def set_LogitUCL(self, v):
    self._LogitUCL = v

class CSFrequencyResults:
  """ Class for CS Frequency results
      which includs a list of CSRow objects
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

class CSTablesResults:
  """ Class for CS Tables results
  """
  def __init__(self):
    self._Rows = []
    self._OddsRatio = None
    self._StandardErrorOR = None
    self._LCLOR = None
    self._UCLOR = None
    self._RiskRatio = None
    self._StandardErrorRR = None
    self._LCLRR = None
    self._UCLRR = None
    self._RiskDifference = None
    self._StandardErrorRD = None
    self._LCLRD = None
    self._UCLRD = None
    self._ErrorMessage = ''
  def get_Rows(self):
    return self._Rows
  def set_Rows(self, v):
    self._Rows = v
  def get_OddsRatio(self):
    return self._OddsRatio
  def set_OddsRatio(self, v):
    self._OddsRatio = v
  def get_StandardErrorOR(self):
    return self._StandardErrorOR
  def set_StandardErrorOR(self, v):
    self._StandardErrorOR = v
  def get_LCLOR(self):
    return self._LCLOR
  def set_LCLOR(self, v):
    self._LCLOR = v
  def get_UCLOR(self):
    return self._UCLOR
  def set_UCLOR(self, v):
    self._UCLOR = v
  def get_RiskRatio(self):
    return self._RiskRatio
  def set_RiskRatio(self, v):
    self._RiskRatio = v
  def get_StandardErrorRR(self):
    return self._StandardErrorRR
  def set_StandardErrorRR(self, v):
    self._StandardErrorRR = v
  def get_LCLRR(self):
    return self._LCLRR
  def set_LCLRR(self, v):
    self._LCLRR = v
  def get_UCLRR(self):
    return self._UCLRR
  def set_UCLRR(self, v):
    self._UCLRR = v
  def get_RiskDifference(self):
    return self._RiskDifference
  def set_RiskDifference(self, v):
    self._RiskDifference = v
  def get_StandardErrorRD(self):
    return self._StandardErrorRD
  def set_StandardErrorRD(self, v):
    self._StandardErrorRD = v
  def get_LCLRD(self):
    return self._LCLRD
  def set_LCLRD(self, v):
    self._LCLRD = v
  def get_UCLRD(self):
    return self._UCLRD
  def set_UCLRD(self, v):
    self._UCLRD = v
  def get_ErrorMessage(self):
    return self._ErrorMessage
  def set_ErrorMessage(self, v):
    self._ErrorMessage = v

class CSField:
  """ Class to hold values from rows of the
      analysis dataset
  """
  def __init__(self):
    self._cnFieldLen = None
    self._csFieldEntry = None
    self._FieldLabel = None
    self._FieldEntry = None
    self._FieldReal = None
    self._cbMissing = None
    self._cenumFieldType = None
  def get_cnFieldLen(self):
    return self._cnFieldLen
  def set_cnFieldLen(self, v):
    self._cnFieldLen = v
  def get_csFieldEntry(self):
    return self._csFieldEntry
  def set_csFieldEntry(self, v):
    self._csFieldEntry = v
  def get_FieldReal(self):
    return self._FieldReal
  def set_FieldReal(self, v):
    self._FieldReal = v
  def get_FieldLabel(self):
    return self._FieldLabel
  def set_FieldLabel(self, v):
    self._FieldLabel = v
  def get_FieldEntry(self):
    return self._FieldEntry
  def set_FieldEntry(self, v):
    self._FieldEntry = v
    self._csFieldEntry = v
  def get_cbMissing(self):
    return self._cbMissing
  def set_cbMissing(self, v):
    self._cbMissing = v
  def get_cenumFieldType(self):
    return self._cenumFieldType
  def set_cenumFieldType(self, v):
    self._cenumFieldType = v
  def get_FieldInt(self):
    if str(self._csFieldEntry).isnumeric() and float(self._csFieldEntry) % 1 == 0:
      return int(self._csFieldEntry)
    else:
      return 0.0

class CSDomain:
  """ Class to hold values from rows of the
      analysis dataset
  """
  def __init__(self):
    self._Domain = None
    self._csFieldEntry = None
    self._SumW = None
    self._N = None
    self._NextDom = None
    self._FirstCat = None
  def get_Domain(self):
    return self._Domain
  def set_Domain(self, v):
    self._Domain = v
  def get_csFieldEntry(self):
    return self._csFieldEntry
  def set_csFieldEntry(self, v):
    self._csFieldEntry = v
  def get_SumW(self):
    return self._SumW
  def set_SumW(self, v):
    self._SumW = v
  def get_N(self):
    return self._N
  def set_N(self, v):
    self._N = v
  def get_NextDom(self):
    return self._NextDom
  def set_NextDom(self, v):
    self._NextDom = v
  def get_FirstCat(self):
    return self._FirstCat
  def set_FirstCat(self, v):
    self._FirstCat = v

class CSTotal:
  """ Class to hold accumulated values
      as the analysis continues
  """
  def __init__(self):
    self._Domain = None
    self._Category = None
    self._YE = None
    self._N = None
    self._qha = None
    self._qha2 = None
    self._Sumqha = None
    self._Sumqha2 = None
    self._VarT = None
    self._NextDom = None
    self._NextCat = None
  def get_Domain(self):
    return self._Domain
  def set_Domain(self, v):
    self._Domain = v
  def get_Category(self):
    return self._Category
  def set_Category(self, v):
    self._Category = v
  def get_YE(self):
    return self._YE
  def set_YE(self, v):
    self._YE = v
  def get_N(self):
    return self._N
  def set_N(self, v):
    self._N = v
  def get_qha(self):
    return self._qha
  def set_qha(self, v):
    self._qha = v
  def get_qha2(self):
    return self._qha2
  def set_qha2(self, v):
    self._qha2 = v
  def get_Sumqha(self):
    return self._Sumqha
  def set_Sumqha(self, v):
    self._Sumqha = v
  def get_Sumqha2(self):
    return self._Sumqha2
  def set_Sumqha2(self, v):
    self._Sumqha2 = v
  def get_VarT(self):
    return self._VarT
  def set_VarT(self, v):
    self._VarT = v
  def get_NextDom(self):
    return self._NextDom
  def set_NextDom(self, v):
    self._NextDom = v
  def get_NextCat(self):
    return self._NextCat
  def set_NextCat(self, v):
    self._NextCat = v

class CSCategory:
  """ Class to hold current values of PSU
      and strata
  """
  def __init__(self):
    self._Category = None
    self._NextCat = None
    self._FirstDom = None
  def get_Category(self):
    return self._Category
  def set_Category(self, v):
    self._Category = v
  def get_NextCat(self):
    return self._NextCat
  def set_NextCat(self, v):
    self._NextCat = v
  def get_FirstDom(self):
    return self._FirstDom
  def set_FirstDom(self, v):
    self._FirstDom = v

class TablesRow:
  """ Class to hold CS Tables row results
  """
  def __init__(self):
    self._Cells = None
    self._RowColPercent = None
  def get_Cells(self):
    return self._Cells
  def set_Cells(self, v):
    self._Cells = v
  def get_RowColPercent(self):
    return self._RowColPercent
  def set_RowColPercent(self, v):
    self._RowColPercent = v

class MeansRow:
  """ Class to hold CS Means row results
  """
  def __init__(self):
    self._Label = None
    self._Count = None
    self._Mean = None
    self._SteErr = None
    self._LCL = None
    self._UCL = None
    self._Min = None
    self._Max = None
  def get_Label(self):
    return self._Label
  def set_Label(self, v):
    self._Label = v
  def get_Count(self):
    return self._Count
  def set_Count(self, v):
    self._Count = v
  def get_Mean(self):
    return self._Mean
  def set_Mean(self, v):
    self._Mean = v
  def get_StdErr(self):
    return self._StdErr
  def set_StdErr(self, v):
    self._StdErr = v
  def get_LCL(self):
    return self._LCL
  def set_LCL(self, v):
    self._LCL = v
  def get_UCL(self):
    return self._UCL
  def set_UCL(self, v):
    self._UCL = v
  def get_Min(self):
    return self._Min
  def set_Min(self, v):
    self._Min = v
  def get_Max(self):
    return self._Max
  def set_Max(self, v):
    self._Max = v

class CSMeansResults:
  """ Class to hold a list of CS Means row results
  """
  def __init__(self):
    self._Rows = []
    self._ErrorMessage = None
  def get_Rows(self):
    return self._Rows
  def set_Rows(self, v):
    self._Rows = v
  def get_ErrorMessage(self):
    return self._ErrorMessage
  def set_ErrorMessage(self, v):
    self._ErrorMessage = v
