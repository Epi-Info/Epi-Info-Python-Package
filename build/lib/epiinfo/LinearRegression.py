from scipy.stats import t as tdist
import math
import time
import itertools
from .randata import randata
from .RegressionUtilities import *

class LinearRegression:
  """ 
      
      It returns a list of dictionaries.

     Author: John Copeland
  """

  def __init__(self):
    self.mstrWeightVar = None
    self.mboolIntercept = True

  def CreateWorkingTable(self, depvar, indepvars, localLOD):
    workingDict = {}
    unsortedList = []
    for r in localLOD:
        rowlist = [r[depvar], r[indepvars[0]]]
        for v in indepvars[1:]:
            rowlist.append(r[v])
        rowlist.append(1)
        unsortedList.append(rowlist)
    sortedList = sorted(unsortedList, key=lambda x: (x[0:]))
    return sortedList

    for r in localLOD:
        if float(r[pNumericVariable]) in workingDict:
            workingDict[float(r[pNumericVariable])] += 1.0
            continue
        workingDict[float(r[pNumericVariable])] = 1.0                                
    for k in workingDict:
        unsortedList.append([float(k), workingDict[float(k)]])
    sortedList = sorted(unsortedList, key=lambda x: (x[0], x[1]))
    return sortedList

  def PfromT(self, T, df):
    PfromT = 0.0
    P = 0.0
    c = 0.0
    s = 0.0
    b = 0.0
    a = 0.0
    G1 = 0.3183098862
    MaxInt = 1000
    if T < 0:
        T *= -1.0
    if df < MaxInt:
        ddf = df
        a = T / ddf ** 0.5
        b = ddf / (ddf + T ** 2)
        i = ddf % 2
        s = 1
        c = 1
        F = 2 + i
        while F <= ddf - 2:
            c = c * b * (F - 1) / F
            s += c
            F += 2
        if i <= 0:
            P = 0.5 - a * b ** 0.5 * s / 2
        else:
            P = 0.5 - (a * b * s + math.atan(a)) * G1
        if P < 0:
            P = 0
        if P > 1:
            P = 1
        PfromT = P
    else:
        PfromT = self.PfromZ(T)
    return PfromT

  def PfromF(self, F, df1, df2):
    PfromF = 0.0
    xx = 0.0
    pp = 0.0
    qq = 0.0
    index = False
    ACU = 0.000000001
    if F == 0.0:
        return float('nan')
    if F < 0 or df1 < 1 or df2 < 1:
        return float('nan')
    if df1 == 1:
        return self.PfromT(F ** 0.5, df2) * 2
    x = df1 * F / (df2 + df1 * F)
    P = df1 / 2
    q = df2 / 2
    psq = P + q
    cx = 1 - x
    if P >= x * psq:
        xx  = x
        pp = P
        qq = q
        index = False
    else:
        xx = cd
        cx = x
        pp = q
        qq = P
        index = True
    term = 1
    ai = 1
    b = 1
    ns = qq + cx * psq
    rx = xx / cx
    term1 = 1
    temp = qq - ai
    if ns == 0:
        rx = xx
    term = term / (pp + ai) * temp * rx
    absterm = term
    if term < 0:
        absterm = -term
    while absterm <= term1:
        b += term
        temp = term
        if term < 0:
            temp = -term
        term1 = temp
        if temp > ACU or temp > ACU * b:
            ai += 1
            ns -= 1
            if ns >= 0:
                temp = qq - ai
                if ns == 0:
                    rx = xx
                term = term / (pp + ai) * temp * rx
                absterm = term
                if term < 0:
                    absterm = -term
                continue
            temp = psq
            psq += 1
            term = term / (pp + ai) * temp * rx
            absterm = term
            if term < 0:
                absterm = -term
    beta = self.algama(P) + self.algama(q) - self.algama(P + q)
    temp = (pp * math.log(xx) + (qq - 1) * math.log(cx) - beta) - math.log(pp)
    if temp > -70:
        b *= math.exp(temp)
    else:
        b *= 0
    if index:
        b = 1 - b
    PfromF = 1 - b
    return PfromF

  def doRegression(self, inputVariableList, dataTable):
    lod = dataTable
    res = []
    DataArray = self.CreateWorkingTable(inputVariableList['dependvar'],
                                           inputVariableList['exposureVariables'],
                                           lod)
    NumRows = len(DataArray)
    NumColumns = len(DataArray[0])
    lintWRows = 0
    lintweight = 0
    ldblweight = False
    lintrowCount = 0
    ldblMagic = 0.0
    if self.mstrWeightVar is not None:
        lintweight = 1
    for i in range(NumRows):
        if lintweight == 1:
            lintWRows += 1
            lintrowCount += DataArray[i][1]
        else:
            lintWRows += 1
            lintrowCount = lintWRows

    fvalue, covb, probf, stdb, y, x, xx, invxx, xy, tx, B, yhat, resid, indx = [], [], [], [], [], [], [], [], [], [], [], [], [], []
    d = 0.0
    mdblToler = 0.00001
    for yhatzero in range(NumRows):
        yhat.append([0])
    
    dlist = [0.0]
    for i in range(NumRows):
        if lintweight == 1:
            y_k = [DataArray[i][0] * DataArray[i][1] ** 0.5]
            y.append(y_k)
            x_k_j = [0]
            for j in range(1 + lintweight, NumColumns):
                x_k_j_num = DataArray[i][j] * DataArray[i][1] ** 0.5
                x_k_j.append(x_k_j_num)
            x.append(x_k_j)
        else:
            y_k = [DataArray[i][0]]
            y.append(y_k)
            x_k_j = []
            for j in range(0 + lintweight, NumColumns - 1):
                x_k_j_num = DataArray[i][j + 1]
                x_k_j.append(x_k_j_num)
            x.append(x_k_j)
    Matrix1 = EIMatrix()
    for r in range(len(x[0])):
        tx.append([0] * len(x))
        xx.append([0] * len(x[0]))
        xy.append([0] * 1)
        B.append([0] * 1)
    Matrix1.trans(len(x), len(x[0]), x, tx)
    Matrix1.mul(len(tx), len(tx[0]), len(x), len(x[0]), tx, x, xx)
    Matrix1.mul(len(tx), len(tx[0]), len(y), len(y[0]), tx, y, xy)
    for xxr in xx:
        invxxr = []
        for v in xxr:
            invxxr.append(v)
        invxx.append(invxxr)
    indx = [0] * len(x[0])
    Matrix1.ludcmp(invxx, len(invxx[0]), indx, dlist)
    d = 1
    for i in range(len(invxx[0])):
        d *= invxx[i][i]
        if d == 0:
            return [["ERROR"], ['Colinear Data']]
    if abs(d) < mdblToler:
        return [['ERROR'], ['Matrix Toleracne Exceeded']]
    Matrix1.inv(xx, invxx)
    Matrix1.mul(len(invxx), len(invxx[0]), len(xy), len(xy[0]), invxx, xy, B)
    Matrix1.mul(len(x), len(x[0]), len(B), len(B[0]), x, B, yhat)
    sse = 0
    meanY = 0
    for i in range(lintWRows):
        if lintweight > 0:
            ldblMagnc = DataArray[i][1] ** 0.5
            resid.append(y[i][0] - yhat[i][0])
            meanY += y[i][0] * ldblMagic
            sse += resid[i] ** 2
        else:
            resid.append(y[i][0] - yhat[i][0])
            meanY += y[i][0]
            sse += resid[i] ** 2
    meanY /= lintrowCount
    ssy = 0
    for i in range(lintWRows):
        if lintweight > 0:
            if DataArray[i][1] != 0:
                ssy += ((y[i][0] * DataArray[i][1] ** -0.5 - meanY)) ** 2 * DataArray[i][1]
        else:
            ssy += (y[i][0] - meanY * abs(int(self.mboolIntercept))) ** 2
    r2 = (ssy - sse) / ssy
    df = lintrowCount - (NumColumns - lintweight + int(self.mboolIntercept) * 0 - 1)
    ra2 = 1 - int(lintrowCount + int(self.mboolIntercept)) * sse / (df * ssy)
    mse = sse / df
    tScore = 1.0
    while self.PfromT(tScore, df) > 0.025:
        tScore += 0.000001
    ftest = (ssy - sse) / int(NumColumns - lintweight - 1 - int(self.mboolIntercept))
    ftest /= mse
    rsme = mse ** 0.5
    for i in range(NumColumns - lintweight - 1):
        covbi = []
        for j in range(NumColumns - lintweight - 1):
            covbi.append(invxx[i][j] * mse)
        covb.append(covbi)
        stdb.append((abs(covb[i][i])) ** 0.5)
        fvalue.append((B[i][0] / stdb[i]) ** 2)
        probf.append(self.PfromF(fvalue[i], 1, df))
    
    for i in range(len(B)):
        association = {}
        if i == len(B) - 1:
            association['variable'] = 'CONSTANT'
        else:
            association['variable'] = inputVariableList['exposureVariables'][i]
        association['beta'] = B[i][0]
        association['lcl'] = B[i][0] - tScore * stdb[i]
        association['ucl'] = B[i][0] + tScore * stdb[i]
        association['stderror'] = stdb[i]
        association['ftest'] = fvalue[i]
        association['pvalue'] = probf[i]
        res.append(association)
    statsdict = {'r2' : r2}
    statsdict['regressionDF'] = lintrowCount - int(self.mboolIntercept) - df
    statsdict['sumOfSquares'] = ssy - sse
    statsdict['meanSquare'] = (ssy - sse) / (NumColumns - lintweight - 1 - int(self.mboolIntercept))
    statsdict['fStatistic'] = ftest
    statsdict['residualsDF'] = df
    statsdict['residualsSS'] = sse
    statsdict['residualsMS'] = mse
    statsdict['totalDF'] = lintrowCount - int(self.mboolIntercept)
    statsdict['totalSS'] = ssy
    res.append(statsdict)
    
    if NumColumns - int(self.mboolIntercept) - lintweight == 2 and r2 >= 0:
        corrdict = {}
        rankArray = [None] * len(y)
        for i in range(len(y)):
            if lintweight == 0:
                rankArray[i] = [DataArray[i][1], 0, DataArray[i][0], float(i + 1), float(1)]
            else:
                rankArray[i] = [DataArray[i][2], 0, DataArray[i][0], float(i + 1), DataArray[i][1]]
        ties, rankAvg, i0, i1 = 0, 0.0, 0, 0
        for i in range(len(y)):
            i0 = i + 1
            ties = 1
            rankAvg = rankArray[i][3]
            if i0 < len(y):
                while rankArray[i][2] == rankArray[i0][2]:
                    ties += 1
                    rankAvg += rankArray[i0][3]
                    i0 += 1
                    if i0 >= len(y):
                        break
            rankAvg /= float(ties)
            if ties > 1:
                for i1 in range(i, i + ties):
                    rankArray[i1][3] = rankAvg
            i += ties - i
        sortedList = sorted(rankArray, key=lambda x: (x[0]))
        rankArray.clear()
        for r in sortedList:
            rankArray.append(r)
        for i in range(len(y)):
            rankArray[i][1] = float(i + 1)
        for i in range(len(y)):
            i0 = i + 1
            ties = 1
            rankAvg = rankArray[i][1]
            if i0 < len(y):
                while rankArray[i][0] == rankArray[i0][0]:
                    ties += 1
                    rankAvg += rankArray[i0][1]
                    i0 += 1
                    if i0 >= len(y):
                        break
            rankAvg /= float(ties)
            if ties > 1:
                for i1 in range(i, i + ties):
                    rankArray[i1][1] = rankAvg
            i += ties - i
        if lintweight == 0:
            pearsonCoefficient = r2 ** 0.5
            if B[0][0] < 0:
                pearsonCoefficiend *= -1
            yRankMean, xRankMean, sumWeights = 0.0, 0.0, 0.0
            for i in range(len(y)):
                yRankMean += rankArray[i][3] * rankArray[i][4]
                xRankMean += rankArray[i][1] * rankArray[i][4]
                sumWeights += rankArray[i][4]
            yRankMean /= sumWeights
            xRankMean /= sumWeights
            numerator, denominatorA, denominatorB = 0.0, 0.0, 0.0
            for i in range(len(y)):
                numerator += rankArray[i][4] * (rankArray[i][1] - xRankMean) * (rankArray[i][3] - yRankMean)
                denominatorA += rankArray[i][4] * (rankArray[i][1] - xRankMean) ** 2.0
                denominatorB += rankArray[i][4] * (rankArray[i][3] - yRankMean) ** 2.0
            spearmanCoefficient = numerator / (denominatorA * denominatorB) ** 0.5
            spearmanCoefficientT = (len(y) - 2) ** 0.5 * (spearmanCoefficient ** 2.0 / (1.0 - spearmanCoefficient ** 2.0)) ** 0.5
            spearmanCoefficientTP = self.PfromT(spearmanCoefficientT, len(y) - 2) * 2
            yMean, xMean = 0.0, 0.0
            sumWeights = 0.0
            for i in range(len(y)):
                yMean += DataArray[i][0]
                xMean += DataArray[i][1]
                sumWeights += 1.0
            yMean /= sumWeights
            xMean /= sumWeights
            numerator, denominatorA, denominatorB = 0.0, 0.0, 0.0
            for i in range(len(y)):
                numerator += (DataArray[i][1] - xMean) * (DataArray[i][0] - yMean)
                denominatorA += (DataArray[i][1] - xMean) ** 2.0
                denominatorB += (DataArray[i][0] - yMean) ** 2.0
            pearsonCoefficient = numerator / (denominatorA * denominatorB) ** 0.5
            pearsonCoefficientT = (len(y) - 2) ** 0.5 * (pearsonCoefficient ** 2.0 / (1.0 - pearsonCoefficient ** 2.0)) ** 0.5
            pearsonCoefficientTP = self.PfromT(pearsonCoefficientT, len(y) - 2) * 2
            corrdict['pearsonCoefficient'] = pearsonCoefficient
            corrdict['pearsonCoefficientT'] = pearsonCoefficientT
            corrdict['pearsonCoefficientTP'] = pearsonCoefficientTP
            corrdict['spearmanCoefficient'] = spearmanCoefficient
            corrdict['spearmanCoefficientT'] = spearmanCoefficientT
            corrdict['spearmanCoefficientTP'] = spearmanCoefficientTP
        else:
            print('Weighted correlations not coded yet.')
            print('See EILinearRegression.vb, line 437.')
        res.append(corrdict)
    return res