import math

class Means:
  """ 
      
      It returns a list of dictionaries.

     Author: John Copeland
  """

  def __init__(self):
    self.weightVariable = None
    self.strataVarNameList = []
    self.strataVarValueList = []

  def CreateWorkingTable(self, pNumericVariable, pCrosstabVariable, localLOD):
    workingDict = {}
    unsortedList = []
    if pCrosstabVariable:
        finalList = []
        xTabValues = []
        for r in localLOD:
            if r[pCrosstabVariable] in xTabValues:
                continue
            xTabValues.append(r[pCrosstabVariable])
        xTabValues.sort()
        for xv in xTabValues:
            unsortedList = []
            workingDict = {}
            for r in localLOD:
                if r[pCrosstabVariable] == xv:
                    if float(r[pNumericVariable]) in workingDict:
                        workingDict[float(r[pNumericVariable])] += 1.0
                        continue
                    workingDict[float(r[pNumericVariable])] = 1.0
            for k in workingDict:
                unsortedList.append([float(k), workingDict[float(k)]])
            sortedList = sorted(unsortedList, key=lambda x: (x[0], x[1]))
            finalList.append({str(xv) : sortedList})
        return finalList

    for r in localLOD:
        if float(r[pNumericVariable]) in workingDict:
            workingDict[float(r[pNumericVariable])] += 1.0
            continue
        workingDict[float(r[pNumericVariable])] = 1.0                                
    for k in workingDict:
        unsortedList.append([float(k), workingDict[float(k)]])
    sortedList = sorted(unsortedList, key=lambda x: (x[0], x[1]))
    return sortedList

  def CalculateSSBetween(self, grandMean, freqs, avgs):
    retval = 0.0
    for x in range(0, len(freqs)):
        retval += freqs[x] * (avgs[x] - grandMean) ** 2
    return retval

  def CalculateSSWithin(self, freqs, vars):
    retval = 0.0
    for x in range(0, len(freqs)):
        retval += (freqs[x] -1) * vars[x]
    return retval

  def PfromX2(self, x, df):
    pi = 3.14159265359
    m = 0
    l = 0
    j = 0
    k = 0
    rr = 0.0
    ii = 0
    PfromX2 = 0
    if x < 0.000000001 or df < 1:
        return 1
    rr = 1.0
    ii = df
    while ii >= 2:
        rr *= ii
        ii -= 2
    k = math.exp(int((df + 1) * 0.5) * math.log(abs(x)) - x * 0.5) / rr
    if k < 0.00001:
        return 0
    if int(df * 0.5) == df * 0.5:
        j = 1
    else:
        j = (2 / x / pi) ** 0.5
    l = 1
    m = 1
    if not math.isnan(x) and not math.isinf(x):
        while True:
            df += 2
            m = m * x / df
            l += m
            if m < 0.00000001:
                break
    PfromX2 = 1 - j * k * l
    return PfromX2

  def CalculateChiSquare(self, dfWithin, pooledVariance, freqs, vars):
    denominator = 0
    result = 0
    for j in range(0, len(freqs)):
        if freqs[j] - 1 != 0 and vars[j] != 0:
            denominator += 1.0 / (freqs[j] - 1)
            result += (freqs[j] - 1) * math.log(vars[j])
        else:
            denominator = 0
            result = 0
    denominator = 1.0 + (1.0 / (3.0 * (len(freqs) - 1))) * (denominator - 1.0 / dfWithin)
    if denominator != 0 and pooledVariance != 0:
        result = (1.0 / denominator) * (dfWithin * math.log(pooledVariance) - result)
    return result

  def algama(self, s):
    algama = 0.0
    Z = 0.0
    F = 0.0
    x = s
    if x < 0:
        return algama
    if x < 7:
        F = 1
        Z = x - 1
        while True:
            Z += 1
            if Z < 7:
                x = Z
                F *= Z
            if Z >= 7:
                break
        x += 1
        F = -math.log(F)
    Z = 1 / x ** 2
    algama = F + (x - 0.5) * math.log(x) - x + 0.918938533204673 + (((-1 / 1680 * Z + 1 / 1260) * Z - 1 / 360) * Z + 1 / 12) / x
    return algama

  def PfromZ(self, Z):
    PfromZ = 0.0
    P = 0.0
    y = 0.0
    x = 0.0
    LTONE = 7
    UTZERO = 12
    CON = 1.28
    x = Z
    if Z < 0.0:
        x = -Z
    if x > UTZERO:
        if Z < 0.0:
            return 1
        else:
            return 0
    y = Z ** 2 / 2
    if x > CON:
        P = x - 0.151679116635 + 5.29330324926 / (x + 4.8385912808 - 15.1508972451 / (x + 0.742380924027 + 30.789933034 / (x + 3.99019417011)))
        P = x + 0.000398064794 + 1.986158381364 / P
        P = x - 0.000000038052 + 1.00000615302 / P
        P = 0.398942280385 * math.Exp(-y) / P
    else:
        P = y / (y + 5.75885480458 - 29.8213557808 / (y + 2.624331121679 + 48.6959930692 / (y + 5.92885724438)))
        P = 0.398942280444 - 0.399903438504 * P
        P = 0.5 - x * P
    if Z < 0:
        PfromZ = 1 - P
    else:
        PfromZ = P
    return PfromZ

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

  def TfromP(self, pdblProportion, pintDF):
    TfromP = 0.0
    if pintDF <= 0:
        return float('nan')
    TOLERANCE = 0.00000001
    ldblProp = 0.0
    ldblLeftX = 0.0
    ldblLeftT = 0.0
    ldblRightX = 0.0
    ldblRightT = 0.0
    ldblX = 0.0
    ldblT = 0.0
    i = 0
    
    if pdblProportion < 0 or pdblProportion > 1:
        return 0
    elif pdblProportion < 0.5:
        ldblProp = pdblProportion
    else:
        ldblProp = 1.0 - pdblProportion
    ldblProp = ldblProp / 2
    ldblLeftX = 0.0
    ldblLeftT = 0.5
    ldblRightX = 1.0
    while True:
        ldblRightX = 10. * ldblRightX
        ldblRightT = self.PfromT(ldblRightX, pintDF)
        if ldblRightT < ldblProp:
            break
    while True:
        ldblX = (ldblRightX + ldblLeftX) / 2
        ldblT = self.PfromT(ldblX, pintDF)
        if ldblT < ldblProp:
            ldblRightT = ldblT
            ldblRightX = ldblX
        else:
            ldblLeftT = ldblT
            ldblLeftX = ldblX
        if abs(ldblT - ldblProp) < TOLERANCE:
            break
    TfromP = ldblX
    return TfromP

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

  def fillFrequencies(self, DT, horizontalFrequencies, allLocalFrequencies):
    DT0 = []
    for dt0 in DT[0]:
        DT0.append(dt0)
    for dtr in DT[1:]:
        RowTotal = 0
        values = []
        c = 0
        for C in DT0[2:]:
            values.append(dtr[c + 2])
            RowTotal += dtr[c + 2]
            c += 1
        allLocalFrequencies.append(values)
        horizontalFrequencies.append(RowTotal)

  def buildDT(self, groups):
    DT = [['__Values__', '__Count__']]
    DT1 = []
    DTdict = {}
    for d in groups:
        for k in d:
            DT[0].append(str(k))
    for i in range(0, len(groups)):
        dict_i = dict(groups[i])
        for k in dict_i:
            list_k = dict_i[k]
            for r in list_k:
                if r[0] in DTdict:
                    DTdict[r[0]][0] += r[1]
                    continue
                DTdict[r[0]] = [r[1]]
    for k in DTdict:
        for c in DT[0][2:]:
            DTdict[k].append(0)
    for i in range(0, len(groups)):
        dict_i = dict(groups[i])
        for k in dict_i:
            kindex = list(DT[0]).index(str(k))
            list_k = dict_i[k]
            for r in list_k:
                DTdict[r[0]][kindex - 1] = r[1]
    for k in DTdict:
        klist = [k]
        for v in DTdict[k]:
            klist.append(v)
        DT1.append(klist)
    DT[1:] = sorted(DT1, key=lambda x: (x[0]))
    return DT

  def CalculateKruskalWallisH(self, freqHorizontal, freqVertical, allLocalFreqs, recordCount, DT):
    DT0 = []
    for dt0 in DT[0]:
        DT0.append(dt0)
    cf = 0.0
    avgr = 0.0
    greaterSize = 0
    if len(freqHorizontal) > len(freqVertical):
        greaterSize = len(freqHorizontal)
    else:
        greaterSize = len(freqVertical)
    sr = [0.0] * greaterSize
    adj = 0.0
    H = 0.0
    for i in range(len(allLocalFreqs)):
        totalHFreq = freqHorizontal[i]
        cf += totalHFreq
        avgr = cf - (totalHFreq - 1) / 2.0
        for l in range(len(allLocalFreqs[0])):
            sr[l] += allLocalFreqs[i][l] * avgr
        adj += totalHFreq * (totalHFreq ** 2 - 1)
    for i in range(len(freqVertical)):
        totalVFreq = freqVertical[i]
        if totalVFreq != 0:
            H += sr[i] * sr[i] / totalVFreq
    H = H * 12 / (recordCount * (recordCount + 1)) - 3 * (recordCount + 1)
    H = H / (1 - adj / (recordCount ** 3 - recordCount))
    return H

  def Execute_CrossTab(self, groups):
    res = []
    for g in groups:
        for k in g:
            m = self.Execute_Means(g[k])
            m[0]['crosstabVariable'] = str(k)
            res.append(m[0])
    
    crosstabs = len(res)
    grandSum = 0.0
    grandObservations = 0.0
    unweightedObservations = 0.0
    averagesList = []
    observationsList = []
    unweightedObservationsList = []
    variancesList = []
    for r in res:
        grandSum += r['total']
        unweightedObservations += r['obs']
        grandObservations += r['obs']
        averagesList.append(r['mean'])
        observationsList.append(r['obs'])
        unweightedObservationsList.append(r['obs'])
        variancesList.append(r['variance'])
    grandMean = grandSum / grandObservations
    ssBetween = self.CalculateSSBetween(grandMean, observationsList, averagesList)
    dfBetween = crosstabs - 1
    msBetween = ssBetween / dfBetween
    ssWithin = self.CalculateSSWithin(observationsList, variancesList)
    dfWithin = grandObservations - crosstabs
    dfError = unweightedObservations - crosstabs
    msWithin = ssWithin / dfWithin
    fStatistic = msBetween / msWithin
    anovaPValue = self.PfromF(fStatistic, dfBetween, dfError)
    chiSquare = self.CalculateChiSquare(dfWithin, msWithin, observationsList, variancesList)
    bartlettPValue = self.PfromX2(chiSquare, dfBetween)
    
    DT = self.buildDT(groups)
    horizontalFrequencies = []
    allLocalFrequencies = []
    verticalFrequencies = []
    self.fillFrequencies(DT, horizontalFrequencies, allLocalFrequencies)
    for o in observationsList:
        verticalFrequencies.append(o)
    kruskalWallisH = self.CalculateKruskalWallisH(horizontalFrequencies, verticalFrequencies, allLocalFrequencies, grandObservations, DT)
    kruskalPValue = self.PfromX2(kruskalWallisH, dfBetween)
    
    if len(groups) == 2:
        ttestdict = {}
        SatterthwaiteDF = ((variancesList[0] / observationsList[0] + variancesList[1] / observationsList[1]) ** 2) / (1.0 / (unweightedObservationsList[0] - 1.0) * ((variancesList[0] / observationsList[0]) ** 2) + 1.0 / (unweightedObservationsList[1] - 1.0) * ((variancesList[1] / observationsList[1]) ** 2))
        SEu = (variancesList[0] / observationsList[0] + variancesList[1] / observationsList[1]) ** 0.5
        SatterthwaiteDF = SEu ** 4 / (1.0 / (unweightedObservationsList[0] - 1.0) * (variancesList[0] / observationsList[0]) ** 2 + 1.0 / (unweightedObservationsList[1] - 1.0) * (variancesList[1] / observationsList[1]) ** 2)
        meansDiff = averagesList[0] - averagesList[1]
        stdDevDiff = (((unweightedObservationsList[0] - 1) * variancesList[0] + (unweightedObservationsList[1] - 1) * variancesList[1]) / (unweightedObservationsList[0] + unweightedObservationsList[1] - 2)) ** 0.5
        df = int(unweightedObservationsList[0]) + int(unweightedObservationsList[1]) - 2
        shortDF = df
        tProbability = 0.05
        intervalLength = self.TfromP(tProbability, shortDF) * stdDevDiff * (1 / observationsList[0] + 1 / observationsList[1]) ** 0.5
        tStatistic = meansDiff / (stdDevDiff * (1.0 / observationsList[0] + 1.0 / observationsList[1]) ** 0.5)
        pEqual = 2.0 * self.PfromT(tStatistic, int(df))
        tStatisticUnequal = meansDiff / SEu
        pUnequalLower = 2.0 * self.PfromT(tStatisticUnequal, int(math.ceil(SatterthwaiteDF)))
        pUnequalUpper = 2.0 * self.PfromT(tStatisticUnequal, int(math.floor(SatterthwaiteDF)))
        pUneqal = pUnequalLower + (SatterthwaiteDF - math.floor(SatterthwaiteDF)) * (pUnequalUpper - pUnequalLower)
        shortDFCeiling = int(math.ceil(SatterthwaiteDF))
        shortDFFloor = int(math.floor(SatterthwaiteDF))
        unEqualIntervalTLower = self.TfromP(tProbability, shortDFCeiling)
        unEqualIntervalTUpper = self.TfromP(tProbability, shortDFFloor)
        unEqualIntervalT = float(unEqualIntervalTLower) + (SatterthwaiteDF - math.floor(SatterthwaiteDF)) * float(unEqualIntervalTUpper - unEqualIntervalTLower)
        
        equalLCLMean = meansDiff - intervalLength
        equalUCLMean = meansDiff + intervalLength
        unequalLCLMean = meansDiff - SEu * unEqualIntervalT
        unequalUCLMean = meansDiff + SEu * unEqualIntervalT
        
        ttestdict['meansDiffPooled'] = meansDiff
        ttestdict['lclPooled'] = equalLCLMean
        ttestdict['uclPooled'] = equalUCLMean
        ttestdict['stdDevDiff'] = stdDevDiff
        ttestdict['meansDiffSatterthwaite'] = meansDiff
        ttestdict['lclSatterthwaite'] = unequalLCLMean
        ttestdict['uclSatterthwaite'] = unequalUCLMean
        ttestdict['pooledDF'] = df
        ttestdict['pooledT'] = tStatistic
        ttestdict['pooledPT'] = pEqual
        ttestdict['SatterthwaiteDF'] = SatterthwaiteDF
        ttestdict['SatterthwaiteT'] = tStatisticUnequal
        ttestdict['SatterthwaitePT'] = pUneqal
        res.append(ttestdict)
    
    anovadict = {}
    anovadict['ssBetween'] = ssBetween
    anovadict['dfBetween'] = dfBetween
    anovadict['msBetween'] = msBetween
    anovadict['fStatistic'] = fStatistic
    anovadict['ssWithin'] = ssWithin
    anovadict['dfWithin'] = dfWithin
    anovadict['msWithin'] = msWithin
    anovadict['ssTotal'] = ssBetween + ssWithin
    anovadict['dfTotal'] = dfBetween + dfWithin
    anovadict['anovaPValue'] = anovaPValue
    anovadict['bartlettChiSquare'] = chiSquare
    anovadict['bartlettPValue'] = bartlettPValue
    anovadict['kruskalWallisH'] = kruskalWallisH
    anovadict['kruskalWallisDF'] = dfBetween
    anovadict['kruskalPValue'] = kruskalPValue
    res.append(anovadict)
    return res

  def Execute_Means(self, ROWS):
    res = []
    obs = 0
    Mode = 0.0
    Total = 0
    Min = 0.0
    Max = 0.0
    Median = 0.0
    Q25 = 0.0
    Q74 = 0.0
    
    AccumulatedTotal = 0
    
    i = 0
    for R in ROWS:
        obs += R[1]
    Total = obs
    obs = 0
    Sum = 0
    Sum_Sqr = 0
    variance = 0.0
    std_ev = 0.0
    for R in ROWS:
        temp = R[0]
        currentCount = R[1]
        obs += currentCount
        Sum += temp * currentCount
    mean = Sum / obs
    modeCount = 0
    for R in ROWS:
        Sum_Sqr += ((R[0] - mean) * (R[0] - mean)) * R[1]
        if i == 0:
            Min = R[0]
        Max = R[0]
        currentCount = R[1]
        if currentCount > modeCount:
            modeCount = currentCount
            Mode = R[0]
        i += R[1]
        
    # Quantile section
    m25 = Total * 0.25
    m50 = Total * 0.5
    m75 = Total * 0.75
    resultObs = 0
    oldValue = 0
    for SortedRowIndex in range(0, len(ROWS)):
        R = ROWS[SortedRowIndex]
        temp = R[0]
        currentCount = R[1]
        p1 = 0.0
        p2 = 0.0
        if m25 > resultObs and m25 <= resultObs + currentCount:
            if m25 == resultObs + currentCount:
                if SortedRowIndex + 1 < len(ROWS):
                    oldValue = ROWS[SortedRowIndex + 1][0]
                    if temp is not None and oldValue is not None:
                        p1 = temp
                        p2 = oldValue
                        Q25 = (p1 + p2) / 2.0
            else:
                Q25 = temp
        if m50 > resultObs and m50 <= resultObs + currentCount:
            if m50 == resultObs + currentCount:
                if SortedRowIndex + 1 < len(ROWS):
                    oldValue = ROWS[SortedRowIndex + 1][0]
                    if temp is not None and oldValue is not None:
                        p1 = temp
                        p2 = oldValue
                        Q50 = (p1 + p2) / 2.0
            else:
                Q50 = temp
        if m75 > resultObs and m75 <= resultObs + currentCount:
            if m75 == resultObs + currentCount:
                if SortedRowIndex + 1 < len(ROWS):
                    oldValue = ROWS[SortedRowIndex + 1][0]
                    if temp is not None and oldValue is not None:
                        p1 = temp
                        p2 = oldValue
                        Q75 = (p1 + p2) / 2.0
            else:
                Q75 = temp
        oldValue = temp
        resultObs += currentCount
    variance = Sum_Sqr / (i - 1)
    std_dev = variance ** 0.5
    statsdict = {'obs': Total, 'total' : Sum, 'mean' : mean, 'variance' : variance, 'std_dev': std_dev}
    statsdict['min'] = Min
    statsdict['q25'] = Q25
    statsdict['q50'] = Q50
    statsdict['q75'] = Q75
    statsdict['max'] = Max
    statsdict['mode'] = Mode
    res.append(statsdict)
    return res

  def Run(self, cols, ulod):
    lod = ulod
    res = []
    xTabVar = None
    if 'crosstabVariable' in cols:
        xTabVar = cols['crosstabVariable']
    workingTable = self.CreateWorkingTable(cols['meanVariable'],
                                           xTabVar,
                                           lod)
    if xTabVar:
        res = self.Execute_CrossTab(workingTable)
    else:
        res = self.Execute_Means(workingTable)
    return res