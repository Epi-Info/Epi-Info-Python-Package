class ResultsFormatting:
  def __init__(self):
    self.cssForTable = '''<style type="text/css">
    th
    {
    border-style:none;
    border-width:0px;
    padding-left:5px;
    padding-right:5px;
    padding-bottom:0px;
    background-color:#318CE7;
    color:#FFFFFF;
    text-align: center;
    }
    td
    {
    border-style:none;
    border-width:0px;
    padding-left:5px;
    padding-right:5px;
    padding-top:0px;
    padding-bottom:0px;
    background-color:#FFFFFF;
    text-align: center;
    }
</style>'''

  def single2X2(self, ivdict, taResults):
    htmlstring = self.cssForTable + '\n<table>'
    htmlstring += '\n<tr><th colspan = ' + str(len(taResults['VariableValues'][0][0]) + 2) + '>'
    htmlstring += ivdict['outcomeVariable'] + '</th></tr>'
    htmlstring += '\n<tr><th>' + str(ivdict['exposureVariables'][0]) + '</th>'
    for vv in taResults['VariableValues'][0][0]:
        htmlstring += '<th>' + str(vv) + '</th>'
    htmlstring += '<th>Total</th></tr>'
    i = 0
    for r in taResults['Tables'][0]:
        htmlstring += '\n<tr><th>' + str(taResults['VariableValues'][0][1][i]) + '</th>'
        j = 0
        for c in r:
            htmlstring += '<td>' + str(c) + '<br>'
            htmlstring += str(round(100 * c / taResults['RowTotals'][0][i], 2)) + '%<br>'
            htmlstring += str(round(100 * c / taResults['ColumnTotals'][0][j], 2)) + '%' + '</td>'
            j += 1
        htmlstring += '<td>' + str(taResults['RowTotals'][0][i]) + '<br>100%<br>'
        htmlstring += str(round(100 * taResults['RowTotals'][0][i] / sum(taResults['RowTotals'][0]), 2)) + '%' + '</tr>'
        i += 1
    htmlstring += '\n<tr><th>Total</th>'
    for ct in taResults['ColumnTotals'][0]:
        htmlstring += '<td>' + str(ct) + '<br>' + str(round(100 * ct / sum(taResults['ColumnTotals'][0]), 2)) + '%<br>100%' + '</td>'
    htmlstring += '<td>' + str(sum(taResults['ColumnTotals'][0])) + '<br>100%<br>100%</td></tr>\n</table>'
    statsdict = dict(taResults['Statistics'][0])
    htmlstring += '\n<table>'
    htmlstring += '\n<tr><td colspan=4><b>Single Table Analysis</b></td></tr>'
    htmlstring += '\n<tr><td></td><td></td><td colspan=2>95% Confidence Interval</td></tr>'
    htmlstring += '\n<tr><td></td><td>Point<br>Estimate</td><td>Lower</td><td>Upper</td></tr>'
    htmlstring += '\n<tr><td>PARAMETERS: Odds-based</td><td colspan=3></td></td>'
    htmlstring += '\n<tr><td>Odds Ratio (cross product)</td><td>' + str(round(statsdict['OR'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['ORLL'], 4)) + '</td>' + '<td>' + str(round(statsdict['ORUL'], 4)) + ' (T)</td></tr>'
    htmlstring += '\n<tr><td>Odds Ratio (MLE)</td><td>' + str(round(statsdict['MidPOR'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['MidPORLL'], 4)) + '</td>' + '<td>'
    htmlstring += str(round(statsdict['MidPORUL'], 4)) + ' (M)</td></tr>'
    htmlstring += '\n<tr><td></td><td></td>'
    htmlstring += '<td>' + str(round(statsdict['FisherORLL'], 4)) + '</td>' + '<td>'
    htmlstring += str(round(statsdict['FisherORUL'], 4)) + ' (F)</td></tr>'
    htmlstring += '\n<tr><td>PARAMETERS: Risk-based</td><td colspan=3></td></td>'
    htmlstring += '\n<tr><td>Risk Ratio (RR)</td><td>' + str(round(statsdict['RiskRatio'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['RiskRatioLL'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['RiskRatioUL'], 4)) + ' (T)</td></tr>'
    htmlstring += '\n<tr><td>Risk Difference (RD%)</td><td>' + str(round(statsdict['RiskDifference'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['RiskDifferenceLL'], 4)) + '</td>'
    htmlstring += '<td>' + str(round(statsdict['RiskDifferenceUL'], 4)) + ' (T)</td></tr>'
    htmlstring += '\n<tr><td colspan=4>(T=Taylor series; C=Cornfield; M=Mid-P; F=Fisher Exact)</td></tr>'
    htmlstring += '\n<tr><td>STATISTICAL TESTS</td><td>Chi-square</td><td>1-tailed p</td><td>2-tailed p</td></tr>'
    htmlstring += '\n<tr><td>Chi-square - uncorrected</td><td>' + str(round(statsdict['UncorrectedX2'], 4)) + '</td>'
    htmlstring += '<td></td>' + '<td>' + str('{:.8f}'.format(statsdict['UncorrectedX2P'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Chi-square - Mantel-Haenszel</td><td>' + str(round(statsdict['MHX2'], 4)) + '</td>'
    htmlstring += '<td></td>' + '<td>' + str('{:.8f}'.format(statsdict['MHX2P'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Chi-square - corrected (Yates)</td><td>' + str(round(statsdict['CorrectedX2'], 4)) + '</td>'
    htmlstring += '<td></td>' + '<td>' + str('{:.8f}'.format(statsdict['CorrectedX2P'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Mid-p exact</td><td></td><td>' + str('{:.8f}'.format(statsdict['MidPExact1Tail'])) + '</td><td></td></tr>'
    htmlstring += '\n<tr><td>Fisher exact</td><td></td><td>' + str('{:.8f}'.format(statsdict['FisherExact1Tail'])) + '</td>'
    htmlstring += '<td>' + str('{:.8f}'.format(statsdict['FisherExact2Tail'])) + '</td></tr>'
    htmlstring += '\n</table>'
    
    return htmlstring

  def strat2X2Summary(self, taSummary):
    htmlstring = '\n<table>'
    htmlstring += '\n<tr><td colspan=4><b>SUMMARY INFORMATION</b></td></tr>'
    htmlstring += '\n<tr><td></td><td>Point</td><td colspan=2>95% Confidence Interval</td></tr>'
    htmlstring += '\n<tr><td>Parameters</td><td>Estimate</td><td>Lower</td><td>Upper</td></tr>'
    htmlstring += '\n<tr><td>Odds Ratio Estimates</td><td colspan=3></td></tr>'
    htmlstring += '\n<tr><td>Adjusted OR (MH)</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedOddsRatio'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedOddsRatioMHLL'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedOddsRatioMHUL'])) + ' (R)</td></tr>'
    htmlstring += '\n<tr><td>Adjusted OR (MLE)</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['mleOR'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['ExactORLL'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['ExactORUL'])) + ' (F)</td></tr>'
    htmlstring += '\n<tr><td>Risk Ratio (RR)</td><td colspan=3></td></tr>'
    htmlstring += '\n<tr><td>Adjusted RR (MH)</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedRR'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedRRMHLL'])) + '</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['computedRRMHUL'])) + ' (T)</td></tr>'
    htmlstring += '\n<tr><td colspan=4>(T=Taylor series; R=RGB; M=Exact mid-P; F=Fisher exact)</td></tr>'
    htmlstring += '\n</table>'
    htmlstring += '\n<table>'
    htmlstring += '\n<tr><td>STATISTICAL TESTS (overall association)</td><td>Chi-square</td><td>1-tailed p</td><td>2-tailed p</td></tr>'
    htmlstring += '\n<tr><td>MH Chi-square - uncorrected</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['uncorrecedChiSquare'])) + '</td>'
    htmlstring += '<td></td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['uncorrectedChiSquareP'])) + '</td></tr>'
    htmlstring += '\n<tr><td>MH Chi-square - corrected</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['corrChisq'])) + '</td>'
    htmlstring += '<td></td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['corrChisqP'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Homogeneity Tests</td><td>Chi-square</td><td>1-tailed p</td><td>2-tailed p</td></tr>'
    htmlstring += '\n<tr><td>Breslow-Day-Tarone test for Odds Ratio</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdtOR'])) + '</td>'
    htmlstring += '<td></td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdtORP'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Breslow-Day test for Odds Ratio</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdOR'])) + '</td>'
    htmlstring += '<td></td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdORP'])) + '</td></tr>'
    htmlstring += '\n<tr><td>Breslow-Day test for Risk Ratio</td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdRR'])) + '</td>'
    htmlstring += '<td></td>'
    htmlstring += '<td>' + str('{:.4f}'.format(taSummary['bdRRP'])) + '</td></tr>'
    htmlstring += '\n</table>'
    return htmlstring

  def meansWithXTab(self, mResults):
    tstats = ''
    m25 = ''
    Pg73 = self.cssForTable + '\n<table>'
    Pg73 += '\n<tr><th colspan=8>Descriptive Statistics for Each Value of Crosstab Variable</th></tr>'
    Pg73 += '\n<tr><td></td><td></td><td>Obs</td><td>Total</td><td>Mean</td><td>Variance</td><td>Std Dev</td><td></td></tr>'
    secondsection = ''
    for m in mResults:
        if 'crosstabVariable' in m:
            Pg73 += '\n<tr><td></td>'
            Pg73 += '<td>' + str(m['crosstabVariable']) + '</td>'
            Pg73 += '<td>' + str(int(m['obs'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(m['total'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(m['mean'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(m['variance'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(m['std_dev'])) + '</td>'
            Pg73 += '<td></td></tr>'
            secondsection += '\n<tr><td>' + str(m['crosstabVariable']) + '</td>'
            secondsection += '<td></td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['min'])) + '</td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['q25'])) + '</td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['q50'])) + '</td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['q75'])) + '</td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['max'])) + '</td>'
            secondsection += '<td>' + str('{:.4f}'.format(m['mode'])) + '</td></tr>'
        if 'meansDiffPooled' in m:
            tstats += '\n<table>'
            tstats += '\n<tr><td></td><td></td><td><b>T-Test</b></td><td></td><td></td><td></td></tr>'
            tstats += '\n<tr><td></td><td>Method</td><td>Mean</td><td colspan=2>95% CL Mean</td><td>Std Dev</td></tr>'
            tstats += '\n<tr><td>Diff (Group 1 - Group 2)</td><td>Pooled</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['meansDiffPooled'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['lclPooled'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['uclPooled'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['stdDevDiff'])) + '</td></tr>'
            tstats += '\n<tr><td>Diff (Group 1 - Group 2)</td><td>Satterthwaite</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['meansDiffSatterthwaite'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['lclSatterthwaite'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['uclSatterthwaite'])) + '</td>'
            tstats += '<td></td></tr>'
            tstats += '\n<tr><td colspan=6></td></tr>'
            tstats += '\n<tr><td>Method</td><td>Variances</td><td>DF</td><td>t Value</td><td>Pr > |t|</td><td></td></tr>'
            tstats += '\n<tr><td>Pooled</td><td>Equal</td>'
            tstats += '<td>' + str(int(m['pooledDF'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['pooledT'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['pooledPT'])) + '</td></tr>'
            tstats += '\n<tr><td>Satterthwaite</td><td>Unequal</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['SatterthwaiteDF'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['SatterthwaiteT'])) + '</td>'
            tstats += '<td>' + str('{:.4f}'.format(m['SatterthwaitePT'])) + '</td></tr>'
            tstats += '\n</table>'
        if 'ssBetween' in m:
            m25 += '\n<table>'
            m25 += '\n<tr><td colspan=5><b>ANOVA, a Parametric Test for Inequality of Population Means</b></td></tr>'
            m25 += '\n<tr><td></td><td colspan=3>(For normally distributed data only)</td><td></td></tr>'
            m25 += '\n<tr><td>Variation</td><td>SS</td><td>df</td><td>MS</td><td>F statistic</td></tr>'
            m25 += '\n<tr><td>Between</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['ssBetween'])) + '</td>'
            m25 += '<td>' + str(int(m['dfBetween'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['msBetween'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['fStatistic'])) + '</td></tr>'
            m25 += '\n<tr><td>Within</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['ssWithin'])) + '</td>'
            m25 += '<td>' + str(int(m['dfWithin'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['msWithin'])) + '</td>'
            m25 += '<td></td></tr>'
            m25 += '\n<tr><td>Total</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['ssTotal'])) + '</td>'
            m25 += '<td>' + str(int(m['dfTotal'])) + '</td>'
            m25 += '<td></td>'
            m25 += '<td></td></tr>'
            m25 += '\n<tr><td colspan=3>P-value = ' + str('{:.5f}'.format(m['anovaPValue'])) + '</td><td></td><td></td></tr>'
            m25 += '\n</table>'
        if 'bartlettChiSquare' in m:
            m25 += '\n<table>'
            m25 += '\n<tr><td><b>Bartlett\'s Test for Inequality of Population Variances</b></td></tr>'
            m25 += '\n<tr><td>Bartlett\'s chi square = ' + str('{:.4f}'.format(m['bartlettChiSquare']))
            m25 += ' df = ' + str(int(m['dfBetween'])) + ' P value = ' + str('{:.5f}'.format(m['bartlettPValue'])) + '</td></tr>'
            m25 += '\n</table>\n<table>'
            m25 += '\n<tr><td>A small p-value (e.g., less than 0.05) suggests that the variances are not homogeneous '
            m25 += 'and that the ANOVA may not be appropriate.</td></tr>'
            m25 += '\n</table>'
        if 'kruskalWallisH' in m:
            m25 += '\n<table>'
            m25 += '\n<tr><td colspan=2><b>Mann-Whitney/Wilcoxon Two-Sample Test (Kruskal-Wallis test for two groups)</b></td></tr>'
            m25 += '\n<tr><td>Kruskal-Wallis H (equivalent to Chi square)</td>'
            m25 += '<td>' + str('{:.4f}'.format(m['kruskalWallisH'])) + '</td></tr>'
            m25 += '\n<tr><td>Degrees of freedom</td><td>' + str(int(m['kruskalWallisDF'])) + '</td></tr>'
            m25 += '\n<tr><td>Degrees of freedom</td><td>' + str('{:.5f}'.format(m['kruskalPValue'])) + '</td></tr>'
            m25 += '\n</table>'
    Pg73 += '\n<tr><td></td><td></td><td>Minimum</td><td>25%</td><td>Median</td><td>75%</td><td>Maximum</td><td>Mode</td></tr>'
    Pg73 += secondsection
    Pg73 += '\n</table>'
    Pg73 += tstats + m25
    return Pg73

  def linearRegression(self, lrResults):
    r2 = ''
    m25 = ''
    Pg73 = self.cssForTable + '\n<table>'
    Pg73 += '\n<tr><td><b>Linear Regression</b></td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><th>Variable</th><th>Coefficient</th><th>95% LCL</th><th>95% UCL</th>'
    Pg73 += '<th>Std Error</th><th>F-test</th><th>P-Value</th></tr>'
    for lr in lrResults:
        if 'variable' in lr:
            Pg73 += '\n<tr><th>' + str(lr['variable']) + '</th>'
            Pg73 += '<td>' + str('{:.4f}'.format(lr['beta'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(lr['lcl'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(lr['ucl'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(lr['stderror'])) + '</td>'
            Pg73 += '<td>' + str('{:.4f}'.format(lr['ftest'])) + '</td>'
            Pg73 += '<td>' + str('{:.6f}'.format(lr['pvalue'])) + '</td></tr>'
        if 'r2' in lr:
            r2 += '\n<table><tr><td><b>Correlation Coefficient: r^2</b></td><td>' + str('{:.4f}'.format(lr['r2'])) + '</td></tr></table>'
            m25 += '\n<table>'
            m25 += '\n<tr><th>Source</th><th>df</th><th>Sum of Squares</th><th>Mean Square</th><th>F-statistic</th></tr>'
            m25 += '\n<tr><th>Regression</th>'
            m25 += '<td>' + str('{:.0f}'.format(lr['regressionDF'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['sumOfSquares'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['meanSquare'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['fStatistic'])) + '</td></tr>'
            m25 += '\n<tr><th>Residuals</th>'
            m25 += '<td>' + str('{:.0f}'.format(lr['residualsDF'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['residualsSS'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['residualsMS'])) + '</td>'
            m25 += '<td></td></tr>'
            m25 += '\n<tr><th>Total</th>'
            m25 += '<td>' + str('{:.0f}'.format(lr['totalDF'])) + '</td>'
            m25 += '<td>' + str('{:.4f}'.format(lr['totalSS'])) + '</td>'
            m25 += '<td></td>'
            m25 += '<td></td></tr>'
            m25 += '\n</table>'
        if 'pearsonCoefficient' in lr:
            m25 += '\n<table>'
            m25 += '\n<tr><th>Pearson\'s Coefficient</th><th>T-Test Value</th><th>P-Value</th></tr>'
            m25 += '\n<tr><td>' + str('{:.4f}'.format(lr['pearsonCoefficient'])) + '</td>'
            m25 += '\n<td>' + str('{:.4f}'.format(lr['pearsonCoefficientT'])) + '</td>'
            m25 += '\n<td>' + str('{:.4f}'.format(lr['pearsonCoefficientTP'])) + '</td></tr>'
            m25 += '\n</table>'
            m25 += '\n<table>'
            m25 += '\n<tr><th>Spearman\'s Coefficient</th><th>T-Test Value</th><th>P-Value</th></tr>'
            m25 += '\n<tr><td>' + str('{:.4f}'.format(lr['spearmanCoefficient'])) + '</td>'
            m25 += '\n<td>' + str('{:.4f}'.format(lr['spearmanCoefficientT'])) + '</td>'
            m25 += '\n<td>' + str('{:.4f}'.format(lr['spearmanCoefficientTP'])) + '</td></tr>'
            m25 += '\n</table>'
    Pg73 += '\n</table>'
    Pg73 += r2 + m25
    return Pg73

  def csMeans(self, ivdict, cmResults):
    Pg73 = self.cssForTable + '\n<table>'
    Pg73 += '\n<tr><th colspan=5>' + ivdict['mainvar'] + '</th><th></th><th></th><th></th></tr>'
    Pg73 += '\n<tr><th rowspan=2></th><th rowspan=2>Count</th><th rowspan=2>Mean</th><th rowspan=2>Std Error</th><th rowspan=1 colspan=2>Confidence Limits</th><th rowspan=2>Minimum</th><th rowspan=2>Maximum</th></tr>'
    Pg73 += '\n<tr><th>Lower</th><th>Upper</th></tr>'
    for r in cmResults.get_Rows():
        Pg73 += '\n<tr><td>' + str(r.get_Label()) + '</td><td>' + str('{:.0f}'.format(r.get_Count())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_Mean())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_StdErr())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_LCL())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_UCL())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_Min())) + '</td>'
        Pg73 += '<td>' + str('{:.4f}'.format(r.get_Max())) + '</td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td colspan=2><b>Sample Design Included</b></td></tr>'
    Pg73 += '\n<tr><td>Weight Variable:</td><td>' + ivdict['weightvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>PSU Variable:</td><td>' + ivdict['psuvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>Stratification Variable:</td><td>' + ivdict['stratavar'] + '</td></tr>'
    Pg73 += '\n</table>'
    return Pg73

  def csFreq(self, ivdict, cfResults):
    Pg73 = self.cssForTable + '\n<table>'
    Pg73 += '\n<tr><th>' + ivdict['outcome_variable'] + '</th><th>Total</th>'
    rawtotal = 0
    weightedtotal = 0
    designeffect = 0
    for r in cfResults.get_Rows():
        rawtotal += r.get_Count()
        weightedtotal += r.get_WeightedCount()
        designeffect = r.get_DesignEffect()
        Pg73 += '\n<tr><td><b>' + str(r.get_Value()) + '</b></td><td>' + str(int(r.get_Count())) + ' (' + str('{:.1f}'.format(r.get_WeightedCount())) + ')</td></tr>'
        Pg73 += '\n<tr><td>Row %</td><td>' + str('{:.3f}'.format(r.get_RowPercent())) + '</td></tr>'
        Pg73 += '\n<tr><td>Col %</td><td>' + str('{:.3f}'.format(r.get_ColPercent())) + '</td></tr>'
        Pg73 += '\n<tr><td>SE %</td><td>' + str('{:.3f}'.format(r.get_SE())) + '</td></tr>'
        Pg73 += '\n<tr><td>LCL %</td><td>' + str('{:.3f}'.format(r.get_LCL())) + '</td></tr>'
        Pg73 += '\n<tr><td>UCL %</td><td>' + str('{:.3f}'.format(r.get_UCL())) + '</td></tr>'
    Pg73 += '\n<tr><td><b>Total</b></td><td>' + str(int(rawtotal)) + ' (' + str('{:.1f}'.format(weightedtotal)) + ')</td></tr>'
    Pg73 += '\n<tr><td>Design Effect</td><td>' + str('{:.2f}'.format(designeffect)) + '</td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td colspan=2><b>Sample Design Included</b></td></tr>'
    Pg73 += '\n<tr><td>Weight Variable:</td><td>' + ivdict['weightvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>PSU Variable:</td><td>' + ivdict['psuvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>Stratification Variable:</td><td>' + ivdict['stratavar'] + '</td></tr>'
    Pg73 += '\n</table>'
    return Pg73

  def csTables(self, ivdict, ctResults):
    outcomes = []
    exposures = []
    exposuresp = []
    for r in ctResults.get_Rows():
        if [r.get_Exposure()] in exposures:
            pass
        else:
            exposures.append([r.get_Exposure()])
            exposuresp.append(r.get_Exposure())
        if r.get_Outcome() in outcomes:
            pass
        else:
            outcomes.append(r.get_Outcome())
    for r in ctResults.get_Rows():
        ind = exposuresp.index(r.get_Exposure())
        rowlist = exposures[ind]
        rowlist.append([r.get_Count(), r.get_RowPercent(), r.get_ColPercent(), r.get_SE(), r.get_LCL(), r.get_UCL(), r.get_DesignEffect()])
    Pg73 = self.cssForTable + '\n<table>'
    Pg73 += '\n<tr><th>' + ivdict['exposure_variable'] + '</th><th colspan=' + str(len(outcomes)) + '>' + ivdict['outcome_variable'] + '</th><th></th></tr>'
    Pg73 += '\n<tr><td></td>'
    for oc in outcomes:
        Pg73 += '<td><b>' + str(oc) + '</b></td>'
    Pg73 += '<td><b>Total</b></td></tr>'
    for r in exposures:
        Pg73 += '\n<tr><td><b>' + str(r[0]) + '</b></td>'
        rowtotal = 0
        for i in range(1, len(outcomes) + 1):
            rowtotal += r[i][0]
            Pg73 += '<td>' + str(r[i][0]) + '</td>'
        Pg73 += '<td>' + str(rowtotal) + '</td></tr>'
        Pg73 += '\n<tr><td>Row %</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][1])) + '</td>'
        Pg73 += '<td>100.000</td></tr>'
        Pg73 += '\n<tr><td>Col %</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][2])) + '</td>'
        Pg73 += '<td></td></tr>'
        Pg73 += '\n<tr><td>SE %</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][3])) + '</td>'
        Pg73 += '<td></td></tr>'
        Pg73 += '\n<tr><td>LCL %</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][4])) + '</td>'
        Pg73 += '<td></td></tr>'
        Pg73 += '\n<tr><td>UCL %</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][5])) + '</td>'
        Pg73 += '<td></td></tr>'
        Pg73 += '\n<tr><td>Design Effect</td>'
        for i in range(1, len(outcomes) + 1):
            Pg73 += '<td>' + str('{:.3f}'.format(r[i][6])) + '</td>'
        Pg73 += '<td></td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td><b>COMPLEX SAMPLE DESIGN ANALYSIS OF 2 X 2 TABLE</b></td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td>Odds Ratio (OR)</td><td>' + str('{:.3f}'.format(ctResults.get_OddsRatio())) + '</td></tr>'
    Pg73 += '\n<tr><td>Standard Error</td><td>' + str('{:.3f}'.format(ctResults.get_StandardErrorOR())) + '</td></tr>'
    Pg73 += '\n<tr><td>95% Conf. Limits</td><td>(' + str('{:.3f}'.format(ctResults.get_LCLOR())) + ', ' + str('{:.3f}'.format(ctResults.get_UCLOR())) + ')</td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td>Risk Ratio (RR)</td><td>' + str('{:.3f}'.format(ctResults.get_RiskRatio())) + '</td></tr>'
    Pg73 += '\n<tr><td>Standard Error</td><td>' + str('{:.3f}'.format(ctResults.get_StandardErrorRR())) + '</td></tr>'
    Pg73 += '\n<tr><td>95% Conf. Limits</td><td>(' + str('{:.3f}'.format(ctResults.get_LCLRR())) + ', ' + str('{:.3f}'.format(ctResults.get_UCLRR())) + ')</td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\nRR = (Risk of ' + ivdict['outcome_variable'] + '=' + str(outcomes[0]) + ' if ' + ivdict['exposure_variable'] + '=' + str(exposuresp[0]) + ') / '
    Pg73 += '(Risk of ' + ivdict['outcome_variable'] + '=' + str(outcomes[0]) + ' if ' + ivdict['exposure_variable'] + '=' + str(exposuresp[1]) + ')<br><br>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td>Risk Difference (RD)</td><td>' + str('{:.3f}'.format(ctResults.get_RiskDifference())) + '</td></tr>'
    Pg73 += '\n<tr><td>Standard Error</td><td>' + str('{:.3f}'.format(ctResults.get_StandardErrorRD())) + '</td></tr>'
    Pg73 += '\n<tr><td>95% Conf. Limits</td><td>(' + str('{:.3f}'.format(ctResults.get_LCLRD())) + ', ' + str('{:.3f}'.format(ctResults.get_UCLRD())) + ')</td></tr>'
    Pg73 += '\n</table>'
    Pg73 += '\nRD = (Risk of ' + ivdict['outcome_variable'] + '=' + str(outcomes[0]) + ' if ' + ivdict['exposure_variable'] + '=' + str(exposuresp[0]) + ') - '
    Pg73 += '(Risk of ' + ivdict['outcome_variable'] + '=' + str(outcomes[0]) + ' if ' + ivdict['exposure_variable'] + '=' + str(exposuresp[1]) + ')<br><br>'
    Pg73 += '\n<table>'
    Pg73 += '\n<tr><td colspan=2><b>Sample Design Included</b></td></tr>'
    Pg73 += '\n<tr><td>Weight Variable:</td><td>' + ivdict['weightvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>PSU Variable:</td><td>' + ivdict['psuvar'] + '</td></tr>'
    Pg73 += '\n<tr><td>Stratification Variable:</td><td>' + ivdict['stratavar'] + '</td></tr>'
    Pg73 += '\n</table>'
    return Pg73