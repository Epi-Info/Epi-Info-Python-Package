import math

class Frequencies:
  """ The Frequencies class also does what you'd expect.
      Its Run function takes two arguments: a list of dicts,
      which is the dataset for which to compute frequencies;
      and a list of columns for which to compute frequencies.
      The list of column names can be a list of strings that
      are the names of keys in the list of dicts, or it can
      contain a single '*' to compute frequencies for all
      columns.
      
      It returns a string that is an HTML table.

     Author: John Copeland
  """

  def __init__(self):
    self.weightVariable = None

  def WILSON(self, a, m1, z, lowerl, upperl):
    LOWER = lowerl[0]
    Upper = upperl[0]
    
    p = a / m1
    q = 1 - p
    c1 = 2 * m1
    c2 = z * z
    c3 = p + c2 / c1
    c4 = c3 - 1
    c5 = z * ((p * (1 - p) + c2 / (4 * m1)) / m1) ** 0.5
    c6 = 1 + c2 / m1
    c7 = 1 / c6
    c8 = z * (c2 + (2 - 1 / m1) + 4 * p * (m1 * q - 1)) ** 0.5
    LOWER = (c3 - c5) * c7
    Upper = (c3 + c5) * c7
    
    lowerl.clear()
    upperl.clear()
    lowerl.append(LOWER)
    upperl.append(Upper)

  def Sub1(self, a, b, r8, r0, passes, a1l, d3l):
    a1 = a1l[0]
    d3 = d3l[0]
    
    if passes == 1:
        r2 = (a + 1) * 2.0
        r1 = b * 2.0
    else:
        r2 = a * 2.0
        r1 = (b + 1) * 2.0
    r3 = r8
    r7 = 2
    absr8 = r8
    if r8 < 0:
        absr8 = -r8
    r4 = (r2 * 0.5) * math.log(absr8)
    i = (r1 - 2) * 0.5
    if i == 0:
        d1 = math.exp(r4)
    else:
        r3 = 1 - r3
        a1 = r3 * (r2 * 0.5) * 10 ** -30
        r5 = a1 + 10 ** -30
        i -= 1
        while i > 0:
            r2 += 2
            r7 += 2
            a1 *= (r2 * r3) / r7
            r5 += a1
            i -= 1
        absr5 = r5
        if r5 < 0:
            absr5 = -r5
        d1 = math.exp(math.log(absr5) + math.log(10 ** 30) + r4)
    if passes == 1:
        d2 = 1 - d1
    else:
        d2 = d1
    d3 = r0 + d2
    a1l.clear()
    d3l.clear()
    a1l.append(a1)
    d3l.append(d3)

  def ExactCI(self, a, m1, p3, lowerl, upperl):
    lowerCI = lowerl[0]
    UpperCI = upperl[0]
    
    x = 0.01745506493
    a1 = 0
    d3 = 0
    b = m1 - a
    passes = 0
    flag2 = 0
    flag3 = 1
    flag4 = 0
    if a == 0:
        flag3 = 0
        tmp = a
        a = b
        b = tmp
    else:
        if a / m1 < 0.5:
            flag4 = 1
            tmp = a
            a = b
            b = tmp
    s1 = a / m1
    s3 = s1 * (1 - s1)
    r0 = (p3 / 100 - 1) * 0.5
    if b == 0:
        flag2 = 1
    while True:
        p = (s3 / m1) ** 0.5
        if p == 0:
            p = 1 / a - x
        if passes == 1:
            p = -1.0 * p
        r8 = s1 - p
        r9 = r8
        a1l = [a1]
        d3l = [d3]
        self.Sub1(a, b, r8, r0, passes, a1l, d3l)
        a1 = a1l[0]
        d3 = d3l[0]
        r6 = d3
        r8 -= x
        while True:
            a1l = [a1]
            d3l = [d3]
            self.Sub1(a, b, r8, r0, passes, a1l, d3l)
            a1 = a1l[0]
            d3 = d3l[0]
            d2 = d3
            d1 = r8
            d3 = ((r9 - r8) / (r6 - d3)) * d2
            r8 -= d3
            r9 = d1
            r6 = d2
            if -0.0000000001 < d3 / r8 < 0.0000000001:
                break
        if flag4 == 1:
            r8 = 1 - r8
        if flag3 == 1:
            flag3 = 0
        else:
            if flag2 == 1:
                flag2 = 0
                r8 = 1 - r8
                lowerCI = 0
            UpperCI = r8
            if r8 < 0:
                UpperCI = -r8
            if flag4 == 1:
                flag4 = 0
                tmp = lowerCI
                lowerCI = UpperCI
                UpperCI = tmp
            break
        lowerCI = r8
        if r8 < 0:
            lowerCI = -r8
        if flag2 == 1:
            flag2 = 0
            UpperCI = 0
            break
        passes = 1
    lowerl.clear()
    upperl.clear()
    lowerl.append(lowerCI)
    upperl.append(UpperCI)

  def GetConfLimit(self, value, frequency, count):
    lower = 0.0
    upper = 0.0
    if frequency == count:
        lower = 1.0
        upper = 1.0
        if count < 300:
            lowerl = [0.0]
            upperl = [1.0]
            self.ExactCI(frequency, count, 95.0, lowerl, upperl)
            lower = lowerl[0]
            upper = upperl[0]
    else:
        lowerl = [lower]
        upperl = [upper]
        if count < 300:
            self.ExactCI(frequency, count, 95.0, lowerl, upperl)
        else:
            self.WILSON(frequency, count, 95.0, lowerl, upperl)
        lower = lowerl[0]
        upper = upperl[0]
    return [value, lower, upper]

  def Run(self, ulod, cols):
    lod = ulod
    hstring = '''<style type="text/css">
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
    padding-bottom:0px;
    background-color:#FFFFFF;
    text-align: center;
    }
    </style>'''
    cls = []
    tot = 0
    if len(cols) == 0:
        return '<h4>No columns specified</h4>'
    if len(lod) == 0:
        return '<h4>Data set is empty</h4>'
    if cols[0] == '*':
        cols = []
        for k in lod[0]:
            cols.append(k)
    for k in cols:
        lod = sorted(ulod, key=lambda d: d[k])
        tot = 1
        hstring += '\n<table>'
        hstring += '\n<tr><th colspan="4">' + k + '</th></tr>'
        hstring += '\n<tr><th>' + 'Value' + '</th>'
        hstring += '<th>' + 'Frequency' + '</th>'
        hstring += '<th>' + 'Percent' + '</th>'
        hstring += '<th>' + 'Cum. Percent' + '</th></tr>'
        weightval = 1
        if self.weightVariable:
            weightval = lod[0][self.weightVariable]
        valsandfreqs = {str(lod[0][k]) : weightval}
        tot = weightval
        for r in lod[1:]:
            weightval = 1
            if self.weightVariable:
                weightval = r[self.weightVariable]
            tot += weightval
            if str(r[k]) in valsandfreqs:
                valsandfreqs[str(r[k])] += weightval
                continue
            valsandfreqs[str(r[k])] = weightval
        cpct = 0.0
        cls.clear()
        for v in valsandfreqs:
            cpct += (valsandfreqs[v] / tot) * 100
            hstring += '\n<tr><td>' + str(v) + '</th><td>' + str(valsandfreqs[v]) + '</td>'
            hstring += '<td>' + str(round((valsandfreqs[v] / tot) * 100, 2)) + '%' + '</td>'
            hstring += '<td>' + str(round(cpct, 2)) + '%' + '</td></tr>'
            cls.append(self.GetConfLimit(str(v), valsandfreqs[v], tot))
        hstring += '\n<tr><td>' + 'Total' + '</th><td>' + str(tot) + '</td><td>' + '100.00%' + '</td>'
        hstring += '<td>' + '100.00%' + '</td></tr>'
        if tot < 300:
            hstring += '\n<tr><td colspan="3"><b>Exact 95% Conf Limits</b></td></tr>'
        else:
            hstring += '\n<tr><td colspan="3"><b>Wilson 95% Conf Limits</b></td></tr>'
        for cl in cls:
            hstring += '\n<tr><td>' + str(cl[0]) + '</td>'
            hstring += '<td>' + str(round(cl[1] * 100, 2)) + '%' + '</td>'
            hstring += '<td>' + str(round(cl[2] * 100, 2)) + '%' + '</td></tr>'
        hstring += '\n</table>'
    return hstring