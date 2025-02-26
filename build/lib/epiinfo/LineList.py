class LineList:
  """ The LineList class does what you'd expect.
      Its Run function takes two arguments: a list of dicts,
      which is the dataset from which to make a line list;
      and a list of columns to list. The list of column
      names can be a list of strings that are the names of
      keys in the list of dicts, or it can contain a single
      '*' to list all columns.
      
      It returns a string that is an HTML table.

     Author: John Copeland
  """

  def Run(self, lod, cols):
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
    </style>
    <table>
    <tr>'''
    if len(cols) == 0:
        return '<h4>No columns specified</h4>'
    if len(lod) == 0:
        return '<h4>Data set is empty</h4>'
    if cols[0] == '*':
        cols = []
        for k in lod[0]:
            cols.append(k)
    for k in cols:
        hstring += '<th>' + k + '</th>'
    hstring += '</tr>'
    for r in lod:
        hstring += '\n<tr>'
        for k in cols:
            if r[k] is None:
                hstring += '<td></td>'
                continue
            hstring += '<td>' + str(r[k]) + '</td>'
        hstring += '</tr>'
    hstring += '\n</table>'
    return hstring