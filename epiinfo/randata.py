import io
import csv
import time

class randata:
    """ The randata class is a dataset as a list of dictionaries
        and allows the indexing of dictionary values for fast
        queries, as well as having other functions for efficient
        data management and analysis.

        Pronounced 'ron data'. Named for the star Ran, or Epsilon
        Eridani.

        It is included here mainly for its countdistinct
        function, which the CS routines use to count PSU
        and strata values.

       Author: John Copeland
    """
    def __init__(self, value=None):
        """ Initializes a randata object with or without a list
            of dictionaries.
            Parameters:
              value (list of dictionaries): the object's dataset
            Returns:
              randata
        """
        self._Dictlist = []
        self._Indexes = {}
        self._Dlindexes = set([])
        if value:
            if type(value) == list and len(value) > 0 and type(value[0]) == dict:
                self._Dictlist = value
                self._Dlindexes = set(range(len(self._Dictlist)))
        
    # list of dictionaries
    def get_Dictlist(self):                   
        """ Returns the Dictlist property, the list of dictionaries.
            Parameters:    
              none
            Returns:                                    
              list
        """
        return self._Dictlist
    def set_Dictlist(self, value):                   
        """ Sets the object's dataset.
            Parameters:    
              value (list of dictionaries)
            Returns:                                    
              none
        """
        self._Dictlist = value
        self._Dlindexes = set(range(len(self._Dictlist)))
    Dictlist = property(get_Dictlist, set_Dictlist)
    
    # set of Dictlist indexes
    def get_Dlindexes(self):                   
        """ Returns the object's set of indexes.
            Parameters:    
              none
            Returns:                                    
              set
        """
        return self._Dlindexes
    def set_Dlindexes(self, value):                   
        """ Initializes a randata object with or without a list
            of dictionaries.
            Parameters:    
              value (set): the object's set of indexes.
            Returns:                                    
              none
        """
        self._Dlindexes = value
    Dlindexes = property(get_Dlindexes, set_Dlindexes)
        
    # dictionary of Indexes
    def get_Indexes(self):
        """ Returns the object's dictionary of indexes.
            Parameters:    
              none
            Returns:                                    
              dict
        """
        return self._Indexes
    def set_Indexes(self, value):
        """ Returns without doing anything.
            Parameters:    
              value (dict)
            Returns:                                    
              none
        """
        return
    Indexes = property(get_Indexes, set_Indexes)
    
    def list(self):
        """ Returns the Dictlist property, the list of dictionaries.
            Parameters:    
              none
            Returns:                                    
              list
        """
        return self._Dictlist
    
    def deduplicate(self, columns):
        """ Returns a list of dictionaries that is a subset
            of the Dictlist property with duplicate dictionaries
            on the columns in the parameter removed.
            Parameters: 
              columns (str or list)
            Returns:                                    
              list
        """
        if type(columns) == str:
            columns = [columns]
        rl = []
        tupes = []
        for row in self._Dictlist:
            tupe = []
            for column in columns:
                tupe.append(row[column])
            if tupe not in tupes:
                tupes.append(tupe)
                rl.append(row)
        return rl
    
    def countdistinct(self, columns):
        """ Returns the number of dictionaries in the Dictlist
            after deduplicating on the columns in the parameter.
            Parameters: 
              columns (str or list)
            Returns:                                    
              int
        """
        if type(columns) != list:
            columns = [columns]
        tupes = []
        for row in self._Dictlist:
            tupe = set([row[c] for c in columns])
            if tupe not in tupes:
                tupes.append(tupe)
        return len(tupes)
    
    def index(self, column):
        """ Creates indexes for the key or keys in the column
            variable by making a dictionary of values and sets
            of list indexes where the key has the value.
            Parameters: 
              column (str or list)
            Returns:                                    
              none
        """
        if len(self._Dictlist) == 0:
            print("No data to index")
            return
        if column in self._Indexes:
            del self._Indexes[column]
        self._Indexes[column] = {}
        rownum = -1
        if type(self._Dictlist[0][column]) == list:
            for row in self._Dictlist:
                rownum += 1
                values = row[column]
                for value in values:
                    if value in self._Indexes[column]:
                        self._Indexes[column][value].add(rownum)
                    else:
                        self._Indexes[column][value] = set([rownum])
            return
        for row in self._Dictlist:
            rownum += 1
            value = row[column]
            if value in self._Indexes[column]:
                self._Indexes[column][value].add(rownum)
            else:
                self._Indexes[column][value] = set([rownum])
    
    def contains(self, col, value=None, cs=False):
        """ Returns a list of dictionaries from Dictlist whose
            values for the indicated key contain the indicated
            value (substring)
            Parameters: 
              col (str): the column name or key to be checked
              value (str): the substring to check
            Returns:                                    
              list of dictionaries
        """
        if col is None or value is None:
            return ['USAGE: randata.contains(str columnName, str substring[, bool CaseSensitive=False])']
        if len(col) == 0 or len(value) == 0:
            return ['USAGE: randata.contains(str columnName, str substring[, bool CaseSensitive=False])']
        if not (type(col) == str and type(value) == str):
            return ['USAGE: randata.contains(str columnName, str substring[, bool CaseSensitive=False])']
        if col not in self._Indexes:
            return ['Column ' + col + ' not indexed']
        
        index = self._Indexes[col]
        allrows = set()
        if cs:
            for c in index:
                if value.lower() in c.lower():
                    allrows = allrows.union(index[c])
            return [self._Dictlist[row] for row in allrows]
        for c in index:
            if value in c:
                allrows = allrows.union(index[c])
        return [self._Dictlist[row] for row in allrows]
    
    def query(self, where, values=None):
        """ Returns a list of dictionaries from Dictlist whose
            values pass the given where condition.
            Parameters: 
              where (str)
              value (str)
            Returns:                                    
              list of dictionaries
        """
        if where is None:
            return
        if len(where) == 0:
            return
        
        if type(where) == str and values is None:
            if ' != ' in where:
                splitter = ' != '
                indexname = where.split(splitter)[0]
                value = where.split(splitter)[1]
                if indexname not in self._Indexes:
                    return ['Column ' + indexname + ' not indexed']
                index = self._Indexes[indexname]
                if value not in index:
                    return self._Dictlist
                rows = self.Dlindexes - index[value]
                return [self._Dictlist[row] for row in rows]
            elif ' == ' in where:
                splitter = ' == '
                indexname = where.split(splitter)[0]
                value = where.split(splitter)[1]
                if indexname not in self._Indexes:
                    return []
                index = self._Indexes[indexname]
                if value not in index:
                    return []
                rows = index[value]
                return [self._Dictlist[row] for row in rows]
        
        elif values is not None:
            if type(where) == str and type(values) == str:
                indexname = where
                value = values
                if indexname not in self._Indexes:
                    print(indexname, 'not indexed')
                    return []
                index = self._Indexes[indexname]
                if value not in index:
                    return []
                rows = index[value]
                return [self._Dictlist[row] for row in rows]
            elif type(where) == list and type(values) == str:
                value = values
                allrows = set()
                for indexname in where:
                    if indexname not in self._Indexes:
                        print(indexname, 'not indexed')
                        continue
                    index = self._Indexes[indexname]
                    if value not in index:
                        continue
                    rows = index[value]
                    allrows = allrows.union(rows)
                return [self._Dictlist[row] for row in allrows]
            elif type(where) == str and type(values) == list:
                indexname = where
                allrows = set()
                for value in values:
                    if indexname not in self._Indexes:
                        print(indexname, 'not indexed')
                        return []
                    index = self._Indexes[indexname]
                    if value not in index:
                        continue
                    rows = index[value]
                    allrows = allrows.union(rows)
                return [self._Dictlist[row] for row in allrows]
            elif type(where) == list and type(values) == list:
                allrows = set()
                for indexname in where:
                    for value in values:
                        if indexname not in self._Indexes:
                            print(indexname, 'not indexed')
                            continue
                        index = self._Indexes[indexname]
                        if value not in index:
                            continue
                        rows = index[value]
                        allrows = allrows.union(rows)
                return [self._Dictlist[row] for row in allrows]
    
    def intcolumns(self, columns):
        """ Converts str values to int in Dictlist
            Parameters: 
              columns (str or list)
            Returns:                                    
              none
        """
        if type(columns) == str:
            columns = [columns]
        for row in self._Dictlist:
            for k in columns:
                row[k] = int(float(row[k]))
    
    def floatcolumns(self, columns):
        """ Converts str values to float in Dictlist
            Parameters: 
              columns (str or list)
            Returns:                                    
              none
        """
        if type(columns) == str:
            columns = [columns]
        for row in self._Dictlist:
            for k in columns:
                row[k] = float(row[k])

def isruddynumeric(val):
    """ Checks if a value is numeric
        even if it is a floating point
        number.
        Parameters:
          val(str)
        Returns:
          bool
    """
    try:
        float(val)
        return True
    except ValueError:
        return False

def csvToRandata(pathandfile):
    """ Reads a CSV file and returns a randata object.
        Parameters:
          pathandfile (str): the path and CSV file name
        Returns:
          randata
    """
    dictlist = []
    strings = {}
    floats = {}
    ints = {}
    changetoint = []
    changetofloat = []
    with open(pathandfile, 'r', errors='replace') as f:
        rv = csv.DictReader(f, skipinitialspace=True)
        for row in rv:
            datarow = {}
            for k, v in row.items():
                if v == None:
                    datarow[k] = None
                elif len(v) == 0:
                    datarow[k] = None
                elif v == 'None':
                    datarow[k] = None
                else:
                    datarow[k] = v
                    if isruddynumeric(v):
                        if int(float(v)) == float(v):
                            ints[k] = True
                        else:
                            floats[k] = True
                    else:
                        strings[k] = True
            dictlist.append(datarow)
        for c in ints:
            if c not in floats and c not in strings:
                changetoint.append(c)
        for c in floats:
            if c not in strings:
                changetofloat.append(c)
        rd = randata(dictlist)
        if len(changetoint) > 0:
            rd.intcolumns(changetoint)
        if len(changetofloat) > 0:
            rd.floatcolumns(changetofloat)
    return rd
