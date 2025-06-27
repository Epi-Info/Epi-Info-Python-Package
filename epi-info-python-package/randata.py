import io
import os
import csv
import json
import ijson
import sqlite3
import ast
import time
from Crypto.Cipher import AES
from Crypto.Protocol.KDF import PBKDF2
import binascii
import base64
import xml.etree.ElementTree as ET
from .EncryptionDecryptionKeys import EncryptionDecryptionKeys

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
        self._Path = None 
        self._TableName = None 
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
        self._Indexes = value
    Indexes = property(get_Indexes, set_Indexes)
        
    # Path to data storage
    def get_Path(self):
        """ Returns the object's path to storage files
            Parameters:    
              none
            Returns:                                    
              str
        """
        return self._Path
    def set_Path(self, value):
        """ Sets the object's path to storage files
            Parameters:    
              value (str)
            Returns:                                    
              none
        """
        self._Path = value
    Path = property(get_Path, set_Path)
        
    # Name of data table for storage in JSON format
    def get_TableName(self):
        """ Returns the object's name for associated JSON file.
            Parameters:    
              none
            Returns:                                    
              str
        """
        return self._TableName
    def set_TableName(self, value):
        """ Sets the object's name for associated JSON file.
            Parameters:    
              value (str)
            Returns:                                    
              none
        """
        self._TableName = value
    TableName = property(get_TableName, set_TableName)
    
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
        if not cs:
            for c in index:
                if c is None:
                    continue
                if value.lower() in c.lower():
                    allrows = allrows.union(index[c])
            return [self._Dictlist[row] for row in allrows]
        for c in index:
            if c is None:
                continue
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
            else:
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
    
    def rowsFromSource(self, setofrows):
        """ Returns a list of dictionaries from JSON source file whose
            values pass the given where condition.
            Parameters: 
              where (str)
              value (str)
            Returns:                                    
              list of dictionaries
        """
        rl = []
        lensetofrows = len(setofrows)
        with open(self.get_Path() + '/' + self.get_TableName() + '.json', 'r') as fin:
            i = 0
            lis = 0
            for li in ijson.items(fin, "item"):
                if i in setofrows:
                    rl.append(li)
                    lis += 1
                    if lis == lensetofrows:
                        break
                i += 1
        return rl
    
    def allRowsFromSource(self):
        """ Returns a list of dictionaries from JSON source file.
            Parameters: 
              where (str)
              value (str)
            Returns:                                    
              list of dictionaries
        """
        with open(self.get_Path() + '/' + self.get_TableName() + '.json', 'r') as fin:
            return json.load(fin)
    
    def querySource(self, where, values=None):
        """ Returns a list of dictionaries from JSON source file whose
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
        if self.get_Path() is None:
            return
        if self.get_TableName() is None:
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
                    return self.allRowsFromSource()
                rows = self.Dlindexes - index[value]
                return self.rowsFromSource(rows)
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
                return self.rowsFromSource(rows)
            else:
                indexname = where
                value = values
                if indexname not in self._Indexes:
                    print(indexname, 'not indexed')
                    return []
                index = self._Indexes[indexname]
                if value not in index:
                    return []
                rows = index[value]
                return self.rowsFromSource(rows)
        
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
                return self.rowsFromSource(rows)
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
                return self.rowsFromSource(allrows)
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
                return self.rowsFromSource(allrows)
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
                return self.rowsFromSource(allrows)
    
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
                if row[k] is not None:
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
                if row[k] is not None:
                    row[k] = float(row[k])
    
    def GetSQLiteDataTypes(self):
        """ Stores the randata object as a sqlite3 table.
            Parameters: 
              force (bool): optional, to force re-creation of sqlite3 table
            Returns:                                    
              none
        """
        datatypes = {}
        sqliteConnection = sqlite3.connect(self.get_Path() + '/' + self.get_TableName() + '.db')
        c = sqliteConnection.cursor()
        c.execute("select sql from sqlite_master where type = 'table' and name = '" + self.get_TableName() + "';")
        columns = list(next(zip(*c.description)))
        lod = []
        for row in c.fetchall():
          rd = {}
          for i in range(len(columns)):
            rd[columns[i]] = row[i]
          lod.append(rd)
        c.close()
        sqliteConnection.close()
        for pair in lod[0]['sql'].split('(')[1][:-1].split(', '):
          datatypes[pair.split(' ')[0]] = pair.split(' ')[1]
        return datatypes
    
    def AppendSQLiteTable(self):
        """ Stores the randata object as a sqlite3 table.
            Parameters: 
              force (bool): optional, to force re-creation of sqlite3 table
            Returns:                                    
              none
        """
        if self.get_Path() == None or self.get_TableName() == None:
          print('Path and Table Name must be set to perform this task')
          return
        if os.path.isdir(self.get_Path()) == False:
          print(self.get_Path(), 'does not exist.')
          return
        sqliteConnection = sqlite3.connect(self.get_Path() + '/' + self.get_TableName() + '.db')
        c = sqliteConnection.cursor()
        c.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='" + self.get_TableName() + "';")
        tableExists = c.fetchone()[0] > 0
        if tableExists == False:
          print(self.get_TableName(), 'does not exist.')
          return
          c.close()
          sqliteConnection.close()
        datatypes = self.GetSQLiteDataTypes()
        for row in self.get_Dictlist():
          insertStatement = 'insert into ' + self.get_TableName() + ' values('
          cols = 0
          for col in row:
            if cols > 0:
              insertStatement += ', '
            cols += 1
            if row[col] == None:
              insertStatement += 'NULL'
              continue
            if datatypes[col] == 'TEXT':
              insertStatement += "'"
            insertStatement += str(row[col])
            if datatypes[col] == 'TEXT':
              insertStatement += "'"
          insertStatement += ');'
          c.execute(insertStatement)
        sqliteConnection.commit()
        c.close()
        sqliteConnection.close()
    
    def CreateSQLiteTable(self, force=False):
        """ Stores the randata object as a sqlite3 table.
            Parameters: 
              force (bool): optional, to force re-creation of sqlite3 table
            Returns:                                    
              none
        """
        if self.get_Path() == None or self.get_TableName() == None:
          print('Path and Table Name must be set to perform this task')
          return
        if os.path.isdir(self.get_Path()) == False:
          print(self.get_Path(), 'does not exist.')
          return
        sqliteConnection = sqlite3.connect(self.get_Path() + '/' + self.get_TableName() + '.db')
        c = sqliteConnection.cursor()
        c.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='" + self.get_TableName() + "';")
        tableExists = c.fetchone()[0] > 0
        if tableExists == True and force == False:
          print(self.get_TableName(), 'already exists.')
          print('Include True parameter to replace the table.')
          c.close()
          sqliteConnection.close()
          return
        if force == True and tableExists == True:
          dropStatement = "drop table " + self.get_TableName()
          c.execute(dropStatement)
        datatypes = {}
        for col in self.get_Dictlist()[0]:
          datatypes[col] = 'NULL'
        for row in self.get_Dictlist():
          for col in row:
            if type(row[col]) is bytes:
              datatypes[col] = 'BLOB'
            elif type(row[col]) is str:
              datatypes[col] = 'TEXT'
            elif type(row[col]) is float and datatypes[col] != str:
              datatypes[col] = 'REAL'
            elif type(row[col]) is int and datatypes[col] != float and datatypes[col] != str:
              datatypes[col] = 'INTEGER'
        createStatement = "create table " + self.get_TableName() + "("
        cols = 0
        for col in datatypes:
          if cols > 0:
            createStatement += ', '
          createStatement += col + ' ' + datatypes[col]
          cols += 1
        createStatement += ');'
        c.execute(createStatement)
        for index in self.get_Indexes():
          createIndex = 'create index ' + self.get_TableName() + '_' + index + ' on ' + self.get_TableName() + '(' + index + ');'
          c.execute(createIndex)
        for row in self.get_Dictlist():
          insertStatement = 'insert into ' + self.get_TableName() + ' values('
          cols = 0
          for col in row:
            if cols > 0:
              insertStatement += ', '
            cols += 1
            if row[col] == None:
              insertStatement += 'NULL'
              continue
            if datatypes[col] == 'TEXT':
              insertStatement += "'"
            insertStatement += str(row[col])
            if datatypes[col] == 'TEXT':
              insertStatement += "'"
          insertStatement += ');'
          c.execute(insertStatement)
        sqliteConnection.commit()
        c.close()
        sqliteConnection.close()
    
    def queryDB(self, where, values=None):
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
        sqliteConnection = sqlite3.connect(self.get_Path() + '/' + self.get_TableName() + '.db')
        c = sqliteConnection.cursor()
        selectStatement = 'select * from ' + self.get_TableName() + ' where '
        
        if type(where) == str and values is None:
            if ' != ' in where:
                splitter = ' != '
                indexname = where.split(splitter)[0]
                value = where.split(splitter)[1]
                selectStatement += indexname + ' <> '
                cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                selectStatement += value
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
            elif ' == ' in where:
                splitter = ' == '
                indexname = where.split(splitter)[0]
                value = where.split(splitter)[1]
                selectStatement += indexname + ' = '
                cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                selectStatement += value
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
            elif ' IS NOT NULL' in where.upper():
                splitter = ' is not '
                indexname = where.lower().split(splitter)[0]
                selectStatement += indexname + ' is not null'
            elif ' IS NULL' in where.upper():
                splitter = ' is '
                indexname = where.lower().split(splitter)[0]
                selectStatement += indexname + ' is null '
            else:
                return []
        
        elif values is not None:
            if type(where) == str and type(values) == str:
                indexname = where
                value = values
                selectStatement += indexname + ' = '
                cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                selectStatement += value
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
            elif type(where) == list and type(values) == str:
                indexname = where[0]
                value = values
                selectStatement += indexname + ' = '
                cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                selectStatement += value
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                if len(where) > 1:
                  for indexname in where[1:]:
                    selectStatement += ' or ' + indexname + ' = '
                    cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
                    selectStatement += value
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
            elif type(where) == str and type(values) == list:
                indexname = where
                value = values[0]
                selectStatement += indexname + ' = '
                cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                selectStatement += value
                if cfetchall00 == 'TEXT':
                  selectStatement += "'"
                if len(values) > 1:
                  for value in values[1:]:
                    selectStatement += ' or ' + indexname + ' = '
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
                    selectStatement += value
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
            elif type(where) == list and type(values) == list:
                useOR = False
                for indexname in where:
                  cfetchall00 = self.GetSQLiteDataTypes()[indexname]
                  for value in values:
                    if useOR:
                      selectStatement += ' or '
                    selectStatement += indexname + ' = '
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
                    selectStatement += value
                    if cfetchall00 == 'TEXT':
                      selectStatement += "'"
                    useOR = True
        selectStatement += ';'
        c.execute(selectStatement)
        returnedcolumns = list(next(zip(*c.description)))
        returnedlod = []
        for queryrow in c.fetchall():
          rowdict = {}
          for returnedcolumnnumber in range(len(returnedcolumns)):
            rowdict[returnedcolumns[returnedcolumnnumber]] = queryrow[returnedcolumnnumber]
          returnedlod.append(rowdict)
        c.close()
        sqliteConnection.close()
        return returnedlod
    
    def CreateTable(self):
        """ Stores the randata object as a JSON file
            and supporting files.
            Parameters: 
              none
            Returns:                                    
              none
        """
        if self.get_Path() == None or self.get_TableName() == None:
          print('Path and Table Name must be set to perform this task')
          return
        if os.path.isdir(self.get_Path()) == False:
          print(self.get_Path(), 'does not exist.')
          return
        if os.path.isfile(self.get_Path() + '/' + self.get_TableName() + '.json') == True:
          print(self.get_TableName(), 'already exists.')
          return
        with open(self.get_Path() + '/' + self.get_TableName() + '.json', 'w') as fout:
          fout.write(json.dumps(self.get_Dictlist()))
        with open(self.get_Path() + '/' + self.get_TableName() + '_Indexes.json', 'w') as fout:
          fout.write(str(self.get_Indexes()))
        with open(self.get_Path() + '/' + self.get_TableName() + '_Dlindexes.json', 'w') as fout:
          fout.write(str(self.get_Dlindexes()))
    
    def SetSource(self, path, tablename):
        """ Points the object to a source JSON file and
            loads associated indexes to memory.
            Parameters: 
              path (str)
              tablename (str)
            Returns:                                    
              none
        """
        if os.path.isdir(path) == False:
          print(path, 'does not exist.')
          return
        if os.path.isfile(path + '/' + tablename + '.json') == False:
          print(tablename, 'does not exist.')
          return
        self.set_Path(path)
        self.set_TableName(tablename)
        indexes = ''
        with open(self.get_Path() + '/' + self.get_TableName() + '_Dlindexes.json', 'r') as fin:
          self.set_Dlindexes(str(fin.read()))
        if os.path.isfile(self.get_Path() + '/' + self.get_TableName() + '_Indexes.json'):
          with open(self.get_Path() + '/' + self.get_TableName() + '_Indexes.json', 'r') as fin:
            indexes = fin.read()
        else:
          print('no index file')
        self.set_Indexes(ast.literal_eval(indexes))

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

def syncToRandata(pathandfile, initVector, passwordSalt, pwd):
    """ Reads an Epi Info sync file (encrypted XML) and
        returns a randata object.
        Parameters:
          pathandfile (str): the path and sync file name
          initVector (str): the Init Vector used to encrypt
            the data (obtained from the mobile app)
          passwordSalt (str): the Salt used to encrypt
            the data (obtained from the mobile app)
          pwd (str): the password used to encrypt the data
        Returns:
          randata
    """
    lod = []
    with open (pathandfile, "r") as myfile:
        data=myfile.readlines()
    try:
        keyData = PBKDF2(pwd.encode("utf-8"), binascii.unhexlify(passwordSalt), 16, count=1000)
        keyArray = binascii.unhexlify(initVector)
        encryptedData = base64.standard_b64decode(data[0].encode("utf-8"))
        cipher = AES.new(keyData, AES.MODE_CBC, keyArray)
        a = cipher.decrypt(encryptedData)
        b = a[:-ord(a[len(a)-1:])]
        c = b.decode('utf-8')
    except Exception as e:
        print('Could not decrypt. Possible incorrect initVector, salt, or password. All are case-sensitive.')
        lod.append({'Error' : 'Could not decrypt. Possible incorrect initVector, salt, or password. All are case-sensitive.'})
        return randata(lod)
    try:
        root = ET.fromstring(c)
        for child in root:
            d = {}
            for k in child.attrib:
                if k.lower() == 'surveyresponseid':
                    d['GlobalRecordId'] = child.attrib[k]
                elif k.lower() == 'fkey':
                    d['FKey'] = child.attrib[k]
                else:
                    d[k] = child.attrib[k]
            for gcs in child:
                for gc in gcs:
                    for k in gc.attrib:
                        d[gc.attrib[k]] = gc.text
                        if gc.text == 'true':
                            d[gc.attrib[k]] = True
                        elif gc.text == 'false':
                            d[gc.attrib[k]] = False
            lod.append(d)
    except Exception as e:
        print('Could not decrypt. Possible incorrect initVector. initVector is case-sensitive.')
        lod.append({'Error' : 'Could not decrypt. Possible incorrect initVector. initVector is case-sensitive.'})
        return randata(lod)
    return randata(lod)

def jsonToRandata(pathandfile):
    """ Reads a JSON file and returns a randata object.
        Parameters:
          pathandfile (str): the path and JSON file name
        Returns:
          randata
    """
    with open(pathandfile, 'r') as fin:
        dictlist = json.loads(fin.read())
    return randata(dictlist)

def sqliteToRandata(pathanddb, tablename, rd, datadict):
    """ Reads a(n) SQLite table and populates a randata object
          and a data dictionary dict.
        Parameters:
          pathanddb (str): the path and SQLite .db file
          tablename (str): the name of the SQLite table
          rd    (randata): the randata object to be built
          datadict (dict): initially {} dict of column names and SQLite data types
        Returns:
          nothing
    """
    # Begin the dictlist and establish the SQLite connection
    dictlist = rd.get_Dictlist()
    sqlconn = sqlite3.connect(pathanddb)
    sqlconn.row_factory = sqlite3.Row
    cursor = sqlconn.cursor()

    # Query the table
    sqlstatement = "SELECT * FROM " + tablename
    qrslt = cursor.execute(sqlstatement).fetchall()

    # Loop the query result and add dicts to dictlist
    for row in qrslt:
        d = {}
        for col in row.keys():
            d[col] = row[col]
        dictlist.append(d)

    # Get the column names and types from the master table
    # Add them to datadict
    cursor.execute("select sql from sqlite_master where type = 'table' and name = '" + tablename + "';")
    columns = list(next(zip(*cursor.description)))
    lod = []
    for row in cursor.fetchall():
      rd = {}
      for i in range(len(columns)):
        rd[columns[i]] = row[i]
      lod.append(rd)
    for pair in lod[0]['sql'].split('(')[1][:-1].split(', '):
      if len(pair.split(' ')) > 1:
        datadict[pair.split(' ')[0]] = pair.split(' ')[1]
      else:
        datadict[pair.split(' ')[0]] = 'TBD'

    # Unlock the database
    cursor.close()
    sqlconn.close()


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
