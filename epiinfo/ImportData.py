from .randata import randata, csvToRandata, jsonToRandata, syncToRandata, sqliteToRandata

def eisync(pathandfile, initVector, passwordSalt, pwd):
    """ Reads an Epi Info sync file (encrypted XML) and
        returns a list of dictionaries.
        Parameters:
          pathandfile (str): the path and sync file name
          initVector (str): the Init Vector used to encrypt
            the data (obtained from the mobile app)
          passwordSalt (str): the Salt used to encrypt
            the data (obtained from the mobile app)
          pwd (str): the password used to encrypt the data
        Returns:
          list of dictionaries
    """
    return syncToRandata(pathandfile, initVector, passwordSalt, pwd).get_Dictlist()

def eicsv(pathandfile):
    """ Reads a CSV file and returns a list of dictionaries.
        Parameters:
          pathandfile (str): the path and CSV file name
        Returns:
          list of dictionaries
    """
    return csvToRandata(pathandfile).get_Dictlist()

def eijson(pathandfile):
    """ Reads a JSON file and returns a list of dictionaries.
        Parameters:
          pathandfile (str): the path and JSON file name
        Returns:
          list of dictionaries
    """
    return jsonToRandata(pathandfile).get_Dictlist()

def eisqlite(pathanddb, tablename, datadict):
    """ Reads a(n) SQLite table and populates a list of dict object
          and a data dictionary dict.
        Parameters:
          pathanddb (str): the path and SQLite .db file
          tablename (str): the name of the SQLite table
          dl       (list): the list of dict to be built
          datadict (dict): initially {} dict of column names and SQLite data types
        Returns:
          list of dictionaries
    """
    rd = randata([])
    sqliteToRandata(pathanddb, tablename, rd, datadict)
    return rd.get_Dictlist()
