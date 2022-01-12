from .randata import randata, csvToRandata, jsonToRandata, syncToRandata

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
