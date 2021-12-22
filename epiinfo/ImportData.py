from .randata import randata, csvToRandata

def eicsv(pathandfile):
    """ Reads a CSV file and returns a list of dictionaries.
        Parameters:
          pathandfile (str): the path and CSV file name
        Returns:
          list of dictionaries
    """
    return csvToRandata(pathandfile).get_Dictlist()
