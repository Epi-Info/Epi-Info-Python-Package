import sqlite3

class MakeSQLite:
    def makeTable(self, dbfile, ds, dsname, sqldatatypes):
        sqlconn = sqlite3.connect(dbfile)
        cursor = sqlconn.cursor()
        cursor.execute("drop table if exists " + dsname)
        createstatement = 'create table ' + dsname + '('
        for k in ds[0]:
            if k in sqldatatypes:
                createstatement += k + " " + sqldatatypes[k] + ", "
            else:
                createstatement += k + ", "
        createstatement = createstatement[:-2] + ');'
        cursor.execute(createstatement)
        for s in ds:
            instmt = 'insert into ' + dsname + ' values('
            for k in s:
                if str(s[k]) == 'None':
                    instmt += "NULL, "
                elif sqldatatypes[k] == 'TEXT':
                    instmt += "'" + str(s[k]) + "', "
                else:
                    instmt += str(s[k]) + ", "
            instmt = instmt[:-2] + ');'
            cursor.execute(instmt)
        sqlconn.commit()
        res = cursor.fetchall()
        cursor.close()
        sqlconn.close()
        return res
