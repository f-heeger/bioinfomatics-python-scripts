from sqlite3 import connect
import atexit

class SqliteCache(object):
    
    def __init__(self, filePath=None, indict=None, table="key2value", 
                 key="key", value="value"):
        if filePath is None:
            self.filePath = self.__class__.__name__+".db"
        else:
            self.filePath = filePath
        
        self.conf = {"table": table, "key": key, "value": value}
        
        #connect to file with sqlite 
        # and create the necessary db if it does not exist
        self.conn = connect(self.filePath)
#        self.conn = connect(":memory:")
        c = self.conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE "
                  "type='table' AND name='%(table)s'" % self.conf)
        if c.fetchone() is None:
            c.execute("CREATE TABLE %(table)s "
                      "(%(key)s TEXT PRIMARY KEY, %(value)s TEXT)" % self.conf)
        self.conn.commit()
        if not indict is None:
            c.executemany("INSERT INTO %(table)s VALUES (?,?)" % self.conf,
                          indict.items())
            self.conn.commit()
        atexit.register(self.close)
        
    def close(self):
        self.conn.close()
    
    def __setitem__(self, key, value):
        c = self.conn.cursor()
        c.execute("INSERT INTO %(table)s VALUES (?,?)" % self.conf,
                  (str(key), str(value)))
        self.conn.commit()
        
    def __getitem__(self, key):
        c = self.conn.cursor()
        c.execute("SELECT %(value)s FROM %(table)s WHERE %(key)s=?" % self.conf, 
                  (str(key),))
        res = c.fetchone()
        if res is None:
            raise KeyError()
        return res[0]
        
    def __delitem(self, key):
        c = self.conn.cursor()
        c.exceute("DELETE FROM %(table)s WHERE %(key)s=?" % self.conf, key)
        self.conn.commit()
        
        
    def keys(self):
        c = self.conn.cursor()
        c.excecute("SELECT %(key)s FROM %(table)s" % self.conf)
        return c.fetchall()
        
    def values(self):
        c = self.conn.cursor()
        c.excecute("SELECT %(value)s FROM %(table)s" % self.conf)
        return c.fetchall()
        
    def items(self):
        return zip(self.keys(), self.values())

class SqliteListCache(object):
    """NOTE: If used in MulticachedDict, take care that the cache is updated,
    if an item is added to the list (eg. with append) it will not automatically
    be added to the cache.
    """
    def __init__(self, filePath=None, indict=None, table="key2value", 
                 key="key", value="value"):
        if filePath is None:
            self.filePath = self.__class__.__name__+".db"
        else:
            self.filePath = filePath
        
        self.conf = {"table": table, "key": key, "value": value}
        
        #connect to file with sqlite 
        # and create the necessary db if it does not exist
        self.conn = connect(self.filePath)
#        self.conn = connect(":memory:")
        c = self.conn.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE "
                  "type='table' AND name='%(table)s'" % self.conf)
        if c.fetchone() is None:
            c.execute("CREATE TABLE %(table)s "
                      "(%(key)s TEXT, %(value)s TEXT)" % self.conf)
        self.conn.commit()
        if not indict is None:
            for key, valueList in indict:
                param = zip([key]*len(valueList), valueList)
                c.executemany("INSERT INTO %(table)s VALUES (?,?)" % self.conf,
                              param)
            self.conn.commit()
        atexit.register(self.close)
        
    def close(self):
        self.conn.close()
    
    def __setitem__(self, key, valueList):
        c = self.conn.cursor()
        param = zip([key]*len(valueList), valueList)
        c.executemany("INSERT INTO %(table)s VALUES (?,?)" % self.conf, param)
        self.conn.commit()
        
    def __getitem__(self, key):
        c = self.conn.cursor()
        c.execute("SELECT %(value)s FROM %(table)s WHERE %(key)s=?" % self.conf, 
                  (str(key),))
        #fetchall returns tuple although they have only one element
        res = [t[0] for t in c.fetchall()]
        if len(res) == 0:
            raise KeyError()
        return res
        
    def __delitem__(self, key):
        c = self.conn.cursor()
        c.exceute("DELETE FROM %(table)s WHERE %(key)s=?" % self.conf, key)
        self.conn.commit()
        
    def keys(self):
        c = self.conn.cursor()
        c.excecute("SELECT %(key)s FROM %(table)s" % self.conf)
        return set(c.fetchall())
        
    def values(self):
        c = self.conn.cursor()
        r = []
        for key in self.keys():
            r.apend(self[key])
        return r
        
    def items(self):
        return zip(self.keys(), self.values())

class MultiCachedDict(object):
    def __init__(self, indict=None, cacheList=[]):
        if indict is None:
            self.ramDict = dict()
        else:
            self.ramDict = dict(indict)
        self.cacheList = cacheList
        
    def __setitem__(self, key, value):
        self.ramDict[key] = value
        for cache in self.cacheList:
            try:
                cache[key] = value
            except NotWritableError:
                pass
        
    def __getitem__(self, key):
        try:
            #try first to find key in the normal python dict we hold in ram
            return self.ramDict[key]
        except KeyError:
            #try all caches one after another
            value = None
            for lvl, cache in enumerate(self.cacheList):
                try:
                    value = cache[key]
                    # set key <-> value in python dict in ram
                    self.ramDict[key] = value
                    # set key <-> value in all caches with a higher level
                    # if they are writable
                    for i in range(lvl):
                        try:
                            self.cacheList[i][key] = value
                        except NotWritableError:
                            pass
                    break
                except KeyError:
                    pass
            if value is None:
                #if key was not found in any of the caches raise exception
                raise KeyError("'%s' was not found in dict "
                               "and in none of the caches"
                               % str(key))
            
            return value
        
    
    def __delitem__(self, key):
        try:
            del self.ramDict[key]
        except KeyError:
            for cache in self.cacheList:
                try:
                    del cache[key]
                except NotWritableError:
                    pass #TODO: should I propagate this up and not pass it silently?
                except KeyError:
                    pass
        
#    def __contains__(self, key):
#        raise NotImplementedError()


class PersistantDict(MultiCachedDict):
    def __init__(self, indict=None, dbpath=None):
        dbdict = SqliteCache(dbpath, indict)
        MultiCachedDict.__init__(self, None, [dbdict])
        
    def keys(self):
        return self.cacheList[0].keys()
        
    def values(self):
        return self.cacheList[0].values()
    
    def items(self):
        return self.cacheList[0].items()
        
class NotWritableError(AttributeError):
    def __init__(self, m):
        AttributeError.__init__(self, m)
