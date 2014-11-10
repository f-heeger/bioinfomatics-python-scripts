from sqlite3 import connect

class SqliteCache(object):
    writable = True # this is a writable dict
    
    def __init__(self, filePath=None):
        if filePath is None:
            self.filePath = self.__class__.__name__+".db"
        else:
            self.filePath = filePath
        
        #connect to file with sqlite 
        # and create the necessary db if it does not exist
        self.conn = connect(self.filePath)
        c = self.conn.cursor()
        c.execute("CREATE TABLE IF NOT EXISTS key2value (key TEXT PRIMARY KEY, value TEXT);")
        self.conn.commit()
        
    def __setitem__(key, value):
        c = self.conn.cursor()
        c.execute("INSERT INTO key2value VALUES (%s,%s)" % (key, value))
        self.conn.commit()
        
    def __getitem__(key):
        c = self.conn.cursor()
        c.execute("SELECT value FROM key2value WHERE key=%s" % key)
        res = c.fetchone()
        return res[0]

    def __del__(self):
        self.conn.close()

class MultiCachedDict(object):
    def __init__(self, indict=None, cacheList=[]):
        self.ramDict = dict(indict)
        self.cacheList = cacheList
        
    def __setitem__(key, value):
        self.ramDict[key] = value
        for cache in self.cacheList:
            if cache.writable:
                cahche[key] = value
        
    def __getitem__(key):
        try:
            #try first to find key in the normal python dict we hold in ram
            return ramDict[key]
        except KeyError:
            #try all caches one after another
            for cache in self.cacheList:
                try:
                    return cache[key]
                except KeyError:
                    pass
            #if key was not found in any of the caches raise exception
            raise KeyError("'%s' was not found in dict "
                           "and in none of the caches"
                           % str(key))
        
    
    def __delitem__(key):
        raise NotImplemented()
        
    def __contains__(key):
        raise NotImplemented()
