from pyopentree import tnrs_match_names, tnrs_contexts, gol_node_info

from MultiCachedDict import MultiCachedDict, SqliteCache, NotWritableError

class OttName2IdMap(dict):
    """Dictionary class to map names to OTT IDs
    
    Uses OTT web service to get information. Returns list of found IDs 
    for the name. A context can be provided to narrow down the search to a 
    certain taxonomic subtree (e.g. Fungi). 
    Retieve a list of valid contexts via the getContexts method
    """
    def __init__(self, indict=None, context=None):
        self.context = context
        if indict is None:
            dict.__init__(self)
        else:
            dict.__init__(self, indict)
        
    def __getitem__(self, name):
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            self[name] = self.ottCall(name)
            return dict.__getitem__(self, name)
            
    def ottCall(self, name):
        """Call the web service to get the ottId for a name.
        
        Returns a set of ottIds.
        """
        response = tnrs_match_names([name], 
                                    context_name=self.context, 
                                    do_approximate_matching=False)
        if len(response["results"]) == 0:
            raise KeyError("No OTT entry found for name %s" % name)
        assert len(response["results"]) == 1
#        assert len(response["results"][0]["matches"]) == 1
        matches = response["results"][0]["matches"]
        if len(matches)==1:
            return set([matches[0]["ot:ottId"]])
        if any([not m["is_synonym"] for m in matches]):
            return set([m["ot:ottId"] for m in matches if not m["is_synonym"]])
        return set([m["ot:ottId"] for m in matches])

    @classmethod
    def getContexts(self):
        """Get a list of a all valid contexts for open tree name search
        
        """
        return [element for l in tnrs_contexts().values() for element in l]

    @property
    def context(self):
        return self.__context
        
    @context.setter
    def context(self, value):
        if not value is None and value not in self.getContexts():
            raise ValueError("Not a valid context: %s. Use classmethod "
                             "getContext for list of valif contexts." % value)
        self.__context = value

class CachedOttName2IdMap(MultiCachedDict):
    """Same as OttName2IdMap, but with a sqlite cache.
    
    """
    def __init__(self, dbpath, context=None, tablename="ottName2Id", 
                 keyname="name", valuename="ottId"):
        
        dbMap = SqliteCache(dbpath, table=tablename, key=keyname, 
                            value=valuename)
        ottMap = OttName2IdMap(context=context)
        
        MultiCachedDict.__init__(self, None, [dbMap, ottMap])

class OttId2otherTaxonomyMap(dict):
    """Dictionary class to map OTT IDs to their source taxonomy IDs
    
    Returns a dictonary of the form {source_database: ID}, where the source 
    databases are: ncbi, if, gbif and irmng
    """
    
    def __getitem__(self, name):
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            self[name] = self.ottCall(name)
            return dict.__getitem__(self, name)
            
    def ottCall(self, ottId):
        """Call the web service to get the source taxonomy info for a ottId.
        
        """
        import pdb; pdb.set_trace()
        response = gol_node_info(ott_id=ottId)
        return dict([entry.split(":") for entry in \
                     response["tax_source"].split(",")])

class OttId2NcbiTaxIdMap(dict):
    """Dictionary class to map OTT IDs to Ncbi taxonomy IDs
    
    Raises KeyError if OTT does not have a Ncbi source entry for this node.
    """
    
    def __getitem__(self, name):
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            self[name] = self.ottCall(name)
            return dict.__getitem__(self, name)
            
    def ottCall(self, ottId):
        """Call the web service to get the source taxonomy info for a ottId.
        
        """
        response = pyopentree.gol_node_info(ott_id=ottId)
        taxonomies = dict([entry.split(":") for entry in \
                          response["tax_source"].split(",")])
        try:
            return taxonomies["ncbi"]
        except KeyError:
            raise KeyError("No NCBI source information was "
                           "found for ott node with ID %s" % ottId)

class CachedOttId2NcbiTaxIdMap(MultiCachedDict):
    """Same as OttId2NcbiTaxIdMap, but with a SQLite cache.
    
    """
    def __init__(self, dbpath, context=None, tablename="ott2ncbi", 
                 keyname="ottId", valuename="ncbiTaxId"):
        
        dbMap = SqliteCache(dbpath, table=tablename, key=keyname, 
                            value=valuename)
        ottMap = OttId2NcbiTaxIdMap()
        
        MultiCachedDict.__init__(self, None, [dbMap, ottMap])
        
    def initWithFlatFile(filepath, gzip=False):
        """
        Initialize the database cache with the content of a ott flat file.
        
        """
        
        if not gzip:
            tOpen = open
        else:
            tOpen = gzip.open
        
        mapDict = {}
        for line in tOpen(filepath):
            arr = line.strip().split("\t|\t")
            if not arr[0] or not arr[4]:
                continue
            source = dict([entry.split(":") for entry in source.split(",")])
            try:
                mapDict[arr[0]] = source["ncbi"]
            except KeyError:
                pass
        self.cacheList[0] = SqliteCache(filePath=self.cacheList[0].filepath, 
                                        indict=mapDict,
                                        **self.cacheList[0].conf)

