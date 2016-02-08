import csv
import time
import sys
from warnings import warn
try:
    from urllib2 import urlopen
except ImportError:
    #this might be due to beeing python3
    from urllib.request import urlopen

from MultiCachedDict import MultiCachedDict, SqliteCache, SqliteListCache, NotWritableError

class CacheNotUsedError(Exception):
    """Exception that is raised if the user tries to do cache operations 
    (save, load) on a Map where the cache is not active"""

class KeggMap(dict):
    """Base class for a dictonary relying on the Kegg Rest API for mapping
    
    This is an abstract base class for a dictonary that uses the Kegg Rest API 
    to represent relations between certain Kegg governed objects.
    Key requests will be looked up from a dictonary. If the key is not available
    a request to Kegg will be send via the Kegg Rest API. The result of 
    this lookup will be stored in the local dictonary. The detaileds of the 
    request have to be implemeneted by inheriting classes. Only one request per 
    second will be send to not overload Kegg servers.
    The local dictonary can be saved to the hard drive as csv file and loaded
    when an object is copnstructed.
    """
    def __init__(self, indict={}, cachePath=None, retry=0, useCache=True):
        """Constuctor for mapping object.
        
        A dictonary of already known assignments can be passed via the indict
        parameter. The cachePath can is the path were the save and load fucntion
        will search for a csv file. retry gives the number of times a request 
        should be send to Kegg  again before giving up (0 means only send 
        once).
        """
        dict.__init__(self, indict)
        self.useCache = useCache
        if useCache:
            if not cachePath:
                #if no cache path was given default to the class name
                cachePath = self.__class__.__name__+".csv"
            self.cachePath = cachePath
            try:
                self.load()
            except IOError:
                pass
                #this happens if no cache was saved previously
        self.retry = retry
        self.lastReq = 0
       
        
    def __getitem__(self, key):
        if key not in self:
            for tries in range(self.retry+1):
                #only allow one request in 1 sec
                if time.time() - self.lastReq < 1:
                    time.sleep(1)
                try:
                    self.lastReq = time.time()
                    handle = self.requestFunction(key)
                    response = handle.read().decode("utf-8")
                except Exception as e:
                    if tries == self.retry:
                        raise KeyError("Problem with key: '%s'. "
                                       "It raised an error: %s" 
                                       % (key, str(e)))
                    else:
                        sys.stderr.write("Problem with key: '%s'. "
                                         "It raised an error: %s\n "
                                         "Will try again...\n" % (key, str(e)) )
                else:
                    #if no exception occured no retry is needed
                    break
            self.readResponse(response, key)
        return dict.__getitem__(self, key)
            
    def requestFunction(key):
        raise NotImplemented("This base class does not implement any request. "
                             "Use a more specific subclass.")
            
    def load(self):
        """load csv cache"""
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "r") as inFile:
            for row in csv.reader(inFile, delimiter=",", quotechar="\""):
                if row[1] == "None":
                    self[row[0]] = None
                else:
                    self[row[0]] = row[1]
    
    def save(self):
        """save csv cache"""
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "w") as out:
            writer = csv.writer(out, delimiter=",", quotechar="\"", 
                                quoting=csv.QUOTE_MINIMAL)
            for key, value in self.items():
                if value is None:
                    writer.writerow([key, "None"])
                else:
                    writer.writerow([key, value])

    def __del__(self):
        if self.useCache:
            self.save()
            
class KeggSetMap(KeggMap):
    """Kegg map object that can have multiple values for one key.
    """

    def save(self):
        """modifed save function to deal with the sets None values"""
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for key, valueList in self.items():
            if valueList is None:
                tab.append([key, "None"])
            for value in valueList:
                tab.append([key, value])
        with open(self.cachePath, "w") as out:
            writer = csv.writer(out, delimiter=",", quotechar="\"", 
                                quoting=csv.QUOTE_MINIMAL)
            for row in tab:
                writer.writerow(row)
    
    def load(self):
        """modifed load function to deal with the sets and None values"""
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "r") as csvFile:
            for row in csv.reader(csvFile):
                key, value = row
                if value == "None":
                    self[key] = None
                elif key not in self:
                    self[key] = set([value])
                else:
                    self[key].add(value)
            
class NcbiGiToKeggMap(KeggMap):
    """Maps Ncbi protein GIs to KEGG gene IDs via the Kegg Rest API.
    
    If a Gi is not found in KEGG the return value is None.
    """
    baseUrl = "http://rest.kegg.jp/conv/genes"
    
    def requestFunction(self, gi):
        return urlopen("%s/ncbi-gi:%s" % (self.baseUrl,str(gi)))
    
    def readResponse(self, resp, key):
        if len(resp.strip()) == 0:
            self[key] = None
        elif "\n" in resp.strip():
            raise ValueError("KEGG respsonse contains more than one line.")
        else:
            giStr, keggStr = resp.strip().split()
            gi = giStr.split(":")[1]
            self[key] = keggStr
            
class KeggGeneToPathwayMap(KeggSetMap):
    """Maps KEGG gene IDs to KEGG pathway IDs vis the KEGG Rest API.
    
    The return value is a set of the IDs of all pathways the gene is part of.
    Pathway IDs are given without the "path" prefix. The key has to be a KEGG
    gene ID.
    """
    baseUrl = "http://rest.kegg.jp/link/pathway"
                
    def requestFunction(self, keggGene):
        return urlopen("%s/%s" % (self.baseUrl, keggGene))
    
    def readResponse(self, resp, key):
        pathes = set([])
        for line in resp.split("\n"):
            if len(line.strip()) == 0:
                continue #skip empty lines
            gene, pathStr = line.strip().split()
            path = pathStr.split(":")[1]
            pathes.add(path)
        self[key] = pathes

class KeggKoToPathwayMap(KeggSetMap):
    """Maps KEGG onthology groups (KO) to KEGG pathway IDs vis the KEGG Rest 
    API.
    
    The return value is a set of the IDs of all (ko) pathways the KO is part of.
    Pathway IDs are given without the "path" prefix. The key has to be a KEGG
    KO ID.
    """
    baseUrl = "http://rest.kegg.jp/link/pathway"
                
    def requestFunction(self, keggKo):
        return urlopen("%s/%s" % (self.baseUrl, keggKo))
    
    def readResponse(self, resp, key):
        pathes = set([])
        for line in resp.split("\n"):
            if len(line.strip()) == 0:
                continue #skip empty lines
            ko, pathStr = line.strip().split()
            path = pathStr.split(":")[1]
            if path.startswith("ko"):
                pathes.add(path)
        self[key] = pathes

class CachedKeggKoToPathwayMap(MultiCachedDict):
    def __init__(self, dbPath):
        kegg = KeggKoToPathwayMap(useCache=False)
        database = SqliteListCache(filePath=dbPath, indict=None, 
                               table="keggKo2pathway", key="ko", 
                               value="pathway")
        MultiCachedDict.__init__(self, None, [database, kegg])

class KeggKoToEnzymeMap(KeggSetMap):
    """Maps KEGG onthology groups (KO) to enzyme EC numbers vis the KEGG Rest 
    API.
    
    The return value is a set of the EC numbers the KO is linked to.
    EC numbers are given without the "ec" prefix. The key has to be a KEGG
    KO ID.
    """
    
    baseUrl = "http://rest.kegg.jp/link/enzyme"
    
    def requestFunction(self, keggKo):
        return urlopen("%s/%s" % (self.baseUrl, keggKo))
    
    def readResponse(self, resp, key):
        ecNumbers = set([])
        for line in resp.split("\n"):
            if len(line.strip()) == 0:
                continue #skip empty lines
            ko, ecStr = line.strip().split()
            ec = ecStr.split(":")[1]
            ecNumbers.add(ec)
        self[key] = ecNumbers
        
class CachedKeggKoToEnzymeMap(MultiCachedDict):
    def __init__(self, dbPath):
        kegg = KeggKoToEnzymeMap(useCache=False)
        database = SqliteListCache(filePath=dbPath, indict=None, 
                               table="keggKo2enzyme", key="ko", 
                               value="enzyme")
        MultiCachedDict.__init__(self, None, [database, kegg])

class KeggPathwayIdToNameMap(KeggMap):
    """Maps KEGG pathway IDs to the pathway name via the KEGG Rest API.
    
    """ 
    
    baseUrl = "http://rest.kegg.jp/get"
    
    def requestFunction(self, pathwayId):
        return urlopen("%s/%s" % (self.baseUrl, pathwayId))
    
    def readResponse(self, resp, key):
        if len(resp.strip()) == 0:
            self[key] = None
        else:
            for line in resp.split("\n"):
                if line[0] in [" ","/"]:
                    #continued field from line before or end of entry
                    continue
                field, value = line.split(None,1)
                if field == "PATHWAY_MAP":
                    pId, name = value.split(None,1)
                    self[key] = name
                    break

class CachedKeggPathwayIdToNameMap(MultiCachedDict):
    def __init__(self, dbPath):
        kegg = KeggPathwayIdToNameMap(useCache=False)
        database = SqliteListCache(filePath=dbPath, indict=None, 
                               table="keggpathwayId2name", key="pathId", 
                               value="pathName")
        MultiCachedDict.__init__(self, None, [database, kegg])

class KeggReactionIdToEcMap(KeggSetMap):
    """Mapp KEGG reaction ID to the involved enzyms EC numbers via the KEGG 
    Rest API.
    
    Returns a set of strings of Enzym Comission (EC) numbers. Key must be a 
    KEGG reaction ID (R[0-9]{5}) or None if no enzmyes are found.
    """
    
    baseUrl = "http://rest.kegg.jp/get"
    
    def requestFunction(self, reactionId):
        return urlopen("%s/%s" % (self.baseUrl, reactionId))
        
    def readResponse(self, resp, key):
        enzymes = []
        for line in resp.split("\n"):
            if line[:6] == "ENZYME":
                for part in line[6:].strip().split(" "):
                    e = part.strip()
                    if len(e) > 0:
                        enzymes.append(e)
                break
        if len(enzymes) == 0:
            self[key] = None
        else:
            self[key] = set(enzymes)
        
class KeggProteinToKoMap(KeggMap):
    """Mapp KEGG protein ID to the Kegg Orthologous group (KO) via the KEGG 
    Rest API.
    
    If the protein is not part of a KO None is returned
    """
    
    baseUrl = "http://rest.kegg.jp/link/ko"
                
    def requestFunction(self, keggProtein):
        return urlopen("%s/%s" % (self.baseUrl, keggProtein))
    
    def readResponse(self, resp, key):
        if len(resp.strip()) == 0:
            #no link found
            self[key] = None
        elif "\n" in resp.strip():
            protStr, koStr = resp.strip().split("\n")[0].split("\t")
            self[key] = koStr
            warn("KEGG respsonse for %s contains more than one line." % key, 
                 AmbiguityWarning)
        else:
            protStr, koStr = resp.strip().split("\t")
            self[key] = koStr
            
class CachedKeggProteinToKoMap(MultiCachedDict):
    def __init__(self, dbPath):
        kegg = KeggProteinToKoMap(useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, 
                               table="keggProt2ko", key="protId", 
                               value="ko")
        MultiCachedDict.__init__(self, None, [database, kegg])
        
        
class AmbiguityWarning(UserWarning):
    def __init__(self, message):
        UserWarning.__init__(self, message)
