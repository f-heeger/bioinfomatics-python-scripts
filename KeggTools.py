import csv
import time
try:
    from urllib2 import urlopen
except ImportError:
    #this might be due to beeing python3
    from urllib.request import urlopen

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
                    response = handle.read()
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
        with open(self.cachePath, "rb") as inFile:
            for row in csv.reader(inFile, delimiter=",", quotechar="\""):
                if row[1] == "None":
                    self[row[0]] = None
                else:
                    self[row[0]] = row[1]
    
    def save(self):
        """save csv cache"""
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "wb") as out:
            writer = csv.writer(out, delimiter=",", quotechar="\"", 
                                quoting=csv.QUOTE_MINIMAL)
            for key, value in self.items():
                writer.writerow([str(key), str(value)])

    def __del__(self):
        if self.useCache:
            self.save()
            
            
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
            
class KeggGeneToPathwayMap(KeggMap):
    """Maps KEGG gene IDs to KEGG pathway IDs vis the KEGG Rest API.
    
    The return value is a set of the IDs of all pathways the gene is part of.
    Pathway IDs are given without the "path" prefix.
    """
    baseUrl = "http://rest.kegg.jp/link/pathway"
    
    def save(self):
        """modifed save function to deal with the sets"""
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for gene, pathList in self.items():
            for path in pathList:
                tab.append([gene, path])
        with open(self.cachePath, "wb") as out:
            for row in tab:
                out.write(",".join([str(field) for field in row])+"\n")
    
    def load(self):
        """modifed load function to deal with the sets"""
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "rb")):
            gene, path = row
            if gene not in self:
                self[gene] = set()
            self[gene].add(path)
            
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
            
