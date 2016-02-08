import urllib.request
import urllib.parse
import re
import gzip
import time
import sys

#based on Uniprot retrieval example on 
# http://www.uniprot.org/help/programmatic_access

from MultiCachedDict import MultiCachedDict, SqliteCache, SqliteListCache, NotWritableError

class CachedUniprotIdMap(MultiCachedDict):
    def __init__(self, dbpath, source="ACC", target="P_REFSEQ_AC", retry=0, 
                 delay=1, contact=None, returnNone=False, tablename=None, 
                 keyname=None, valuename=None):
        
        self.source = source
        self.target = target
        
        dbMap = SqliteCache(dbpath, table=tablename, key=keyname, 
                            value=valuename)
        uniprotMap = UniprotIdMap(source, target, retry, delay, contact, 
                                  returnNone)
        
        MultiCachedDict.__init__(self, None, [dbMap, uniprotMap])
        
    def initWithUniprotFlatFile(filepath, gzip=True):
        """
        Initialize the database cache with the content of a uniprot flat file.
        
        This will only work if source and target ID type a part of the 
        flat file. If the file is not gzipped this has to be specified.
        """
        
        if not gzip:
            tOpen = open
        else:
            tOpen = gzip.open
        
        name2pos = {"ACC" : 0, "ID" : 1, "P_ENTREZGENEID" : 2, 
                    "P_REFSEQ_AC" : 3, "P_GI" : 4, "PDB_ID" : 5, "NF100" : 7,
                    "NF90" : 8, "NF50" : 9, "UPARC" : 10, "PIR" : 11, 
                    "MIM_ID" : 13, "UNIGENE_ID" : 14, "EMBL_ID" : 16, 
                    "EMBL" : 17, "ENSEMBL_ID" : 18, "ENSEMBL_TRS_ID" : 19,
                    "ENSEMBL_PRO_ID" : 20, }
        if self.target not in name2pos:
            raise ValueError("Map target '%s' is not available from flat file."
                             % self.target)
        targetPos = name2pos[target]
        if self.source not in name2pos:
            raise ValueError("Map source '%s' is not available from flat file."
                             % self.source)
        sourcePos = name2pos[source]
        
        
        mapDict = {}
        for line in tOpen(filepath):
            arr = line.strip().split("\t")
            if not arr[targetPos] or not arr[sourcePos]:
                continue
            mapDict[arr[sourcePos]] = arr[targetPos]
        self.cacheList[0] = SqliteCache(filePath=self.cacheList[0].filepath, 
                                        indict=mapDict,
                                        **self.cacheList[0].conf)

class UniprotIdMap(object):

    """Base class for a ID dictionary based on uniprot web service
    """
    def __init__(self, source="ACC", target="P_REFSEQ_AC", retry=0, 
                 delay=1, contact=None, returnNone=False):
        """Constructor for ID dictionary
        
           source and target need to be abbreviations according to definition
           of the uniprot web service.
           retry and delay define retry behavior. On failure to connect to 
           uniprot the object will make retry attempts with delay seconds 
           between them.
           contact should be an e-mail where uniprot can contact you for 
           debugging.
        """
        self.source = source
        self.target = target
        self.retry = retry
        self.delay = delay
        self.contact = contact
        self.returnNone = returnNone
        
        self.response = None
        
        
    def __getitem__(self, key):
        for tries in range(self.retry+1):
            try:
                self.response = self.requestFunction(key)
                break
            except Exception as e:
                if tries == self.retry:
                    raise KeyError("Problem with key: '%s'. "
                                   "It raised an error: %s" 
                                   % (key, str(e)))
                else:
                    time.sleep(self.delay)
                    sys.stderr.write("Problem with key: '%s'. "
                                     "It raised an error: %s\n "
                                     "Will try again...\n" % (key, str(e)) )
            
        return self.readResponse(key)
            
    def __setitem__(self, key, value):
        raise NotWritableError("Values can not be set, because this "
                               "dictionary is based on Uniprot web service.")
    
    def __delitem__(self, key):
        raise NotWritableError("Values can not be deleted, because this "
                               "dictionary is based on Uniprot web service.")
    
    def requestFunction(self, key):
        url = "http://www.uniprot.org/mapping/"

        params = {
        "from": self.source,
        "to": self.target,
        "format": "tab",
        "query": key
        }

        data = urllib.parse.urlencode(params)
        request = urllib.request.Request(url, data.encode("ascii"))
        if self.contact:
            request.add_header(b'User-Agent', 
                               b'Python '+self.contact.encode("ascii"))
        return urllib.request.urlopen(request)

    def readResponse(self, key):
        lines = self.response.read().decode("utf-8").strip().split("\n")
        if len(lines) < 2:
            if self.returnNone:
                return None
            raise KeyError("Uniprot web service did not return information for "
                           "key: '%s'" % key)
        if len(lines) > 2:
            raise ValueError("Response from Uniprot has more than two lines. "
                             "Response was :\n%s" % "\n".join(lines))
        # line 0 is the header
        key, value = lines[1].strip().split()
        self.response = None
        return value
        



class UniprotInterface(object):
    """Interface class to retrieve UniProt entries and parse them. 
    
    Has functions to read GO, KEGG ans Cazy information.
    """
    def __init__(self, retry=0, delay=1, contact=None):
        """
        """
        self.retry = retry
        self.delay = delay
        self.contact = contact
        
        self.response = None
        self.uid = None
        
        
    def getData(self, uid):
        self.response = None
        self.uid = uid
        for tries in range(self.retry+1):
            try:
                self.response = self.requestFunction()
                break
            except Exception as e:
                if tries == self.retry:
                    raise ValueError("Problem with ID: '%s'. "
                                     "It raised an error: %s" 
                                      % (self.uid, str(e)))
                else:
                    time.sleep(delay)
                    sys.stderr.write("Problem with ID: '%s'. "
                                     "It raised an error: %s\n "
                                     "Will try again...\n" 
                                     % (self.uid, str(e)) )
            
    
    def clear(self):
        self.response = None
        self.uid = None
    
    def requestFunction(self):
        url = "http://www.uniprot.org/uniprot/%s.txt" % self.uid

        request = urllib.request.Request(url)
        if self.contact:
            request.add_header('User-Agent', 'Python %s' % self.contact)
        return urllib.request.urlopen(request)

    def readCazyInfo(self):
        if self.response is None:
            raise ValueError("No response available.")        
        try:
            line = next(self.response).decode("utf-8")
            #skip lines until one that starts with "DR    CAZy"
            while not line.startswith("DR   CAZy"):
                line = next(self.response).decode("utf-8")
        except StopIteration:
            # if the while loop run all the way no cazy info was found
            raise KeyError("Uniprot web service did not return CAZy information"
                           " for ID: '%s'" % self.uid)
        _, cId, cName = line.strip(".\n").split(";")
        return {"id": cId.strip(), "name": cName.strip()}
    
    def readGoInfo(self):
        if self.response is None:
            raise ValueError("No response available.")
        goInfo = []
        try:
            line = next(self.response).decode("utf-8")
            #skip lines until the first that starts with "DR"
            while line[:2] != "DR":
                line = next(self.response).decode("utf-8")
            #match lines while they start with DR
            while line[:2] == "DR":
                m=re.match("DR   GO; (?P<goId>GO:\d{7}); "
                           "(?P<aspect>[FPC]):(?P<function>[^;]+); "
                           "(?P<evidenceCode>[A-Z]+):(?P<evidenceSource>[^.]+)", 
                           line)
                if not m is None:
                    goInfo.append(m.groupdict())
                line = next(self.response).decode("utf-8")
        except StopIteration:
            pass
        if len(goInfo) == 0:
            raise ValueError("Uniprot web service did not return GO information"
                           " for ID: '%s'" % self.uid)
        return goInfo
        
    def readKeggInfo(self):
        if self.response is None:
            raise ValueError("No response available.")
        keggInfo = []
        try:
            line = next(self.response).decode("utf-8")
            #skip lines until the first that starts with "DR"
            while line[:2] != "DR":
                line = next(self.response).decode("utf-8")
            
            #match lines while they start with DR
            while line[:2] == "DR":
                if line[5:9] == "KEGG":
                    assert len(line.split(";")) == 3
                    keggInfo.append(line.split(";")[1].strip())
                line = next(self.response).decode("utf-8")
        except StopIteration:
            pass
        if len(keggInfo) == 0:
            raise ValueError("Uniprot web service did not return KEGG "
                             "information for ID: '%s'" % self.uid)
        return keggInfo

class UniprotToKeggMap(dict):
    """Map a uniprot ID to the Kegg gene IDs via the Uniprot REST API"""
    def __init__(self, indict={}, contact=None):
        dict.__init__(self, indict)
        self.interface = UniprotInterface(contact=contact)
        
    def __getitem__(self, key):
        if key not in self:
            try:
                self.interface.getData(key)
            except ValueError:
                raise KeyError("Key %s was not in the dict and Uniprot did not "
                               "return any information for it." % key)
            try:
                self[key] = set(self.interface.readKeggInfo())
            except ValueError as e:
                self[key] = set([])
        return dict.__getitem__(self, key)

class CachedUniprotToKeggMap(MultiCachedDict):
    def __init__(self, dbPath, email):
        uniprot = UniprotToKeggMap(contact=email)
        database = SqliteListCache(filePath=dbPath, indict=None, 
                                   table="uniprot2kegg", key="uId", 
                                   value="keggGene")
        MultiCachedDict.__init__(self, None, [database, uniprot])
        
class UniprotToGoMap(dict):
    """Map a uniprot ID to theGO IDs via the Uniprot REST API"""
    def __init__(self, indict={}, contact=None):
        dict.__init__(self, indict)
        self.interface = UniprotInterface(contact=contact)
        
    def __getitem__(self, key):
        if key not in self:
            try:
                self.interface.getData(key)
            except ValueError:
                raise KeyError("Key %s was not in the dict and Uniprot did not "
                               "return any information for it." % key)
            try:
                self[key] = set([entry["goId"] for entry in self.interface.readGoInfo()])
            except ValueError as e:
                self[key] = set([])
        return dict.__getitem__(self, key)
        
class CachedUniprotToGoMap(MultiCachedDict):
    def __init__(self, dbPath, email):
        uniprot = UniprotToGoMap(contact=email)
        database = SqliteListCache(filePath=dbPath, indict=None, 
                                   table="uniprot2go", key="uId", 
                                   value="go")
        MultiCachedDict.__init__(self, None, [database, uniprot])
