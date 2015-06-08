import csv
import sys
import time
import re

from suds.client import Client as SoapClient
from suds import plugin

from MultiCachedDict import MultiCachedDict, SqliteCache

class CacheNotUsedError(Exception):
    """Exception that is raised if the user tries to do cache operations 
    (save, load) on a Map where the cache is not active"""

class NcbiSoapMap(dict):
    """Base class for a dictonary relying in NCBI soap interface for assignments.
    
    This is an abstract base class for a dictonary that uses the NCBI soap 
    interface to represent relations between certain NCBI governed objects like
    genes, taxonomic objects or genomes.
    Key requests will be looked up from a dictonary. If the key is not available
    a request to NCBI will be send via the SOAP interface. The result of this 
    lookup will be stored in the local dictonary. The detaileds of the request
    have to be implemeneted by inheriting calsses.
    The local dictonary can be saved to the hard drive as csv file and loaded
    when an object is copnstructed.
    """
    def __init__(self, indict={}, cachePath=None, retry=0, useCache=True):
        """Constuctor for mapping object.
        
        A dictonary of already known assignments can be passed via the indict
        parameter. The cachePath can is the path were the save and load fucntion
        will search for a csv file. retry gives the number of times a request 
        should be send to NCBI after again before giving up (0 means only send 
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
       
        
    def __getitem__(self, key):
        if key not in self:
            for tries in range(self.retry+1):
                try:
                    response = self.requestFunction(key)
                    break
                except Exception as e:
                    if tries == self.retry:
                        raise KeyError("Problem with key: '%s'. "
                                       "It raised an error: %s" 
                                       % (key, str(e)))
                    else:
                        time.sleep(1)
                        sys.stderr.write("Problem with key: '%s'. "
                                         "It raised an error: %s\n "
                                         "Will try again...\n" % (key, str(e)) )
            self.readResponse(response, key)
        return dict.__getitem__(self, key)
            
    def requestFunction(key):
        raise NotImplemented("This base class does not implement any request. "
                             "Use a more specific subclass.")
            
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "rb")):
            self[row[0]] = row[1]
    
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "wb") as out:
            for key, value in self.items():
                out.write("%s,%s\n" % (str(key), str(value)))

    def __del__(self):
        if self.useCache:
            self.save()

class NcbiSoapBatchMap(NcbiSoapMap):
    """Mapping object that will take a list of key and return a list of values 
    instead of a single key value pair. Requests to NCBI will be done in batch 
    mode with all keys that are not in the cache."""
    
    def __getitem__(self, keyList):
        for key in keyList:
            ncbiReqList = [key for key in keyList if key not in self]
            
            for tries in range(self.retry+1):
                try:
                    response = self.requestFunction(",".join([str(k) for k 
                                                    in ncbiReqList]))
                    break
                except Exception as e:
                    if tries == self.retry:
                        raise ValueError("Problem with key: '%s'. "
                                         "It raised an error: %s" 
                                         % (key, str(e)))
                    else:
                        time.sleep(1)
                        sys.stderr.write("Problem with key: '%s'. "
                                         "It raised an error: %s\n "
                                         "Will try again...\n" % (key, str(e)) )
            self.readResponse(response, keyList)
        return [dict.__getitem__(self, key) for key in keyList]

class SpeciesName2TaxId(NcbiSoapMap):
    """Map scientific names of species to their NCBI taxonomy DB ID.
    
    """
    wsdlUrl="http://www.ncbi.nlm.nih.gov/soap/v2.0/eutils.wsdl"
    
    def requestFunction(self, key):
        cl = SoapClient(self.wsdlUrl)
        return cl.service.run_eSearch("taxonomy", key)
        
    def readResponse(self, resp, key):
        if int(resp.Count) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. It got multiple answers.")
        self[key] = resp.IdList[0][0]

class LineageMap(NcbiSoapMap):
    """Map NCBI taxonomy IDs to full lineage information from NCBI taxonomy.
    
    Returns a list of tuples of the form (<Rank>, <Taxonomy Node ID>, 
    <Taxonomy Node Name>).
    """
    wsdlUrl = "http://www.ncbi.nlm.nih.gov/soap/v2.0/efetch_taxon.wsdl"
    
    def requestFunction(self, key):
        cl = SoapClient(self.wsdlUrl)
        return cl.service.run_eFetch(str(key))
    
    def readResponse(self, resp, key):
        m = []
        for r in resp.TaxaSet[0][0]["LineageEx"][0]:
            m.append((r["Rank"], r["TaxId"], r["ScientificName"]))
        self[key] = m
        
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for tax, m in self.items():
            for rank, lTax, lName in m:
                tab.append([tax, rank, lTax, lName])
        with open(self.cachePath, "wb") as out:
            for row in tab:
                out.write(",".join([str(field) for field in row])+"\n")
    
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "rb")):
            tax, rank, lTax, lName = row
            if tax not in self:
                self[tax] = []
            self[tax].append((rank, lTax, lName))

class SingleLevelLineageMap(LineageMap):
    """Map NCBI taxonomy ID to a certain taxonomic level.
    
    
    Returns a list of tuples of the form (<Rank>, <Taxonomy Node ID>, 
    <Taxonomy Node Name>).
    """
    
    def __init__(self, level, indict={}, cachePath=None, retry=0, useCache=True):
        NcbiSoapMap.__init__(self, indict, cachePath, retry, useCache)
        self.level = level
        
    def readResponse(self, resp, key):
        m = []
        for r in resp.TaxaSet[0][0]["LineageEx"][0]:
            if r["Rank"] == self.level:
                m.append((r["Rank"], r["TaxId"], r["ScientificName"]))
        self[key] = m

class TaxonomyParentMap(NcbiSoapMap):
    """Map NCBI taxonomy IDs to the ID of its parent node in NCBI taxonomy.
    """
    wsdlUrl = "http://www.ncbi.nlm.nih.gov/soap/v2.0/efetch_taxon.wsdl"
    
    def requestFunction(self, key):
        cl = SoapClient(self.wsdlUrl)
        return cl.service.run_eFetch(str(key))
    
    def readResponse(self, resp, key):
        self[key] = resp.TaxaSet[0][0]["ParentTaxId"]

class NuclId2TaxIdMap(NcbiSoapMap):
    """Map NCBI nucleotide GIs to the taxonomy ID of the species they come from.
    
    """
    wsdlUrl = "http://www.ncbi.nlm.nih.gov/soap/v2.0/efetch_seq.wsdl"
    
    def requestFunction(self, key):
        cl = SoapClient(self.wsdlUrl)
        return cl.service.run_eFetch("nucleotide", str(key))
    
    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response did not contain taxonomy information.")
        if len(resp) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        
        #if there is only one feature for some reason 
        # the feature list is not actually a list
        #here is a work around:
        try:
            if type(resp["GBSet"][0]["GBSeq_feature-table"][0]) != type([]):
                feature_list = [resp["GBSet"][0]["GBSeq_feature-table"][0]]
            else:
                feature_list = resp["GBSet"][0]["GBSeq_feature-table"][0]
            for feature in feature_list:
                if feature["GBFeature_key"] == "source":
                    for qual in feature["GBFeature_quals"][0]:
                        if qual["GBQualifier_name"] == "db_xref":
                            match = re.match("taxon:(\d+)", 
                                             qual["GBQualifier_value"])
                            if match:
                                self[key] = match.group(1)
                                return True
        except Exception as e:
            sys.stderr.write(str(resp))
            raise KeyError("'%s' is not in the dictonary. "
                           "Reading NCBI response caused an exception."
                           "NCBI response was:\n%s" % (key, str(resp)))
        raise KeyError("'%s' is not in the dictonary. "
                       "NCBI response did not contain taxonomy inforamtion. "
                       "NCBI response was:\n%s" % (key, str(resp)[:100]))

class NuclId2SpeciesNameMap(NcbiSoapMap):
    """Map NCBI nucleotide GIs to the scientific name of the species they 
    come from.
    
    """
    
    wsdlUrl = "http://www.ncbi.nlm.nih.gov/soap/v2.0/efetch_seq.wsdl"
    
    def requestFunction(self, key):
        cl = SoapClient(self.wsdlUrl)
        return cl.service.run_eFetch("nucleotide", str(key))
    
    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response did not contain taxonomy information.")
        if len(resp) > 1:
            self[key] = None
            import pdb; pdb.set_trace()
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        organism = resp[0][0].GBSeq_organism
        if "Unknown" in organism:
            raise ValueError("Source of sequence is given as 'unknown' in NCBI")
        self[key] = organism
        

class CachedNuclId2TaxIdMap(MultiCachedDict):
    def __init__(self, dbPath):
        ncbi = NuclId2TaxIdMap(useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, table="gi2tax", 
                               key="gi", value="tax")
        MultiCachedDict.__init__(self, None, [database, ncbi])
        
class CachedTaxonomyParentMap(MultiCachedDict):
    def __init__(self, dbPath):
        ncbi = TaxonomyParentMap(useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, table="tax2parent", 
                               key="tax", value="parent")
        MultiCachedDict.__init__(self, None, [database, ncbi])

class NcbiTaxonomyTree(object):
    """Representation of the NCBI taxonoy tree.
    
    Can be querried for information on the tree with a taxonomy ID."""
    def __init__(self, cachePath=None):
        """Constructor for NCBI taxonomy tree class.
        
        If cache path is gicen and not None the database at that path will be 
        used as persistent cache.
        """
        if cachePath is None:
            self._parent = TaxonomyParentMap(useCache=False)
            self.cached = False
        else:
            self._parent = CachedTaxonomyParentMap(cachePath)
            self.cached = True
        
    def initialize(self, dmpPath):
        """Initialize the tree object with data from a NCBI taxonomy dump file.
        
        This only works if the cache is enabled.
        """
        if not self.cached:
            raise NotImplementedError("This tree is not using a cache. "
                                      "It can not be initialized.")
        with open(dmpPath) as dmp:
            for line in dmp:
                tax, parent, _ = line.split("\t|\t", 2)
                self._parent[tax] = parent
    
    def parent(self, taxId):
        """
        Return NCBI taxonomy ID of the parent node of the node given by the 
        NCBI taxonomy ID."""
        return self._parent[taxId]
        
    def lineage(self, taxId):
        """
        Return the list of NCBI taxonomy IDs of the ancestors (in the tree 
        sense) of the node given by the NCBI taxonomy ID. the list includes the 
        root node (1) as the first element."""
        t = taxId
        lineage = [str(t)]
        while int(t) != 1: # 1 == root
            t = str(self._parent[t])
            lineage.append(t)
        return lineage[::-1]
    
    def lca(self, taxId1, taxId2):
        """
        Return the NCBI taxonomy ID of the lowest common ancestor in the NCBI
        taxonomy tree of the two given taxonomy IDs. Note that there is 
        always at least on common ancestor (the root node)."""
        l1 = self.lineage(taxId1)
        l2 = self.lineage(taxId2)
        i = 0
        while i < len(l1) and i < len(l2) and l1[i] == l2[i]:
            i+=1
        return l1[i-1]

if __name__ == "__main__":
    import logging
    import re

    class DebugPlugin(plugin.MessagePlugin):
        def received(self, context):
            import pdb; pdb.set_trace()

    class InitDebugPlugin(plugin.InitPlugin):
        def initialized(self, context):
            import pdb; pdb.set_trace()

    def lineageTest(a,b):
        for rowA, rowB in zip(a, b):
            for colA, colB in zip(rowA, rowB):
                if colA != colB:
                    return False
        return True

#    logging.basicConfig(level=logging.INFO)
#    logging.getLogger('suds.client').setLevel(logging.DEBUG)
    
    sciNam2IdTestSet = {"Homo sapiens": "9606",
                        "Mus musculus": "10090",
                        "Human immunodeficiency virus 1": "11676",
                        "Escherichia coli": "562",
                        "Clavariopsis aquatica": "253306",
                        }
                        
    nucl2taxNameTestSet = {297822720: "Arabidopsis lyrata subsp. lyrata",
                           237825467: "Clavariopsis aquatica",
                           21694053: "Homo sapiens",
                           34559759: "Rupicapra rupicapra",
                           666359714: "Ursus maritimus",
                           }
                       
    nucl2taxIdTestSet = {297822720: "81972",
                         237825467: "253306",
                         21694053: "9606",
                         34559759: "34869",
                         666359714: "29073",
                         }
                         
    taxParentTestSet = {666681: "1055487",
                        4890: "451864",
                        63221: "9606",
                        9606: "9605",
                        153057: "327045",
                        } 
                        
                         
    
    def test(mapObj, data):
        for key, value in data.items():
            sys.stdout.write("%s\t" % key)
            try:
                returned = str(mapObj[key])
            except Exception as e:
                import traceback
                sys.stderr.write(traceback.format_exc())
                sys.stdout.write("\tFailed. Exception raised: \"%s\". "
                                 "See log file for more information\n" % str(e))
                continue
            sys.stdout.write(returned)
            if returned == value:
                sys.stdout.write("\t\033[92mOK\033[0m\n")
            else:
                sys.stdout.write("\t\033[91mFailed. Unexpected value: %s\033[0m\n" % tId)
            time.sleep(1)
    
    print("Running Ncbi Soap Tool tests")
    
    print("testing cache usage option")
    baseMap = NcbiSoapMap(useCache=False)
    try:
        baseMap.save()
    except CacheNotUsedError:
        print("Deactivation of cache is working as expected")
    
    print("testing scientific name to taxonomy ID map:")
    sciName2taxId = SpeciesName2TaxId()
    print("Successfully build mapping object")
    test(sciName2taxId, sciNam2IdTestSet)

    
    print("testing nuclear ID to scientific name (taxonomy) map:")
    nuclId2taxName = NuclId2SpeciesNameMap()
    print("Successfully build mapping object")
    test(nuclId2taxName, nucl2taxNameTestSet)
    
    print("testing nuclear ID to taxonomy ID map:")
    nuclId2taxId = NuclId2TaxIdMap()
    print("Successfully build mapping object")
    test(nuclId2taxId, nucl2taxIdTestSet)
    
    print("testing lineage map")
    lineageMap = LineageMap()
    print("Successfully build mapping object")
    sys.stdout.write("polar bear <-> brown bear\t")
    brownBear = lineageMap[9644]
    polarBear = lineageMap[29073]
    if lineageTest(brownBear, polarBear):
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write("Failed. Lineages should be the same\n")
    sys.stdout.write("E. coli <-> E. alberti\t")
    ecoli = lineageMap[562]
    ealbertii = lineageMap[208962]
    if lineageTest(ecoli, ealbertii):
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write("Failed. Lineages sould be the same\n")
    sys.stdout.write("mouse <-> rat\t")
    mouse = lineageMap[10090]
    rat = lineageMap[10116]
    if lineageTest(mouse, rat):
        sys.stdout.write("Failed. Lineages should be different\n")
    else:
        sys.stdout.write("OK\n")
    print("testing Cached gi to tax id mapping")
    cGi2tax = CachedNuclId2TaxIdMap("/tmp/testDb.db")
    print("Succsessfully build mapping object")
    #TODO write test for functionality
    print("testing taxonomy parent map")
    taxParent = TaxonomyParentMap()
    print("Succsessfully build mapping object")
    test(taxParent, taxParentTestSet)
    print("testing save/load functions")
    sciName2taxId.cachePath = "/tmp/testSave.csv"
    sciName2taxId.save()
    print("saved successful")
    sciName2taxId = SpeciesName2TaxId(cachePath="/tmp/testSave.csv")
    test(sciName2taxId, sciNam2IdTestSet)
    print("loaded successful")
    print("done testing")
