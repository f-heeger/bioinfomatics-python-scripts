import csv
import sys
import time
import re
import atexit

from suds.client import Client as SoapClient
from suds import plugin

from Bio import Entrez
from Bio.Entrez.Parser import ValidationError

from MultiCachedDict import MultiCachedDict, SqliteCache

class CacheNotUsedError(Exception):
    """Exception that is raised if the user tries to do cache operations 
    (save, load) on a Map where the cache is not active"""

class NcbiMap(dict):
    """Base class for a dictonary relying on the BioPython NCBI interface for 
    assignments.
    
    This is an abstract base class for a dictonary that uses the BioPython NCBI 
    interface to represent relations between certain NCBI governed objects like
    genes, taxonomic objects or genomes.
    Key requests will be looked up from a dictonary. If the key is not available
    a request to NCBI will be send via the BioPython interface. The result of 
    this lookup will be stored in the local dictonary. The detaileds of the 
    request have to be implemeneted by inheriting calsses.
    The local dictonary can be saved to the hard drive as csv file and loaded
    when an object is copnstructed.
    """
    def __init__(self, email, indict={}, cachePath=None, retry=0, useCache=True):
        """Constuctor for mapping object.
        
        A dictonary of already known assignments can be passed via the indict
        parameter. The cachePath can is the path were the save and load fucntion
        will search for a csv file. retry gives the number of times a request 
        should be send to NCBI after again before giving up (0 means only send 
        once).
        """
        dict.__init__(self, indict)
        Entrez.email = email
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
            atexit.register(self.save)
        self.retry = retry
       
        
    def __getitem__(self, key):
        if key not in self:
            for tries in range(self.retry+1):
                try:
                    handle = self.requestFunction(key)
                    response = Entrez.read(handle, validate=False)
                except ValidationError as e:
                    sys.stderr.write("Problem while parsing response for key '%s'. "
                                     "This might be due to outdatetd DTD files in "
                                     "Biopython. Try updating from http://www.ncbi."
                                     "nlm.nih.gov/data_specs/dtd/" % key)
                    raise e
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
                else:
                    #if no exception occured no retry is needed
                    break
            self.readResponse(response, key)
        return dict.__getitem__(self, key)
            
    def requestFunction(key):
        raise NotImplemented("This base class does not implement any request. "
                             "Use a more specific subclass.")
            
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "rb") as inFile:
            for row in csv.reader(inFile, delimiter=",", quotechar="\""):
                self[row[0]] = row[1]
    
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        with open(self.cachePath, "wb") as out:
            writer = csv.writer(out, delimiter=",", quotechar="\"", 
                                quoting=csv.QUOTE_MINIMAL)
            for key, value in self.items():
                writer.writerow([str(key), str(value)])


#TODO: rwrite for non soap version
#class NcbiSoapBatchMap(NcbiMap):
#    """Mapping object that will take a list of key and return a list of values 
#    instead of a single key value pair. Requests to NCBI will be done in batch 
#    mode with all keys that are not in the cache."""
#    
#    def __getitem__(self, keyList):
#        for key in keyList:
#            ncbiReqList = [key for key in keyList if key not in self]
#            
#            for tries in range(self.retry+1):
#                try:
#                    response = self.requestFunction(",".join([str(k) for k 
#                                                    in ncbiReqList]))
#                    break
#                except Exception as e:
#                    if tries == self.retry:
#                        raise ValueError("Problem with key: '%s'. "
#                                         "It raised an error: %s" 
#                                         % (key, str(e)))
#                    else:
#                        time.sleep(1)
#                        sys.stderr.write("Problem with key: '%s'. "
#                                         "It raised an error: %s\n "
#                                         "Will try again...\n" % (key, str(e)) )
#            self.readResponse(response, keyList)
#        return [dict.__getitem__(self, key) for key in keyList]

class SpeciesName2TaxId(NcbiMap):
    """Map scientific names of species to their NCBI taxonomy DB ID.
    
    """
    
    def requestFunction(self, key):
        return Entrez.esearch(db="taxonomy", term=key)
        
    def readResponse(self, resp, key):
        if int(resp["Count"]) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. It got multiple answers.")
        self[key] = resp["IdList"][0]

class TaxonomyNodeName2TaxId(NcbiMap):
    """Map names of NCBI taxonomy DB nodes to a list of IDs of nodes with this 
    name.
    
    One name can be mapped to several nodes (IDs), so the return value is always
    a set.
    """
    
    def requestFunction(self, key):
        return Entrez.esearch(db="taxonomy", term=key)
        
    def readResponse(self, resp, key):
        if int(resp["Count"]) == 0:
            raise KeyError()
        self[key] = set(resp["IdList"])
        
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for name, taxList in self.items():
            for tax in taxList:
                tab.append([name, tax])
        with open(self.cachePath, "wb") as out:
            for row in tab:
                out.write(",".join([str(field) for field in row])+"\n")
    
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "r")):
            name, tax = row
            if tax not in self:
                self[tax] = set()
            self[tax].add(name)

class LineageMap(NcbiMap):
    """Map NCBI taxonomy IDs to full lineage information from NCBI taxonomy.
    
    Returns a list of tuples of the form (<Rank>, <Taxonomy Node ID>, 
    <Taxonomy Node Name>).
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="taxonomy", id=str(key))
    
    def readResponse(self, resp, key):
        if len(resp) == 0:
            raise KeyError("No result for lineage of species '%s'" % key)
        m = []
        if "LineageEx" not in resp[0]:
            if resp[0]["ScientificName"] != "root":
                raise ValueError("Wired NCBi reponse without lineage info: %s" 
                                 % str(resp))
        else:
            for r in resp[0]["LineageEx"]:
                m.append((r["Rank"], r["TaxId"], r["ScientificName"]))
        m.append((resp[0]["Rank"], resp[0]["TaxId"], resp[0]["ScientificName"]))
        self[key] = m
        
    def save(self):
        if not self.useCache:
            raise CacheNotUsedError()
        tab = []
        for tax, m in self.items():
            for rank, lTax, lName in m:
                tab.append([tax, rank, lTax, '"%s"' % lName])
        with open(self.cachePath, "w") as out:
            for row in tab:
                out.write(",".join([str(field) for field in row])+"\n")
    
    def load(self):
        if not self.useCache:
            raise CacheNotUsedError()
        for row in csv.reader(open(self.cachePath, "r")):
            tax, rank, lTax, lName = row
            if tax not in self:
                self[tax] = []
            self[tax].append((rank, lTax, lName))

class SingleLevelLineageMap(LineageMap):
    """Map NCBI taxonomy ID to a certain taxonomic level.
    
    
    Returns a list of tuples of the form (<Rank>, <Taxonomy Node ID>, 
    <Taxonomy Node Name>).
    """
    
    def __init__(self, email, level, indict={}, cachePath=None, retry=0, useCache=True):
        NcbiMap.__init__(self, email, indict, cachePath, retry, useCache)
        self.level = level
        
    def readResponse(self, resp, key):
        if len(resp) == 0:
            raise KeyError("No result for lineage of species '%s'" % key)
        m = []
        if resp[0]["Rank"] == self.level:
            m.append((resp[0]["Rank"], resp[0]["TaxId"], resp[0]["ScientificName"]))
        for r in resp[0]["LineageEx"]:
            if r["Rank"] == self.level:
                m.append((r["Rank"], r["TaxId"], r["ScientificName"]))
                break
        self[key] = m

class TaxonomyParentMap(NcbiMap):
    """Map NCBI taxonomy IDs to the ID of its parent node in NCBI taxonomy.
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="taxonomy", id=str(key))
    
    def readResponse(self, resp, key):
        try:
            self[key] = resp[0]["ParentTaxId"]
        except TypeError:
            self[key] = -1

class NuclId2TaxIdMap(NcbiMap):
    """Map NCBI nucleotide GIs to the taxonomy ID of the species they come from.
    
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="nucleotide", id=str(key), retmode="xml")
    
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
            feature_list = resp[0]["GBSeq_feature-table"]
            for feature in feature_list:
                if feature["GBFeature_key"] == "source":
                    for qual in feature["GBFeature_quals"]:
                        if qual["GBQualifier_name"] == "db_xref":
                            match = re.match("taxon:(\d+)", 
                                             qual["GBQualifier_value"])
                            if match:
                                self[key] = match.group(1)
                                return
        except Exception as e:
            sys.stderr.write(str(resp))
            raise KeyError("'%s' is not in the dictonary. "
                           "Reading NCBI response caused an exception."
                           "NCBI response was:\n%s" % (key, str(resp)))
        raise KeyError("'%s' is not in the dictonary. "
                       "NCBI response did not contain taxonomy inforamtion. "
                       "NCBI response was:\n%s" % (key, str(resp)[:100]))


class TaxIdToSpeciesNameMap(NcbiMap):
    """Map NCBI taxonomy IDs to the scientific name of the species
    
    """
    def requestFunction(self, key):
        return Entrez.efetch(db="taxonomy", id=str(key), retmode="xml")

    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response was empty." % key)
        if len(resp) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        self[key] = resp[0]['ScientificName']

class NuclId2SpeciesNameMap(NcbiMap):
    """Map NCBI nucleotide GIs to the scientific name of the species they 
    come from.
    
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="nucleotide", id=str(key), retmode="xml")
    
    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response did not contain taxonomy information.")
        if len(resp) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        organism = resp[0]["GBSeq_organism"]
        if "Unknown" in organism:
            raise ValueError("Source of sequence is given as 'unknown' in NCBI")
        self[key] = organism
        
class ProtId2ProtNameMap(NcbiMap):
    """Map NCBI protein GIs to the title of the protein.
    
    """
    
    def requestFunction(self, key):
        return Entrez.efetch(db="protein", id=str(key), retmode="xml")
    
    def readResponse(self, resp, key):
        if len(resp) < 1:
            raise KeyError("'%s' is not in the dictionary. "
                           "NCBI response did not contain taxonomy information.")
        if len(resp) > 1:
            self[key] = None
            raise ValueError("Problem with key: %s. "
                             "It got multiple answers." % key)
        name = resp[0]["GBSeq_definition"]
        self[key] = name
        

#The following Cached* classes are versions of the before defined maps that are
# multi-level-cached (in RAM and in a SQLite DB) the actual map class is only
# used if the mapping can not be found in the caches

class CachedProtId2ProtNameMap(MultiCachedDict):
    def __init__(self, dbPath, email):
        ncbi = ProtId2ProtNameMap(email, useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, table="gi2name", 
                               key="gi", value="name")
        MultiCachedDict.__init__(self, None, [database, ncbi])
        
    def initialize(self, csvPath, delimiter=",", quotechar="\""):
        """Initialize the DB cache with a csv file.
        
        """
        with open(csvPath) as csvFile:
            for row in csv.reader(csvFile, delimiter=delimiter, 
                                  quotechar=quotechar):
                self[row[0]] = row[1]

class CachedNuclId2TaxIdMap(MultiCachedDict):
    def __init__(self, dbPath, email):
        ncbi = NuclId2TaxIdMap(email, useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, table="gi2tax", 
                               key="gi", value="tax")
        MultiCachedDict.__init__(self, None, [database, ncbi])
        
class CachedTaxonomyParentMap(MultiCachedDict):
    def __init__(self, dbPath, email):
        ncbi = TaxonomyParentMap(email, useCache=False)
        database = SqliteCache(filePath=dbPath, indict=None, table="tax2parent", 
                               key="tax", value="parent")
        MultiCachedDict.__init__(self, None, [database, ncbi])

class NcbiTaxonomyTree(object):
    """Representation of the NCBI taxonoy tree.
    
    Can be querried for information on the tree with a taxonomy ID."""
    def __init__(self, email, cachePath=None):
        """Constructor for NCBI taxonomy tree class.
        
        If cache path is given and not None the database at that path will be 
        used as persistent cache.
        """
        if cachePath is None:
            self._parent = TaxonomyParentMap(email, useCache=False)
            self.cached = False
        else:
            self._parent = CachedTaxonomyParentMap(cachePath, email)
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
    
    def lca(self, taxList):
        """
        Return the NCBI taxonomy ID of the lowest common ancestor in the NCBI
        taxonomy tree of a list of given taxonomy IDs. Note that there is 
        always at least on common ancestor (the root node)."""
        lin = [self.lineage(tax) for tax in set(taxList)]
        lvl = 0
        while all([lvl < len(l) for l in lin]) \
              and all([lin[0][lvl] == lin[i][lvl] for i in range(1, len(lin))]):
            lvl+=1
        return lin[0][lvl-1]
        
    def lcn(self, taxList):
        """Return the NCBI taxonomy ID of the lowest common node in the 
        NCBI taxonomy tree of a list of given taxonomy IDs.
        
        This is very similar to the LCA but for nodes that are part of a path
        from the root to one of them the lowest node will be choosen."""
        
        lin = [self.lineage(tax) for tax in set(taxList)]
        lvl = 0
        while lvl < max([len(l) for l in lin]):
            active = [l for l in lin if len(l) > lvl]
            if not all([active[0][lvl] == l[lvl] for l in active[1:]]):
                break
            lvl += 1
        return active[0][lvl-1]


