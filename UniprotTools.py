import urllib
import urllib2
import re

#based on Uniprot retrieval example on 
# http://www.uniprot.org/help/programmatic_access

from MultiCachedDict import MultiCachedDict, SqliteCache, NotWritableError

class CachedUniprotIdMap(MultiCachedDict):
    def __init__(self, dbpath, source="ACC", target="P_REFSEQ_AC", retry=0, 
                 delay=1, contact=None, returnNone=False):
        
        dbMap = SqliteCache(dbpath)
        uniprotMap = UniprotIdMap(source, target, retry, delay, contact, 
                                  returnNone)
        
        MultiCachedDict.__init__(self, None, [dbMap, uniprotMap])

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
                    time.sleep(delay)
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

        data = urllib.urlencode(params)
        request = urllib2.Request(url, data)
        if self.contact:
            request.add_header('User-Agent', 'Python %s' % self.contact)
        return urllib2.urlopen(request)

    def readResponse(self, key):
        lines = self.response.readlines()
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

    """
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

        request = urllib2.Request(url)
        if self.contact:
            request.add_header('User-Agent', 'Python %s' % self.contact)
        return urllib2.urlopen(request)

    def readCazyInfo(self):
        if self.response is None:
            raise ValueError("No response available.")
        lines = self.response.readlines()
        if len(lines) < 2:
            raise KeyError("Uniprot web service did not return information for "
                           "ID: '%s'" % self.uid)
        l=0
        #skip lines until one that starts with "DR    CAZy"
        while l < len(lines) and not lines[l].startswith("DR   CAZy"):
            l+=1
        if l >= len(lines):
            # if the while loop run all the way ne cazy info was found
            raise KeyError("Uniprot web service did not return CAZy information"
                           " for ID: '%s'" % self.uid)
        _, cId, cName = lines[l].strip(".\n").split(";")
        return {"id": cId, "name": cName}
    
    def readGoInfo(self):
        if self.response is None:
            raise ValueError("No response available.")
        lines = self.response.readlines()
        if len(lines) < 2:
            raise KeyError("Uniprot web service did not return information for "
                           "ID: '%s'" % self.uid)
        l=0
        #skip lines until the first that starts with "DR"
        while lines[l][:2] != "DR":
            l+=1
        goInfo = []
        #match lines while they start with DR
        while lines[l][:2] == "DR":
            m=re.match("DR   GO; (?P<goId>GO:\d{7}); "
                       "(?P<aspect>[FPC]):(?P<function>[^;]+); "
                       "(?P<evidenceCode>[A-Z]+):(?P<evidenceSource>[^.]+)", 
                       lines[l])
            if not m is None:
                goInfo.append(m.groupdict())
            l+=1
        return goInfo
        
    
        
if __name__ == "__main__":
    import sys
    goTestData = [("N4UXJ0",[{'function': 'ATP binding', 
                              'evidenceCode': 'IEA', 
                              'aspect': 'F', 
                              'goId': 'GO:0005524', 
                              'evidenceSource': 'InterPro'}, 
                             {'function': 'protein kinase activity', 
                              'evidenceCode': 'IEA', 
                              'aspect': 'F', 
                              'goId': 'GO:0004672', 
                              'evidenceSource': 'InterPro'}
                             ]
                  ),
                ]
    cazyTestData = [("A2QYE5", {"id": "GH28", 
                                "name": "Glycoside Hydrolase Family 28"}),
                    ("A0A077W732", None),
                    ("Q39Z37", {"id": "GT2", 
                                "name": "Glycosyltransferase Family 2"}),
                   ]

    print("Testing Uniprot ID map")
    umap = UniprotIdMap("P_GI","ACC", retry=3, contact="fheeger@mi.fu-berlin.de")
    print(umap["477524024"])

    print("Testing Uniprot Interface")
    uniprot = UniprotInterface(contact="fheeger@mi.fu-berlin.de")
    print("Interface object build successful")
#    print("Testing go info")
#    for uid, result in goTestData:
#        uniprot.getData(uid)
#        print uniprot.readGoInfo()
    print("Testing CAZy info")
    for uid, result in cazyTestData:
        uniprot.getData(uid)
        try:
            print uniprot.readCazyInfo()
        except KeyError as e:
            print e
