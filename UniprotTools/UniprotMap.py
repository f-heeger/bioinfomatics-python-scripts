#based on uniprot retirval sample on 
# http://www.uniprot.org/help/programmatic_access

import urllib
import urllib2

from MultiCachedDict import *


class UniprotIdMap(object):
    writable = False #this is not a writable dict

    """Base class for a dictonary based on uniprot
    """
    def __init__(self, source="ACC", target="P_REFSEQ_AC", retry=0, 
                 delay=1, contact=None):
        self.source = source
        self.target = target
        self.retry = retry
        self.delay = delay
        self.contact = contact
        
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
            
        return self.readResponse()
            
    def requestFunction(self, key):
        url = 'http://www.uniprot.org/mapping/'

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

    def readResponse(self):
        lines = self.response.readlines()
        if len(lines) > 2:
            raise ValueError("Response from Uniprot has more than teo lines.")
        # line 0 is the header
        key, value = lines[1].strip().split()
        self.response = None
        return value
        
if __name__ == "__main__":
    import sys
    umap = UniprotIdMap("P_GI","ACC", retry=3, contact="fheeger@mi.fu-berlin.de")
    umap["477524024"]
