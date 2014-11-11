import urllib
import urllib2
import re

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

    def readGoInfo(self):
        if self.response is None:
            raise ValueError("No response available.")
        lines = self.response.readlines()
        if len(lines) < 2:
            raise KeyError("Uniprot web service did not return information for "
                           "ID: '%s'" % self.uid)
        l=0
        #skip lines until the first that starts with DR
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

    print("Testing Uniprot Interface")
    uniprot = UniprotInterface(contact="fheeger@mi.fu-berlin.de")
    print("Interface object build successful")
    print("Testing go info")
    for uid, result in goTestData:
        uniprot.getData(uid)
        print uniprot.readGoInfo()
