from urllib2 import urlopen
import json

class EolName2IdMap(dict):
    """Dictionary class to map names to EOL IDs
    
    Uses EOL web service 'search' to get information. Returns list of found IDs 
    for the name.
    """
    urlBase = "http://eol.org/api/search/1.0.json?q=%s&page=%i&"
    
    def __init__(self, indict={}):
        """Initalize with config: return exact results and filter by nothing"""
        self.config = {"exact": "true", "filter_by_taxon_concept_id": "false",
                       "filter_by_hierarchy_entry_id": "false",
                       "filter_by_string": ""}
        dict.__init__(self, indict)
        
    @property
    def urlTmpl(self):
        return self.urlBase + "&".join(["%s=%s" % (key, str(value)) \
                                        for key, value in self.config.items()])
    
    def __getitem__(self, name):
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            self[name] = self.eolCall(name)
            return dict.__getitem__(self, name)
            
    def eolCall(self, name, page=1):
        """Call the web service to get information on a name.
        
        Will call itself recursively if necessary to get the content of all 
        available pages."""
        print self.urlTmpl % (name.replace(" ","+"), page)
        resp = urlopen(self.urlTmpl % (name.replace(" ","+"), page))
        respData = json.loads(resp.read())
        if respData["totalResults"] == 0:
            raise KeyError()
        result = [r["id"] for r in respData["results"]]
        if "next" in respData:
            return result + self.eolCall(name, page+1)
        return result
        
        
class EolInterface(object):
    """Interface class to read information about an EOL page.
    
    Data is loaded once through EOL web service 'page' and can be queried for 
    information with specialized functions or accessed directly.
    """
    
    urlBase = "http://eol.org/api/pages/1.0/%i.json?"
    def __init__(self):
        self.config = {"images": 0, "vieos": 0, "sounds": 0,"maps": 0, 
                       "text": 2, "iucn": "false", "subjects": "overview",
                       "licenceses": "all", "details": "false", 
                       "common_names": "false", "synonyms": "false",
                       "refereneces": "false", "vetted": 0}
        self.data = None
    
    @property
    def urlTmpl(self):
        return self.urlBase + "&".join(["%s=%s" % (key, str(value)) \
                                        for key, value in self.config.items()])
    
    def getData(self, eolId):
        """Get data set of a certain EOL entry identified by its ID"""
        self.eolId = eolId
        print self.urlTmpl % self.eolId
        self.data = json.loads(urlopen(self.urlTmpl % self.eolId).read())
        
    def getTaxonId(self, source):
        """From the loaded data get the taxonomy ID from a certain source"""
        if self.data is None:
            raise ValueError("No data loaded")
        for taxon in self.data["taxonConcepts"]:
            if taxon["nameAccordingTo"] == source:
                return taxon["identifier"]
        raise KeyError("No information available from '%s'" % source)
        
if __name__ == "__main__":
    name2id = EolName2IdMap()
    e = name2id["Candida albicans"][0]
    print "eolID:", e
    i = EolInterface()
    i.getData(e)
    t = i.getTaxonId("Index Fungorum")
    print "Index Fugorum ID:", t
    
