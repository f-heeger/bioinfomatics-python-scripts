import re
from collections import deque

class GoTerm(object):
    """Class to hold the OBO 1.2 entry for a GO term.
    
    The following information that can be in an OBO term will not be hold here:
    comment, intersection_of, union_of, disjoint_from, created_by, creation_date
    """
    def __init__(self, gId=None, namespace=None, name=None, definition=None, 
                 obsolete=False):
        self.gId = gId
        self.namespace = namespace
        self.name = name
        self.definition = definition
        self.obsolte = obsolete
        self.xref = []
        self.subset = []
        self.synonym = []
        self.is_a = set([])
        self.relationship={"part_of" : set([]),
                           "has_part" : set([]),
                           "regulates" : set([]),
                           "positively_regulates" : set([]),
                           "negatively_regulates" : set([]),
                           "occurs_in" : set([]),
                           "happens_during" : set([]),
                           "ends_during" : set([]),
        }
        self.consider = set([])
        self.replaced_by = []
        self.alt_id = set([])
        
    def __repr__(self):
        return 'GoTerm("%s", "%s", "%s")' % (self.gId, self.namespace, self.name)

class GoNode(object):
    """Class to represent the DAG node of a GO Term.
    
    This is basically an additional layer around GoTerm objects to hold direct 
    references to the parents' GoNodes. Only the "is_a" and "part_of" realtions
    in GO are considered to be child-parent relationships for the sake of this 
    tree. They are stored seperatly, but can be retived together with the parent
    property.
    """
    def __init__(self, term):
        self.term = term
        self.is_a_parents = []
        self.part_of_parents = []
        
    @property
    def parents(self):
        return self.is_a_parents + self.part_of_parents
        
    def __repr__(self):
        return "GoNode(%s)" % str(self.term)
        
    def ancestors(self):
        """Function retrieve all ancestors of this Node, ie. its parents and
        their parents and so on. 
        """
        aTerms = set([])
        toDo = deque()
        toDo.append(self)
        while len(toDo) > 0:
            tNode = toDo.pop()
            aTerms.add(tNode.term.gId)
            if len(tNode.parents)==0:
                print("root node: %s" % str(tNode))
            for parent in tNode.parents:
                if parent.term.gId not in aTerms:
                    toDo.append(parent)
        return aTerms

def parseGoObo(inStream, ignore_obsolete=False, log=None):
    """Parse a GO OBO 1.2 file from inStream. This is not a fully featured OBO 
    parser!
    
    If ignore_obsolete is True all obsolete terms will be ignored. If log is 
    given, it has to be a file like object (with a write method) and verbose 
    status messages will be written to it.
    """
    term = None
    count = 0
    for line in inStream:
        if len(line.strip()) == 0 or line[0] == "!":
            continue #skip empty lines and comment lines
        if line[0] == "[":
            if not term is None:
                count += 1
                yield term
            if line.strip() == "[Term]":
                #entering a new term
                term = GoTerm()
            else:
                #entering any other block
                term = None
        elif not term is None:
            if "!" in line:
                data, comment = line.strip().split("!", 1)
            else:
                data = line
            key, value = data.strip().split(":", 1)
            if key == "id":
                term.gId = value.strip()
            elif key == "namespace":
                term.namespace = value.strip()
            elif key == "name":
                term.name = value.strip()
            elif key == "def":
                defRegex = "\"(.*)\" +\[(.+)\] *({.*})?"
                definition, refs, trailing = re.match(defRegex, 
                                                      value.strip()
                                                      ).groups()
                term.definition = definition
                if log:
                    log.write("Ignoring trailing modifiers:\n%s\n" 
                               % trailing)
            elif key == "is_obsolete" and value.strip() == "true":
                term.obsolete = True
                if ignore_obsolete:
                    if log:
                        log.write("Ignoring obsolete term:\n%s\n" % str(term))
                    term = None
            elif key == "xref":
                term.xref.append(value.strip())
            elif key == "subset":
                term.subset.append(value.strip())
            elif key == "synonym":
                term.synonym.append(value.strip())
            elif key == "is_a":
                term.is_a.add(value.strip())
            elif key == "relationship":
                relType, target = value.strip().split()
                term.relationship[relType].add(target)
            elif key == "consider":
                term.consider = value.strip()
            elif key == "replaced_by":
                term.replaced_by = value.strip()
            elif key == "alt_id":
                term.alt_id.add(value.strip())
            else:
                if log:
                    log.write("Ignored key-value pair: %s - %s\n" % (key, value))
    if not term is None:
        count += 1
        yield term
    if log:
        log.write("Parsed %i GO terms.\n" % count)

def parseGoOboDict(inStream, ignore_obsolete=False, log=None):
    """Returns a dictonary with GO term IDs as keys and GoTerm objects of the 
    corresponding GO terms.
    
    Uses the parseGoObo function. Parameters are passed through to it and 
    have the same meaning.
    """
    goDict = {}
    for g in parseGoObo(inStream, ignore_obsolete, log):
        goDict[g.gId] = g
        for a in g.alt_id:
            goDict[a] = g
    return goDict
    
def createGoTree(goDict):
    """Creates a dictonary of correctly linked GoNodes to represent the GO DAG.
    
    The dictonary allows you to retrieve the GoNode according to a GO ID and
    travers the GO DAG bottom up along the parnets attribute.
    """
    goTree = {}
    for gId, term in goDict.items():
        if gId != term.gId:
            #this is an entrie because of an alt_id -> do not create a node
            continue
        goTree[term.gId] = GoNode(term)
        for a in term.alt_id:
            goTree[a] = goTree[term.gId]
    
    for tNode in goTree.values():
        for i in tNode.term.is_a:
            tNode.is_a_parents.append(goTree[i])
        for i in tNode.term.relationship["part_of"]:
            tNode.part_of_parents.append(goTree[i])
    return goTree

if __name__ == "__main__":
    import sys
    from urllib2 import urlopen
    d=parseGoOboDict(urlopen("http://purl.obolibrary.org/obo/go.obo"), True, open("gopareser.log","w"))
    t=createGoTree(d)
