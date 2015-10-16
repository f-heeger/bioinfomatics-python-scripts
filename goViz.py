import sys
import urllib2
from itertools import chain
from collections import deque
import subprocess

from goOboParser import _parseGOOBO

localOboFile="/home/heeger/data/go/go.obo"


class GoNode(object):
    def __init__(self, gId, name):
        self.gId = gId
        self.name = name
        self.pvalue = None
        self.count = None
        self.size = None
        self.parents = []
        self.nameWidth=30
        
    def __repr__(self):
        return "GoNode(%s, %s)" % (self.gId, self.name)
        
    @property
    def wrappedName(self):
        if len(self.name) < self.nameWidth:
            return self.name
        lines = []
        for word in self.name.split():
            if len(lines) != 0 and len(lines[-1] + " " + word) < self.nameWidth:
                lines[-1] += " "+word
            else:
                lines.append(word)
        return "\n".join(lines)
        
    @property
    def gNumber(self):
        return self.gId[3:]

def readGoObo(localOboPath=None):
    goTree = {}
    if localOboPath is None:
        oboFile = urllib2.urlopen("http://purl.obolibrary.org/obo/go.obo")
    else:
        oboFile = open(localOboFile)

    terms = list(_parseGOOBO(oboFile))

    for term in terms:
        if "is_obsolete" in term and term["is_obsolete"] == "true":
            continue #skip obsolete
        tId = term["id"]
        goTree[tId] = GoNode(tId, term["name"])
        #goNamespace[tId] = term["namespace"]
    for term in terms:
        if "is_a" not in term:
            continue #root nodes
        if type(term["is_a"]) != type([]):
            #workaround for wired parsing behaviour
            parents = [term["is_a"]]
        else:
            parents = term["is_a"]
        for entry in parents:
            goTree[term["id"]].parents.append(goTree[entry.split(" ! ")[0]])
    return goTree


def readEnrichmentFile(path, goTree, sigLvl=0.05):
    sig = []
    enrichmentFile = open(path)
    next(enrichmentFile) #skip header
    for line in enrichmentFile:
        arr = line.split("\t")
        gId = arr[0]
        goTree[gId].pvalue = float(arr[1])
        goTree[gId].count = int(arr[4])
        goTree[gId].size = int(arr[5])
        if goTree[gId].pvalue < sigLvl:
            sig.append(gId)
    return sig

def createDot(goTree, sig):
    plotGo = {}
    for go in sig:
        plotGo[go] = None
        toDo = deque()
        toDo.append(go)
        while len(toDo) > 0:
            tGo = toDo.pop()
            try:
                for parent in goTree[tGo].parents:
                    if parent.gId not in plotGo:
                        toDo.append(parent.gId)
            except KeyError:
                sys.stderr.write("root node: %s %s\n" % (tGo, goTree[tGo].name))
                pass # reached tree root
            plotGo[tGo] = None


    maxValue = max([int(g.count) for g in goTree.values() if not g.count is None])
    scale = 255./maxValue

    dot = ["digraph go_graph{"]
    for go in plotGo.keys():
        par = {"id": "\"%s\"" % go[3:]}
        if go in sig:
            par["color"] = "red"
            par["penwidth"] = "4"
        if goTree[go].pvalue is None:
            par["label"] = "\"{g.gId}\n{g.wrappedName}\"".format(g=goTree[go])
        else:
            par["label"] = "\"{g.gId}\n{g.wrappedName}\n{g.count}/{g.size} p={g.pvalue:.6f}\"".format(g=goTree[go])
            col = 255-int(goTree[go].count*scale)
            colStr = hex(col)[2:].upper().rjust(2,"0")
            par["style"] = "filled"
            par["fillcolor"] = "\"#FF{0}{0}\"".format(colStr)
        dot.append("    %s [%s];" % (go[3:], ", ".join(["%s=%s" % i for i in par.items()])))
        try:
            for parent in goTree[go].parents:
                if parent.gId in plotGo:
                    dot.append("     %s -> %s;" % (go[3:], parent.gId[3:]))
        except KeyError:
            sys.stderr.write("root node: %s %s\n" % (go, gotree[go].name))
    dot.append("}")
    return dot

goTree = readGoObo(localOboFile)
sig = readEnrichmentFile(sys.argv[1], goTree)
dotList = createDot(goTree, sig)

dotCmd = ["dot", "-Tsvg"]
dotProc = subprocess.Popen(dotCmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
graphSvg, err = dotProc.communicate("\n".join(dotList))
sys.stderr.write(str(err))
while dotProc.returncode is None:
    out, err = dotProc.communicate()
    sys.stderr.write(str(err))
    graphSvg.append(out)
if dotProc.returncode != 0:
    raise(subprocess.CalledProcessError("Exitcode of dot not 0"))

tab = "<table id=\"go_table\">\n"
tab += "<thead><tr><th>GO term ID</th><th>GO term</th><th>differential expressed genes</th><th>p-value</th></tr></thead>\n"
tab += "<tbody>\n"
for go in goTree.values():
    if go.pvalue is None:
        continue
    tab+="<tr><td id=\"tab_{g.gNumber}\">{g.gId}</td><td>{g.name}</td><td>{g.count}/{g.size}</td><td>{g.pvalue}</td></tr>\n".format(g=go)
tab += "</tbody>\n</table>"


print("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" 
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>GO Term Enrichment</title>
    
    <!-- DataTables CSS -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.9/css/jquery.dataTables.css">
    <!-- jQuery -->
    <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-1.11.3.min.js"></script>
    <!-- DataTables -->
    <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.9/js/jquery.dataTables.js"></script>
    <script type="text/javascript" charset="utf8">
        $(document).ready( function () {
            $('#go_table').DataTable({
                "scrollY": "400px",
                "scrollCollapse": true,
                "paging": false
            });
            $('svg').attr("width", "100%%");
            $('svg').attr("height", "100%%");
            document.getElementsByTagName('svg')[0].setAttribute("viewBox", "0 0 1000 600");
            
            $("g.node").each(function () {$(this).attr("onclick", "searchInTable(this.id)")});
            $("#go_table tr").each(function () {$(this).attr("onclick", "searchInGraph($(this).find('td[id]').attr('id'))")});
        } );
        function scrollSvg(elem, dist, dir) {
            coord = elem.getAttribute("viewBox").split(" ")
            minX = parseFloat(coord[0]);
            minY = parseFloat(coord[1]);
            width = parseFloat(coord[2]);
            height = parseFloat(coord[3]);
            if (dir=="up") {
                minY -= dist;
            } else if (dir=="down") {
                minY += dist;
            } else if (dir=="left") {
                minX -= dist;
            } else if (dir=="right") {
                minX += dist;
            } else {
                alert("Unknown direction:"+dir);
            }
            elem.setAttribute("viewBox", [minX, minY, width, height].join(" "));
        }
        
        function zoomSvg(elem, amount) {
            var coord = elem.getAttribute("viewBox").split(" ")
            var minX = parseFloat(coord[0]);
            var minY = parseFloat(coord[1]);
            var width = parseFloat(coord[2]);
            var height = parseFloat(coord[3]);
            width = width/amount;
            height = height/amount;
           
            elem.setAttribute("viewBox", [minX, minY, width, height].join(" "));
        }
        
        function searchInTable(g) {
            $("#go_table .selected").toggleClass('selected');
            $("tr").has("#tab_"+g).toggleClass('selected');
            $('.dataTables_scrollBody').animate({
                scrollTop:  $("tr").has("#tab_"+g).prop('offsetTop') - $('.dataTables_scrollBody').height()/2
            }, 400);
        }
        
        function searchInGraph(g) {
            var x = $("#"+g.slice(4)+" ellipse").attr("cx");
            var y = $("#"+g.slice(4)+" ellipse").attr("cy");
            
            var svg=document.getElementsByTagName('svg')[0]
            
            var coord = svg.getAttribute("viewBox").split(" ")
            var width = parseFloat(coord[2]);
            var height = parseFloat(coord[3]);
            svg.setAttribute("viewBox", [x-width/2, y-height/2, width, height].join(" "));
        }
    </script>
  </head>
  <body>
    <div style="width:49%%;float:left;">
        %(svg)s
        <div>
            <button type="button", onclick="scrollSvg(document.getElementsByTagName('svg')[0], 100, 'left')">left</button>
            <button type="button", onclick="scrollSvg(document.getElementsByTagName('svg')[0], 100, 'right')">right</button>
            <button type="button", onclick="scrollSvg(document.getElementsByTagName('svg')[0], 100, 'up')">up</button>
            <button type="button", onclick="scrollSvg(document.getElementsByTagName('svg')[0], 100, 'down')">down</button>
            <button type="button", onclick="zoomSvg(document.getElementsByTagName('svg')[0], 1.2)">+</button>
            <button type="button", onclick="zoomSvg(document.getElementsByTagName('svg')[0], 0.8)">-</button>
        </div>
    </div>
    <div style="width:49%%;float:right;">
    %(table)s
    </div>
  </body>
</html>""" % {"svg": graphSvg, "table": tab}
)

