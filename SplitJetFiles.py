import yoda
import itertools 
import os, optparse
import re
from sys import stdout


op = optparse.OptionParser()

##add args to name plot so we have meaningful legends later
##add axis labels by default


opts, args = op.parse_args()
name = os.path.splitext(args[0])

aodict = yoda.core.read(args[0], True)

yfiledict = {}

partlabels = ["Unlabeled", "Photon", "Gluon", "Bottom", "Charm", "Light",
        "Electron", "Muon", "Tau"]
jetalgs = ["MaxPt", "Reclustered", "Akt", "Kt"]
jettypes = ["Inclusive", "Jet0", "Jet1", "Jet2", "Jet3"]

labels = ["%s_%s_%s_" % x for x in itertools.product(jettypes,
    partlabels, jetalgs)]
label = re.compile("|".join(labels))

for k, ao in aodict.iteritems():
    m = label.search(k)
    if m:
        lab = m.group(0)
    else:
        continue

    ao.path = label.sub("", ao.path)
    if lab in yfiledict:
        yfiledict[lab].append(ao)
    else:
        yfiledict[lab] = [ao]

    continue

# loop over rho, aos
for lab, aos in yfiledict.iteritems():
    yoda.write(aos, "%s_%s.yoda" % (args[0].replace(".yoda", ""),
        lab.strip('_')))
