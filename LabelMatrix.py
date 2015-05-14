from __future__ import print_function
import yoda
import itertools 
import os, optparse
import re
from sys import stdout

binning = list(enumerate([
    ("unlabeled", [0]),
    ("$g$", [21]),
    ("$q$", [1, 2, 3]),
    ("$\\bar q$", [-1, -2, -3]),
    ("$c$", [4]),
    ("$\\bar c$", [-4]),
    ("$b$", [5]),
    ("$\\bar b$", [-5]),
    ("$\\gamma$", [22]),
    ("$e^-$", [11]),
    ("$e^+$", [-11]),
    ("$\\mu^-$", [13]),
    ("$\\mu^+$", [-13]),
    ("$\\tau$", [15])
    ]))

binningAnnotation = '\t'.join(map(
        lambda (n, (lab, _)): "%02.2f\t\"%s\"" % ((n+0.5), lab),
        binning))


def histBinToMatrixBin(b):
    mx = b.xMid
    my = b.yMid
    xbin = -1
    ybin = -1
    
    for binNum, (l, pids) in binning:
        for pid in pids:
            if (xbin < 0) and (pid-0.5) < mx and mx < (pid+0.5):
                xbin = binNum

            if (ybin < 0) and (pid-0.5) < my and my < (pid+0.5):
                ybin = binNum

            continue

            # we've already found both bin numbers...
            if xbin >= 0 and ybin >= 0:
                break


        continue

    return (xbin, ybin)


def histToMatrix(hobj):
    hMatrix = yoda.Histo2D(14, 0, 14, 14, 0, 14,
            hobj.path + "_matrix", "jet label matrix")

    hMatrix.setAnnotation("XCustomMajorTicks", binningAnnotation)
    hMatrix.setAnnotation("YCustomMajorTicks", binningAnnotation)
    hMatrix.setAnnotation("PlotTickLabels", "1")


    for b in hobj.bins:
        xbin, ybin = histBinToMatrixBin(b)
        hMatrix.fill(xbin, ybin, b.sumW)
        continue

    # TODO
    # bin labels

    return hMatrix

def main():
    op = optparse.OptionParser()

    ##add args to name plot so we have meaningful legends later
    ##add axis labels by default


    opts, args = op.parse_args()
    name = os.path.splitext(args[0])

    aodict = yoda.core.read(args[0], True)

    yfiledict = {}

    jetcombos = ["%sLabVs%sLab" % x for x in itertools.combinations(["Akt", "Kt", "Reclustered", "MaxPt"], 2)]
    jettypes = ["Inclusive", "Jet0", "Jet1", "Jet2", "Jet3"]

    labels = ["%s_%s" % x for x in itertools.product(jettypes,
        jetcombos)]

    label = re.compile("|".join(labels))

    aos = []
    for k, ao in aodict.iteritems():
        m = label.search(k)
        if m:
            lab = m.group(0)
        else:
            continue

        print("found label " + lab); stdout.flush()

        aos.append(histToMatrix(ao))
        continue

    yoda.write(aos, args[1])

    return 0

if __name__ == '__main__':
    main()
