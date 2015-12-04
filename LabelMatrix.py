from __future__ import print_function
import yoda
import itertools 
import os, optparse
import re
from sys import stdout

binning = list(enumerate([
    ("none", [0]),
    ("$g$", [21]),
    ("$q$", [-3, -2, -1, 1, 2, 3]),
    ("$c$", [-4, 4]),
    ("$b$", [-5, 5]),
    ("$\\gamma$", [22]),
    ("$e^-$", [11]),
    ("$e^+$", [-11]),
    ("$\\mu^-$", [13]),
    ("$\\mu^+$", [-13]),
    ("$\\tau$", [15])
    ]))


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


def histToMatrix(hobj, nbins):
    hMatrix = yoda.Histo2D(nbins, 0, nbins, nbins, 0, nbins,
            hobj.path + "_matrix", "jet label matrix")

    binningAnnotation = '\t'.join(map(
            lambda (n, (lab, _)): "%02.2f\t%s" % ((n+0.5), lab),
            binning[:nbins]))

    hMatrix.setAnnotation("XCustomMajorTicks", binningAnnotation)
    hMatrix.setAnnotation("YCustomMajorTicks", binningAnnotation)
    hMatrix.setAnnotation("PlotTickLabels", "1")
    hMatrix.setAnnotation("PlotXMajorTicks", "0")
    hMatrix.setAnnotation("ZLabel", "arbitrary")
    hMatrix.setAnnotation("XLabel", hobj.annotation("XLabel"))
    hMatrix.setAnnotation("YLabel", hobj.annotation("YLabel"))
    hMatrix.setAnnotation("ZCustomMajorTicks", "0.5\t$ $")


    hobj.normalize(1.0)
    for b in hobj.bins:
        xbin, ybin = histBinToMatrixBin(b)
        hMatrix.fill(xbin, ybin, b.sumW)
        continue

    return hMatrix


def roundText(n):
    if n > 0.1:
        return "%0.1f" % n
    else:
        return "-"


def matrixToLatex(m):
    labels = m.annotation("XCustomMajorTicks").split('\t')[1::2]
    xlab = m.annotation("XLabel")
    ylab = m.annotation("YLabel")

    l = len(labels)

    valuestable = [ map(lambda b: roundText(100*b.volume),
                m.bins[l*i:l*(i+1)] ) for i in range(len(labels))]


    fstrow = "    \\multicolumn{2}{|c|}{ } & \\multicolumn{%s}{c|}{%s} \\\\" % (l, ylab)
    secrow = "    \\multicolumn{2}{|c|}{ } & %s \\\\" % (" & ".join(labels))
    trdrow = "    \\multirow{%s}{*}{\\rotatebox{90}{%s}} & %s & %s \\\\" \
            % (l, xlab, labels[0], " & ".join(valuestable[0]) )

    restrows = ["     & %s & %s \\\\" % (labels[i], " & ".join(valuestable[i]) ) \
            for i in range(1, len(labels)) ]

    rows = [fstrow, secrow, trdrow] + restrows

    rows.insert(0, "    \\hline")
    rows.insert(2, "    \\cline{3-%s}" % (l+2))
    rows.insert(4, "    \\hline")
    rows.append("    \\hline")
    rowstext = '\n'.join(rows)


    tex = """\
\\begin{table}[t]
  \\centering
  \\begin{tabular}{%s}
%s
  \\end{tabular}
\\end{table}""" % ("|c|c|" + 'c'*l + "|", rowstext)

    return tex


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
    texs = []
    for k, ao in aodict.iteritems():
        m = label.search(k)
        if m:
            lab = m.group(0)
        else:
            continue

        print("found label " + lab); stdout.flush()

        m = histToMatrix(ao, 6)
        aos.append(m)
        texs.append(m.path + "\n" + matrixToLatex(m))
        continue

    yoda.write(aos, args[1])

    open(args[2], 'w').write("\n\n".join(texs))

    return 0

if __name__ == '__main__':
    main()
