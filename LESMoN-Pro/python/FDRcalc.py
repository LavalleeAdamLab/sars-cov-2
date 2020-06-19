#!/usr/bin/env python3

# functions for network analysis of sars-cov2 interactome
import os
import random
import pandas as pd
import numpy as np

def getfdrs(realfile, shuffile, outfile):

    # read in tsv files
    f = open(realfile, "r")
    realdf = pd.read_csv(f, sep="\t")
    f.close()
    f = open(shuffile, "r")
    shufdf = pd.read_csv(f, sep="\t")
    f.close()

    pylist = []
    shufpylist = []

    # create xlist using adaptive iteration so more thresholds are generated as values get lower
    # xlist = list(np.arange(0.1, 0, -0.0001))
    xlist = []
    thresh = float(1)

    while thresh > 0:
        if thresh > 0.1:
            thresh = thresh - 0.05
        elif thresh > 0.05:
            thresh = thresh - 0.01
        elif thresh > 0.001:
            thresh = thresh - 0.0001
        else:
            thresh = thresh - 0.000001
        xlist.append(thresh)


    for i in xlist:
        pylist.append(len(list(filter(lambda x: x < i, realdf["P-value"]))))
        shufpylist.append(len(list(filter(lambda x: x < i, shufdf["P-value"]))))

    # scratch thoughts:
    # pylist[i] = for threshold in xlist[i], the number of significant hits
    # shufpylist = "" for shuffled
    # xlist = pval thresholds

    # normalize based on total number of motifs we got
    # TODO unhardcode this you lazy sob
    normpylist = [x / 48128026 for x in pylist]
    normshufpylist = [x / 30787618 for x in shufpylist]


    # create file of actual FDRs at all thresholds
    with open(outfile, "w+") as fout:
        fout.write("threshold\tshuffledPositives\trealPositives\tshuffledNormalizedPositives\t" +
                   "realNormalizedPositives\tfdr\n")
        min = 1  # for monotonic transformation
        for i in range(len(normpylist)):
            try:
                print("cutoff: " + str(xlist[i]))
                print("Shuffled positives: " + str(shufpylist[i]))
                print("Real positives: " + str(pylist[i]))
                print("Shufnorm positives: " + str(normshufpylist[i]))
                print("Realnorm positives: " + str(normpylist[i]))
                fdr = normshufpylist[i] / normpylist[i]

                # monotonic transformation
                if fdr < min:
                    min = fdr
                    monofdr = min
                else:
                    monofdr = min

                print("FDR: " + str(fdr))
                print("Monotonic FDR: " + str(monofdr) + "\n")

                fout.write(str(xlist[i]) + "\t" + str(shufpylist[i]) + "\t" + str(pylist[i]) + "\t" +
                           str(normshufpylist[i]) + "\t" + str(normpylist[i]) + "\t" + str(fdr) + "\t" +
                           str(monofdr) + "\n")
            except ZeroDivisionError:
                print("end of divisible values")
                break

if __name__ == "__main__":
    # this is hard coded for our files, switch them up or import functions as necessary for your workflow
    getfdrs("krogan_MotifAnnotations_real_s10X5.txt", "krogan_MotifAnnotations_shuffled_s10X5.txt", "covFDRsFULL.tsv")
