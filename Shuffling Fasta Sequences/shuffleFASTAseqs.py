#!/usr/bin/env python3

# functions for network analysis of sars-cov2 interactome
import os
import random
import pandas as pd
import numpy as np

# utility fn for splitting string into chunks (ie windows of size 10)
def chunkstring(string, length):
    return list((string[0+i:length+i] for i in range(0, len(string), length)))

# utility fn for converting list to string
def list2string(s):
    # initialize an empty string
    str1 = ""
    # traverse in the string
    for ele in s:
        str1 += ele
        # return string
    return str1


# takes FASTA file of Krogan et al. proteins and returns shuffled fasta
def shuffasta(infile, outfile, numshufs=10000):
    try:
        os.remove(outfile)  # delete old file so we start fresh
    except FileNotFoundError as err:
        pass

    f2 = open(outfile, "w")  # open outfile for writing

    # open file
    with open(infile, "r") as f:
        lineslist = f.readlines()

        # iterate through lines
        for i in lineslist:
            # look for header lines
            if i.startswith(">"):  # new header, write windows to file and
                f2.write(i)

            else:
                splitline = chunkstring(i, 10)

                # iterate through windows and shuffle them
                shufline = []
                endflag = False
                for w in splitline:
                    # check if line contains endline character
                    if "\n" in w:
                        window = w[:-1]
                        endflag = True
                    else:
                        window = w
                    if len(window) != 0:
                        for n in range(numshufs):
                            window = list(window)
                            # get two random indices from window
                            r1 = random.randint(0, len(window)-1)
                            r2 = random.randint(0, len(window)-1)
                            # get those characters
                            c1 = window[r1]
                            c2 = window[r2]
                            # swap them
                            window[r1] = c2
                            window[r2] = c1
                            window = list2string(window)
                    else: window = ""
                    shufline.append(window)
                    if endflag: shufline.append("\n")
                shufline = list2string(shufline)
                f2.write(shufline)

if __name__ == "__main__":
    # this is hard coded for our files, switch them up or import functions as necessary for your workflow
    shuffasta("Krogan_Protein_database_REDUCED.fasta", "Krogan_Protein_database_REDUCED_SHUFFLED.fasta")