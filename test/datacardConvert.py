#!/usr/bin/env python

#import sys, os
#sys.path.insert(0,f"{os.environ['HOME']}/git/combine2pyhf" )
#sys.path.insert(0,f"{os.environ['HOME']}/git/" )

import json
from optparse import OptionParser
from HiggsAnalysis.CombinedLimit.DatacardConverter import *

def main():
    parser = OptionParser(usage="usage: %prog [options] datacard.txt -o output \nrun with --help to get list of options")
    parser.add_option("-o", "--out", default="./output", type="string", help="output file")
    parser.add_option("--bbl", action="store_true", help="use Barlow-Beeston lite approach for statistical uncertainties")
    parser.add_option("--normshape", action="store_true", help="split shape uncertainties into pure shape and normalization components")
    parser.add_option("--prune", action="store_true", help="remove shape systematics with no effect")
    parser.add_option("--neg", action="store_true", help="remove negative predictions in bins")
    (options, args) = parser.parse_args()
    if options.out.endswith ( ".json" ):
        options.out = options.out.replace(".json","")

    if len(args) == 0:
        parser.print_usage()
        exit(1)
                
    opts = type("opts", (object,), dict(bin=True, noJMax=False, stat=False, nuisancesToExclude=[], allowNoSignal=True, allowNoBackground=True))
        
    file = open(args[0], "r")
    convertCard(args[0], file, opts, options.out, options.bbl, options.normshape, options.prune, options.neg)

if __name__ == "__main__":
    main()
