#!/usr/bin/env python

import json
from optparse import OptionParser
from HiggsAnalysis.CombinedLimit.DatacardConverter import *

parser = OptionParser(usage="usage: %prog [options] datacard.txt -o output \nrun with --help to get list of options")
parser.add_option("-o", "--out", default="output", type="string", help="output file")
parser.add_option("--bbl", action="store_true", help="use Barlow-Beeston lite approach for statistical uncertainties")
(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_usage()
    exit(1)
            
opts = type("opts", (object,), dict(bin=True, noJMax=False, stat=False, nuisancesToExclude=[], allowNoSignal=True, allowNoBackground=True))
    
file = open(args[0], "r")
convertCard(args[0], file, opts, options.out, options.bbl)
