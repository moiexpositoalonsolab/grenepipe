#!/usr/bin/python

import sys
import os
import os.path
from tabulate import tabulate
from operator import itemgetter

# Usage: Either call with no arguments, in which case the current directory is scanned,
# or provide a single argument that is the run dir or the benchmark dir within a run dir to scan.

# Get the directory of benchmarks to scan.
benchdir = "."
if len(sys.argv) >= 2:
    benchdir = sys.argv[1]
if os.path.isdir( os.path.join( benchdir, "benchmarks" )):
    benchdir = os.path.join( benchdir, "benchmarks" )
print "Summarizing", benchdir, "\n"

# Scan the given directory for all `.bench.log` files, which is what our rules produce.
entries = list()
for dirpath, dirnames, filenames in os.walk( benchdir ):
    for filename in [f for f in filenames if f.endswith(".bench.log")]:
        path = os.path.join(dirpath, filename)
        name = path[2:-4]
        # print name

        with open(path) as fp:
            # Check first line. Hardcoded to the current output of Snakemake benchmark.
            # If this fails, we cannot parse the file anyway, and have to adapt the script.
            line = fp.readline()
            if line != "s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load\n":
                print( "Job " + name + " invalid!" )
                break

            # Parse all remaining lines, and add them to the list.
            line = fp.readline()
            while line:
                entry = line.split("\t")
                entry.insert( 0, name )
                entries.append(entry)
                line = fp.readline()

# Print a table of all job benchmarks.
entries = sorted( entries, key=itemgetter(0) )
headers = ["job", "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load"]
print( tabulate( entries, headers=headers ))
