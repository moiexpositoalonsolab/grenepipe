#!/usr/bin/python2

import pandas as pd
# from tabulate import tabulate
import os, sys
import gzip

# Usage: Call with no arguments to load `samples.tsv` from the current dir, or provide
# a single argument to a `samples.tsv` file or directory containing such a file.

# Get the samples file to parse.
smpfile = sys.argv[1]
print "Summarizing", smpfile, "\n"

# Change to the dir of the sample file, so that we can find sample files relative to it.
os.chdir(os.path.dirname(smpfile))

# Read samples and units table
samples = pd.read_csv(smpfile, sep='\t', dtype=str)
# samples = pd.read_csv(smpfile, sep='\t', dtype=str).set_index(["sample", "unit"], drop=False)
# samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index

# Add the file sizes as
samples["fq1_size"] = [
    os.path.getsize(str(f))
    if os.path.isfile(str(f))
    else os.path.getsize("../" + str(f))
    for f in samples["fq1"]
]
samples["fq2_size"] = [
    os.path.getsize(str(f))
    if os.path.isfile(str(f))
    else (
        os.path.getsize("../" + str(f)) if os.path.isfile("../" + str(f)) else 0
    ) for f in samples["fq2"]
]

def count_fastq_gzip_lines(fastqfile):
    count = 0
    with gzip.open(fastqfile,'rt') as f:
        for line in f:
            count += 1
    return int(count / 4)

samples["fq1_seqs"] = [
    count_fastq_gzip_lines(str(f))
    if os.path.isfile(str(f))
    else count_fastq_gzip_lines("../" + str(f))
    for f in samples["fq1"]
]
samples["fq2_seqs"] = [
    count_fastq_gzip_lines(str(f))
    if os.path.isfile(str(f))
    else (
        count_fastq_gzip_lines("../" + str(f)) if os.path.isfile("../" + str(f)) else 0
    ) for f in samples["fq2"]
]

print samples
# print( tabulate( samples ))
