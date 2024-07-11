#!/usr/bin/python3

import pandas as pd

# from tabulate import tabulate
import os, sys
import gzip

# Get the samples file to parse.
smpfile = snakemake.params[0]
# smpfile = sys.argv[1]

# sys.stdout = open(snakemake.output[0], 'w')
# print ("Summarizing", smpfile, "\n")

# Read samples and units table
samples = pd.read_csv(smpfile, sep="\t", dtype=str)
# samples = pd.read_csv(smpfile, sep='\t', dtype=str).set_index(["sample", "unit"], drop=False)
# samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index

# Change to the dir of the sample file, so that we can find sample files relative to it.
olddir = os.getcwd()
os.chdir(os.path.dirname(smpfile) if os.path.dirname(smpfile) else ".")

# Add the file sizes as
samples["fq1_size"] = [
    (
        os.path.getsize(str(f))
        if os.path.isfile(str(f))
        else (os.path.getsize("../" + str(f)) if os.path.isfile("../" + str(f)) else 0)
    )
    for f in samples["fq1"]
]
samples["fq2_size"] = [
    (
        os.path.getsize(str(f))
        if os.path.isfile(str(f))
        else (os.path.getsize("../" + str(f)) if os.path.isfile("../" + str(f)) else 0)
    )
    for f in samples["fq2"]
]


def count_fastq_gzip_lines(fastqfile):
    count = 0
    with gzip.open(fastqfile, "rt") as f:
        for line in f:
            count += 1
    return int(count / 4)


samples["fq1_seqs"] = [
    (
        count_fastq_gzip_lines(str(f))
        if os.path.isfile(str(f))
        else count_fastq_gzip_lines("../" + str(f))
    )
    for f in samples["fq1"]
]
samples["fq2_seqs"] = [
    (
        count_fastq_gzip_lines(str(f))
        if os.path.isfile(str(f))
        else (count_fastq_gzip_lines("../" + str(f)) if os.path.isfile("../" + str(f)) else 0)
    )
    for f in samples["fq2"]
]

os.chdir(olddir)
samples.to_csv(snakemake.output[0], index=False)
# print (samples)
# print( tabulate( samples ))
