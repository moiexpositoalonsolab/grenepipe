#!/usr/bin/env python3

import sys, os
import re
from glob import glob
from termcolor import colored

# =================================================================================================
#     Usage and Description
# =================================================================================================

# Generate the samples table to be used as grenepipe input that lists paired end read samples.
# Usage: ./generate-table.py <directory-with-fastq-files> [<samples.tsv>]

# =================================================================================================
#     Settings
# =================================================================================================

# Valid extensions for fastq files that we are looking for.
extensions = [ ".fastq", ".fa", ".fq", ".fastq.gz", ".fa.gz", ".fq.gz" ]

# If no output file name is given as second command line argument, use this one.
outfile = "samples.tsv"

# =================================================================================================
#     Functions
# =================================================================================================

def get_fastq_files(indir):
    # Get all files in the given directory that match our file extensions, and sort them,
    # so that adjacent entries are potential paired files.
    seqfiles = []
    for dirpath, subdirs, filenames in os.walk(indir):
        # print(filenames)
        for f in filenames:
            path = os.path.abspath(os.path.join(dirpath, f))
            if path.endswith( tuple( extensions )):
                seqfiles.append(path)
    seqfiles.sort()
    return seqfiles

def match_pe_files(s1, s2):
    # Find the single position at which the two strings differ,
    # or -1 if they differ in more than one position
    # (or are the same string, but that should not happen in our case).
    pos = -1
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 != c2:
            if pos != -1:
                return -1
            else:
                pos = i

    # In case that they do _not_ differ in exactly one position, we are done.
    if pos == -1:
        return -1

    # Now check whether that position is `1` and `2` respectively,
    # for forward and backward strand, and return the result.
    # Any value > -1 means, yes, those are forward and backward strand.
    # As our strings come in sorted, we only need to check one permutation.
    if s1[pos] == '1' and s2[pos] == '2':
        return pos
    else:
        return -1

def get_sample_name(tup):
    # Expects one element of the get_mates function output list,
    # and cleans up the beginning and end of the file names to get a clean sample name.
    return re.split( '/', tup[0][:tup[2]] )[-1].strip( '-_R' )

def get_mates(seqfiles):
    # Walk through that list and find pairs of reads by checking for single character differences,
    # where one of the file names is `R1`, and the other is `R2`.
    # print(seqfiles)
    result = []
    while len(seqfiles):
        seq1 = seqfiles.pop(0)
        # print(seq)
        found_mate = False
        for j in range(len(seqfiles)):
            seq2 = seqfiles[j]
            match_pos = match_pe_files(seq1, seq2)
            if match_pos > -1:
                # print(seq1, seqfiles[j])
                # Print the matches, with the difference in red, so that the user
                # can check that the pairing is correct.
                print()
                print("Match:")
                print(seq1[0:match_pos] + colored(seq1[match_pos], "red") + seq1[match_pos+1:])
                print(seq2[0:match_pos] + colored(seq2[match_pos], "red") + seq2[match_pos+1:])
                print("--> " + get_sample_name(( seq1, seq2, match_pos )))
                result.append(( seq1, seq2, match_pos ))
                seqfiles.pop(j)
                found_mate = True
                break
        if not found_mate:
            print("Did not find mate for " + seq1)
    return result

def yes_or_no(question):
    while True:
        reply = str(input(question+' (y/n): ')).lower().strip()
        if len(reply):
            if reply[0] == 'y':
                return True
            if reply[0] == 'n':
                return False

# =================================================================================================
#     Main function that actually does stuff
# =================================================================================================

if __name__ == "__main__":
    # Need at least input dir, and maybe output file name as command line args.
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        sys.exit("Need directory path as input, and optionally output file name")
    indir = sys.argv[1]
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

    # Output file check.
    if os.path.exists(outfile):
        if not yes_or_no("Output file `" + outfile + "` exists. Do you want to overwrite?"):
            sys.exit()

    print("We will now process the directoy and look for paired end mates.")
    print("Please use the provided output to check that the mates are correctly assigned to each other,")
    print("and that the extracted sample names make sense!")
    print("(If not, this code here has to be adapted accordingly.)")

    # Get the read mate pairs...
    seqfiles = get_fastq_files( indir )
    mates = get_mates( seqfiles )

    # ... and write them to a table.
    # Potential future expansion: If the same sample name occurs, make this different units
    # (1..n) of that sample name, instead of complaining here with the error message.
    samplenames = []
    with open( outfile, 'w' ) as out:
        out.write("sample\tunit\tplatform\tfq1\tfq2\n")
        for tup in mates:
            name = get_sample_name( tup )
            if name in samplenames:
                print(colored( "Sample name " + name + " occurs multiple times!", "red" ))
            samplenames.append(name)
            out.write(name + "\t1\t-\t" + tup[0] + "\t" + tup[1] + "\n")
