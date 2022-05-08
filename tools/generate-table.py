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
    # Split into dir and file name. If the dirs already differ, we are done.
    # We are doing this check so that we can rely on having the same dir below.
    # We also add the trailing slash to the dir, so that we have the full path length correctly.
    s1s = os.path.split(s1)
    s2s = os.path.split(s2)
    if s1s[0] != s2s[0] or len(s1s[1]) != len(s2s[1]):
        return -1;
    dirn = str(s1s[0])
    if not dirn.endswith('/'):
        dirn += '/'

    # Find the single position at which the two file basenames differ,
    # or -1 if they differ in more than one position
    # (or are the same string, but that should not happen in our case).
    pos = -1
    for i, (c1, c2) in enumerate(zip(s1s[1], s2s[1])):
        if c1 != c2:
            if pos != -1:
                return -1
            else:
                pos = len(dirn) + i

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

def get_mates(seqfiles):
    # Walk through that list and find pairs of reads by checking for single character differences,
    # where one of the file names is `R1`, and the other is `R2`.
    # print(seqfiles)
    mates = []
    while len(seqfiles):
        seq1 = seqfiles.pop(0)
        # print(seq)
        found_mate = False
        for j in range(len(seqfiles)):
            seq2 = seqfiles[j]
            match_pos = match_pe_files(seq1, seq2)
            if match_pos > -1:
                # print(seq1, seqfiles[j])
                mates.append(( seq1, seq2, match_pos ))
                seqfiles.pop(j)
                found_mate = True
                break
        if not found_mate:
            mates.append(( seq1, "", -1 ))
    return mates

def get_sample_name(tup):
    # Expects one element of the get_mates function output list,
    # and cleans up the beginning and end of the file names to get a clean sample name.
    if tup[2] >=0:
        return re.split( '/', tup[0][:tup[2]] )[-1].strip( '-_R' )
    else:
        # Remove directory, and try to find the extension to remove
        # (which should always succeed, as we only picked those files in the first place).
        n = re.split( '/', tup[0] )[-1]
        for s in extensions:
            if n.endswith(s):
                return n[:-len(s)]
        return n

def write_table(mates, outfile):
    # If the same sample name occurs, make this different units (1..n) of that sample name.
    samplenames = {}
    mate_cnt = 0
    with open( outfile, 'w' ) as out:
        # Write the table header. We currently also write the `platform` column, but as do not know
        # what the platform is, we just provide a dummy '-' here. Users can then edit this as needed.
        # We could make this an argument of the script, but let's keep it simple for now.
        out.write("sample\tunit\tplatform\tfq1\tfq2\n")
        for tup in mates:

            # Get the name and unit for the sample.
            # If samples with the same name appear, we give them increasing unit numbers.
            name = get_sample_name( tup )
            if name not in samplenames.keys():
                samplenames[name] = 1
            else:
                samplenames[name] += 1
            unit = samplenames[name]

            # If this is a mate pair, print its two files. If it is single, print the file.
            if tup[2] >= 0:
                # Print the matches, with the difference in red, so that the user
                # can check that the pairing is correct.
                match_pos = tup[2]
                print()
                print(colored("Match (paired end):", "green"))
                print(tup[0][0:match_pos] + colored(tup[0][match_pos], "red") + tup[0][match_pos+1:])
                print(tup[1][0:match_pos] + colored(tup[1][match_pos], "red") + tup[1][match_pos+1:])
            else:
                print()
                print(colored("No match (single end):", "yellow"))
                print(tup[0])
            print(colored("Name: " + name + ", unit: " + str(unit), "blue"))

            # Now print the result table.
            out.write(name + "\t" + str(unit) + "\t-\t" + tup[0] + "\t" + tup[1] + "\n")
            if tup[2] >= 0:
                mate_cnt += 1

    # Summary for the user
    print()
    print("Summary:")
    print(colored("Found " + str(mate_cnt) + " mates.", "green"))
    print(colored("Found " + str(len(mates) - mate_cnt) + " singles.", "yellow"))

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

    print("We will now process the directory and look for paired end mates.")
    print("Please use the provided output to check that the mates are correctly assigned to each other,")
    print("and that the extracted sample names make sense!")
    print("(If not, this script here has to be adapted accordingly.)")

    # Get the read mate pairs and write them to a table.
    seqfiles = get_fastq_files( indir )
    mates = get_mates( seqfiles )
    write_table( mates, outfile )
