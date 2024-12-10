#!/usr/bin/env python3

import sys, os
import re
import gzip
import argparse
from glob import glob
from collections import namedtuple
from termcolor import colored

# =================================================================================================
#     Usage and Description
# =================================================================================================

# Generate the samples table to be used as grenepipe input that lists paired end read samples.
# Usage: ./generate-table.py [--single] <directory-with-fastq-files> [<samples.tsv>]

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

# We want nicer data access than just addressing matches by tuple indices.
# We store the file names of the two matches, as well as the position where the match of the '1'/'2'
# char match occurrs. While finding match candidates, we also need the indices in the input list.
Match = namedtuple("Match", ['seq1', 'seq2', 'pos', 'idx1', 'idx2'])

def get_fastq_files(indir, shallow):
    # Get all files in the given directory that match our file extensions, and sort them,
    # so that adjacent entries are potential paired files.
    seqfiles = []
    if shallow:
        for f in os.listdir(indir):
            path = os.path.abspath(os.path.join(indir, f))
            if path.endswith( tuple( extensions )):
                seqfiles.append(path)
    else:
        for dirpath, subdirs, filenames in os.walk(indir):
            # print(filenames)
            for f in filenames:
                path = os.path.abspath(os.path.join(dirpath, f))
                if path.endswith( tuple( extensions )):
                    seqfiles.append(path)
    seqfiles.sort()
    return seqfiles

def is_interleaved_fastq(fq):
    with open(fq, 'rb') as test_f:
        is_gzipped = ( test_f.read(2) == b'\x1f\x8b' )

    # We check line 1 and 5 (names of the first and second sequence),
    # and see if they contain " 1:" and " 2:" respectively, which is a strong indicator
    # for interleaved fastq files (at least the ones we have seen so far).
    with ( gzip.open(fq,'r') if is_gzipped else open(fq) ) as f:
        cnt = 1
        for line in f:
            if cnt == 1:
                line1 = str(line)
            if cnt == 5:
                line5 = str(line)
            if cnt > 8:
                break
            cnt += 1
    if cnt <= 8:
        return False
    return " 1:" in line1 and " 2:" in line5

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
    # or -1 if they differ in more than one position (except 'R'/'F')
    # (or are the same string, but that should not happen in our case).
    pos = -1
    for i, (c1, c2) in enumerate(zip(s1s[1], s2s[1])):
        if c1 != c2:
            if c1.upper() in ['R', 'F'] and c2.upper() in ['R']:
                continue
            elif pos != -1:
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
        candidates = []
        for j in range(len(seqfiles)):
            seq2 = seqfiles[j]
            match_pos = match_pe_files(seq1, seq2)
            if match_pos > -1:
                # print(seq1, seqfiles[j])
                # We found a candidate. Add it, as well as the index of the second match file.
                match = Match(seq1=seq1, seq2=seq2, pos=match_pos, idx1=0, idx2=j)
                candidates.append(match)
        if len(candidates) == 0:
            # Did not find any matches. Single end read.
            match = Match(seq1=seq1, seq2="", pos=-1, idx1=0, idx2=None)
            mates.append(match)
        elif len(candidates) == 1:
            # Found one match, which we now use. Add to results,
            # and remove the second match file from the inputs.
            mates.append(candidates[0])
            seqfiles.pop(candidates[0].idx2)
        else:
            # Found multiple matches. We look for one that has "R" or "r" in front of the match
            # position, which often denotes the forward/backward. If that does not work, we fail,
            # as in that case, we cannot determine the right match from the file naming scheme.
            cand_idx = -1
            multiple_r = False
            for i in range(len(candidates)):
                cand = candidates[i]
                if cand.pos == 0:
                    # If the match position is 0 (that is, the first char is '1'/'2' in the files),
                    # these cannot have an 'R' in front of them.
                    continue
                if cand.seq1[ cand.pos - 1 ].upper() in ['R', 'F']:
                    assert cand.seq2[ cand.pos - 1 ].upper() in ['R']
                    # We found an 'R'. If we already found another, that is also a fail.
                    if cand_idx > -1:
                        multiple_r = True
                    cand_idx = i
            if cand_idx == -1 or multiple_r:
                # Here, we either found no 'R', or multiple ones. Either way, we fail.
                print(colored("Found multiple conflicting match pair candidates.", "red"))
                print("File " + seq1 + " matches with:")
                for cand in candidates:
                    print("  -  " + cand.seq2)
                print(
                    "We cannot automatically determine match pairs with this. Please consider "
                    "renaming your files so that the intended match pairs have the letter 'R'/'F' "
                    "in front of the two numbers, like 'R1'/'F1' and 'R2', case independent."
                )
                sys.exit()
            # Here, we have determined a pair that works. Add it to the results, and remove
            # the second sequence from the inputs.
            mates.append(candidates[cand_idx])
            seqfiles.pop(candidates[cand_idx].idx2)

    return mates

def get_singles(seqfiles):
    # Walk through the list of fastq files, but treat all of them as single end.
    singles = []
    for seq in seqfiles:
        match = Match(seq1=seq, seq2="", pos=-1, idx1=0, idx2=None)
        singles.append(match)

    return singles

def get_sample_name(match):
    # Expects one element of the get_mates function output list,
    # and cleans up the beginning and end of the file names to get a clean sample name.
    if match.pos >= 0:
        # If the char immediately before the 1/2 is an R, we clean this.
        # Then, we also clean trailing dashes, dots, and underscores.
        # We do not do this in one step, as that might truncate sample names ending in 'r'.
        name = re.split( '/', match.seq1[:match.pos] )[-1]
        if name[-1:].upper() in ['R', 'F']:
            name = name[:-1]
        return name.rstrip( '-_.' )
    else:
        # Remove directory, and try to find the extension to remove
        # (which should always succeed, as we only picked those files in the first place).
        n = re.split( '/', match.seq1 )[-1]
        for s in extensions:
            if n.endswith(s):
                return n[:-len(s)]
        return n

# Helper to check if a file path contains weird characters.
def valid_filepath( fn ):
    # Only accept alnum, underscore, dash, dot, and slashes.
    clean = fn.replace('_', '').replace('-', '').replace('.', '').replace('/', '').replace('\\', '')
    return clean.isalnum() and clean.isascii()

def write_table(mates, outfile):
    # If the same sample name occurs, make this different units (1..n) of that sample name.
    samplenames = {}
    mate_cnt = 0
    invalid_samples = 0
    invalid_files = 0
    interleaved_files = 0
    with open( outfile, 'w' ) as out:
        # Write the table header. We currently also write the `platform` column, but as do not know
        # what the platform is, we just provide a dummy '-' here. Users can then edit this as needed.
        # We could make this an argument of the script, but let's keep it simple for now.
        out.write("sample\tunit\tplatform\tfq1\tfq2\n")
        for match in mates:

            # Get the name and unit for the sample.
            # If samples with the same name appear, we give them increasing unit numbers.
            name = get_sample_name( match )
            if name not in samplenames.keys():
                samplenames[name] = 1
            else:
                samplenames[name] += 1
            unit = samplenames[name]

            # If this is a mate pair, print its two files. If it is single, print the file.
            if match.pos >= 0:
                # Print the matches, with the difference in red, so that the user
                # can check that the pairing is correct.
                match_pos = match.pos
                print()
                print(colored("Match (paired end):", "green"))
                print(
                    match.seq1[0:match_pos] +
                    colored(match.seq1[match_pos], "red") +
                    match.seq1[match_pos+1:]
                )
                if is_interleaved_fastq(match.seq1):
                    print(colored( "The file is likely an interleaved fastq file.", "red"))
                    interleaved_files += 1
                print(
                    match.seq2[0:match_pos] +
                    colored(match.seq2[match_pos], "red") +
                    match.seq2[match_pos+1:]
                )
                if is_interleaved_fastq(match.seq2):
                    print(colored("The file is likely an interleaved fastq file.", "red"))
                    interleaved_files += 1
            else:
                print()
                print(colored("No match (single end):", "yellow"))
                print(match.seq1)
                if is_interleaved_fastq(match.seq1):
                    print(colored("The file is likely an interleaved fastq file.", "red"))
                    interleaved_files += 1
            print(colored("Name: " + name + ", unit: " + str(unit), "blue"))

            # Now print the result table.
            out.write(name + "\t" + str(unit) + "\t-\t" + match.seq1 + "\t" + match.seq2 + "\n")
            if match.pos >= 0:
                mate_cnt += 1

            # Count the still invalid names. Can only happen without --clean
            if not valid_filepath( name ):
                invalid_samples += 1
            if not valid_filepath( match.seq1 ):
                invalid_files += 1
            if len(match.seq2) > 0 and not valid_filepath( match.seq2 ):
                invalid_files += 1

    # Summary for the user
    print()
    print("Summary:")
    print(colored("Found " + str(mate_cnt) + " mates.", "green"))
    print(colored("Found " + str(len(mates) - mate_cnt) + " singles.", "yellow"))

    # Warn about invalid file names
    if invalid_samples > 0:
        print(colored(
            "Out of these, " + str(invalid_samples) + " sample names contain invalid characters. "
            "This might cause trouble when using grenepipe with these files. "
            "Please consider to use the `copy-samples.py` script with the `--clean` option "
            "to fix these sample names.", "red"
        ))
    if invalid_files > 0:
        print(colored(
            "Out of these, " + str(invalid_files) + " fastq files contain invalid characters. "
            "This might cause trouble when using grenepipe with these files. "
            "Please consider to use the `copy-samples.py` script with the `--clean` option "
            "to fix these file names.", "red"
        ))
    if interleaved_files > 0:
        print(colored(
            "Out of these, " + str(interleaved_files) + " fastq files are likely interleaved "
            "and should not be used with grenepipe, unless you know what you are doing and set "
            "the appropriate params for all tools. Use for example https://github.com/lh3/seqtk "
            "to de-interleave them, and then run this script again.", "red"
        ))

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

    # Set up the command line interface and parse the user provided options
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", help="Input directory to recursively scan for fastq files")
    parser.add_argument(
        "table", nargs='?', default=outfile,
        help="Output table name, optional, defaults to " + outfile
    )
    parser.add_argument(
        "--single",
        action='store_true',
        help="If set, all fastq files are considered to be single-end, " +
        "and no attempt to match their names is performed"
        )
    parser.add_argument(
        "--shallow",
        action='store_true',
        help="If set, only the given directory is searched for files, and not its sub-directories"
        )
    args = parser.parse_args()

    # Output file check.
    if os.path.exists(args.table):
        if not yes_or_no("Output file `" + args.table + "` exists. Do you want to overwrite?"):
            sys.exit()

    # Some output for the user
    if args.single:
        print("We will now process the directory and look for fastq files,")
        print("treating them as single-ended read files.")
    else:
        print("We will now process the directory and look for paired end mates.")
        print("Please use the provided output to check that the mates are correctly")
        print("assigned to each other, and that the extracted sample names make sense!")
        print("(If not, this script here has to be adapted accordingly.)")

    # Get the read mate pairs and write them to a table.
    seqfiles = get_fastq_files( args.dir, args.shallow )
    if args.single:
        mates = get_singles( seqfiles )
    else:
        mates = get_mates( seqfiles )
    write_table( mates, args.table )
