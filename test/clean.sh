#!/bin/bash

# Change to the parent directory of where this script is located, which is the main grenepipe
# directory, so that we can call this script from anywhere, and run everything from there.
cd `dirname ${0}`/..

# Clean all files produced by the run script.
rm -rf ./test/out-*
rm -rf ./test/reference
rm -f  ./test/samples.tsv
rm -f  ./test/samples-pe.tsv
