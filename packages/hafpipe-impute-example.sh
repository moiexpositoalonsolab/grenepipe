#!/usr/bin/env bash

# Example script to showcase how to make a custom SNP impute script for Task 2 of our HAFpipe
# integration in grenepipe. The script takes one argument - the original SNP table computed
# by Task 1 on HAFpipe, and is expected to procude the imputed table, with the same file name,
# appended by `.my_method` for some method name of how you want to name your imputation.
# This `my_method` name then needs to match the `impmethod` setting in the grenepipe config file,
# so that when running grenepipe, we know which imputed table files to use.
#
# If you are actually running a more elaborate method for this, of course your script can
# process the table as needed. You can for example read parameters from a file, and use them
# for the method, or make a call to some other tools here, as needed.
#
# As this is just an example, we do not actually impute anything, but simply make a symlink
# to the origina table here, thus pretending that we did something useful here.

cd `dirname ${1}`
f=`basename ${1}`
ln -s ${f} ${f}.my_method
