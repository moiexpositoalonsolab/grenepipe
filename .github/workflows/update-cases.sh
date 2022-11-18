#!/bin/bash

# Helper script to automatically update the test cases in the CI workflow.

# Change to main grenepipe directory, so that we can run this from everywhere.
cd `dirname ${0}`/../..

# Get all test cases.
BASE=".github/workflows"
ls test/cases/*.txt | \
    sed "s/.txt//" | \
    sed "s?test/cases/??" | \
    sed "/^_/d" | \
    sed "s/^/          - /" \
    > ${BASE}/cases.txt

# Replace them in the CI yaml.
# Adapted from here: https://unix.stackexchange.com/a/212011
ed ${BASE}/ci.yaml <<EOF
/^ *# CASES_BEGIN
+,/^ *# CASES_END/-1d
-r ${BASE}/cases.txt
w
q
EOF

rm ${BASE}/cases.txt
