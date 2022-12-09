#!/bin/bash

# Go to the example dir
cd "$(dirname "$0")"

# Get the directory, so that we can use it for replacing paths
DIR=`pwd -P $(dirname "$0")`

# Copy the config file
cp ../config.yaml ./config.yaml

# Replace paths in config file
sed -i.bak -e "s?/path/to/data?${DIR}?g" config.yaml
sed -i.bak -e "s?genome.fa?TAIR10_chr_all.fa.gz?g" config.yaml
rm config.yaml.bak

# Replace paths in samples table
sed -i.bak -e "s?/path/to/grenepipe/example?${DIR}?g" samples.tsv
rm samples.tsv.bak
