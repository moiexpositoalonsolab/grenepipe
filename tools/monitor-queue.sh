#!/bin/bash

while true; do
    date "+%F %T"
    squeue -u `whoami`
    echo
    sleep 5
done
