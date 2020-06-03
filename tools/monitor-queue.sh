#!/bin/bash

while true; do
    squeue -u `whoami`
    sleep 5
done
