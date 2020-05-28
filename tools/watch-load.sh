#!/bin/bash

#cnt=0
while true; do 

    # Print processor load and 10 most cpu consuming processes
    uptime >> watch_uptime.log
    (echo "%CPU %MEM ARGS $(date)" && ps -e -o pcpu,pmem,args --sort=pcpu | cut -d" " -f1-5 | tail | tac -) >> watch_ps.log
    
    # Every 5sec, print 10 most cpu consuming processes
    #if ! ((cnt % 5)); then
        #(echo "%CPU %MEM ARGS $(date)" && ps -e -o pcpu,pmem,args --sort=pcpu | cut -d" " -f1-5 | tail | tac -) >> ps.log
    #fi  
    
    sleep 5
    #sleep 1
    #cnt=$((cnt+1))
done

#while true; do uptime >> uptime.log; sleep 1; done
#while true; do (echo "%CPU %MEM ARGS $(date)" && ps -e -o pcpu,pmem,args --sort=pcpu | cut -d" " -f1-5 | tail | tac -) >> ps.log; sleep 5; done


