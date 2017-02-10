#!/bin/bash

for i in {1..30} 
do
    echo "Running TestB$i ..."
    time ./run example1.cfg "testB$i"
    echo "Completed." 
    echo ""
done
