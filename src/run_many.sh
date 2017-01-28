#!/bin/bash

for i in {1..30} 
do
    echo "Running Test$i ..."
    time ./run example1.cfg "test$i"
    echo "Completed." 
    echo ""
done
