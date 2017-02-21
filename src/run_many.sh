#!/bin/bash
TAG=$1
NUM=$2
for ((i=1; i<=$NUM; i++)) 
do
    echo "Running $TAG$i ..."
    time ./run example1.cfg "$TAG$i"
    echo "Completed." 
    echo ""
done
