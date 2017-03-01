#!/bin/bash
CONFIG=$1
TAG=$2
NUM=$3
for ((i=1; i<=$NUM; i++)) 
do
    echo "Running $TAG $i ..."
    time ./run "$CONFIG" "$TAG-$i"
    echo "Completed." 
    echo ""
done
