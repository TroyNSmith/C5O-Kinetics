#!/bin/bash

cd $1

if [ ! -f "submit.sh" ]; then
    echo "No submit.sh found, please check $DIR..."
    exit 1
fi

for word in $(sq --me); do
    if [[ $word == *"$1"* ]]; then
        echo "Job is already submitted and/or running, please wait for it to finish..."
        exit 1
    fi
done

if [ -f "CCSDT/CCSDT.log" ]; then
    echo "CCSDT ended abruptly or already completed, please delete the log file if it needs to be resubmitted..."
    exit 1
else
    echo "Starting $1 optimization..."
    bash submit.sh
fi