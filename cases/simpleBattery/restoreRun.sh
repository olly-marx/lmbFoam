#!/bin/sh

# bash to restore data from a run
# takes all files from a data directory and puts them back in the case
# use args to give the script the directory name

echo "Restoring data from $1"
mv $1/* .
