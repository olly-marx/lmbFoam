# bash script to store the run data in a directory
# Usage: ./storeRun.sh <directory name>
#!/bin/sh

dir=$1

rm -rf data/$dir
mkdir -p data/$dir
mv ./[0-9]* data/$dir

# move the log files to the data directory
mv ./log.* ./data/$dir/

# copy other case files to the data directory
cp -r ./system ./data/$dir/
cp -r ./constant ./data/$dir/
cp -r ./save ./data/$dir/

