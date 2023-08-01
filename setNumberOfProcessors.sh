#Small bash script to run 2 sed commands
#!/bin/bash
n=$1

find ./run/ -type f -exec grep -q 'numberOfSubdomains\| -np [0-9]\{1,2\} ' {} \; -exec sed -i -e 's/numberOfSubdomains .*/numberOfSubdomains '$n';/g' -e 's/-np . /-np '$n' /g' -e 's/-np .. /-np '$n' /g' {} \;
