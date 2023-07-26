#Small bash script to run 2 sed commands
#!/bin/bash

n=$1

find ./run/ -type f -exec sed -i -e 's/numberOfSubdomains .*/numberOfSubdomains '$n';/g' {} \;
find ./run/ -type f -exec sed -i -e 's/-np . /-np '$n' /g' {} \;
find ./run/ -type f -exec sed -i -e 's/-np .. /-np '$n' /g' {} \;
