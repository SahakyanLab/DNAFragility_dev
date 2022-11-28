# !/bin/bash

check_to_clean=$(/usr/bin/awk -F, -v "row=1" 'NR==row { print $16; exit }' ../../data/org_file.csv)
if [ ! -z "${check_to_clean}" ]
then
    /usr/bin/cut -d ',' -f 1-15 ../../data/org_file.csv > ../../data/temp.csv
    /bin/mv ../../data/temp.csv ../../data/org_file.csv
fi