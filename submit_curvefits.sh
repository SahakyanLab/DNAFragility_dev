# !/bin/bash

RUN_SCRIPT="TRUE"

# args for scripts
cores=1
nr_of_lines=$(wc -l < Raw_data/org_file.csv)

for auto_fit in "TRUE" "FALSE"
do
    echo "Auto-fitting Gaussian curves = ${auto_fit}"

    if [ "${auto_fit}" == "TRUE" ]
    then
        check_to_clean=$(awk -F, -v "row=1" 'NR==row { print $15; exit }' Raw_data/org_file.csv)
        if [ ! -z "${check_to_clean}" ]
        then
            cut -d ',' -f 1-14 Raw_data/org_file.csv > Raw_data/temp.csv
            mv Raw_data/temp.csv Raw_data/org_file.csv
        fi
    fi

    if [ "${RUN_SCRIPT}" == "TRUE" ]
    then
        for ((i=2; i<=$nr_of_lines; i++))
        # for i in 52
        do
            dsb_map=$(awk -F, -v "row=$i" 'NR==row { print $11; exit }' Raw_data/org_file.csv)
            if [ "${dsb_map}" == "TRUE" ]
            then
                breakpoint_type=$(awk -F, -v "row=$i" 'NR==row { print $1; exit }' Raw_data/org_file.csv)
                exp=$(awk -F, -v "row=$i" 'NR==row { print $2; exit }' Raw_data/org_file.csv)
                to_process=$(awk -F, -v "row=$i" 'NR==row { print $4; exit }' Raw_data/org_file.csv)
                alignment_strands=$(awk -F, -v "row=$i" 'NR==row { print $12; exit }' Raw_data/org_file.csv)
                ref_path=$(awk -F, -v col1=$breakpoint_path -v col2=$exp '$1==col1 || $2==col2 { print $3; exit }' Raw_data/org_file.csv)
                category=$(awk -F, -v "row=$i" 'NR==row { print $10; exit }' Raw_data/org_file.csv)

                echo "Processing ${breakpoint_type}/${exp}..."

                # calculate sequence driven effects near breakpoints
                cd ./01_RMSD/scripts/
                bash script.sh $breakpoint_type $exp $ref_path $cores 1 "FALSE" $category $auto_fit
                cd ../../
            fi
        done
    else 
        echo "Script did not run. Set RUN_SCRIPT=TRUE to run."
    fi

    echo "Processing heatmap tracts..."
    cd ./02_Correlations/scripts/
    bash script_tracts.sh "z-score" $auto_fit "1"
    cd ../../

    if [ "${auto_fit}" == "TRUE" ]
    then
        cd ./01_RMSD/scripts/
        bash script_rmsdexploration.sh
        cd ../../
    fi
done