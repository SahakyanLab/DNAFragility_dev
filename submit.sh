# !/bin/bash

RUN_SCRIPT="TRUE"

# args for scripts
cores=1
interval=5000000
chromosome=1
upper_limit=1
control="FALSE"
nr_of_lines=$(wc -l < Raw_data/org_file.csv)

if [ "${RUN_SCRIPT}" == "TRUE" ]
then
    # for ((i=2; i<=$nr_of_lines; i++))
    # for ((i=2; i<=26; i++))
    # for ((i=84; i<=108; i++))
    # for ((i=32; i<=35; i++))
    # do
        i=28
        # i=27
        # i=37
        # i=112
        # i=117
        breakpoint_type=$(awk -F, -v "row=$i" 'NR==row { print $1; exit }' Raw_data/org_file.csv)
        exp=$(awk -F, -v "row=$i" 'NR==row { print $2; exit }' Raw_data/org_file.csv)
        to_process=$(awk -F, -v "row=$i" 'NR==row { print $4; exit }' Raw_data/org_file.csv)
        dsb_map=$(awk -F, -v "row=$i" 'NR==row { print $11; exit }' Raw_data/org_file.csv)
        alignment_strands=$(awk -F, -v "row=$i" 'NR==row { print $12; exit }' Raw_data/org_file.csv)
        ref_path=$(awk -F, -v col1=$breakpoint_path -v col2=$exp '$1==col1 || $2==col2 { print $3; exit }' Raw_data/org_file.csv)
        category=$(awk -F, -v "row=$i" 'NR==row { print $10; exit }' Raw_data/org_file.csv)

        echo "Processing ${breakpoint_type}/${exp}..."
        mkdir -p Raw_data/${breakpoint_type}/${exp}/{breakpoint_positions,kmertone}

        if [ "${to_process}" == "TRUE" ] && [ "${dsb_map}" == "TRUE" ]
        then
            start_idx=$(awk -F, -v "row=$i" 'NR==row { print $6; exit }' Raw_data/org_file.csv)
            end_idx=$(awk -F, -v "row=$i" 'NR==row { print $7; exit }' Raw_data/org_file.csv)
            file_format=$(awk -F, -v "row=$i" 'NR==row { print $9; exit }' Raw_data/org_file.csv)

            path_to_bp_files="Raw_data/${breakpoint_type}/${exp}/breakpoint_positions"
            nr_of_files=$((ls $path_to_bp_files/*.csv | wc -l) 2> /dev/null)

            if [[ ! -f $path_to_bp_files/*.csv ]]
            then
                # if [ $nr_of_files -lt 22 ]
                if [ $nr_of_files -lt 24 ]
                then
                    # Process files 
                    cd ./00_Preprocessing/scripts/
                    bash script.sh $breakpoint_type $exp $file_format $start_idx $end_idx $alignment_strands
                    cd ../../
                fi
            fi
        else
            # align raw sequencing reads
            cd ./00_ReadsAlignment/scripts/
            # bash script.sh $breakpoint_type $exp $ref_path $cores $interval $alignment_strands
            cd ../../
        fi

        # calculate sequence driven effects near breakpoints
        cd ./01_RMSD/scripts/
        # bash script.sh $breakpoint_type $exp $ref_path $cores $chromosome $control $category
        cd ../../

        # perform kmertone enrichment/depletion analysis
        first_idx=$(awk -F, -v "row=$i" 'NR==row { print $5; exit }' Raw_data/org_file.csv)

        cd ./02_Correlations/scripts/
        bash script.sh $breakpoint_type $exp $ref_path $cores $upper_limit $first_idx
        cd ../../
    # done
else 
    echo "Script did not run. Set RUN_SCRIPT=TRUE to run."
fi