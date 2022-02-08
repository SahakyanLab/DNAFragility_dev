# !/bin/bash

RUN_SCRIPT="TRUE"

# args for scripts
breakpoint_type="03-Ancient_DNA"
exp="Vindija_Neandertal"
cores=1
interval=5000000
chromosome=1
upper_limit=1
alignment_strands="both"
exp_num=1
control="FALSE"

# match the breakpoint type with the ref seq
declare -A refpath_map=(
    [00-Ultrasonication/Simons_exp]=Simons_exp
    [04-Enzymatic/exp]=Simons_exp
    [03-Ancient_DNA/Chagyrskaya_Neandertal]=Simons_exp
    [03-Ancient_DNA/Vindija_Neandertal]=Simons_exp

    [01-Nebulization/1000_Genomes_exp]=1000_Genomes_exp  # need to re-run btw...
    [03-Ancient_DNA/Altai_Neandertal]=1000_Genomes_exp
    [03-Ancient_DNA/Denisovan_Genome]=1000_Genomes_exp

    [01-Nebulization/Pilot/1000_Genomes_exp]=1000_Genomes_Pilot # need to re-run btw...
)

for key in ${!refpath_map[@]}
do 
    if [ "${breakpoint_type}/${exp}" == $key ]
    then
        ref_path=${refpath_map[${breakpoint_type}/${exp}]}
    fi
done

if [ "${RUN_SCRIPT}" == "TRUE" ]
then
    echo "Running script..."

    # align raw sequencing reads
    cd ./00_ReadsAlignment/scripts/
    bash script.sh $breakpoint_type $exp $ref_path $cores $interval $alignment_strands
    cd ../../

    # calculate sequence driven effects near breakpoints
    cd ./01_RMSD/scripts/
    bash script.sh $breakpoint_type $exp $ref_path $cores $exp_num $chromosome $control
    cd ../../

    # perform kmertone enrichment/depletion analysis
    cd ./02_Correlations/scripts/
    bash script.sh $breakpoint_type $exp $ref_path $cores $exp_num $upper_limit
    cd ../../
else 
    echo "Script did not run. Set RUN_SCRIPT=TRUE to run."
fi