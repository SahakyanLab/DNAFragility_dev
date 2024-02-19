#!/usr/bin/bash

# directory placeholders
mkdir -p ./data/00_Breakage/Ultrasonication/{Simons_exp_1,Simons_exp_12,Simons_exp_15,Simons_exp_18,Simons_exp_20,Simons_exp_23,Simons_exp_3,Simons_exp_6,Simons_exp_9,Simons_exp_10,Simons_exp_13,Simons_exp_16,Simons_exp_19,Simons_exp_21,Simons_exp_24,Simons_exp_4,Simons_exp_7,Simons_exp_11,Simons_exp_14,Simons_exp_17,Simons_exp_2,Simons_exp_22,Simons_exp_25,Simons_exp_5,Simons_exp_8}
mkdir -p ./data/00_Breakage/Nebulization/{1000_Genomes_exp_2,1000_Genomes_exp_3}
mkdir -p ./data/00_Breakage/Ancient_DNA/{Altai_Neandertal,Chagyrskaya_Neandertal,Denisovan_Genome,Ust_Ishim,Vindija_Neandertal}
mkdir -p ./data/00_Breakage/cfDNA/{Bladder_cancer,Esophageal_cancer,Kidney_cancer_Clear_cell,Pancreatic_cancer_Ductal_adenocarcinoma,Testicular_cancer_Seminomatous,Breast_cancer_Invasive_infiltratingductal,Head_and_neck_cancer,Liver_cancer_Hepatocellular_carcinoma,Prostate_cancer,Ulcerative_Colitis,Colorectal_cancer,Healthy,Lung_cancer_Adenocarcinoma,Skin_cancer_Melanoma,Uterine_cancer,Crohns_disease,Healthy_person,Ovarian_cancer,Systemic_Lupus_Erythematosus}
mkdir -p ./data/00_Breakage/Enzymatic/{AsiSI_AID-DlvA_cells,EcoRV_HeLa_cells,Nt_BbvCI_K562_cells/NT,Nt_BbvCI_K562_cells/Sap,Twist_library}
mkdir -p ./data/00_Breakage/Recombination_rates/
mkdir -p ./data/00_Breakage/Apoptoseq/Apoptotic_DSBs/HL-60-leukemic/{AHH001,AHH002}
mkdir -p ./data/00_Breakage/ENDseq/WRN_loss_DSBs/{HCT116_APH-ATR-inhib,HCT116_siWRN,HCT116_shWRN,KM12_APH-ATR-inhib,KM12_siWRN,KM12_shWRN}
mkdir -p ./data/00_Breakage/SSB/{G5,Imatinib,SN,Talazoparib}/{6h,48h}
mkdir -p ./data/00_Breakage/BrITL/MDA-MB-231/{APH-ATR-inhib,DMSO}
mkdir -p ./data/00_Breakage/DSBCapture/NHEK_DSBs
mkdir -p ./data/00_Breakage/SARseq/iNeuron
mkdir -p ./data/00_Breakage/CCseq/RPE-1_Top2-linked_DNA_DSBs/{WG13-WG21_RPE1_WT_G1_UT,WG14-WG22_RPE1_WT_G1_VP16,WG17-WG25A_RPE1_T2B_G1_UT,WG18-WG26_RPE1_T2B_G1_VP16}
mkdir -p ./data/00_Breakage/sBLISS/{TK6_Top2_mediated_DSBs,CD34_Top2_mediated_DSBs}/{DMSO_rep1,DMSO_rep2,ETO_rep1,ETO_rep2}
mkdir -p ./data/00_Breakage/sBLISS/K562_Top2_mediated_DSBs/{ETO,DMSO}
mkdir -p ./data/00_Breakage/sBLISS/Colibactin_Ecoli_induced_DSBs/{ETO,DMSO,Caco-2_etoposide_rep1,Caco-2_etoposide_rep2,Caco-2_etoposide_rep3,Caco-2_etoposide_rep4,Caco-2_pks_minus_Ecoli_rep1,Caco-2_pks_minus_Ecoli_rep2,Caco-2_pks_minus_Ecoli_rep3,Caco-2_pks_minus_Ecoli_rep4,Caco-2_pks_plus_Ecoli_rep1,Caco-2_pks_plus_Ecoli_rep2}
mkdir -p ./data/04_QM_parameters/6mer/Raw

# install relevant software
cd ./setup/
bash install_packages.sh

# download natural decay files
for folder in ./data/00_Breakage/Ancient_DNA/*
do
    cd "$folder"
    echo "Downloading files in $folder"
    bash $sh_file
    cd ../../../../
done

# process transcription factors
cd ./data/02_JASPAR/
Rscript jaspar_download.R
cd ../../