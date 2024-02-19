# DNAFragility_dev

Research and Development Code for Ultrasound Sonication Induced DNA Breakpoints in Sequencing

## Setup

Clone the project:

```
git clone https://github.com/SahakyanLab/DNAFragility_dev.git
```

Please follow the instructions below on how to acquire the public datasets, setup the directory stucture, and software necessary to run all the studies from the publication.  At the end of this `README` file, you can find two separate bash script commands that runs the majority of the setup and runs the calculations sequentially. 

## 1. Software requirements
The resource-demanding computations were performed on a single NVIDIA RTX A6000 GPU with 40GB RAM. The developed workflows and analyses employed the [R programming language 4.3.2](https://www.r-project.org/).

Please run the below script to install the latest versions of the R packages necessary to perform the calculations and analyses. 

```bash
bash ./setup/install_packages.sh
```

Please also download and install the below software.

### Edlib
* Please clone the repo from [this link](https://github.com/Martinsos/edlib) (Edlib >= 1.2.7). Place the [edlib.h](https://github.com/Martinsos/edlib/tree/master/edlib/include) and [edlib.cpp](https://github.com/Martinsos/edlib/tree/master/edlib/src) into [02_Alignment/lib/edlib](https://github.com/SahakyanLab/SonicBreaks_dev/tree/master/02_Alignment/lib/edlib) folder.

### kseq.h

Please download the kseq.h file from [this link](https://github.com/attractivechaos/klib). Place this into [03_FitCurves/lib](https://github.com/SahakyanLab/SonicBreaks_dev/tree/master/03_FitCurves/lib) folder.

### phmap.hpp via gtl
* Please clone the repo from [this link](https://github.com/greg7mdp/gtl). Place the contents of gtl into [03_FitCurves/lib](https://github.com/SahakyanLab/SonicBreaks_dev/tree/master/03_FitCurves/lib) folder.

### `factoextra`

Please note, if you are using Ubuntu, you may have trouble installing the [factoextra](https://github.com/kassambara/factoextra) R package. However, the below steps has worked for us. 

1. sudo apt-get update
2. sudo apt-get install r-cran-car
3. install.packages("factoextra")

This will be needed for the [04_KmericAnalysis](https://github.com/SahakyanLab/SonicBreaks_dev/tree/master/04_KmericAnalysis) folder.

### `Cytoscape` desktop app

Please download the `Cytoscape` desktop application for your operating system from [their website](https://cytoscape.org/download.html). This will be needed for the [04_KmericAnalysis](https://github.com/SahakyanLab/SonicBreaks_dev/tree/master/04_KmericAnalysis) folder.

You also need to install `clustermaker2`, `enrichmentmap`, and `autoannotate` apps within `Cytoscape`.

## 2. Public files to download
### Ultrasonication DNA breakages

25 ultrasonication-induced breakage datasets were retrieved from  NCBI’s [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra) Run Selector for the Bioproject [PRJEB9586](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB9586) with the SRAs as outlined below.

For each dataset, we used the first run. On the SRA website, we clicked `Alignment`, removed any range from `Set Range`, and selected the whole `FASTA` run. To download the chromosome-separated files, we selected the numeric chromosome number in `Choose Reference`. We repeated this for chromosomes `1` to `22`. Then, we click `File` to download the whole `FASTA` run.

* [ERX1104497 (Simons_exp_1)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025627&display=alignment)
* [ERX1104508 (Simons_exp_2)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025638&display=alignment)
* [ERX1104516 (Simons_exp_3)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025646&display=alignment)
* [ERX1104519 (Simons_exp_4)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025649&display=alignment)
* [ERX1104521 (Simons_exp_5)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025651&display=alignment)
* [ERX1104529 (Simons_exp_6)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025659&display=alignment)
* [ERX1097982 (Simons_exp_7)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019042&display=alignment)
* [ERX1097983 (Simons_exp_8)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019043&display=alignment)
* [ERX1097984 (Simons_exp_9)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019044&display=alignment)
* [ERX1097986 (Simons_exp_10)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019046&display=alignment)
* [ERX1097988 (Simons_exp_11)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019048&display=alignment)
* [ERX1097992 (Simons_exp_12)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019052&display=alignment)
* [ERX1097994 (Simons_exp_13)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019054&display=alignment)
* [ERX1097998 (Simons_exp_14)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019058&display=alignment)
* [ERX1098018 (Simons_exp_15)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019078&display=alignment)
* [ERX1098019 (Simons_exp_16)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1019079&display=alignment)
* [ERX1104476 (Simons_exp_17)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025606&display=alignment)
* [ERX1104481 (Simons_exp_18)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025611&display=alignment)
* [ERX1104486 (Simons_exp_19)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025616&display=alignment)
* [ERX1104491 (Simons_exp_20)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025621&display=alignment)
* [ERX1104498 (Simons_exp_21)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025628&display=alignment)
* [ERX1104512 (Simons_exp_22)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025642&display=alignment)
* [ERX1104523 (Simons_exp_23)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025653&display=alignment)
* [ERX1104527 (Simons_exp_24)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025657&display=alignment)
* [ERX1104532 (Simons_exp_25)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR1025662&display=alignment)

Place each of the 22 chromosome-separated FASTA files into the respective folders as indicated inside the brackets. All chromosome fasta file has to have the name structure `chr$num.fasta.gz` where `$num` is any number between `1` and `22`. The downloaded files will automatically be compressed as a `.gz` file. Do not uncompress this as we do it automatically when processing the files.

### Nebulization DNA breakages

Two nebulization-induced DNA breakage datasets were obtained from the below data accession codes from NCBI's NCBI’s [Sequence Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra) Run Selector for the Bioproject [PRJNA33851](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA33851).

* [NA18794 (1000_Genomes_exp_2)](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=ERR753353&display=alignment)
* [NA18795 (1000_Genomes_exp_3)](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=ERR257979&display=alignment)

Place each of the 22 chromosome-separated FASTA files into the respective folders as indicated inside the brackets. Please follow the same downloading process as from [the Ultrasonication DNA breakages section](#ultrasonication-dna-breakages)

### Naturay decay DNA breakages

Five natural decay and fossilisation datasets were retrieved from the Max Planck Institute for Evolutionary Anthropology. These files will be automatically downloaded in the appropriate folders when you run the bash script at the end of this `README` file.

```bash
bash run_all_setup_files.sh
```

### Cell-free DNA breakages

16 datasets of cell-free DNA fragments coming from individual human peripheral blood plasma were obtained from NCBI’s SRA Run Selector for the Bioproject [PRJNA291063](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA291063) with the SRAs as outlined below.

* [SRX1120757 (Healthy)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2129993&display=alignment)
* [SRX1120768 (Breast_cancer_Invasive_infiltratingductal)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130004&display=alignment)
* [SRX1120766 (Ovarian_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130002&display=alignment)
* [SRX1120767 (Skin_cancer_Melanoma)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130003&display=alignment)
* [SRX1120769 (Lung_cancer_Adenocarcinoma)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130005&display=alignment)
* [SRX1120771 (Uterine_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130007&display=alignment)
* [SRX1120774 (Colorectal_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130010&display=alignment)
* [SRX1120776 (Prostate_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130012&display=alignment)
* [SRX1120777 (Head_and_neck_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130013&display=alignment)
* [SRX1120780 (Liver_cancer_Hepatocellular_carcinoma)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130016&display=alignment)
* [SRX1120779 (Bladder_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130015&display=alignment)
* [SRX1120781 (Kidney_cancer_Clear_cell)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130017&display=alignment)
* [SRX1120782 (Testicular_cancer_Seminomatous)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130018&display=alignment)
* [SRX1120784 (Pancreatic_cancer_Ductal_adenocarcinoma)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130020&display=alignment)
* [SRX1120793 (Esophageal_cancer)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2130029&display=alignment)
* [SRX1120758 (Crohns_disease)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2129994&display=alignment)
* [SRX1120760 (Ulcerative_Colitis)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2129996&display=alignment)
* [SRX1120762 (Systemic_Lupus_Erythematosus)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2129998&display=alignment)

Place each of the 22 chromosome-separated FASTA files into the respective folders as indicated inside the brackets. Please follow the same downloading process as from [the Ultrasonication DNA breakages section](#ultrasonication-dna-breakages)

### Physiological DNA breakages

46 physiological breakage datasets were retrieved from various tissues and cell lines.

* [GSE115623 (BrITL/MDA-MB-231/APH-ATR-inhib)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115623/suppl/GSE115623%5FATRi%5FAPH%5F9h.narrowPeak.gz)
* [GSE115623 (BrITL/MDA-MB-231/DMSO)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115623/suppl/GSE115623%5FDMSO%5F9h.narrowPeak.gz)
* [GSE78172 (DSBCapture/NHEK_DSBs)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78172/suppl/GSE78172%5FBreak%5Fprimary%5F20ug.N1.subset.NOdups.q005%5Fpeaks.narrowPeak.gz)
* [GSE167257 SARseq/iNeuron](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167257/suppl/GSE167257%5FSARseq%5FiNeuron%5FOverlapRep123.peaks.bed.gz)
* [GSE136943 (CCseq/RPE-1_Top2-linked_DNA_DSBs/WG13-WG21_RPE1_WT_G1_UT)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136943/suppl/GSE136943%5FWG13%2DWG21%5Fhg19%5FRPE1%5FWT%5FG1%5FUT.FullMap.txt.gz)
* [GSE136943 (CCseq/RPE-1_Top2-linked_DNA_DSBs/WG14-WG22_RPE1_WT_G1_VP16)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136943/suppl/GSE136943%5FWG14%2DWG22%5Fhg19%5FRPE1%5FWT%5FG1%5FVP16.FullMap.txt.gz)
* [GSE136943 (CCseq/RPE-1_Top2-linked_DNA_DSBs/WG17-WG25A_RPE1_T2B_G1_UT)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136943/suppl/GSE136943%5FWG17%2DWG25A%5Fhg19%5FRPE1%5FT2B%5FG1%5FUT.FullMap.txt.gz)
* [GSE136943 (CCseq/RPE-1_Top2-linked_DNA_DSBs/WG18-WG26_RPE1_T2B_G1_VP16)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136943/suppl/GSE136943%5FWG18%2DWG26%5Fhg19%5FRPE1%5FT2B%5FG1%5FVP16.FullMap.txt.gz)
* [GSE121742 (sBLISS/CD34_Top2_mediated_DSBs/DMSO_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3687nnn/GSM3687236/suppl/GSM3687236%5FCD34%5FDMSO%5Frep1.bed.gz)
* [GSE121742 (sBLISS/CD34_Top2_mediated_DSBs/DMSO_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3687nnn/GSM3687238/suppl/GSM3687238%5FCD34%5FDMSO%5Frep2.bed.gz)
* [GSE121742 (sBLISS/CD34_Top2_mediated_DSBs/ETO_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3687nnn/GSM3687237/suppl/GSM3687237%5FCD34%5FETO%5Frep1.bed.gz)
* [GSE121742 (sBLISS/CD34_Top2_mediated_DSBs/ETO_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3687nnn/GSM3687239/suppl/GSM3687239%5FCD34%5FETO%5Frep2.bed.gz)
* [GSE121742 (sBLISS/K562_Top2_mediated_DSBs/ETO)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444989/suppl/GSM3444989%5FK562%5FETO.bed.gz)
* [GSE121742 (sBLISS/K562_Top2_mediated_DSBs/DMSO)](ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444988/suppl/GSM3444988%5FK562%5FDMSO.bed.gz)
* [GSE121742 (sBLISS/TK6_Top2_mediated_DSBs/DMSO_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444984/suppl/GSM3444984%5FTK6%5FDMSO%5Frep1.bed.gz)
* [GSE121742 (sBLISS/TK6_Top2_mediated_DSBs/DMSO_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444986/suppl/GSM3444986%5FTK6%5FDMSO%5Frep2.bed.gz)
* [GSE121742 (sBLISS/TK6_Top2_mediated_DSBs/ETO_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444985/suppl/GSM3444985%5FTK6%5FETO%5Frep1.bed.gz)
* [GSE121742 (sBLISS/TK6_Top2_mediated_DSBs/ETO_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444987/suppl/GSM3444987%5FTK6%5FETO%5Frep2.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_etoposide_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322006/suppl/GSM4322006%5FBB76%5FCaco2%5FETO%5F1B%5FGTCGTATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_etoposide_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322010/suppl/GSM4322010%5FBB78%5FCaco2%5FETO%5F2B%5FGTCGTATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_etoposide_rep3)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322014/suppl/GSM4322014%5FBB82%5FCaco2%2DETO%2D3B%5FGTCGTATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_etoposide_rep4)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322022/suppl/GSM4322022%5FBB191%5FCaco2ETO%5FCTAACGAG.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_minus_Ecoli_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322008/suppl/GSM4322008%5FBB77%5FCaco2%5FEcoliMut%5F1D%5FTGATGATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_minus_Ecoli_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322012/suppl/GSM4322012%5FBB79%5FCaco2%5FEcoliMut%5F2D%5FTGATGATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_minus_Ecoli_rep3)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322016/suppl/GSM4322016%5FBB83%5FCaco2%2DEcoliMut%2D3D%5FTGATGATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_minus_Ecoli_rep4)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322024/suppl/GSM4322024%5FBB192%5FCaco2coliWT%5FGCTTGTCA.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_plus_Ecoli_rep1)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322007/suppl/GSM4322007%5FBB77%5FCaco2%5FEcoliWT%5F1C%5FACGACATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_plus_Ecoli_rep2)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322011/suppl/GSM4322011%5FBB79%5FCaco2%5FEcoliWT%5F2C%5FACGACATC.bed.gz)
* [GSE145594 (sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_pks_plus_Ecoli_rep4)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4322nnn/GSM4322023/suppl/GSM4322023%5FBB192%5FCaco2coliMut%5FAGCCATCA.bed.gz)

For the GEO accession code GSE139011, download the [tar file](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139011&format=file). For each of `G5`, `Imatinib`, `SN`, and `Talazoparib`, unpack and extract the two biological replicates for each of `6h` and `48h` treatment. Place them into their folder structures as follows:

* [data/00_Breakage/SSB/G5/6h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/G5/6h) and [data/00_Breakage/SSB/G5/48h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/G5/48h)
* [data/00_Breakage/SSB/Imatinib/6h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/Imatinib/6h) and [data/00_Breakage/SSB/Imatinib/48h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/Imatinib/48h)
* [data/00_Breakage/SSB/SN/6h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/SN/6h) and [data/00_Breakage/SSB/SN/48h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/SN/48h)
* [data/00_Breakage/SSB/Talazoparib/6h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/Talazoparib/6h) and [data/00_Breakage/SSB/Talazoparib/48h](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/Talazoparib/48h)

For each folder, concatenate the files and retain distinct DNA strand break positions. To do this, please run the below file located in [data/00_Breakage/SSB/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/SSB/)

```R
Rscript unpack_files.R
```

For the WRN-associated DNA breakages, download `Source Data Fig. 1` from [this paper](https://doi.org/10.1038/s41586-020-2769-8). In the `1e` subsheet, you will find the below datasets. Copy and paste the `Chromosome`, `Start`, and `Stop` header for each and create a `.csv`. Then, rename `Stop` to `End`. Place them into their folder structures as follows:

* **ENDseq peaks in HCT116 cells (siWRN)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_siWRN](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_siWRN) as `HCT116_siWRN.csv`
* **ENDseq peaks in HCT116 cells (shWRN)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_shWRN](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_shWRN) as `HCT116_shWRN.csv`
* **ENDseq peaks in KM12 cells (siWRN)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_siWRN](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_siWRN) as `KM12_siWRN.csv`
* **ENDseq peaks in KM12 cells (shWRN)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_shWRN](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_shWRN) as `KM12_shWRN.csv`

For the remaining WRN-associated DNA breakages, download `Source Data Fig. 3` from [this paper](https://doi.org/10.1038/s41586-020-2769-8). In the `3c` subsheet, you will find the below datasets. Copy and paste the `Chromosome`, `Start`, and `Stop` header for each and create a `.csv`. Then, rename `Stop` to `End`. Place them into their folder structures as follows:

* **ENDseq peaks with APH+ATRi treatment (HCT116)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_APH-ATR-inhib](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/HCT116_APH-ATR-inhib) as `HCT116_APH-ATR-inhib.csv`
* **ENDseq peaks with APH+ATRi treatment (KM12)** header into [data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_APH-ATR-inhib](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/ENDseq/WRN_loss_DSBs/KM12_APH-ATR-inhib) as `KM12_APH-ATR-inhib.csv`

For apoptotic DNA breakages, download [Table S2 from the supporting information](https://doi.org/10.1371/journal.pone.0026054.s009). The tar file contains the excel file Table_S2-Apoptotic_Breakpoints.xlsx. Create a new `AHH001.csv` file, and copy over the `Chromosome`, `Start`, and `End` columns from the `Ahh001` subsheet. Save this into [data/00_Breakage/Apoptoseq/Apoptotic_DSBs/HL-60-leukemic/AHH001/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Apoptoseq/Apoptotic_DSBs/HL-60-leukemic/AHH001/) folder. Repeat the same for `AHH002.csv` into [data/00_Breakage/Apoptoseq/Apoptotic_DSBs/HL-60-leukemic/AHH002/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Apoptoseq/Apoptotic_DSBs/HL-60-leukemic/AHH002/) folder.

For the human recombination map, clone [this GitHub repo](https://github.com/adimitromanolakis/geneticMap-GRCh37.git) into [data/00_Breakage/Recombination_rates/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Recombination_rates/). Then, change all file names from `genetic_map_GRCh37_chr$num.txt.gz` to `chr$num.txt.gz`, where `$num` is any number from `1` to `22` for the autosomes, and uncompress the files.

### Enzymatic DNA breakages

* [GSE149709 (EcoRV_HeLa_cells)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78172/suppl/GSE78172%5FEco.first%5Fnodup.cov.bed.gz)
* [SRX7808529 (Twist_library/C25cl48-cells)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11188283&display=alignment). For this dataset only, please follow the same downloading process as from [the Ultrasonication DNA breakages section](#ultrasonication-dna-breakages).

For the GEO accession code GSE149709, download the [GSE78172_AsiSI.tar.gz file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78172/suppl/GSE78172%5FAsiSI.tar.gz). Unpack and uncompress the `AsiSI_restriction_sites.single.bed` file and move to [data/00_Breakage/Enzymatic/AsiSI_AID-DIvA_cells/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Enzymatic/AsiSI_AID-DIvA_cells/). Then, inside this folder, run the below R script to process the enzymatic breaks.

```R
Rscript Process.R
```

For the GEO accession code GSE139011, download the [tar file](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139011&format=file). There are 6 files of interest as outlined below.

* **Nt_BbvCI_K562_cells/Sap** Uncompress the `GSM4126200_Nt.BbvCI.Sap.B$num.without.polyA.untailing.bed.gz` files, where `$num` is a number between `1` and `3`. Combine all 3 files with 
* **Nt_BbvCI_K562_cells/NT** For the untreated sample, Uncompress the `GSM4126200_Nt.BbvCI.B$num.without.polyA.untailing.bed.gz` files, where `$num` is a number between `1` and `3`.

Place them into their respective folders, where `.Sap` files go into [data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap/) and the other 3 files go to [data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/NT/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/NT/). When done, go one level up [data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/), and run the below bash script to process the enzymatic breaks.

```bash
bash process_files.sh
```

### Transcription Factor data

We retrieved 247 core-validated vertebrate transcript factor binding sites (TFBS) from the [JASPAR 2024 database](https://jaspar.elixir.no/downloads/).

* Download the [vertebrate dataset from here](https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt)
* Download the [bed files from here](https://jaspar.elixir.no/download/data/2024/bed.tar.gz) as "jaspar_beds"

Unpack and extract the relevant files from above. Place the contents into [data/TFBS/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/JASPAR/) folder. 

Alternatively, run the below Rscript or run it as part the setup file at the end of this `README` file.

```R
Rscript jaspar_download.R
```

### Quantum mechanical k-mer parameters

The quantum mechanical hexameric parameters `denergy.txt.gz` can be downloaded from [DNAkmerQM](https://github.com/SahakyanLab/DNAkmerQM/tree/master/6-mer). Uncompress the file and place it into [data/04_QM_parameters/6mer/Raw/denergy.txt](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/04_QM_parameters/6mer/Raw/denergy.txt) folder.

### Epigenome marks

**ATACseq datasets**
* [ENCFF296ZZB (ATACseq/HCT116)](https://www.encodeproject.org/files/ENCFF296ZZB/@@download/ENCFF296ZZB.bed.gz)
* [ENCFF558BLC (ATACseq/K562)](https://www.encodeproject.org/files/ENCFF558BLC/@@download/ENCFF558BLC.bed.gz)
* [ENCFF821OEF (ATACseq/MCF7)](https://www.encodeproject.org/files/ENCFF821OEF/@@download/ENCFF821OEF.bed.gz)

**FAIREseq datasets**
* [ENCFF001UYM (FAIREseq/HeLa_S3)](https://www.encodeproject.org/files/ENCFF001UYM/@@download/ENCFF001UYM.bed.gz)
* [ENCFF001UYS (FAIREseq/K562)](https://www.encodeproject.org/files/ENCFF001UYS/@@download/ENCFF001UYS.bed.gz)
* [ENCFF001UYW (FAIREseq/MCF7)](https://www.encodeproject.org/files/ENCFF001UYW/@@download/ENCFF001UYW.bed.gz)

**DNaseseq datasets**
* [ENCFF579UXQ (DNaseseq/Caco2)](https://www.encodeproject.org/files/ENCFF579UXQ/@@download/ENCFF579UXQ.bed.gz)
* [ENCFF240LRP (DNaseseq/HCT116)](https://www.encodeproject.org/files/ENCFF240LRP/@@download/ENCFF240LRP.bed.gz)
* [ENCFF024UJN (DNaseseq/HeLa_S3)](https://www.encodeproject.org/files/ENCFF024UJN/@@download/ENCFF024UJN.bed.gz)
* [ENCFF773SFA (DNaseseq/HL60)](https://www.encodeproject.org/files/ENCFF773SFA/@@download/ENCFF773SFA.bed.gz)
* [Hypersensitive clusters (DNaseseq/HyperSensitive_Clusters)](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz)
* [ENCFF274YGF (DNaseseq/K562)](https://www.encodeproject.org/files/ENCFF274YGF/@@download/ENCFF274YGF.bed.gz)
* [ENCFF438LQM (DNaseseq/MCF7)](https://www.encodeproject.org/files/ENCFF438LQM/@@download/ENCFF438LQM.bed.gz)

**Chipseq datasets**
* [ENCFF579UXQ (Chipseq/H3K4me2_NHEK_cells)](https://www.encodeproject.org/files/ENCFF579UXQ/@@download/ENCFF579UXQ.bed.gz)

**Histone**
* [GSM5501177 (Chipseq/Histone/H3K4me2_DMSO_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501177/suppl/GSM5501177%5FDMSO%5FH3K4%5Fme2%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM5501178 (Chipseq/Histone/H3K4me2_LSD1-inhib_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501178/suppl/GSM5501178%5Finib%5F2%5F5%5FH3K4%5Fme2%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM5501179 (Chipseq/Histone/H3K4me3_DMSO_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501179/suppl/GSM5501179%5FDMSO%5FH3K4%5Fme3%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM5501180 (Chipseq/Histone/H3K4me3_LSD1-inhib_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501180/suppl/GSM5501180%5Finib%5F2%5F5%5FH3K4%5Fme3%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM5501175 (Chipseq/Histone/H3K4me1_DMSO_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501175/suppl/GSM5501175%5FDMSO%5FH3K4%5Fme1%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM5501176 (Chipseq/Histone/H3K4me1_LSD1-inhib_GBM_tumor_initiating_cell_line)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5501nnn/GSM5501176/suppl/GSM5501176%5Finib%5F2%5F5%5FH3K4%5Fme1%5Fvs%5Finput%5Fpeaks.txt.gz)
* [GSM1035424 (Chipseq/SUMO1/Proliferative_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035424/suppl/GSM1035424%5Fprolif%5FSUMO1%5Frep1%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035433 (Chipseq/SUMO1/Ras-induced-senescent_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035433/suppl/GSM1035433%5Fras%5FSUMO1%5Frep1%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035426 (Chipseq/SUMO2/Proliferative_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035426/suppl/GSM1035426%5Fprolif%5FSUMO2%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035435 (Chipseq/SUMO2/Ras-induced-senescent_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035435/suppl/GSM1035435%5Fras%5FSUMO2%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035427 (Chipseq/Ubc9/Proliferative_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035427/suppl/GSM1035427%5Fprolif%5FUbc9%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035436 (Chipseq/Ubc9/Ras-induced-senescent_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035436/suppl/GSM1035436%5Fras%5FUbc9%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035441 (Chipseq/PIASy/Proliferative_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035441/suppl/GSM1035441%5Fprolif%5FPIASy%5FChIPSeq%5Fpeaks.bed.gz)
* [GSM1035442 (Chipseq/PIASy/Ras-induced-senescent_fibro_WI38_cells)](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1035nnn/GSM1035442/suppl/GSM1035442%5Fras%5FPIASy%5FChIPSeq%5Fpeaks.bed.gz)

Download the two NhekH3k4me3 replicas from the GEO accession number GSM945175 [rep1](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM945nnn/GSM945175/suppl/GSM945175%5Fhg19%5FwgEncodeUwHistoneNhekH3k4me3StdPkRep1.narrowPeak.gz) and [rep2](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM945nnn/GSM945175/suppl/GSM945175%5Fhg19%5FwgEncodeUwHistoneNhekH3k4me3StdPkRep2.narrowPeak.gz). Concatenate the two files, retain unique ranges, and keep one file in [Chipseq/Histone/H3K4me3_NHEK_cells](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/01_Epigenome/Chipseq/Histone/H3K4me3_NHEK_cells) folder.

Download the two NhekH3k36me3 replicas from the GEO accession number GSM945175 [rep1](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM945nnn/GSM945174/suppl/GSM945174%5Fhg19%5FwgEncodeUwHistoneNhekH3k36me3StdPkRep1.narrowPeak.gz) and [rep2](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM945nnn/GSM945174/suppl/GSM945174%5Fhg19%5FwgEncodeUwHistoneNhekH3k36me3StdPkRep2.narrowPeak.gz). Concatenate the two files, retain unique ranges, and keep one file in [Chipseq/Histone/H3K36me3_NHEK_cells](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/01_Epigenome/Chipseq/Histone/H3K36me3_NHEK_cells) folder.

The following 9 histone marks were obtained from GEO accession [GSE29611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29611) with individual files already saved in the appropriate file structure, as also used in the [DSBCapture study](https://doi.org/10.1038/nmeth.3960), including: H2AZ_NHEK_cells, EZH2_NHEK_cells, H3K27ac_NHEK_cells, H3K4me2_NHEK_cells, H3K27me3_NHEK_cells, H3K79me2_NHEK_cells, H3K4me1_NHEK_cells, H3K9me3_NHEK_cells, POL2B.

### ENCODE Blacklist regions

Please download the ENCODE Blacklist regions for [hg19 from here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz) and [for hg38 from here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz). Uncompress them, and name them `hg19.bed` and `hg38.bed`, respectively.

Place the contents into [data/blacklists/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/blacklists/).

### Hi-C annotations

We retrieved the hg19-annotated Hi-C subcompartment data of the K562 and HeLa cell lines [from here](https://cmu.app.box.com/s/n4jh3utmitzl88264s8bzsfcjhqnhaa0/folder/86847304302).

Place the contents into [data/HiC_annotations/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/HiC_annotations/).

## Liftover files

We processed all datasets in the reference genome version used as per the deposition. For the TFBS, we lifted them over from hg38 to hg19.

* [hg38 to hg19](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz)

Unpack and extract the relevant files. Place the contents into [data/liftover/](https://github.com/SahakyanLab/DNAFragility_dev/tree/master/data/liftover/) folder. 

## Other notes

* All cpp files are interfaced *via* the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) library in R with `omp.h` when possible. Please ensure you have this installed.
* `RcppArmadillo.h` and `RcppEigen.h` are necessary for the feature extraction process. Please ensure you have this installed. By default, will not use it in case you have not installed it.

## Run all setup files

If you wish to run all setups, including all the aforementioned bash scripts, please run the below bash script. 

Please note that many of the calculations were computationally intensive. Most things were run in parallel in smaller batches. However, if you submit the below bash script, it runs all scripts sequentially. This can **take several months** to complete. 
Most tasks take up several tens to hundreds of GBs worth of RAM. The entire study requires between 2-4 TB of hard drive space.

You may need to monitor your memory usage, memory cache, and swap to ensure calculations run smoothly.

```bash
bash run_all_setup_files.sh
```