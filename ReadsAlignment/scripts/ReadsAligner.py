from Bio import SeqIO
import itertools
import regex as re
import pandas as pd
import edlib
import subprocess
import os
import progressbar
import sys
import math

# read arguments from job submission
nr_cur = sys.argv[1]
ind = int(sys.argv[2])
ref_path = sys.argv[3]
reads_path = sys.argv[4]
my_path = sys.argv[5]

# define path and chromosome file
path = "/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/"
chromosome = f"chr{nr_cur}"

# load reference genome
for record in SeqIO.parse(f"{path}data/ref/{ref_path}/{chromosome}.fasta", "fasta"):
    record_ref = record.seq

# obtain all files from dir and sort them
search_folder = f"{path}data/reads/{reads_path}/fasta_chunks/"

# remove hidden files
FILE = f"{path}data/reads/{reads_path}/fasta_chunks/.DS_Store"
if os.path.exists(FILE):
    subprocess.call(f"rm {FILE}", shell = True)
    
if os.path.exists(f"{path}data/reads/{reads_path}/fasta_chunks/._*"):
    subprocess.call(f"{path}data/reads/{reads_path}/fasta_chunks/._*", shell = True)

results = [os.path.join(root, f) for root, dirs, files in os.walk(search_folder) for f in files] 
sorted_results = sorted(results, key = lambda x: os.path.splitext(x)[1])

# split files into 2 batches to reduce memory consumption
first_half = math.ceil(len(sorted_results)/2)
results_batch_1 = sorted_results[:first_half]
results_batch_2 = sorted_results[first_half:]
    
# obtain total two-mer count in reference genome
two_mers = [''.join(x) for x in itertools.product('ACGT', repeat = 2)]
two_mers = {i:0 for i in two_mers}

ref_seq = str(record_ref)
for keys in two_mers:
    string_matches = [m for m in re.finditer(pattern = keys, string = ref_seq)]
    two_mers[keys] = len(string_matches)
    
with open(f'{path}data/two_mers/{reads_path}/{chromosome}_ref.csv', 'w') as f:
    f.write('value,count\n')
    for key, value in two_mers.items():
        f.write(f'{key},{value}\n')
    
# assign correct batch of files to run alignment on
if ind == 1:
    run_alignment = results_batch_1
elif ind == 2:
    run_alignment = results_batch_2
    
# text progress bar
bar = progressbar.ProgressBar(maxval = len(run_alignment)+1, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()
    
# initialise lists
bp_start_pos = [] 
lev_dist = []
iteration = 0

for txt_file in run_alignment:   
    with open(txt_file, mode = 'r') as handle:
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(handle, 'fasta'):
            description = record.description

            if(int(description.split()[0].split("/")[1]) <= 2):
                # Extract individual parts of the FASTA record
                sequence = record.seq
                sequence_revcomp = sequence.reverse_complement()
                both_reads = [sequence, sequence_revcomp]

                # starting position of given read
                genome_pos = int(description.split()[3].split("=")[1])
                start = genome_pos-1
                end = start+100
                ref_seq_fwd = record_ref[start:end]

                # use edit distance to choose the "correct" alignment between read and ref seq
                levenshtein_dist = [edlib.align(str(ref_seq_fwd), str(sequence), 
                                                mode = "HW", task = "path")['editDistance'],
                                    edlib.align(str(ref_seq_fwd), str(sequence_revcomp), 
                                                mode = "HW", task = "path")['editDistance']]
                
                # read with lowest levenshtein distance
                read = both_reads[levenshtein_dist.index(min(levenshtein_dist))]
                
                # align the first 2 nucleotides of read against ref seq
                if levenshtein_dist.index(min(levenshtein_dist)) == 0:
                    # fwd strand
                    if (str(ref_seq_fwd)[0:2] == str(sequence)[0:2]):
                        bp_start_pos.append(start-1)
                        lev_dist.append(min(levenshtein_dist))
                else:
                    # rev strand
                    if (str(ref_seq_fwd)[-2:] == str(sequence_revcomp)[-2:]):
                        bp_start_pos.append(end-1)
                        lev_dist.append(min(levenshtein_dist))

        # update iteration
        iteration += 1
        bar.update(iteration+1)
                
# finish text progress bar
bar.finish()

# save output as data frame
df = pd.DataFrame({
    'bp_start_pos': bp_start_pos,
    'lev_dist':  lev_dist
})
df['freq'] = df.groupby('bp_start_pos')['bp_start_pos'].transform('count')
df = df.drop_duplicates(subset = 'bp_start_pos')

# save positions as txt file
if ind == 1:    
    df.to_csv(f'{path}data/reads/{reads_path}/breakpoint_positions/{chromosome}/alignment_file_1.txt', 
                index = False, 
                sep = ',', 
                mode = 'a')
elif ind == 2:
    df.to_csv(f'{path}data/reads/{reads_path}/breakpoint_positions/{chromosome}/alignment_file_2.txt', 
                index = False, 
                sep = ',', 
                mode = 'a')