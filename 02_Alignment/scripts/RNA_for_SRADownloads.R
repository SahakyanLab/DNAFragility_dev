# nature paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5161557/
# select downloads from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP010710&o=acc_s%3Aa
# for each randomly selected Run, we choose the first Run of that given experiment
# only select the data set if it is deposited as aligned raw paired-end reads
# E.g.: Experiment ERX1097982 -> First Run ERR1019042 -> download all aligned sequencing runs

set.seed(1234)
rng <- sample(x = 400, size = 5, replace = FALSE)
print(rng)

# experiments only from center name: harvard medical school, department of genetics
set.seed(1997)
rng <- sample(x = 116, size = 20, replace = FALSE)
print(rng)