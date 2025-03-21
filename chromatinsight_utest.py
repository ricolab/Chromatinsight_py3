
# Both pip (sudo apt install python3-pip) and 
# pandas (pip install pandas) are needed
import os
import chromatinsight as ci

# Define file paths
myRegionFile = "./utest/input/regionFile_chrX_selection.bed"
myGroupingFile = "./utest/input/grouping.txt"
myOutputFolder = "./utest/output"

# Ensure input files exist
import sys

def check_file_exists(file):
    if not os.path.isfile(file):
        sys.exit(f"Error: Required input file '{file}' not found.")

check_file_exists(myRegionFile)
check_file_exists(myGroupingFile)

        
# Create output directory if it doesn't exist
os.makedirs(myOutputFolder, exist_ok=True)

if not os.path.isdir(myOutputFolder): os.mkdir(myOutputFolder)

# RF_seed = 0 and label_seed = 0 are set to generate a deterministic random behaviour
# Run testPrediction with deterministic behavior
resultObserved = ci.testPrediction(
    groupingFile=myGroupingFile,   # (str) Path to the file defining groupings
    regionFile=myRegionFile,       # (str) Path to the BED file with genomic regions
    histmod="ac",                  # (str) Histone modification type (e.g., "ac" for acetylation)
    chrom="chrX",                  # (str) Chromosome to analyze (e.g., "chrX")
    outputFolder=myOutputFolder,   # (str) Directory where output files will be saved
    output="mono_ac_chrX_real.txt", # (str) Name of the output file
    totRandomStates=11,            # (int) Total number of random states to consider
    randomize=False,               # (bool) If False, runs with real labels (no randomization)
    verbose=True,                  # (bool) If True, prints detailed execution logs
    RF_seed=0                       # (int) Random seed for the Random Forest model (ensures reproducibility)
)

# Run testPrediction with randomization
resultRnd = ci.testPrediction(
    groupingFile=myGroupingFile,   # (str) Path to the file defining groupings
    regionFile=myRegionFile,       # (str) Path to the BED file with genomic regions
    histmod="ac",                  # (str) Histone modification type (e.g., "ac" for acetylation)
    chrom="chrX",                  # (str) Chromosome to analyze (e.g., "chrX")
    outputFolder=myOutputFolder,   # (str) Directory where output files will be saved
    output="mono_ac_chrX_rnd.txt", # (str) Name of the output file
    totRandomStates=11,            # (int) Total number of random states to consider
    randomize=True,                # (bool) If True, shuffles labels for a randomized test
    verbose=True,                  # (bool) If True, prints detailed execution logs
    label_seed=0,                   # (int) Random seed for label shuffling (ensures reproducibility)
    RF_seed=0                       # (int) Random seed for the Random Forest model (ensures reproducibility)
)
