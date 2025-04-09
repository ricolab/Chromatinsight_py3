# Chromatinsight
An application to extract meaningful information from datasets of ChIP-seq experiments.

The code in Python represents the core part of the application, which uses the random forest model in the sklearn library (scikit-learn) to detect differential features between two groups of ChIP-seq experiments binarised using ChromHMM.

This software has been tested to detect the human sex dimorphism between two groups (males and females). The output can be further analysed using the R library associated, chromatinsight.tools https://github.com/ricolab/chromatinsight.tools .

**Input**:
* The output text files of ChromHMM binarisation function BinarizeBed (one for each chromosome) that assigns presence (1) or absence (0) of epigenomic features at intervals of 200bp (see https://ernstlab.biolchem.ucla.edu/software-and-resources/chromhmm ). At least ten samples in each group are needed to provide enough statistical power.
* A bed file with the genomic regions of interest (such as TADs).

**Output**:
* A text file with the degree of dimorphism for each genomic region and each trial of the algorithm.
* A second text file (optional) with the degree of dimorphism of each genomic region after randomising the sample labels (important to display the random behaviour and calculating the FDR).

A **basic unit test** is provided, which also serves as example of how to run the code; it includes files for twenty samples (ten men and ten women). To process it, use the data in ./utest/input and use the code provided in chromatinsight_utest.py .

**Requirements**:\
Python 3.x\
pandas installed, see https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html \
sklearn installed, see https://scikit-learn.org/stable/install.html \
Quick way to install both:\
pip install pandas\
pip install scikit-learn

---

# Dockerfile for Chromatinsight_py3

We’ve created a Dockerfile that allows the application to run in an isolated Docker environment. The Dockerfile includes all necessary dependencies and follows best practices for Docker image creation.

## Usage Instructions

### 1. Build the Docker Image

Clone the repository and build the image:

```bash
git clone https://github.com/ricolab/Chromatinsight_py3.git
cd Chromatinsight_py3
docker build -t chromatinsight:1.0 .
```

### 2. Run the Container

To run the included unit test:

```bash
docker run --rm chromatinsight:1.0
```

### 3. Run Custom Analyses

To run with your own data, mount a volume containing your input files:

```bash
docker run --rm -v /path/to/your/data:/data chromatinsight:1.0 python /app/your_script.py --input /data/input --output /data/output
```

Replace `your_script.py` with the script you want to run, adjusting parameters as needed.

### 4. Enter the Container for Interactive Use

```bash
docker run --rm -it chromatinsight:1.0 /bin/bash
```

## Possible Customizations

1. **Python Version**: If you need a specific version of Python, change the base image (e.g., `python:3.8-slim`).

2. **Additional Dependencies**: Add more libraries to the `pip install` line if your project requires them.

3. **Persistent Volumes**: For long-running analyses, consider mounting a named volume:

```bash
docker run --rm -v chromatinsight_data:/data chromatinsight:1.0
```

4. **Environment Variables**: Use environment variables to configure runtime behavior:

```bash
docker run --rm -e NUM_THREADS=4 chromatinsight:1.0
```

## Notes

- The Dockerfile is optimized to balance image size and functionality.
- Security has been prioritized by running the app as a non-root user.
- For large datasets, use mounted volumes for input/output to avoid image bloat.
- In production environments, consider setting memory and CPU limits appropriately.

This Dockerfile provides a solid foundation for running **Chromatinsight_py3** in any Docker-compatible environment, ensuring consistency and reproducibility for your ChIP-seq data analyses.

You may also find the Docker image in https://hub.docker.com/repository/docker/juanochoteco0/chromatinsight_py3/general
