# WES tumor samples - VCF filtering Pipeline

## Overview

This project provides a pipeline for processing VCF files, specifically for analyzing genetic variations and annotations. 
This project was created for analizing three patients with glioblastoma, but it can be used for any disease and any numbers of patients.
BAM files were analized through strelka pipeline https://github.com/Illumina/strelka/tree/v2.9.x and annotated through Ensembl Variant Effect Predictor (VEP) 

## Project Structure
- Snakemake: file for managing the pipeline.
- data_filtration.ipynb
- dag.png: illustration of this pipeline.
- README.md

## Prerequisites

Make sure you have the following tools installed:

	•	bcftools
	•	Snakemake
	•	Python 3

## Usage
1. Clone the repository:
    ```sh
    git clone https://github.com/DariaKIL/Glioblastoma_mutations.git
    ```
2. If necessary create and activate a virtual environment:
    ```sh
    python -m venv myenv
    source myenv/bin/activate
    ```
3. Navigate to the project directory:
    ```sh
    cd /path/to/Glioblastoma_mutations
    ```
4. Run the Snakemake workflow:
    ```sh
    snakemake 
    ```



