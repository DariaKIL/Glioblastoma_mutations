# WES tumor samples - VCF filtering Pipeline

## Overview

This project provides a pipeline for processing VCF files, specifically for analyzing genetic variations and annotations. 
This project was created for analizing three patients with glioblastoma, but it can be used for any disease and any numbers of patients.
BUM files were analized through strelka pipeline https://github.com/Illumina/strelka/tree/v2.9.x and annotated through Ensembl Variant Effect Predictor (VEP) 

 ```sh
    /usr/local/bin/configureStrelkaGermlineWorkflow.py \
    --bam sample1_sample2_sample3/sample1.dedup.sorted.bam \
    --bam sample1_sample2_sample3/sample2.dedup.sorted.bam \
    --bam sample1_sample2_sample3/sample3.dedup.sorted.bam \
    --referenceFasta /path/to/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --runDir /path/to/dir/strelka2 \
    --exome \
    --ploidy /path/to/dir/gender.vcf.gz \
    --callRegions /path/to/dir/targets38.bed.gz
    ```

## Project Structure
- Snakemake: file for managing the pipeline.
- filtrating.ipynb
- dag.pdf: illustration of this pipeline.
- README.md: This file.

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



