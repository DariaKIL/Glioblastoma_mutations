## VCF filtering pipeline for tumor samples (WES) 

## Overview

This project provides a pipeline for processing VCF files, specifically for analyzing and filtrating genetic variations. Originally created for analyzing three glioblastoma patients, this pipeline can be adapted for any disease and any number of patients. BAM files were analyzed using the [Strelka pipeline](https://github.com/Illumina/strelka/tree/v2.9.x) and annotated using the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html).

## Project Structure

- `Snakefile`: Manages the pipeline with Snakemake.
- `data_filtration.ipynb`: Jupyter notebook for filtering data.
- `dag.png`: Illustration of the pipeline.
- `README.md`: Project documentation.

## Prerequisites

Ensure you have the following tools installed:

- bcftools
- Snakemake
- Python 3
- Jupyter Notebook

## Usage

1. Clone the repository:
    ```sh
    git clone https://github.com/DariaKIL/Glioblastoma_mutations.git
    ```

2. If necessary, create and activate a virtual environment:
    ```sh
    python -m venv myenv
    source myenv/bin/activate
    ```

3. Navigate to the project directory:
    ```sh
    cd /path/to/Glioblastoma_mutations
    ```

4. Make sure your annotated VCF file (if it contains multiple samples) is named according to this type: 'samplename1_samplename2.anno.vcf.' Run the Snakemake workflow to generate tables for each sample:
    ```sh
    snakemake
    ```

5. Open `data_filtration.ipynb` and execute all cells to generate filtered tables for each sample. You can easily change filtering criteria according your task.

6. The filtered tables will contain a number of rows allowing you to manually review and analyze the obtained mutations based on your task, literature data, etc.



