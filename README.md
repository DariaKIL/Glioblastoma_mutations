## VCF filtering pipeline for tumor samples (WES) 

## Overview

This project provides a pipeline for processing VCF files, specifically for analyzing and filtrating genetic variations. Originally created for analyzing three glioblastoma patients, this pipeline can be adapted for any disease and any number of patients. BAM files were analyzed using the [Strelka pipeline](https://github.com/Illumina/strelka/tree/v2.9.x) and annotated using the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html).

## Project Structure

- Snakefile: Snakemake workflow for processing combined VCF files with multiple samples.
- Snakefile_separate: Snakemake workflow for processing separate VCF files per sample using gnomADv4 allele frequencies.
- data_filtration.ipynb: Jupyter notebook for filtering data using gnomADv3 allele frequencies.
- data_filtration_gnomADv4.ipynb: Jupyter notebook for filtering data using gnomADv4 allele frequencies.
- dag.png: Illustration of the pipeline.
- README.md: Project documentation.

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

4. **For combined VCF files:** Make sure your annotated VCF file (if it contains multiple samples) is named according to this type: **samplename1_samplename2.anno.vcf** <br> Run the Snakemake workflow to generate tables for each sample:
    ```sh
    snakemake -s Snakefile
    ```
    **For separate VCF files per sample:** Place your separate annotated VCF files (named samplename.anno.vcf) in a data folder.

Run the Snakemake workflow designed for separate files:
    ```sh
    snakemake -s Snakefile_separate
    ```

5. Filtering:
	•	For filtering based on gnomAD v3, use data_filtration.ipynb.
	•	For filtering based on gnomAD v4, use data_filtration_gnomADv4.ipynb.

Open the appropriate notebook in Jupyter and execute all cells to generate filtered tables for each sample. You can easily adjust the filtering criteria according to your analysis needs.

6. The filtered tables will contain a number of rows allowing you to manually review and analyze the obtained mutations based on your task, literature data, etc.



