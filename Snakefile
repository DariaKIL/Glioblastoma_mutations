import os
import glob

# Define the working directory
work_dir = os.getcwd()
print("work_dir: ", work_dir)

# Find all .anno.vcf.gz files in the data directory
vcf_files = glob.glob(os.path.join(work_dir, 'data/*.anno.vcf.gz'))
print('VCF files:', vcf_files)

# Extract sample names and create a dictionary mapping files to sample names
samples_dict = {}
for vcf in vcf_files:
    samples = os.path.basename(vcf).split('.')[0].split('_')
    for sample in samples:
        samples_dict[sample] = vcf
sample_names = list(samples_dict.keys())
print('Sample names:', sample_names)
print('Samples dictionary:', samples_dict)

rule all:
    input:
        expand(os.path.join(work_dir, 'output/{sample}.tsv'), sample=sample_names) +
        expand(os.path.join(work_dir, 'output/{sample}_vep.tsv'), sample=sample_names)

rule split_vcf:
    input:
        vcf=lambda wildcards: samples_dict[wildcards.sample]
    output:
        os.path.join(work_dir, 'temp/{sample}.vcf')
    shell:
        """
        bcftools view -Ov -s {wildcards.sample} {input.vcf} -o {output}
        """

rule fix_GERP_notation:
    input:
        os.path.join(work_dir, 'temp/{sample}.vcf')
    output:
        os.path.join(work_dir, 'temp/{sample}_modified.vcf')
    shell:
        """
        sed 's/GERP++_RS_rankscore/GERP_RS_rankscore/g' {input} > {output}
        """

rule to_tsv_main:
    input:
        os.path.join(work_dir, 'temp/{sample}_modified.vcf')
    output:
        os.path.join(work_dir, 'output/{sample}.tsv')
    shell:
        """
        bcftools view -i 'FILTER="PASS"' {input} | bcftools query -f '%CHROM %POS %ID %REF %ALT %QUAL [%GT %DP]\n' -s {wildcards.sample} > {output}
        """

rule to_vcf_VEP:
    input:
        os.path.join(work_dir, 'temp/{sample}_modified.vcf')
    output:
        os.path.join(work_dir, 'output/{sample}_vep.tsv')
    shell:
        """
        bcftools +split-vep {input} -f '%Allele %Consequence %IMPACT %SYMBOL %Gene %Feature_type %Feature %BIOTYPE %EXON %INTRON %HGVSc %HGVSp %cDNA_position %CDS_position %Protein_position %Amino_acids %Codons %Existing_variation %DISTANCE %STRAND %FLAGS %VARIANT_CLASS %SYMBOL_SOURCE %HGNC_ID %CANONICAL %MANE_SELECT %MANE_PLUS_CLINICAL %TSL %APPRIS %RefSeq %SOURCE %GENE_PHENO %SIFT %PolyPhen %DOMAINS %HGVS_OFFSET %HGVSg %MAX_AF %MAX_AF_POPS %CLIN_SIG %SOMATIC %PHENO %PUBMED %VAR_SYNONYMS %MOTIF_NAME %MOTIF_POS %HIGH_INF_POS %MOTIF_SCORE_CHANGE %TRANSCRIPTION_FACTORS %CADD_phred %DEOGEN2_pred %FATHMM_pred %GERP_RS_rankscore %LRT_Omega %LRT_pred %MetaSVM_pred %MutationTaster_pred %PROVEAN_pred %Polyphen2_HVAR_pred %PrimateAI_pred %SIFT_pred %phastCons17way_primate_rankscore %ada_score %rf_score %DisGeNET %SpliceAI_cutoff %SpliceAI_pred_DP_AG %SpliceAI_pred_DP_AL %SpliceAI_pred_DP_DG %SpliceAI_pred_DP_DL %SpliceAI_pred_DS_AG %SpliceAI_pred_DS_AL %SpliceAI_pred_DS_DG %SpliceAI_pred_DS_DL %SpliceAI_pred_SYMBOL %Mastermind_MMID3 %Mastermind_counts %ExACpLI %LoF %LoF_filter %LoF_flags %LoF_info %CLINVAR %CLINVAR_CLNSIG %CLINVAR_CLNDN %CLINVAR_CLNSIGCONF %dgWGS %dgWGS_AF %dgWGS_AC %dgWGS_AC_Het %dgWGS_AC_Hom %dgWGS_AC_Hemi %dgWGS_HWE %dgWES %dgWES_AF %dgWES_AC %dgWES_AC_Het %dgWES_AC_Hom %dgWES_AC_Hemi %dgWES_HWE %gnomADv3 %gnomADv3_AF %gnomADv3_AC_raw %gnomADv3_AC_XY %gnomADv3_AC_XX %gnomADv3_nhomalt' -A tab > {output}
        """

