import os
import glob

# Define the working directory
work_dir = os.getcwd()
#print("work_dir: ", work_dir)

# Find all .anno.vcf.gz files in the data directory
vcf_files = glob.glob(os.path.join(work_dir, 'data/*.anno.vcf.gz'))
#print('VCF files:', vcf_files)

# Extract sample names and create a dictionary mapping files to sample names
samples_dict = {}
for vcf in vcf_files:
    samples = os.path.basename(vcf).split('.')[0].split('_')
    for sample in samples:
        samples_dict[sample] = vcf
sample_names = list(samples_dict.keys())
#print('Sample names:', sample_names)
#print('Samples dictionary:', samples_dict)

rule all:
    input:
        expand(os.path.join(work_dir, 'temp/{sample}_combined.tsv'), sample=sample_names)

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
        os.path.join(work_dir, 'temp/{sample}_main.tsv')
    shell:
        """
        bcftools query -f '%CHROM %POS %ID %REF %ALT %QUAL %FILTER [%GT %DP]\n' -s {wildcards.sample} {input} > {output}
        """

rule to_vcf_VEP:
    input:
        os.path.join(work_dir, 'temp/{sample}_modified.vcf')
    output:
        os.path.join(work_dir, 'temp/{sample}_vep.tsv')
    shell:
        """
        bcftools +split-vep {input} -f '%Allele %Consequence %IMPACT %SYMBOL %Gene %Feature_type %Feature %BIOTYPE %EXON %INTRON %HGVSc %HGVSp %cDNA_position %CDS_position %Protein_position %Amino_acids %Codons %Existing_variation %DISTANCE %STRAND %FLAGS %VARIANT_CLASS %SYMBOL_SOURCE %HGNC_ID %CANONICAL %MANE_SELECT %MANE_PLUS_CLINICAL %TSL %APPRIS %RefSeq %SOURCE %GENE_PHENO %SIFT %PolyPhen %DOMAINS %HGVS_OFFSET %HGVSg %MAX_AF %MAX_AF_POPS %CLIN_SIG %SOMATIC %PHENO %PUBMED %VAR_SYNONYMS %MOTIF_NAME %MOTIF_POS %HIGH_INF_POS %MOTIF_SCORE_CHANGE %TRANSCRIPTION_FACTORS %CADD_phred %DEOGEN2_pred %FATHMM_pred %GERP_RS_rankscore %LRT_Omega %LRT_pred %MetaSVM_pred %MutationTaster_pred %PROVEAN_pred %Polyphen2_HVAR_pred %PrimateAI_pred %SIFT_pred %phastCons17way_primate_rankscore %ada_score %rf_score %DisGeNET %SpliceAI_cutoff %SpliceAI_pred_DP_AG %SpliceAI_pred_DP_AL %SpliceAI_pred_DP_DG %SpliceAI_pred_DP_DL %SpliceAI_pred_DS_AG %SpliceAI_pred_DS_AL %SpliceAI_pred_DS_DG %SpliceAI_pred_DS_DL %SpliceAI_pred_SYMBOL %Mastermind_MMID3 %Mastermind_counts %ExACpLI %LoF %LoF_filter %LoF_flags %LoF_info %CLINVAR %CLINVAR_CLNSIG %CLINVAR_CLNDN %CLINVAR_CLNSIGCONF %dgWGS %dgWGS_AF %dgWGS_AC %dgWGS_AC_Het %dgWGS_AC_Hom %dgWGS_AC_Hemi %dgWGS_HWE %dgWES %dgWES_AF %dgWES_AC %dgWES_AC_Het %dgWES_AC_Hom %dgWES_AC_Hemi %dgWES_HWE %gnomADv3 %gnomADv3_AF %gnomADv3_AC_raw %gnomADv3_AC_XY %gnomADv3_AC_XX %gnomADv3_nhomalt' -A tab > {output}
        """

rule header_to_main:
    input:
        os.path.join(work_dir, 'temp/{sample}_main.tsv')
    output:
        os.path.join(work_dir, 'temp/{sample}_main_with_header.tsv')
    shell: 
        """
        awk 'BEGIN{{print "CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\t{wildcards.sample}_GT\\t{wildcards.sample}_DP"}} {{print}}' {input} > {output}
        """

rule header_to_VEP:
    input:
        os.path.join(work_dir, 'temp/{sample}_vep.tsv')
    output:
        os.path.join(work_dir, 'temp/{sample}_vep_with_header.tsv')
    shell:
        """
        awk 'BEGIN{{print "Allele\\tConsequence\\tIMPACT\\tSYMBOL\\tGene\\tFeature_type\\tFeature\\tBIOTYPE\\tEXON\\tINTRON\\tHGVSc\\tHGVSp\\tcDNA_position\\tCDS_position\\tProtein_position\\tAmino_acids\\tCodons\\tExisting_variation\\tDISTANCE\\tSTRAND\\tFLAGS\\tVARIANT_CLASS\\tSYMBOL_SOURCE\\tHGNC_ID\\tCANONICAL\\tMANE_SELECT\\tMANE_PLUS_CLINICAL\\tTSL\\tAPPRIS\\tRefSeq\\tSOURCE\\tGENE_PHENO\\tSIFT\\tPolyPhen\\tDOMAINS\\tHGVS_OFFSET\\tHGVSg\\tMAX_AF\\tMAX_AF_POPS\\tCLIN_SIG\\tSOMATIC\\tPHENO\\tPUBMED\\tVAR_SYNONYMS\\tMOTIF_NAME\\tMOTIF_POS\\tHIGH_INF_POS\\tMOTIF_SCORE_CHANGE\\tTRANSCRIPTION_FACTORS\\tCADD_phred\\tDEOGEN2_pred\\tFATHMM_pred\\tGERP_RS_rankscore\\tLRT_Omega\\tLRT_pred\\tMetaSVM_pred\\tMutationTaster_pred\\tPROVEAN_pred\\tPolyphen2_HVAR_pred\\tPrimateAI_pred\\tSIFT_pred\\tphastCons17way_primate_rankscore\\tada_score\\trf_score\\tDisGeNET\\tSpliceAI_cutoff\\tSpliceAI_pred_DP_AG\\tSpliceAI_pred_DP_AL\\tSpliceAI_pred_DP_DG\\tSpliceAI_pred_DP_DL\\tSpliceAI_pred_DS_AG\\tSpliceAI_pred_DS_AL\\tSpliceAI_pred_DS_DG\\tSpliceAI_pred_DS_DL\\tSpliceAI_pred_SYMBOL\\tMastermind_MMID3\\tMastermind_counts\\tExACpLI\\tLoF\\tLoF_filter\\tLoF_flags\\tLoF_info\\tCLINVAR\\tCLINVAR_CLNSIG\\tCLINVAR_CLNDN\\tCLINVAR_CLNSIGCONF\\tdgWGS\\tdgWGS_AF\\tdgWGS_AC\\tdgWGS_AC_Het\\tdgWGS_AC_Hom\\tdgWGS_AC_Hemi\\tdgWGS_HWE\\tdgWES\\tdgWES_AF\\tdgWES_AC\\tdgWES_AC_Het\\tdgWES_AC_Hom\\tdgWES_AC_Hemi\\tdgWES_HWE\\tgnomADv3\\tgnomADv3_AF\\tgnomADv3_AC_raw\\tgnomADv3_AC_XY\\tgnomADv3_AC_XX\\tgnomADv3_nhomalt"}} {{print}}' {input} > {output}
        """
rule combine_tsv:
    input:
        os.path.join(work_dir, 'temp/{sample}_main_with_header.tsv'),
        os.path.join(work_dir, 'temp/{sample}_vep_with_header.tsv')
    output:
        os.path.join(work_dir, 'temp/{sample}_combined.tsv')
    shell:
        """
        paste {input} | sed 's/ /\t/g' > {output}
        """
