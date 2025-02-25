Downloaded data files and curated files can be found in this lambda folder:
/nfs/ml_lab/projects/Pilot1_PreclinicalHPC/priyanka/Data_Curation_final

CCLE_Multiomics_Data -
-Contains 10 files (Release 22Q2) downloaded from https://depmap.org/portal/download/
-Argonne_combined_rnaseq_data - File containing Argonne IDs of the cell line . Used to generate the mapping file
-CCLE_Multiomics_Data/README.txt - Information on the CCLE multiomics files
-CCLE_Multiomics_Data/sample_info.csv - Sample information on the CCLE multiomics files
Curated_CCLE_Multiomics_files/ - Contains the curated files

# Scripts

1) maps.py - Creates the following files:
        a) DepMap_Argonne_Mapping.csv - Mapping between DepMap ID, Argonne ID, and RRID (Cellosaurus ID)
        b) es_gs_orig.csv - Contains original gene symbols and EntrenzID from CCLE_expression.csv, CCLE_gene_cn.csv, and CCLE_wes_gene_cn.csv
        c) ens_gs_gef_orig.csv - Contains original gene symbols and Ensembl IDs from CCLE_expression_full.csv
2) eg_gs.R - Generate mapping files to map entrenz gene ids, ensembl ids and gene symbols. Creates the following files:
        a) eg_gs_map_mut.csv -  Maps Entrenz ID to gene symbol and ensembl ids (from the CCLE_mutation.csv data file)
        b) eg_gs_map_cn_ge.csv - Maps Entrenz ID to gene symbol and ensembl ids (from es_gs_orig.csv)
        c) gs_eg_gs_map_rrbs.csv - Maps gene symbol to Entrnz id and back to corrected gene symbol. Also maps entrez id to ensembl id (from CCLE_RRBS_TSS_1kb_20180614.txt)
        d) ens_eg_gs_map_gef.csv - Maps Ensembl ID to Entrenz ID to gene symbol (from ens_gs_gef_orig.csv)
3) mutation_curation.py - Code for curating mutation data. Generates curated mutation files -
        a) Mutation_AID_count.csv - File showing the number of mutations in genes for cell-lines. Column headers - Entrenz gene ids(Second row shows Gene symbols, third row - Ensembl ID). Rows - cell-lines
        b) Mutation_AID_binary.csv - Binary file showing the presence of specific genome changes in the cell-lines. Column headers - cell-lines. Rows - Genome-change. Several columns before 'Genome_Change' contains more information on the mutations
        c) Mutation_AID_long_format.csv - Mutation file in the long format
4) mutation_curation_with_isDel.py - Code for curating mutation data.Included isDeleterious and Annotation_Transcripts columns. Also does some processing on the columns to unify metadata for rows having different Annotation_Transcript for the same Genome_Change. Generates:
        a) Mutation_RRID_binary_isDel - .csv and .parquet
        b) Mutation_RRID_count_isDel - .csv 
        c) Mutation_RRID_long_format_isDel - .csv and .parquet
5) Data_curation.py - Code for curating all the other CCLE multiomics data. The files generated are:
        a) CCLE_AID_expression.csv - curated gene expression file. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines
        b) CCLE_AID_expression_scaled.csv - curated gene expression file after scaling the TPM values such that the TPM sum is 1 million
        c) CCLE_AID_expression_full.csv - Curated gene expression file from CCLE_expression_full.csv. Sum of TPM values is ~ 1 million. Features - Ensemble gene IDs. (First row contains Entrenz ID and second row contains gene symbols
        d) CCLE_AID_gene_cn.csv - Curated copy number data. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines
        e) CCLE_AID_gene_cn_discretized.csv - Discretized version of CCLE_AID_gene_cn.csv using:deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        f) CCLE_AID_wes_gene_cn.csv - Curated copy number wes data. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines
        g) CCLE_AID_wes_gene_cn_discretized.csv - Discretized version of CCLE_AID_wes_gene_cn.csv using:deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
        h) CCLE_AID_miRNA_20180525.csv - Curated miRNA data
        i) CCLE_AID_RPPA_20180123.csv - Curated Reverse Phase Protein Array (RPPA) data. Features - Antibiotic name
        j) CCLE_AID_RRBS_TSS_1kb_20180614.csv - Curated DNA Methylation - Reduced Representation Bisulfite Sequencing (RRBS) data. Features - TSS-ID of genes. Second row contains Entrenz gene id, third row contains gene symbols; and third row contains ensembl ids
6) samples.csv -  Mapping file from https://github.com/PNNL-CompBio/candleDataProcessing/tree/main/data

# Directories
1) CCLE_Multiomics_Data - Data from DepMapp 22Q2
        -Contains 10 files downloaded from https://depmap.org/portal/download/
        -Argonne_combined_rnaseq_data - File containing Argonne IDs of the cell line . Used to generate the mapping file
        -README.txt - Information on the CCLE multiomics files
        -sample_info.csv - Sample information on the CCLE multiomics files
2) Curated_CCLE_Multiomics_files - Contains the curated files
3) Maps - Contains the generated maps
        - DepMap_Argonne_Mapping.csv - Map DepMap ID with Argonne ID and RRID
        - eg_gs_map_mut.csv - Map entrenz id with gene symbol and ensembl id for mutation file
        - eg_gs_map_cn_ge.csv - Map entrenz id with gene symbol and ensembl id for copy number and gene expression files
        - es_gs_orig.csv - Contains the tabulated gene symbol and entrenz gene ids extracted from the gene expression and copy number files downloaded from DepMap. This file is input to the R script, eg_gs.R to generate eg_gs_map_cn_ge.csv
        - gs_eg_gs_map_rrbs.csv - Mapping file generated by eg_gs.R for the RRBS data. This file contains the orginal gene symbol, mapped entrenz id, corrected gene symbol, and mapped ensembl id.
        - ens_eg_gs_map_gef.csv Mapping file generated from ens_gs_gef_orig.csv. This file contains original gene symbols, original Ensembl IDs, mapped Entrenz ID, and mapped gene symbols.

