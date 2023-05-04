
Downloaded data files and curated files can be found in this lambda folder:<br>
/lambda_stor/homes/pvasanthakumari/Data_Curation_final<br>
1) CCLE_Multiomics_Data -<br>
        -Contains 10 files (Release 22Q2) downloaded from https://depmap.org/portal/download/<br>
        -Argonne_combined_rnaseq_data - File containing Argonne IDs of the cell line . Used to generate the mapping file<br>
        -README.txt - Information on the CCLE multiomics files<br>
        -sample_info.csv - Sample information on the CCLE multiomics files<br>
2) Curated_CCLE_Multiomics_files - Contained the curated files<br>


## Data curation scripts 
1) maps.py - Creates the following files:<br>
        a) DepMap_Argonne_Mapping.csv - Mapping between DepMap ID, Argonne ID, and RRID (Cellosaurus ID)<br>
        b) es_gs_orig.csv - Contains original gene symbols and EntrenzID from CCLE_expression.csv, CCLE_gene_cn.csv, and CCLE_wes_gene_cn.csv<br>
        c) ens_gs_gef_orig.csv - Contains original gene symbols and Ensembl IDs from CCLE_expression_full.csv<br>
2) eg_gs.R - Generate mapping files to map entrenz gene ids, ensembl ids and gene symbols. Creates the following files:<br>
        a) eg_gs_map_mut.csv -  Maps Entrenz ID to gene symbol and ensembl ids (from the CCLE_mutation.csv data file)<br>
        b) eg_gs_map_cn_ge.csv - Maps Entrenz ID to gene symbol and ensembl ids (from es_gs_orig.csv)<br>
        c) gs_eg_gs_map_rrbs.csv - Maps gene symbol to Entrnz id and back to corrected gene symbol. Also maps entrez id to ensembl id (from CCLE_RRBS_TSS_1kb_20180614.txt)<br>
        d) ens_eg_gs_map_gef.csv - Maps Ensembl ID to Entrenz ID to gene symbol (from ens_gs_gef_orig.csv)<br>
3) mutation_curation.py - Code for curating mutation data. Generates curated mutation files -<br>
        a) Mutation_AID_count.csv - File showing the number of mutations in genes for cell-lines. Column headers - Entrenz gene ids(Second row shows Gene symbols, third row - Ensembl ID). Rows - cell-lines<br>
        b) Mutation_AID_binary2.csv - Binary file showing the presence of specific genome changes in the cell-lines. Column headers - cell-lines. Rows - Genome-change. Several columns before 'Genome_Change' contains more information on the mutations<br>
        c) Mutation_AID_binary2_filter.csv - File b after filtering mutations that cause protein changes.<br>
        d) Mutation_AID_count_filter.csv - File a after filtering mutations that cause protein changes.<br>
4) Data_curation.py - Code for curating all the other CCLE multiomics data. The files generated are:<br>
        a) CCLE_AID_expression.csv - curated gene expression file. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines<br>
        b) CCLE_AID_expression_scaled.csv - curated gene expression file after scaling the TPM values such that the TPM sum is 1 million<br>
        c) CCLE_AID_expression_full.csv - Curated gene expression file from CCLE_expression_full.csv. Sum of TPM values is ~ 1 million. Features - Ensemble gene IDs. (First row contains Entrenz ID and second row contains gene symbols<br>
        d) CCLE_AID_gene_cn.csv - Curated copy number data. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines<br>
        e) CCLE_AID_gene_cn_discretized.csv - Binarized version of CCLE_AID_gene_cn.csv using:deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp <br>
        f) CCLE_AID_wes_gene_cn.csv - Curated copy number wes data. Features - Entrez gene IDs (Hugo gene symbols in the first row; Ensembl id in the second row). Rows - Cell-lines<br>
        g) CCLE_AID_wes_gene_cn_discretized.csv - Binarized version of CCLE_AID_wes_gene_cn.csv using:deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp <br>
        h) CCLE_AID_miRNA_20180525.csv - Curated miRNA data<br>
        i) CCLE_AID_RPPA_20180123.csv - Curated Reverse Phase Protein Array (RPPA) data. Features - Antibiotic name<br>
        j) CCLE_AID_RRBS_TSS_1kb_20180614.csv - Curated DNA Methylation - Reduced Representation Bisulfite Sequencing (RRBS) data. Features - TSS-ID of genes. Second row contains Entrenz gene id, third row contains gene symbols; and third row contains ensembl ids<br>
5) samples.csv -  Mapping file from https://github.com/PNNL-CompBio/candleDataProcessing/tree/main/data <br>

## Directories 
1) Maps - Contains the generated maps<br>
        - DepMap_Argonne_Mapping.csv - Map DepMap ID with Argonne ID and RRID<br>
        - eg_gs_map_mut.csv - Map entrenz id with gene symbol and ensembl id for mutation file<br>
        - eg_gs_map_cn_ge.csv - Map entrenz id with gene symbol and ensembl id for copy number and gene expression files<br>
        - es_gs_orig.csv - Contains the tabulated gene symbol and entrenz gene ids extracted from the gene expression and copy number files downloaded from DepMap. This file is input to the R script, eg_gs.R to generate eg_gs_map_cn_ge.csv<br>
        - gs_eg_gs_map_rrbs.csv - Mapping file generated by eg_gs.R for the RRBS data. This file contains the orginal gene symbol, mapped entrenz id, corrected gene symbol, and mapped ensembl id.<br>
        - ens_eg_gs_map_gef.csv Mapping file generated from ens_gs_gef_orig.csv. This file contains original gene symbols, original Ensembl IDs, mapped Entrenz ID, and mapped gene symbols.<br>
