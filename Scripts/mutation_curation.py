import pandas as pd
import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")

file_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) # Move to one level lower
path = file_path+'/CCLE_Multiomics_Data/CCLE_mutations.csv'
mut = pd.read_csv(path)
mut = mut[mut['Entrez_Gene_Id']!=0].reset_index(drop=True) # Remove rows with entrenz id = 0
mut = mut.drop_duplicates().reset_index(drop=True) # Remove duplicate rows
#Map entrez gene id to corrected gene symbol
eg2gs_map = pd.read_csv(str(file_path+'/Maps/eg_gs_map_mut.csv')).drop(columns={'Unnamed: 0'}) # Map generaed from R script
map = dict(zip(eg2gs_map.old_entrez_id, eg2gs_map.en2gs_GS )) # Mapping dictionary (entrez id to gene symbol)
ens_map = dict(zip(eg2gs_map.old_entrez_id, eg2gs_map.ens_id )) # Mapping dictionary (entrez id to ensembl id)
mut.insert(loc=1, column='Gene_symbol', value=['' for i in range(mut.shape[0])])
mut['Gene_symbol'] = mut['Entrez_Gene_Id'].map(map)
mut.insert(loc=3, column='Ensembl_id', value=['' for i in range(mut.shape[0])])
mut['Ensembl_id'] = mut['Entrez_Gene_Id'].map(ens_map)
mut = mut.drop(columns = 'Hugo_Symbol')
mut = mut[~pd.isna(mut['Gene_symbol'])].reset_index(drop=True) # Dropping rows where entrnz gene id cannot map to a gene symbol

# Check for cases where different genes have same mutations ('Genome_Change')
mut_dup = mut.drop_duplicates(subset = ['Entrez_Gene_Id','Genome_Change'], keep = False ).reset_index(drop=True)
mut_dup = mut_dup[mut_dup.duplicated(subset = 'Genome_Change', keep = False)]
mut_dup = mut_dup.sort_values(by = 'Genome_Change').reset_index(drop=True)
mut_dup = mut_dup[['Gene_symbol', 'Entrez_Gene_Id', 'Ensembl_id', 'Genome_Change']]

## Adding manual corrections to those genes that are mapped to same mutation (genome_change)
# keys - genes to be replaced; values - genes to be replaced with
gene_map = {'IPO4':'REC8', 
            'CFB':'C2',
            'PCDHA13':'PCDHA9',
            'ATP5MF-PTCD1':'PTCD1',
            'NR1D1':'THRA',
            'TRIM39-RPP21':'TRIM39',
            'BIVM-ERCC5':'ERCC5',
            'PPAN-P2RY11':'PPAN',
            'ANKMY1':'DUSP28',
            'POC1B-GALNT4':'GALNT4',
            'UGT1A8':'UGT1A1'}

for key in gene_map.keys():
    gc = mut_dup['Genome_Change'][np.where(mut_dup['Gene_symbol']==key)[0]].values[0]
    ind = mut.loc[(mut['Gene_symbol'] == key) & (mut['Genome_Change']==gc)].index[0]
    mut.loc[ind,'Gene_symbol'] = gene_map[key]
    mut.loc[ind,'Entrez_Gene_Id'] = mut_dup['Entrez_Gene_Id'][np.where(mut_dup['Gene_symbol']==gene_map[key])[0]].values[0]
    mut.loc[ind,'Ensembl_id'] = mut_dup['Ensembl_id'][np.where(mut_dup['Gene_symbol']==gene_map[key])[0]].values[0]

mut = mut.drop_duplicates().reset_index(drop=True) # Remove duplicate rows

# All Genes
en_gs_id = mut[['Entrez_Gene_Id', 'Gene_symbol', 'Ensembl_id']]
#en_gs_id = [en_gs_id.drop(i) for i in np.where(en_gs_id['Entrez_Gene_Id'] == 0)]
en_gs_id = en_gs_id.drop_duplicates() # Unique Rows
#All Cell-lines
cl = mut['DepMap_ID'].unique()

# Check if 1 cell line has more than 1 occurence for 'Genome_Change'. This assertion must be equal to 0
assert(len(mut[mut.duplicated(subset = ['DepMap_ID','Genome_Change'], keep=False)]) ==0)

def mutation_long(savename):
    mut_long = mut.copy()
    mut_long.insert(loc=3, column='Source', value='CCLE')
    mut_long = mut_long.drop(columns=['dbSNP_RS', 'dbSNP_Val_Status', 'Annotation_Transcript', 'cDNA_Change', 'Codon_Change'])
    mut_long = mut_long.iloc[:,:-13] # Drop last 13 columns
    col = mut_long.pop('DepMap_ID')
    mut_long.insert(loc =3, column = col.name, value= col)
    # Add RRID ID
    # Load DepMap mapping file
    # Load DepMap mapping file
    map = pd.read_csv(str(file_path+'/CCLE_Multiomics_Data/sample_info.csv')) # Mapping file downloaded from DepMap
    map_imp_id = pd.read_csv(str(file_path+ '/samples.csv')) # Mapping file from Sara Gosline's git repo
    map_imp_id = map_imp_id[map_imp_id['id_source'] == 'DepMap'].reset_index(drop=True) # filtering id source as DepMap
    ind_imp_id = np.where(pd.isna(map['RRID'])==True)[0] # Indices where RRID = NaN
    #Assign improve ids to map file with missing RRID
    for ind in ind_imp_id:
        ind2 = np.where(map_imp_id['other_id']==map['DepMap_ID'][ind])[0]
        if len(ind2)!=0:
            map['RRID'][ind] = map_imp_id['improve_sample_id'][ind2].values[0]
    map = map.dropna(subset = ['RRID']).reset_index(drop=True) # Remove rows where RRID/Improve_id == NaN
    # Some DepMap IDs in this mapping file is mapped to the same RRIDs. Manual changes to address this issue:
    cl_rm = ['ACH-001189','ACH-001024', 'ACH-002222', 'ACH-002260'] #cell-lines to be removed
    for cell in cl_rm:
        map = map[map.DepMap_ID != cell].reset_index(drop=True)
    
    mut_long.insert(loc=0, column='RRID', value=['' for i in range(mut_long.shape[0])])
    col = mut_long.columns
    rows = []
    for i in range(len(map)):
        id_all = np.where(mut_long['DepMap_ID']==map['DepMap_ID'][i])[0]
        if len(id_all)!=0:
            for id in id_all:
                row_append = mut_long.iloc[id]
                row_append['RRID'] = map['RRID'][i]
                rows.append(list(row_append.values))     
    mut_long_new = pd.DataFrame(rows, columns = col)
    mut_long_new = mut_long_new.drop(columns = ['DepMap_ID'])
    mut_long_new = mut_long_new.drop_duplicates().reset_index(drop=True) # Remove duplicate rows
    mut_long_new = mut_long_new.rename(columns= {'New_Hugo_Symbol': 'Gene_symbol', 'Entrez_Gene_Id': 'Entrez_id', })
    save_path = file_path+'/Curated_CCLE_Multiomics_files/'+savename+'.csv'
    mut_long_new.to_csv(str(save_path),index=False)


#Mutation count data frame - V1
def mutation_count(savename, protein_filter = False):
    mut_cp = mut.copy()
    mut_count = pd.DataFrame(0, index = range(len(cl)),columns = [en_gs_id['Entrez_Gene_Id'].values])
    mut_count.insert(loc=0, column='DepMap_ID', value=list(cl))
    if protein_filter:#Filter mutations with protein change
        mut_cp = mut_cp[~pd.isna(mut_cp.Protein_Change)].reset_index()
    for i in range(len(mut_count)):
        gene_id = np.where(mut_cp['DepMap_ID'] == mut_count['DepMap_ID'].iloc[i]['DepMap_ID'])
        gene = mut_cp['Entrez_Gene_Id'][gene_id[0]].reset_index(drop = True)
        for j in range(len(gene)):
            if gene[j] in mut_count.columns:
                mut_count.loc[i, gene[j]] = mut_count.loc[i, gene[j]] + 1

    # Add RRID ID
    # Load DepMap mapping file
    map = pd.read_csv(str(file_path+'/CCLE_Multiomics_Data/sample_info.csv')) # Mapping file downloaded from DepMap
    map_imp_id = pd.read_csv(str(file_path+ '/samples.csv')) # Mapping file from Sara Gosline's git repo
    map_imp_id = map_imp_id[map_imp_id['id_source'] == 'DepMap'].reset_index(drop=True) # filtering id source as DepMap
    ind_imp_id = np.where(pd.isna(map['RRID'])==True)[0] # Indices where RRID = NaN
    #Assign improve ids to map file with missing RRID
    for ind in ind_imp_id:
        ind2 = np.where(map_imp_id['other_id']==map['DepMap_ID'][ind])[0]
        if len(ind2)!=0:
            map['RRID'][ind] = map_imp_id['improve_sample_id'][ind2].values[0]
    map = map.dropna(subset = ['RRID']).reset_index(drop=True) # Remove rows where RRID/Improve_id == NaN
    # Some DepMap IDs in this mapping file is mapped to the same RRIDs. Manual changes to address this issue:
    cl_rm = ['ACH-001189','ACH-001024', 'ACH-002222', 'ACH-002260'] #cell-lines to be removed
    for cell in cl_rm:
        map = map[map.DepMap_ID != cell].reset_index(drop=True)
    ##
    mut_count.insert(loc=0, column='RRID', value=['' for i in range(mut_count.shape[0])])
    new_row1 = pd.DataFrame([['', '']+list(en_gs_id['Gene_symbol'].values)],columns = mut_count.columns )
    new_row2 = pd.DataFrame([['', '']+list(en_gs_id['Ensembl_id'].values)],columns = mut_count.columns )
    mut_count = pd.concat([new_row1, new_row2, mut_count]).reset_index(drop=True)
    col = mut_count.columns
    rows = [mut_count.iloc[0].values]
    rows.append(list(mut_count.iloc[1].values))
    for i in range(len(map)):
        id = np.where(mut_count['DepMap_ID']==map['DepMap_ID'][i])
        if len(id[0])!=0:
            row_append = mut_count.iloc[id[0][0]]
            row_append['RRID'] = map['RRID'][i]
            rows.append(list(row_append.values))     
    mut_count_new = pd.DataFrame(rows, columns = col)
    mut_count_new = mut_count_new.drop(columns = ['DepMap_ID'])
    save_path = file_path+'/Curated_CCLE_Multiomics_files/'+savename+'.csv'
    mut_count_new = mut_count_new.drop_duplicates().reset_index(drop=True) # Remove duplicate rows
    mut_count_new.to_csv(str(save_path),index=False)


#Mutation  - Binary file
def mutation_binary2(savename, protein_filter = False):
    mut_cp = mut.copy()
    if protein_filter:#Filter mutations with protein change
        mut_cp = mut_cp[~pd.isna(mut_cp.Protein_Change)].reset_index(drop=True)
    mut_v2 = mut_cp.iloc[:,range(15)].copy() # Preserve the first 15 columns
    mut_v2['Protein_Change'] = mut_cp['Protein_Change'] # Add protein change column - users can decide whether to filter basedon protein change or not
    mut_v2['DepMap_ID'] = mut_cp['DepMap_ID'] # Add DepMap column
    mut_v2 = mut_v2.drop(columns = ['dbSNP_RS', 'dbSNP_Val_Status'])
    mut_v2 = mut_v2.astype({'Chromosome': 'str'})
    mut_nodup = mut_v2.iloc[:,range(14)].drop_duplicates() # Drop duplicates 
    mut_nodup_gc = mut_nodup.drop_duplicates(subset = ['Genome_Change']) # Get unique genome changes
    mut_nodup_gc=mut_nodup_gc.reset_index(drop = True) 

    # Generate binary mutation file
    for col in cl:
        mut_nodup_gc[col] = 0
    for i in range(len(cl)):
        id = np.where(mut_v2['DepMap_ID'] == cl[i])[0]
        gc = mut_v2['Genome_Change'][id]
        id2 = mut_nodup_gc[mut_nodup_gc['Genome_Change'].isin(list(gc))].index
        mut_nodup_gc[cl[i]][id2] = 1

    # Map DepMap ID to RRID
    # To preserve depmap cell line ids
    # head = mut_nodup_gc.columns
    # new_row = pd.DataFrame([head],columns = mut_nodup_gc.columns )
    # new_row.iloc[:,:14] = ''
    # mut_nodup_gc = pd.concat([new_row, mut_nodup_gc]).reset_index(drop = True)
    # Load DepMap mapping file
    map = pd.read_csv(str(file_path+'/CCLE_Multiomics_Data/sample_info.csv')) # Mapping file downloaded from DepMap
    map_imp_id = pd.read_csv(str(file_path+ '/samples.csv')) # Mapping file from Sara Gosline's git repo
    map_imp_id = map_imp_id[map_imp_id['id_source'] == 'DepMap'].reset_index(drop=True) # filtering id source as DepMap
    ind_imp_id = np.where(pd.isna(map['RRID'])==True)[0] # Indices where RRID = NaN
    #Assign improve ids to map file with missing RRID
    for ind in ind_imp_id:
        ind2 = np.where(map_imp_id['other_id']==map['DepMap_ID'][ind])[0]
        if len(ind2)!=0:
            map['RRID'][ind] = map_imp_id['improve_sample_id'][ind2].values[0]
    map = map.dropna(subset = ['RRID']).reset_index(drop=True) # Remove rows where RRID/Improve_id == NaN
    # Some DepMap IDs in this mapping file is mapped to the same RRIDs. Manual changes to address this issue:
    cl_rm = ['ACH-001189','ACH-001024', 'ACH-002222', 'ACH-002260'] #cell-lines to be removed
    for cell in cl_rm:
        map = map[map.DepMap_ID != cell].reset_index(drop=True)

    mutbin_file = mut_nodup_gc.iloc[:,range(14)].copy()
    for i in range(len(map)):
        if (map['DepMap_ID'][i] in mut_nodup_gc.columns):
            id = mut_nodup_gc.columns.get_loc(map['DepMap_ID'][i])
            mutbin_file[map['RRID'][i]] = mut_nodup_gc.iloc[:,id]
    save_path = file_path+'/Curated_CCLE_Multiomics_files/'+savename+'.csv'
    mutbin_file = mutbin_file.drop_duplicates().reset_index(drop=True) # Remove duplicate rows
    # Rename some column names
    mutbin_file = mutbin_file.rename(columns= {'New_Hugo_Symbol': 'Gene_symbol', 'Entrez_Gene_Id': 'Entrez_id', })
    #mutbin_file.to_csv(save_path)
    save_path = file_path+'/Curated_CCLE_Multiomics_files/'+savename+'.parquet'
    mutbin_file.columns = mutbin_file.columns.astype(str)
    mutbin_file.to_parquet(save_path)

# Function calls
if __name__ == "__main__":
    #Unfiltered
    print('Processing Mutation_AID_binary')
    mutation_binary2('Mutation_AID_binary', False)
    #print('Processing Mutation_AID_count')
    #mutation_count('Mutation_AID_count', False)
    #print('Processing Mutation_AID_long_format')
    #mutation_long('Mutation_AID_long_format')

    #Filtering based on protein change
    #print('Processing Mutation_AID_binary2_filter')
    #mutation_binary2('Mutation_AID_binary2_filter', True)
    #print('Processing Mutation_AID_count_filter')
    #mutation_count('Mutation_AID_count_filter', True)

