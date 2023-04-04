import pandas as pd
import os
import numpy as np

########### Map DepMapID with Argonne ID #################

file_path = os.path.dirname(os.path.realpath(__file__))
path = str(file_path + '/CCLE_Multiomics_Data/Argonne_combined_rnaseq_data')
anl_file = pd.read_table(path)
cell_lines = []
cell_lines_strp = []
map = {}
map_strp = {}
cnt = 0

#All the cell lines
for ind in range(len(anl_file['Sample'])):
    if(anl_file['Sample'][ind].split('.')[0] == 'GDC' or anl_file['Sample'][ind].split('.')[0] =='NCIPDM'):
        continue
    cl = anl_file['Sample'][ind].split('.')[1]
    strp_cl = (''.join(e for e in cl if e.isalnum())).upper()
    cell_lines.append(cl)
    cell_lines_strp.append(strp_cl)

# Dictionary that maps stripped cell line names to Argonne IDs 
    if strp_cl not in map_strp.keys():
        map_strp[strp_cl] = [anl_file['Sample'][ind]]
        cnt+=1
        continue
    if (anl_file['Sample'][ind] not in map_strp[strp_cl]):
        map_strp[strp_cl].append(anl_file['Sample'][ind])
        cnt+=1

# Dictionary that maps cell line names to Argonne IDs 
    if cl not in map.keys():
        map[cl] = [anl_file['Sample'][ind]]
        
        continue
    if (anl_file['Sample'][ind] not in map[cl]):
        map[cl].append(anl_file['Sample'][ind])

cl_unique = set(cell_lines)
cl_strp_unique = set(cell_lines_strp)

# Argonne ID's are a mix of stripped and unstripped cell-line names. So we need to check that first
# CCLE - stripped names
check = {}
cl_unique = list(cl_unique)
strp_cl_un = []
count = 0
for i in range(len(cl_unique)):
    strp_cl_un.append((''.join(e for e in cl_unique[i] if e.isalnum())).upper())
    if strp_cl_un[i] not in check.keys():
        check[strp_cl_un[i]] = [cl_unique[i]]
        count+=1
    else:
        check[strp_cl_un[i]].append(cl_unique[i])
        count+=1
df =  pd.DataFrame(columns = ['cell_line_strip', 'cell_line', 'Argonne_ID'])
df['cell_line_strip'] = check.keys()
for i in range(len(df)):
    df['cell_line'][i] = check[df['cell_line_strip'][i]]
    df['Argonne_ID'][i] = map_strp[df['cell_line_strip'][i]]


#sample info - CCLE Multiomics data
path = str(file_path +'/CCLE_Multiomics_Data/sample_info.csv')
sample_info = pd.read_csv(path)
#Unstripped - check for common cell-line names with Argonne data
sample_cl = sample_info.iloc[:,1]
sample_cl_un = set(sample_cl)
common = set(cl_unique).intersection(sample_cl_un)
#Stripped - check for common cell-line names with Argonne data
sample_cl_strp = sample_info.iloc[:,2]
sample_cl_strp_un = set(sample_cl_strp)
common_strp = cl_strp_unique.intersection(sample_cl_strp_un)

#mapping file
map_file = pd.DataFrame(columns = ['stripped_cell_line_name', 'cell_line_name', 'DepMap_ID', 'Argonne_ID', 'CCLE_Name', 'RRID','Manual_edit'])
map_file['stripped_cell_line_name'] = sample_info['stripped_cell_line_name']
map_file['cell_line_name'] = sample_info['cell_line_name']
map_file['DepMap_ID'] = sample_info['DepMap_ID']
map_file['CCLE_Name'] = sample_info['CCLE_Name']
map_file['RRID'] = sample_info['RRID']
left_out = []
# Add Argonne ID
for i in range(len(map_file)):
    id = np.where(df['cell_line_strip'] == map_file['stripped_cell_line_name'][i])
    if len(id[0]) !=0:
        map_file['Argonne_ID'][i] = df['Argonne_ID'][id[0][0]]
    if (not pd.isna(sample_info['alias'][i])) and (sample_info['alias'][i] != map_file['stripped_cell_line_name'][i]):
        id1 = np.where(df['cell_line_strip'] == sample_info['alias'][i])
        if len(id1[0]) !=0:
            a_id = map_file['Argonne_ID'][i]
            if type(a_id) is not list:
                a_id = [a_id]
            if pd.isna(a_id[0]):
                map_file['Argonne_ID'][i] =df['Argonne_ID'][id1[0][0]]
            else:
                for k in range(len(df['Argonne_ID'][id1[0][0]])):
                    map_file['Argonne_ID'][i].append(df['Argonne_ID'][id1[0][0]][k])
        left_out.append(map_file['stripped_cell_line_name'][i]) # cell lines not in Argonne ID

## Manual additions
id = np.where(map_file['stripped_cell_line_name'] == '786O')
map_file['Argonne_ID'][id[0][0]].append('NCI60.786-0')
map_file['Argonne_ID'][id[0][0]].append('GDSC.786-0')
map_file['Manual_edit'][id[0][0]] = 'Yes'

id = np.where(map_file['stripped_cell_line_name'] == 'U251MG')
map_file['Argonne_ID'][id[0][0]].append('NCI60.U251')
map_file['Argonne_ID'][id[0][0]].append('GDSC.U251')
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'INA6')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.INA-6-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'SKNO1')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.SKNO-1-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'PLB985')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.PLB-985-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'HNT34')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.HNT-34-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'OCIMY5')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.OCI-MY5-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'KMS18')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.KMS-18-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'OPM1')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.OPM-1-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'


id = np.where(map_file['stripped_cell_line_name'] == 'OCIMY7')
map_file['Argonne_ID'][id[0][0]] = ['CCLE.OCI-MY7-2']
map_file['Manual_edit'][id[0][0]] = 'Yes'

# Remove DepMap IDs corresponding to cell-line names GR-M and CHL-1. Two DepMap IDs corresponded to the same cell-lines
id = np.where(map_file['stripped_cell_line_name']=='GRM')[0]
map_file = map_file.drop(id).reset_index(drop=True)

id = np.where(map_file['stripped_cell_line_name']=='CHL1')[0]
map_file = map_file.drop(id).reset_index(drop=True)

map_file2 = pd.DataFrame(columns = ['stripped_cell_line_name', 'cell_line_name', 'DepMap_ID', 'Argonne_ID','CCLE_Name', 'RRID','Manual_edit'])
count = 0
for i in range(len(map_file)):
    a_id = map_file['Argonne_ID'][i]
    if type(a_id) is not list:
        a_id = [a_id]
    if pd.isna(a_id[0]):
        continue
    for j in range(len(a_id)):
        map_file2 = map_file2.append(map_file.iloc[i,:], ignore_index=True)
        map_file2.reset_index(drop=True, inplace=True)
        map_file2['Argonne_ID'][count] = a_id[j]
        count+=1

# Reorder columns
map_file2 = map_file2[['Argonne_ID','DepMap_ID','RRID','cell_line_name','stripped_cell_line_name','CCLE_Name', 'Manual_edit']]
#Save
map_file2.to_csv(str(file_path+'/Maps/DepMap_Argonne_Mapping.csv'))

############ Map Entrenz ID to gene symbol ###############
files = ['CCLE_wes_gene_cn.csv', 'CCLE_gene_cn.csv', 'CCLE_expression.csv']
genes = []
for file in files:
    path = str(file_path +'/CCLE_Multiomics_Data/'+file)
    read_file = pd.read_csv(path)
    genes.extend(list(read_file.columns[1:]))
genes_unique = list(set(genes))
gs = []
es = []
count = 0
for gene in genes_unique:
    if '(' in gene:
        es.append(gene.split('(', 1)[1].split(')')[0])
        gs.append(gene.split(' ', 1)[0])
    else:
        count = count+1
        es.append(gene)
        gs.append(np.nan)
es_gs = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es})
es_gs2 = es_gs[es_gs.entrenz_id!='nan'].reset_index(drop=True)
es_gs2.to_csv(file_path+'/Maps/es_gs_orig.csv')
#This file is used to call the R function eg2gs. See R script for map generation

#Gene_expression_full file - This file is downloaded from DepMap. It has sum of TPM = 1 million
#Genes in this file are represeneted by their gene symbols and the corresponding Ensembl ID. 
# So we need to map ensembl ID to gene symbols
file = 'CCLE_expression_full.csv'
path = str(file_path +'/CCLE_Multiomics_Data/'+file)
read_file = pd.read_csv(path)
genes = list(read_file.columns[1:])
genes_unique_gef = list(set(genes))
gs = []
es = []
count = 0
for gene in genes_unique_gef:
    if '(' in gene:
        es.append(gene.split('(', 1)[1].split(')')[0])
        gs.append(gene.split(' ', 1)[0])
    else:
        count = count+1
        es.append(gene)
        gs.append(np.nan)
ens_gs = pd.DataFrame({'gene_symbol':gs, 'ensembl_id':es})
ens_gs2 = ens_gs[ens_gs.ensembl_id!='nan'].reset_index(drop=True)
ens_gs2.to_csv(file_path+'/Maps/ens_gs_gef_orig.csv')