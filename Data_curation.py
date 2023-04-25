import pandas as pd
import os
import numpy as np
import warnings
import math
warnings.filterwarnings("ignore")
import logging
import sys
# setup logging
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger('Processing CCLE Multiomics data')
file_path = os.path.dirname(os.path.realpath(__file__))
# Load Argonne-DepMap mapping file
map = pd.read_csv(str(file_path+'/Maps/DepMap_Argonne_Mapping.csv'))
# Load entrz id-gene symbol mapping file for copy number and gene expression files
path = file_path+'/Maps/eg_gs_map_cn_ge.csv'
eg_gs = pd.read_csv(path).drop(columns={'Unnamed: 0'}) 
eg_gs = eg_gs.astype({'old_entrez_id': str})
eg_gs_map = dict(zip(eg_gs.old_entrez_id, eg_gs.en2gs_GS )) # Mapping dictionary
eg_ens_map = dict(zip(eg_gs.old_entrez_id, eg_gs.ens_id )) # Mapping dictionary


#########Copy number data - WES#############
logger.info('Processing Copy number WES data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_wes_gene_cn.csv')
cn_file = pd.read_csv(path)
cn_file=cn_file.rename(columns={'Unnamed: 0':'DepMap_ID'})
# Map Entrenz gene id to new gene symbols and replace old gene symbols
genes = list(cn_file.columns[1:])
gs = []
es = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    es.append(gene.split('(', 1)[1].split(')')[0])
    gs.append(gene.split(' ', 1)[0])
    col_ind.append(cn_file.columns.get_loc(gene))
es_gs = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es, 'col_ind':col_ind})
ind = np.where(es_gs['entrenz_id']=='nan')[0]
cn_file = cn_file.drop(cn_file.columns[ind+1], axis=1) # Drop columns where entrenz_id = nan
genes = list(cn_file.columns[1:])
gs = []
es = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    es.append(gene.split('(', 1)[1].split(')')[0])
    gs.append(gene.split(' ', 1)[0])
    col_ind.append(cn_file.columns.get_loc(gene))
es_gs2 = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es, 'col_ind':col_ind})
#es_gs2 = es_gs[es_gs.entrenz_id!='nan'].reset_index(drop=True)

# Several Entrenz gene ids are repeating. So we need to remove columns corresponding to wrong entrenz id
dup_ind = np.where(es_gs2.duplicated(subset = ['entrenz_id'], keep=False)==True)[0]
es_gs_dup = es_gs2.iloc[dup_ind].reset_index(drop=True)
es_gs_dup['gs_map'] = es_gs_dup['entrenz_id'].map(eg_gs_map)
ind_drop = np.where((es_gs_dup['gene_symbol'] == es_gs_dup['gs_map'])==False)[0]
col_drop = es_gs_dup['col_ind'][ind_drop]
cn_file = cn_file.drop(cn_file.columns[col_drop], axis=1)
es_gs2 = es_gs2.drop(np.where(es_gs2['col_ind'].isin(col_drop)==True)[0]).reset_index(drop=True)

#mapping eg to gs and eg to ens_id
new_gs = es_gs2['entrenz_id'].map(eg_gs_map)
new_ens_id = es_gs2['entrenz_id'].map(eg_ens_map)
es_gs2['New_gene_symbol'] = new_gs
es_gs2['New_ensembl_id'] = new_ens_id
ind = np.where(pd.isna(new_gs)==True)[0] #index to drop columns
cn_file = cn_file.drop(cn_file.columns[ind+1], axis=1) # Drop columns where entrenz_id cannot be mapped to a gene symbol
es_gs2 = es_gs2.dropna(subset=['New_gene_symbol'])
# Add the new gene symbols to the data file
cn_file.columns = ['DepMap_ID']+list(es_gs2['entrenz_id'])
new_row1 = pd.DataFrame([['']+list(es_gs2['New_gene_symbol'])],columns = cn_file.columns )
new_row2 = pd.DataFrame([['']+list(es_gs2['New_ensembl_id'])],columns = cn_file.columns )
cn_file = pd.concat([new_row1, new_row2, cn_file]).reset_index(drop = True)
#Map DepMap ID with Argonne ID
cn_file.insert(loc=0, column='RRID', value=['' for i in range(cn_file.shape[0])])
col = cn_file.columns
rows = [cn_file.iloc[0].values]
rows.append(list(cn_file.iloc[1].values))
for i in range(len(map)):
    id = np.where(cn_file['DepMap_ID']==map['DepMap_ID'][i])
    if len(id[0])!=0:
        row_append = cn_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
cn_file_new = pd.DataFrame(rows, columns = col)
cn_file_new = cn_file_new.drop(columns = ['DepMap_ID'])
cn_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_wes_gene_cn.csv'),index=False)
# Binary copy number
# deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
cn_file_new_b = cn_file_new.drop(columns = ['RRID'])
cn_file_new_b = cn_file_new_b.drop([0,1])
cn_file_new_b[cn_file_new_b <= 0.5210507] = -2
cn_file_new_b[(cn_file_new_b > 0.5210507) &  (cn_file_new_b <= 0.7311832)] = -1
cn_file_new_b[(cn_file_new_b > 0.7311832) &  (cn_file_new_b <= 1.214125)] = 0
cn_file_new_b[(cn_file_new_b > 1.214125) &  (cn_file_new_b <= 1.422233)] = 1
cn_file_new_b[cn_file_new_b > 1.422233] = 2
cn_file_new_b.insert(column = 'RRID', loc=0, value = cn_file_new['RRID'][2:])
cn_file_new_b = pd.concat([cn_file_new[0:2], cn_file_new_b])
cn_file_new_b.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_wes_gene_cn_binary.csv'),index=False)


#########Copy number data - gene#############
logger.info('Processing Copy number gene data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_gene_cn.csv')
cn_file = pd.read_csv(path)
cn_file=cn_file.rename(columns={'Unnamed: 0':'DepMap_ID'})
# Map Entrenz gene id to new gene symbols and replace old gene symbols
genes = list(cn_file.columns[1:])
gs = []
es = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    es.append(gene.split('(', 1)[1].split(')')[0])
    gs.append(gene.split(' ', 1)[0])
    col_ind.append(cn_file.columns.get_loc(gene))
es_gs = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es, 'col_ind':col_ind})
ind = np.where(es_gs['entrenz_id']=='nan')[0]
cn_file = cn_file.drop(cn_file.columns[ind+1], axis=1) # Drop columns where entrenz_id = nan
genes = list(cn_file.columns[1:])
gs = []
es = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    es.append(gene.split('(', 1)[1].split(')')[0])
    gs.append(gene.split(' ', 1)[0])
    col_ind.append(cn_file.columns.get_loc(gene))
es_gs2 = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es, 'col_ind':col_ind})
#es_gs2 = es_gs[es_gs.entrenz_id!='nan'].reset_index(drop=True)

# Several Entrenz gene ids are repeating. So we need to remove columns corresponding to wrong entrenz id
dup_ind = np.where(es_gs2.duplicated(subset = ['entrenz_id'], keep=False)==True)[0]
es_gs_dup = es_gs2.iloc[dup_ind].reset_index(drop=True)
es_gs_dup['gs_map'] = es_gs_dup['entrenz_id'].map(eg_gs_map)
ind_drop = np.where((es_gs_dup['gene_symbol'] == es_gs_dup['gs_map'])==False)[0]
col_drop = es_gs_dup['col_ind'][ind_drop]
cn_file = cn_file.drop(cn_file.columns[col_drop], axis=1)
es_gs2 = es_gs2.drop(np.where(es_gs2['col_ind'].isin(col_drop)==True)[0]).reset_index(drop=True)

#mapping eg to gs
new_gs = es_gs2['entrenz_id'].map(eg_gs_map)
new_ens_id = es_gs2['entrenz_id'].map(eg_ens_map)
es_gs2['New_gene_symbol'] = new_gs
es_gs2['New_ensembl_id'] = new_ens_id
ind = np.where(pd.isna(new_gs)==True)[0] #index to drop columns
cn_file = cn_file.drop(cn_file.columns[ind+1], axis=1) # Drop columns where entrenz_id cannot be mapped to a gene symbol
es_gs2 = es_gs2.dropna(subset=['New_gene_symbol'])
# Add the new gene symbols to the data file
cn_file.columns = ['DepMap_ID']+list(es_gs2['entrenz_id'])
new_row1 = pd.DataFrame([['']+list(es_gs2['New_gene_symbol'])],columns = cn_file.columns )
new_row2 = pd.DataFrame([['']+list(es_gs2['New_ensembl_id'])],columns = cn_file.columns )
cn_file = pd.concat([new_row1, new_row2, cn_file]).reset_index(drop = True)
#Map DepMap ID with Argonne ID
cn_file.insert(loc=0, column='RRID', value=['' for i in range(cn_file.shape[0])])
col = cn_file.columns
rows = [cn_file.iloc[0].values]
rows.append(list(cn_file.iloc[1].values))
for i in range(len(map)):
    id = np.where(cn_file['DepMap_ID']==map['DepMap_ID'][i])
    if len(id[0])!=0:
        row_append = cn_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
cn_file_new = pd.DataFrame(rows, columns = col)
cn_file_new = cn_file_new.drop(columns = ['DepMap_ID'])
cn_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_gene_cn.csv'),index=False)
# Binary copy number
# deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
cn_file_new_b = cn_file_new.drop(columns = ['RRID'])
cn_file_new_b = cn_file_new_b.drop([0,1])
cn_file_new_b[cn_file_new_b <= 0.5210507] = -2
cn_file_new_b[(cn_file_new_b > 0.5210507) &  (cn_file_new_b <= 0.7311832)] = -1
cn_file_new_b[(cn_file_new_b > 0.7311832) &  (cn_file_new_b <= 1.214125)] = 0
cn_file_new_b[(cn_file_new_b > 1.214125) &  (cn_file_new_b <= 1.422233)] = 1
cn_file_new_b[cn_file_new_b > 1.422233] = 2
cn_file_new_b.insert(column = 'RRID', loc=0, value = cn_file_new['RRID'][2:])
cn_file_new_b = pd.concat([cn_file_new[0:2], cn_file_new_b])
cn_file_new_b.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_gene_cn_binary.csv'),index=False)


######### Gene expression #############
logger.info('Processing gene expression data V1')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_expression.csv')
ge_file = pd.read_csv(path)
ge_file=ge_file.rename(columns={'Unnamed: 0':'DepMap_ID'})
# Map Entrenz gene id to new gene symbols and replace old gene symbols
genes = list(ge_file.columns[1:])
gs = []
es = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    es.append(gene.split('(', 1)[1].split(')')[0])
    gs.append(gene.split(' ', 1)[0])
    col_ind.append(ge_file.columns.get_loc(gene))
es_gs = pd.DataFrame({'gene_symbol':gs, 'entrenz_id':es, 'col_ind':col_ind})
ind = np.where(es_gs['entrenz_id']=='nan')[0]
ge_file = ge_file.drop(ge_file.columns[ind+1], axis=1) # Drop columns where entrenz_id = nan
es_gs2 = es_gs[es_gs.entrenz_id!='nan'].reset_index(drop=True)

# Several Entrenz gene ids are repeating. So we need to remove columns corresponding to wrong entrenz id
dup_ind = np.where(es_gs2.duplicated(subset = ['entrenz_id'], keep=False)==True)[0]
es_gs_dup = es_gs2.iloc[dup_ind].reset_index(drop=True)
es_gs_dup['gs_map'] = es_gs_dup['entrenz_id'].map(eg_gs_map)
ind_drop = np.where((es_gs_dup['gene_symbol'] == es_gs_dup['gs_map'])==False)[0]
col_drop = es_gs_dup['col_ind'][ind_drop]
ge_file = ge_file.drop(ge_file.columns[col_drop], axis=1)
es_gs2 = es_gs2.drop(np.where(es_gs2['col_ind'].isin(col_drop)==True)[0]).reset_index(drop=True)

#mapping eg to gs
new_gs = es_gs2['entrenz_id'].map(eg_gs_map)
new_ens_id = es_gs2['entrenz_id'].map(eg_ens_map)
es_gs2['New_gene_symbol'] = new_gs
es_gs2['New_ensembl_id'] = new_ens_id
ind = np.where(pd.isna(new_gs)==True)[0] #index to drop columns
ge_file = ge_file.drop(ge_file.columns[ind+1], axis=1) # Drop columns where entrenz_id cannot be mapped to a gene symbol
es_gs2 = es_gs2.dropna(subset=['New_gene_symbol'])
# Add the new gene symbols to the data file
ge_file.columns = ['DepMap_ID']+list(es_gs2['entrenz_id'])
new_row1 = pd.DataFrame([['']+list(es_gs2['New_gene_symbol'])],columns = ge_file.columns )
new_row2 = pd.DataFrame([['']+list(es_gs2['New_ensembl_id'])],columns = ge_file.columns )
ge_file = pd.concat([new_row1, new_row2, ge_file]).reset_index(drop = True)
#Map DepMap ID with Argonne ID
ge_file.insert(loc=0, column='RRID', value=['' for i in range(ge_file.shape[0])])
col = ge_file.columns
rows = [ge_file.iloc[0].values]
rows.append(list(ge_file.iloc[1].values))
for i in range(len(map)):
    id = np.where(ge_file['DepMap_ID']==map['DepMap_ID'][i])
    if len(id[0])!=0:
        row_append = ge_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
ge_file_new = pd.DataFrame(rows, columns = col)
ge_file_new = ge_file_new.drop(columns = ['DepMap_ID'])
ge_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_expression.csv'),index=False)

############ Gene expression - scaled rows ############
#Since several cell lines have total TPM values significantly less than 1 million, 
#we need to generate another version of gene expression data such that each row/cell-line is scaled by a factor
#that makes the total TPM for the cell-line to be 1 million.
#Check Gene TPM
logger.info('Processing gene expression scaled data')
ge_copy = ge_file_new.copy()
ge_copy.drop(columns = ['RRID'], inplace = True)
ge_copy = ge_copy.drop([0,1]).reset_index(drop=True)
tpm = ge_copy.apply(lambda a:(2**a)-1, axis = 1 )
tpm_sum = tpm.sum(axis = 1)
factor = 1000000/tpm_sum
tpm_factor = tpm.mul(factor, axis=0)
new_ge = tpm_factor.apply(lambda a:np.log2(a+1), axis = 1 )
new_tpm = new_ge.apply(lambda a:(2**a)-1, axis = 1 )
new_tpm_sum = new_tpm.sum(axis = 1)
new_ge.insert(loc = 0, column = 'RRID', value =ge_file_new['RRID'][1:].reset_index(drop=True))
new_row1 = pd.DataFrame([ge_file_new.iloc[0].values],columns = new_ge.columns )
new_row2 = pd.DataFrame([ge_file_new.iloc[1].values],columns = new_ge.columns )
new_ge = pd.concat([new_row1,new_row2, new_ge]).reset_index(drop = True)
new_ge.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_expression_scaled.csv'),index=False)

######### Gene Expression full ############# This gene expression file has the sum of TPM ~1 million
logger.info('Processing gene expression full data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_expression_full.csv')
gef_file = pd.read_csv(path)
gef_file=gef_file.rename(columns={'Unnamed: 0':'DepMap_ID'})
# Map Ensembl gene id to new gene symbols and replace old gene symbols
genes = list(gef_file.columns[1:])
gs = []
ens = []
col_ind = []
for gene in genes: #Separating column header into gene symbols and entrenz ids
    if '(' in gene:
        ens.append(gene.split('(', 1)[1].split(')')[0])
        gs.append(gene.split(' ', 1)[0])
    else:
        ens.append(gene)
        gs.append(np.nan)
    col_ind.append(gef_file.columns.get_loc(gene))
es_gs = pd.DataFrame({'gene_symbol':gs, 'ensembl_id':ens, 'col_ind':col_ind})
ind = np.where(es_gs['ensembl_id']=='nan')[0]
gef_file = gef_file.drop(gef_file.columns[ind+1], axis=1) # Drop columns where entrenz_id = nan
es_gs2 = es_gs[es_gs.ensembl_id!='nan'].reset_index(drop=True)
# Check if Ensembl Ids are repeating
dup_ind = np.where(es_gs2.duplicated(subset = ['ensembl_id'], keep=False)==True)[0] #len(dup_ind) = 0, so there are no duplicate ensembl ids
#mapping ens to gs
# Load Ensembl_ID-Entrz id-Gene symbol mapping file - Gene expression full file
path = file_path+'/Maps/ens_eg_gs_map_gef.csv'
ens_eg_gs = pd.read_csv(path).drop(columns={'Unnamed: 0'}) 
col_drop = ens_eg_gs['orig_ensembl_id'][np.where(pd.isna(ens_eg_gs.entrenzID)==True)[0]] 
ens_eg_gs = ens_eg_gs.dropna(subset = ['entrenzID']).reset_index(drop=True)
ens_eg_gs = ens_eg_gs.astype({'entrenzID': int})
ens_eg_gs = ens_eg_gs.astype({'entrenzID': str})
ens_gs_map = dict(zip(ens_eg_gs.orig_ensembl_id, ens_eg_gs.new_GS )) # Mapping dictionary EnsemblID to new gene symbol
ens_eg_map = dict(zip(ens_eg_gs.orig_ensembl_id, ens_eg_gs.entrenzID )) # Mapping dictionary Ensembl ID to Entrenz ID
new_eg = es_gs2['ensembl_id'].map(ens_eg_map)
new_gs = es_gs2['ensembl_id'].map(ens_gs_map)
es_gs2['New_gene_symbol'] = new_gs
es_gs2['New_entrenzID'] = new_eg
ind = np.where(pd.isna(new_gs)==True)[0] #index to drop columns
gef_file = gef_file.drop(gef_file.columns[ind+1], axis=1) # Drop columns where entrenz_id cannot be mapped to a gene symbol
es_gs2 = es_gs2.dropna(subset = ['New_gene_symbol']).reset_index(drop=True)
### Check for duplicate gene symbols in the mapping file
dup_ind = np.where(es_gs2.duplicated(subset = ['New_gene_symbol'], keep=False)==True)[0]
es_gs_dup = es_gs2.iloc[dup_ind].reset_index(drop=True)
es_gs_dup = es_gs_dup.sort_values(by = ['New_entrenzID']).reset_index(drop=True)
en = []
c_drop = []
special_case = pd.DataFrame(columns = es_gs_dup.columns)
for i in range(len(es_gs_dup)):
    if es_gs_dup['New_entrenzID'][i] in en:
        continue
    en.append(es_gs_dup['New_entrenzID'][i])
    sub = es_gs_dup.iloc[np.where(es_gs_dup['New_entrenzID']==es_gs_dup['New_entrenzID'][i])[0]].reset_index(drop=True)
    id = np.where(sub['gene_symbol'] == es_gs_dup['New_gene_symbol'][i])[0]
    if len(id) == 0 or len(id)>1:
        special_case = pd.concat([special_case, sub]).reset_index(drop = True)
        continue
    id_no = [k for k in range(len(sub)) if k !=id ]   
    #c_drop.extend(np.where(es_gs2['col_ind'].values==sub['col_ind'][id_no].values)[0]+1)     
    c_drop.extend(np.where(es_gs2['col_ind'].isin(sub['col_ind'][id_no].values)==True)[0]+1)     
    #es_gs2 = es_gs2.drop(index = np.where(es_gs2['col_ind'].values==sub['col_ind'][id_no].values)[0])
    #sub = sub.drop(index = id_no)

# Manually editing special cases
c_drop.extend(np.where(es_gs2.ensembl_id == 'ENSG00000277526')[0]+1)
c_drop.extend(np.where(es_gs2.ensembl_id == 'ENSG00000248991')[0]+1)
c_drop.extend(np.where(es_gs2.ensembl_id == 'ENSG00000258325')[0]+1)
c_drop.extend(np.where(es_gs2.ensembl_id == 'ENSG00000234394')[0]+1)

gef_file = gef_file.drop(gef_file.columns[c_drop], axis=1) 
es_gs2 = es_gs2.drop(index = [x - 1 for x in c_drop]).reset_index(drop = True)

###
# Add the new gene symbols to the data file
gef_file.columns = ['DepMap_ID']+list(es_gs2['ensembl_id'])
new_row1 = pd.DataFrame([['']+list(es_gs2['New_entrenzID'])],columns = gef_file.columns )
new_row2 = pd.DataFrame([['']+list(es_gs2['New_gene_symbol'])],columns = gef_file.columns )
gef_file = pd.concat([new_row2, gef_file]).reset_index(drop = True)
gef_file = pd.concat([new_row1, gef_file]).reset_index(drop = True)
#Map DepMap ID with Argonne ID
gef_file.insert(loc=0, column='RRID', value=['' for i in range(gef_file.shape[0])])
col = gef_file.columns
rows = [gef_file.iloc[0].values]
rows.append(list(gef_file.iloc[1].values))
for i in range(len(map)):
    id = np.where(gef_file['DepMap_ID']==map['DepMap_ID'][i])
    if len(id[0])!=0:
        row_append = gef_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
gef_file_new = pd.DataFrame(rows, columns = col)
gef_file_new = gef_file_new.drop(columns = ['DepMap_ID'])
gef_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_expression_full.csv'),index=False)


######## Reverse Phase Protein Array (RPPA) ##########
logger.info('Processing RPPA data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_RPPA_20180123.csv')
rppa_file = pd.read_csv(path)
rppa_file= rppa_file.rename(columns={'Unnamed: 0':'CCLE_Name'})
rppa_file.insert(loc=0, column='RRID', value=['' for i in range(rppa_file.shape[0])])
col = rppa_file.columns
rows = []
for i in range(len(map)):
    id = np.where(rppa_file['CCLE_Name']==map['CCLE_Name'][i])
    if len(id[0])!=0:
        row_append = rppa_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
rppa_file_new = pd.DataFrame(rows, columns = col)
rppa_file_new = rppa_file_new.drop(columns = ['CCLE_Name'])
rppa_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_RPPA_20180123.csv'),index=False)


########### DNA Methylation - Reduced Representation Bisulfite Sequencing (RRBS) ##########
logger.info('Processing RRBS data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_RRBS_TSS_1kb_20180614.txt')
rrbs_file = pd.read_table(path)
#rrbs_file = pd.read_table(path,index_col=0)
#Loading Gene symbol - Entrz id map file
gs_eg_rrbs_map = pd.read_csv(file_path +'/Maps/gs_eg_gs_map_rrbs.csv').drop(columns={'Unnamed: 0'}) 
# Some of the genes in the mapping file contains aliases of genes already present in the feature set.
# This created some duplicate TSS_ids
# We need to remove those genes from the mapping file as well as the original file
dup_genes = ['C9orf47', 'ERICH1-AS1', 'ZASP', 'NARR', 'PALM2', 'AKAP2']
rrbs_file = rrbs_file[~rrbs_file['gene'].isin(dup_genes)].reset_index(drop=True)
gs_eg_rrbs_map = gs_eg_rrbs_map[~gs_eg_rrbs_map['orig_GS'].isin(dup_genes)].reset_index(drop=True)

ind_equal = np.where(rrbs_file['gene']==gs_eg_rrbs_map['new_GS'])[0]
ind_diff = np.setdiff1d(list(range(len(rrbs_file))), ind_equal) 
gene = rrbs_file['gene']
e_id = gs_eg_rrbs_map['new_entrenz_id']
ens_id = gs_eg_rrbs_map['ens_id']
tss_id = rrbs_file['TSS_id']
drop_ind = []
for ind in ind_diff:
    if not pd.isna(gs_eg_rrbs_map['new_GS'][ind]):
        gene[ind] = gs_eg_rrbs_map['new_GS'][ind]
        tss_id[ind] = gs_eg_rrbs_map['new_GS'][ind]+'_'+tss_id[ind].split('_', 1)[1]
    else:
        gene = gene.drop(ind)
        e_id = e_id.drop(ind)
        ens_id = ens_id.drop(ind)
        tss_id = tss_id.drop(ind)
        drop_ind.append(ind)
rrbs_file=rrbs_file.drop(index=drop_ind).reset_index(drop=True)
rrbs_file['TSS_id'] = tss_id.values
rrbs_file['gene'] = gene.values

rrbs_file.drop(rrbs_file.columns[[0,1, 2,3, 4, 5,6]], axis = 1, inplace = True)
cell_lines = rrbs_file.columns
rrbs_file = rrbs_file.transpose().reset_index(drop=True)
rrbs_file.columns = tss_id.values
rrbs_file.insert(loc = 0, column = 'CCLE_Name', value = cell_lines)
new_row = pd.DataFrame([['']+list(gene)],columns = rrbs_file.columns )
rrbs_file = pd.concat([new_row, rrbs_file]).reset_index(drop=True)
new_row1 = pd.DataFrame([['']+list(e_id.astype('int64'))],columns = rrbs_file.columns )
new_row2 = pd.DataFrame([['']+list(ens_id)],columns = rrbs_file.columns )
rrbs_file = pd.concat([new_row1,new_row2, rrbs_file]).reset_index(drop=True)
rrbs_file.insert(loc=0, column='RRID', value=['' for i in range(rrbs_file.shape[0])])
col = rrbs_file.columns
rows = [rrbs_file.iloc[0].values, rrbs_file.iloc[1].values]
for i in range(len(map)):
    id = np.where(rrbs_file['CCLE_Name']==map['CCLE_Name'][i])
    if len(id[0])!=0:
        row_append = rrbs_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
rrbs_file_new = pd.DataFrame(rows, columns = col)
rrbs_file_new = rrbs_file_new.drop(columns = ['CCLE_Name'])
rrbs_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_RRBS_TSS_1kb_20180614.csv'),index=False)

########## miRNA #############
logger.info('Processing miRNA data')
path = str(file_path +'/CCLE_Multiomics_Data/CCLE_miRNA_20180525.gct')
mirna_file = pd.read_table(path)
mirna_file.drop(mirna_file.columns[[0]], axis = 1, inplace = True)
mirna_file = mirna_file.set_index('Description')
mirna_file = mirna_file.transpose()
mirna_file.reset_index(inplace=True)
mirna_file=mirna_file.rename(columns={'index':'CCLE_Name'})
mirna_file.insert(loc=0, column='RRID', value=['' for i in range(mirna_file.shape[0])])
col = mirna_file.columns
rows = []
for i in range(len(map)):
    id = np.where(mirna_file['CCLE_Name']==map['CCLE_Name'][i])
    if len(id[0])!=0:
        row_append = mirna_file.iloc[id[0][0]]
        row_append['RRID'] = map['RRID'][i]
        rows.append(list(row_append.values))     
mirna_file_new = pd.DataFrame(rows, columns = col)
mirna_file_new = mirna_file_new.drop(columns = ['CCLE_Name'])
mirna_file_new.to_csv(str(file_path+'/Curated_CCLE_Multiomics_files/CCLE_AID_miRNA_20180525.csv'),index=False)

logger.info('Data curation completed!')
