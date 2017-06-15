import pandas as pd
import numpy as np

'''
script that takes in firehose methylation beta files ('.meth.by_min_expr_corr.data.txt')
for each cancer type and for each gene, Z score scale methylation beta values
combine all cancer types after scaling
output an additional table that converts scaled values to binary cutoff 
UPDATE: altered script so that it outputs in long rather than wide format
'''
# read in DF with cancer types
CANCERS=pd.read_csv('/nrnb/users/abuckley/methylation/firehose/stddata__2016_01_28/min_exp_corr/cancers', sep='\t', header=None, names=['CANCER'])

# function to log transform and scale methylation beta values
def scale_beta(x):
    transformed_vals=np.log10(x)
    mean=np.mean(transformed_vals)
    std=np.std(transformed_vals)
    scaled_vals=(transformed_vals - mean)/std
    return(scaled_vals)

# read in each cancer type methylation file
# scale values by gene
# append scaled DF to list

DF_list=[]
for i in CANCERS.CANCER:
	file='/nrnb/users/abuckley/methylation/firehose/stddata__2016_01_28/min_exp_corr/' + i + '.methyl.clean.txt'
	DF=pd.read_csv(file, sep='\t')
	mat=DF.iloc[:,1:].apply(scale_beta, axis=1)
	genes=DF.ix[:,0]
	mat['GENE']=genes
	mat_long=pd.melt(mat, id_vars='GENE', var_name="ID", value_name="GENO")
	DF_list.append(mat_long)

# merge all cancer types together after scaling
#FIN=reduce(lambda x, y: pd.merge(x, y, on = 'GENE', how="outer"), DF_list)
FIN=pd.concat(DF_list)
FIN.fillna("NA", inplace=True)

# Write scaled data to TSV
FIN.to_csv('/nrnb/users/abuckley/methylation/firehose/meth_expr/scaled_files/pancan.methylation.Zscore.long.txt', sep='\t', index=False)

# Binarize methylation values based on z-score > 3 cutoff
FIN.GENO[FIN.GENO < 3]=0
FIN.GENO[FIN.GENO > 0]=1
FIN.to_csv('/nrnb/users/abuckley/methylation/firehose/meth_expr/scaled_files/pancan.methylation.Zscore.binary.long.txt', sep='\t', index=False)
