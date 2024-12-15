# merging unphased bedmethyls from bedtools map
# 12-2024 
# Author: Melissa Meredith; UCSC, NIH CARD

import terra_notebook_utils as tnu
import terra_pandas as tp
import os
import io
import gzip
import math
import pandas as pd
import numpy as np
import tarfile
import matplotlib.pyplot as plt
import pysam
import bgzip
from matplotlib.lines import Line2D
import seaborn as sns
from datetime import datetime
from IPython.display import Image, display

from importlib import reload
import import_ipynb
import phased_methylation_merge_functions as pm


COHORTS = ['NABEC','HBCC']

# Methylation tables have phased bedMethyls
nabec_cohort_meth = tp.table_to_dataframe("NABEC_Cohort_Methylation", workspace=WORKSPACE, workspace_namespace=PROJECT)
hbcc_cohort_meth = tp.table_to_dataframe("HBCC_Cohort_Methylation", workspace=WORKSPACE, workspace_namespace=PROJECT)


base_gs_url="gs://fc-secure-578c07ca-f188-44a1-a254-b524a5e73ecf/NABEC_plots"

#-------------
#1kb Whole Genome
#-------------

def load_in_bedtoolsMap_bed(filePath, sample):
    """ bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end',
               'ratio_'+sample, 'avgMod_'+sample],
    usecols=['chrom', 'start', 'end', 'ratio_'+sample, 'avgMod_'+sample])

    return df

region = 'wholeGenome_1Kb_50CpG'
bed_Fields = ['chrom', 'start', 'end'] 

print('starting ',datetime.now())

pm.analyze_regional_calcmeth_beds('hbcc', 
                               hbcc_cohort_meth, 
                               "WG_1kb", 
                               region, 
                               load_in_bedtoolsMap_bed, 
                              bed_Fields, do_cov_work=False, variance_plots=True)

print('done ',datetime.now())



directory_path="WG_1kb"
cohort='hbcc'


pm.calc_var_std(cohort_dfm_min_cov_avgmod[cohort_dfm_min_cov_avgmod.index.get_level_values('chrom').str.len() < 6], 
                     cohort, directory_path, 
                     regional_bed_columns[1], high_var_cutoff=0.1)


cohort_dfm_min_cov_avgmod[cohort_dfm_min_cov_avgmod.index.get_level_values('chrom').str.len() < 6]


directory_path="WG_1kb"
cohort='nabec'

default_ratio_value = 0.9

pm.plot_heatmap(pm.regional_sample_coverage(cohort_dfm_min_cov_avgmod), cohort, 
             directory_path, default_ratio_value)


#--------
#TSS 5 CpG Min
#--------
# nabec_cohort_meth = tp.table_to_dataframe("NABEC_Cohort_Methylation", workspace=WORKSPACE, workspace_namespace=PROJECT)
# hbcc_cohort_meth = tp.table_to_dataframe("HBCC_Cohort_Methylation", workspace=WORKSPACE, workspace_namespace=PROJECT)


region = 'cagePeaks_5cpgs'
cage_peak_Fields = ['chrom', 'start', 'end', 'geneInfo','num','strand'] #,'start2','end2']

print('starting ',datetime.now())

pm.analyze_regional_calcmeth_beds('hbcc', 
                               hbcc_cohort_meth, 
                               "cagePeaks_5cpgs", 
                               region, 
                               pm.load_in_bedtoolsMap_cagePeaksTSS, 
                              cage_peak_Fields, do_cov_work=False, variance_plots=True)

print('done ',datetime.now())



cohort = 'nabec'
directory_path = "cagePeaks_5cpgs"

dfm = pd.read_csv(directory_path+"/nabec_cagePeaks_5cpgs_avgmod.bed", sep="\t")





#----------
#Promoters 2kb
#---------

region = 'Hs_EPDnew_006_hg38_2kb_bed'
promoterFields = ['chrom', 'start', 'end', 'promoter_name',  'number', 'strand']

print('starting ',datetime.now())
region = 'Hs_EPDnew_006_hg38_2kb_bed'
pm.analyze_regional_calcmeth_beds('hbcc', 
                               hbcc_cohort, 
                               "promoters2Kb_bmap3", 
                               region, 
                               pm.load_in_bedtoolsMap_promoters, 
                               promoterFields, do_cov_work=False, variance_plots=True)

print('done ',datetime.now())


!bedtools intersect -a promoters2Kb_bmap3/nabec_promoters2Kb_bmap_avgmod.unfiltered.bed -b promoters2Kb_bmap3/nabec_promoters2Kb_bmap3_avgmod.bed -v | awk -v OFS="/t" '{print $1,$2,$3,$4}' >   promoters2Kb_bmap3/nabec_CoVfiltered_regions.bed


nabec_pro = pd.read_csv("promoters2Kb_bmap3/nabec_promoters2Kb_bmap3_avgmod.bed", 
                          sep="\t")
hbcc_pro = pd.read_csv("promoters2Kb_bmap3/hbcc_promoters2Kb_bmap3_avgmod.bed", 
                          sep="\t")

nabec_pro.drop(['promoter_name',  'number', 'strand'], axis=1, inplace=True)
hbcc_pro.drop(['promoter_name',  'number', 'strand'], axis=1, inplace=True)

for col in nabec_pro.columns[6:]:
    nabec_pro[col] = pd.to_numeric(nabec_pro[col], errors='coerce')

for col in hbcc_pro.columns[6:]:
    hbcc_pro[col] = pd.to_numeric(hbcc_pro[col], errors='coerce')
    

outer_promoters = pd.merge(nabec_pro, hbcc_pro, on=['chrom', 'start', 'end'], how='outer', indicator=True)

# Rows only in nabec
only_nabec_promoters = outer_promoters[outer_promoters['_merge'] == 'left_only']

# Rows only in hbcc
only_hbcc_promoters = outer_promoters[outer_promoters['_merge'] == 'right_only']

both_promoters = outer_promoters[outer_promoters['_merge'] == 'both']

outer_promoters.set_index(['chrom','start','end','_merge'], inplace=True)
both_promoters.set_index(['chrom','start','end','_merge'], inplace=True)

pro_pca_result, pro_pca = perform_pca(both.dropna().reset_index(drop=True).T, n_components=10)

cohorts = ['NABEC']*206 + ['HBCC']*155
pro_pca_result['Category'] = cohorts

plt.figure(figsize=(10, 7))
sns.scatterplot(x='PC1', y='PC2', hue='Category', data=pro_pca_result, palette='viridis')
plt.title('PCA of Promoters')
plt.xlabel(f"PC 1 ({round(explained_variance[0],4)} variance explained)")
plt.ylabel(f"PC 2 ({round(explained_variance[1],4)} variance explained)")
plt.legend()
plt.show()



#------------
CpG Islands
#------------

islandFields = ['chrom','start','end','name','length','cpgNum','gcNum','perCpg','perGc','obsExp']
pm.analyze_regional_calcmeth_beds('nabec', 
                               nabec_cohort, 
                               "cpg_islands_bmap", 
                               "cpg_hg38_bed", 
                               pm.load_in_bedtoolsMap_islands, 
                               islandFields, do_cov_work=False, variance_plots=False)

print('done ',datetime.now())


print('starting ',datetime.now())
pm.analyze_regional_calcmeth_beds('hbcc', 
                               hbcc_cohort, 
                               "cpg_islands_bmap", 
                               "cpg_hg38_bed", 
                               pm.load_in_bedtoolsMap_islands, 
                               islandFields, do_cov_work=False, variance_plots=True)

print('done ',datetime.now())


# Which cpg islands are in both cohorts? 

nabec_cpgI = pd.read_csv("cpg_islands_bmap/nabec_cpg_islands_bmap_avgmod.bed", 
                          sep="\t", dtype="object")
hbcc_cpgI = pd.read_csv("cpg_islands_bmap/hbcc_cpg_islands_bmap_avgmod.bed", 
                          sep="\t", dtype="object")

nabec_cpgI.drop(['gcNum','perCpg','perGc','obsExp'], axis=1, inplace=True)
hbcc_cpgI.drop(['gcNum','perCpg','perGc','obsExp'], axis=1, inplace=True)

for col in nabec_cpgI.columns[6:]:
    nabec_cpgI[col] = pd.to_numeric(nabec_cpgI[col], errors='coerce')

for col in hbcc_cpgI.columns[6:]:
    hbcc_cpgI[col] = pd.to_numeric(hbcc_cpgI[col], errors='coerce')


outer = pd.merge(nabec_cpgI, hbcc_cpgI, on=['chrom', 'start', 'end'], how='outer', indicator=True)

# Rows only in nabec
only_nabec = outer[outer['_merge'] == 'left_only']

# Rows only in hbcc
only_hbcc = outer[outer['_merge'] == 'right_only']

both = outer[outer['_merge'] == 'both']

outer.set_index(['chrom','start','end','name_x','length_x','cpgNum_x','_merge'], inplace=True)

outer.drop(['name_y','length_y','cpgNum_y'], axis=1, inplace=True)

both.set_index(['chrom','start','end','name_x','length_x','cpgNum_x','_merge'], inplace=True)

both.drop(['name_y','length_y','cpgNum_y'], axis=1, inplace=True)



Q1 = outer.quantile(0.25)
Q2 = outer.median()
Q3 = outer.quantile(0.75)
Q4 = outer.max()
quartiles_df = pd.DataFrame({'Q3': Q3, 'Q4': Q4})

nabec_color = '#1f77b4'
hbcc_color = '#ff7f0e'


plt.figure(figsize=(10, 7))
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 6))

sns.scatterplot(x=np.arange(206), y=Q1[:206], color=nabec_color, ax=axs[0]) 
sns.scatterplot(x=np.arange(207,len(Q1)), y=Q1[207:], color=hbcc_color, ax=axs[0])

sns.scatterplot(x=np.arange(206), y=Q3[:206], color=nabec_color) 
sns.scatterplot(x=np.arange(207,len(Q3)), y=Q3[207:], color=hbcc_color, ax=axs[0])

sns.scatterplot(x=np.arange(206), y=Q4[:206], color=nabec_color) 
sns.scatterplot(x=np.arange(207,len(Q4)), y=Q4[207:], color=hbcc_color, ax=axs[0])

plt.xticks([])
plt.xlabel('Principal Components')
plt.ylabel('Methylation Ratio ')
plt.title('CpG Island Q1, Q3 and Q4 Sample Distribution')



#-------
# quartiles all regions
#---------

""" Islands """
Q1 = both.quantile(0.25)
Q2 = both.median()
Q3 = both.quantile(0.75)
Q4 = both.max()
quartiles_df = pd.DataFrame({
    'Q1_NABEC': Q1[:206], 'Q1_HBCC': Q1[207:],
    'Q2_NABEC': Q2[:206], 'Q2_HBCC': Q2[207:],
    'Q3_NABEC': Q3[:206], 'Q3_HBCC': Q3[207:],
    'Q4_NABEC': Q4[:206], 'Q4_HBCC': Q4[207:],
}).reset_index().melt(id_vars='index', var_name='Quartile', value_name='Value')

""" Promoters 2kb """
Q1_pro = both_promoters.quantile(0.25)
Q2_pro = both_promoters.median()
Q3_pro = both_promoters.quantile(0.75)
Q4_pro = both_promoters.max()
quartiles_pro = pd.DataFrame({
    'Q1_NABEC': Q1_pro[:206], 'Q1_HBCC': Q1_pro[207:],
    'Q2_NABEC': Q2_pro[:206], 'Q2_HBCC': Q2_pro[207:],
    'Q3_NABEC': Q3_pro[:206], 'Q3_HBCC': Q3_pro[207:],
    'Q4_NABEC': Q4_pro[:206], 'Q4_HBCC': Q4_pro[207:],
}).reset_index().melt(id_vars='index', var_name='Quartile', value_name='Value')

""" Gene Bodies """ 
Q1_gb = both_gb.quantile(0.25)
Q2_gb = both_gb.median()
Q3_gb = both_gb.quantile(0.75)
Q4_gb = both_gb.max()
quartiles_gb = pd.DataFrame({
    'Q1_NABEC': Q1_gb[:206], 'Q1_HBCC': Q1_gb[207:],
    'Q2_NABEC': Q2_gb[:206], 'Q2_HBCC': Q2_gb[207:],
    'Q3_NABEC': Q3_gb[:206], 'Q3_HBCC': Q3_gb[207:],
    'Q4_NABEC': Q4_gb[:206], 'Q4_HBCC': Q4_gb[207:],
}).reset_index().melt(id_vars='index', var_name='Quartile', value_name='Value')

""" enhancers """
Q1_enhancers = both_enhancers.quantile(0.25)
Q2_enhancers = both_enhancers.median()
Q3_enhancers = both_enhancers.quantile(0.75)
Q4_enhancers = both_enhancers.max()
quartiles_enhancers = pd.DataFrame({
    'Q1_NABEC': Q1_enhancers[:206], 'Q1_HBCC': Q1_enhancers[207:],
    'Q2_NABEC': Q2_enhancers[:206], 'Q2_HBCC': Q2_enhancers[207:],
    'Q3_NABEC': Q3_enhancers[:206], 'Q3_HBCC': Q3_enhancers[207:],
    'Q4_NABEC': Q4_enhancers[:206], 'Q4_HBCC': Q4_enhancers[207:],
}).reset_index().melt(id_vars='index', var_name='Quartile', value_name='Value')

""" TSS """
both_tss
Q1_tss = both_tss.quantile(0.25)
Q2_tss = both_tss.median()
Q3_tss = both_tss.quantile(0.75)
Q4_tss = both_tss.max()
quartiles_tss = pd.DataFrame({
    'Q1_NABEC': Q1_tss[:206], 'Q1_HBCC': Q1_tss[207:],
    'Q2_NABEC': Q2_tss[:206], 'Q2_HBCC': Q2_tss[207:],
    'Q3_NABEC': Q3_tss[:206], 'Q3_HBCC': Q3_tss[207:],
    'Q4_NABEC': Q4_tss[:206], 'Q4_HBCC': Q4_tss[207:],
}).reset_index().melt(id_vars='index', var_name='Quartile', value_name='Value')

nabec_color = '#1f77b4'
hbcc_color = '#ff7f0e'

# custom color palette for each quartile
custom_palette = {
    'Q1_NABEC': nabec_color,
    'Q1_HBCC': hbcc_color,
    'Q2_NABEC': nabec_color,
    'Q2_HBCC': hbcc_color,
    'Q3_NABEC': nabec_color,
    'Q3_HBCC': hbcc_color,
    'Q4_NABEC': nabec_color,
    'Q4_HBCC': hbcc_color
}

regions = ['CpG Island','Promoter','Gene Body','Enhancers','TSS']
cols = len(regions)
fig, axs = plt.subplots(nrows=1, ncols=cols, figsize=(12, 5)) 

# swarm plot
sns.swarmplot(x='Quartile', y='Value', data=quartiles_df, hue='Quartile', palette=custom_palette, size=2, ax=axs[0])

# promoter swarm plot
sns.swarmplot(x='Quartile', y='Value', data=quartiles_pro, hue='Quartile', palette=custom_palette, size=2, ax=axs[1])

# gene body swarm plot
sns.swarmplot(x='Quartile', y='Value', data=quartiles_gb, hue='Quartile', palette=custom_palette, size=2, ax=axs[2])

# enhancer swarm 
sns.swarmplot(x='Quartile', y='Value', data=quartiles_enhancers, hue='Quartile', palette=custom_palette, size=2, ax=axs[3])


for i in np.arange(cols):
    axs[i].set_title(f'{regions[i]}')
    axs[i].set_ylim(-5,105)
    axs[i].set_xlabel('')
    if i==0:
        axs[i].set_ylabel('Methylation Frequency')
    else:
        axs[i].set_ylabel('')
    axs[i].tick_params(axis='x', rotation=90)
    axs[i].grid(True)


# Rotate x-ticks for better readability
# axs[1].set_xticklabels(quartiles_df.index, rotation=45)

plt.tight_layout()
plt.show()

islands_pca_result, islands_pca = perform_pca(outer.dropna().reset_index(drop=True).T, n_components=10)
# islands_pca_result, islands_pca = perform_pca(both.dropna().reset_index(drop=True).T, n_components=10)

explained_variance = islands_pca.explained_variance_ratio_

plt.figure(figsize=(10, 6))
plt.bar(range(1, len(explained_variance[1:10]) + 1), explained_variance[1:10], alpha=0.7, align='center',
        label='Individual explained variance')
# plt.step(range(1, len(explained_variance[:10]) + 1), explained_variance[:10].cumsum(), where='mid',
#          label='Cumulative explained variance')

plt.xlabel('Principal Components')
plt.ylabel('Explained Variance Ratio')
plt.title('Variance Explained by CpG Island Principal Components')
plt.legend(loc='best')
plt.show()


cohorts = ['NABEC']*206 + ['HBCC']*155
islands_pca_result['Category'] = cohorts

plt.figure(figsize=(10, 7))
sns.scatterplot(x='PC1', y='PC2', hue='Category', data=islands_pca_result, palette='viridis')
plt.axis('equal')

min_val = min(islands_pca_result['PC1'].min(), islands_pca_result['PC2'].min())
max_val = max(islands_pca_result['PC1'].max(), islands_pca_result['PC2'].max())
print(min_val, max_val)
plt.xlim(min_val, max_val)
plt.ylim(min_val, max_val)

plt.title('PCA of CpG Islands')
plt.xlabel(f"PC 1 ({round(explained_variance[0],4)} variance explained)")
plt.ylabel(f"PC 2 ({round(explained_variance[1],4)} variance explained)")
plt.legend()
plt.show()


#------------
# Gene Bodies
#------------

geneBody_bed_headers = ['chrom', 'start', 'end', 'name','dot','strand']
print('starting ',datetime.now())
pm.analyze_regional_calcmeth_beds("nabec", 
                                    nabec_cohort, 
                                    "geneBodies", 
                                    'hg38_gencodeV44_genes_bed', 
                                    pm.load_in_bedtoolsMap_geneBodies, 
                                    geneBody_bed_headers, do_cov_work=False, variance_plots=False)
print('done ',datetime.now())

print('starting ',datetime.now())
pm.analyze_regional_calcmeth_beds('hbcc', 
                               hbcc_cohort, 
                               "geneBodies", 
                               "hg38_gencodeV44_genes_bed", 
                               pm.load_in_bedtoolsMap_geneBodies, 
                               geneBody_bed_headers, do_cov_work=False, variance_plots=True)

print('done ',datetime.now())




nabec_gb = pd.read_csv("geneBodies/nabec_geneBodies_avgmod.bed", 
                          sep="\t")
hbcc_gb = pd.read_csv("geneBodies/hbcc_geneBodies_avgmod.bed", 
                          sep="\t")

nabec_gb.drop(['name','dot','strand'], axis=1, inplace=True)
hbcc_gb.drop(['name','dot','strand'], axis=1, inplace=True)

for col in nabec_gb.columns[6:]:
    nabec_gb[col] = pd.to_numeric(nabec_gb[col], errors='coerce')

for col in hbcc_gb.columns[6:]:
    hbcc_gb[col] = pd.to_numeric(hbcc_gb[col], errors='coerce')
  

outer_gb = pd.merge(nabec_gb, hbcc_gb, on=['chrom', 'start', 'end'], how='outer', indicator=True)

# Rows only in nabec
only_nabec_gb = outer_gb[outer_gb['_merge'] == 'left_only']

# Rows only in hbcc
only_hbcc_gb = outer_gb[outer_gb['_merge'] == 'right_only']

both_gb = outer_gb[outer_gb['_merge'] == 'both']


print(both_gb.shape, only_nabec_gb.shape, only_hbcc_gb.shape)

outer_gb.set_index(['chrom','start','end','_merge'], inplace=True)
both_gb.set_index(['chrom','start','end','_merge'], inplace=True)


explained_variance = geneBody_pca.explained_variance_ratio_

plt.figure(figsize=(10, 6))
plt.bar(range(1, len(explained_variance[1:10]) + 1), explained_variance[1:10], alpha=0.7, align='center',
        label='Individual explained variance')


plt.xlabel('Principal Components')
plt.ylabel('Explained Variance Ratio')
plt.title('Variance Explained by Gene Body Principal Components')
plt.legend(loc='best')
plt.show()

plt.figure(figsize=(10, 7))
sns.scatterplot(x='PC1', y='PC2', hue='Category', data=geneBody_pca_result, palette='viridis')
plt.axis('equal')
plt.title('PCA of Gene Body')
plt.xlabel(f"PC 1 ({round(explained_variance[0],4)} variance explained)")
plt.ylabel(f"PC 2 ({round(explained_variance[1],4)} variance explained)")
plt.legend()
plt.show()



from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def perform_pca(dataframe, n_components=None):
    """
    Perform PCA on a DataFrame.

    Parameters:
    - dataframe: pandas DataFrame
      The input DataFrame containing numerical data.
    - n_components: int or None, optional (default=None)
      Number of components to keep. If None, all components are kept.

    Returns:
    - pca_result: pandas DataFrame
      DataFrame containing the principal components.
    - pca: sklearn.decomposition.PCA
      The fitted PCA model.
    """

    # Check if the DataFrame contains numerical data
    if pd.api.types.is_numeric_dtype(dataframe.dtypes):
        raise ValueError("Input DataFrame must contain numerical data.")

    # Standardize the data (mean=0 and variance=1)
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(dataframe)

    # Apply PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(standardized_data)

    # Create a DataFrame with the principal components
    columns = [f"PC{i+1}" for i in range(pca_result.shape[1])]
    pca_result = pd.DataFrame(pca_result, columns=columns, index=dataframe.index)

    return pca_result, pca

