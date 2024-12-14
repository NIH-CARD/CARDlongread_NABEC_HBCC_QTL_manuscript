# Genome Window Methylation and Age Regressions
# 12-2024 
# Author: Melissa Meredith; UCSC, NIH CARD

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import math
from numpy.polynomial import Polynomial
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import seaborn as sns


def regression_covs(data_df, cohort_region, date):
    """ regression function """

    region_lreg = {}

    for column in data_df.columns:
        if column not in ['Age','Gender_male','GiB_of_uBAM','PMI','Region']:
            clean_data_df = data_df[['Age', 'PMI', 'Gender_male', column]].dropna()

            X = clean_data_df[['Age', 'PMI', 'Gender_male']]  # I variables
            y = clean_data_df[[column]]  # D variable must be 2D array  
            # if there are no rows, skip it
            if X.shape[0] == 0:
                continue
        
            # add a constant term (intercept) to the independent variable
            X = sm.add_constant(X)
    
            # fit the regression model
            model = sm.OLS(np.asarray(y), np.asarray(X)).fit()
            
            # get p-values  
            region_lreg[column] = { 'r2': model.rsquared,
                                      'intercept':model.params[0],
                                      'Age_Coeff': model.params[1], 
                                      'PMI_Coeff': model.params[2],
                                      'Gender_coeff':model.params[3],
                                      'P_Value_Age': model.pvalues[1],
                                      'P_Value_PMI': model.pvalues[2],
                                      'P_Value_Gender_Female': model.pvalues[3],
                                      
                                     }

    # store the LR output as dataframe, separate parameters, and write to csv/bed
    data_cov_lr_df = pd.DataFrame(region_lreg).T
    
    data_cov_lr_df.to_csv(f'{cohort_region}_olsRcov_Age_PMI_Gender_{date}.csv', sep='\t', header=False,index=False)

    return data_cov_lr_df


# get Metadata for both cohorts
# NABEC metadata
nabec_meta = pd.read_csv('NABEC_Metadata.tsv',  delimiter='\t')
nabec_meta.set_index(['entity:NABEC_meta_cohort_id'], inplace=True)

# HBCC metadata filterd
hbcc_meta = pd.read_csv('HBCC_covs.csv')
hbcc_meta['HBCC_Cohort_id'] = hbcc_meta['SampleId'] + "_FTX"
hbcc_meta.drop(["Ethnicity", "Race", "SampleId"], axis=1, inplace=True)
hbcc_meta.set_index(["HBCC_Cohort_id"], inplace=True)
hbcc_meta.rename(columns={"Sex at Birth": "Gender", "AgeDeath":"Age"}, inplace=True)

# drop samples that are not FTX 
hbcc_meta_FTX = hbcc_meta.loc[hbcc_meta['Region']=="PFC"]
not_FTX_samples = hbcc_meta.loc[hbcc_meta['Region']!="PFC"].index.tolist()

# one-hot encode gender 
nabec_meta = pd.get_dummies(nabec_meta, columns=['Gender'], drop_first=True)  
nabec_meta['Gender_male'] = nabec_meta['Gender_male'].astype(int)

hbcc_meta_FTX['Gender'] = hbcc_meta_FTX['Gender'].str.lower()
hbcc_meta_FTX = pd.get_dummies(hbcc_meta_FTX, columns=['Gender'], drop_first=True)


# reformat the group beds, methylation was averaged over these windows on Terra using bedtools map
""" NABEC """
nabec_wg1kb_bed = pd.read_csv('nabec_WG_1kb_avgmod.bed',  delimiter='\t', na_values=['.'])

dataCols = ['chrom', 'start','end']
wg1kb_nabec_metadata = nabec_wg1kb_bed[dataCols]
nabec_wg1kb_bed['position'] = nabec_wg1kb_bed['chrom'] + '_' + nabec_wg1kb_bed['start'].astype(str) + '_' + nabec_wg1kb_bed['end'].astype(str)
nabec_wg1kb_bed.drop(['chrom', 'start','end'], axis=1, inplace=True)
nabec_wg1kb_bed.set_index(['position',], inplace=True)

nabec_wg1kb_bed.rename(columns=lambda x: x.replace('avgMod_', '') if isinstance(x, str) and x.startswith('avgMod_') else x, inplace=True)

"""  HBCC"""
hbcc_wg1kb_bed = pd.read_csv('hbcc_WG_1kb_avgmod.bed',  delimiter='\t', na_values=['.'])

dataCols = ['chrom', 'start','end']
wg1kb_hbcc_metadata = hbcc_wg1kb_bed[dataCols]
hbcc_wg1kb_bed['position'] = hbcc_wg1kb_bed['chrom'] + '_' + hbcc_wg1kb_bed['start'].astype(str) + '_' + hbcc_wg1kb_bed['end'].astype(str)
hbcc_wg1kb_bed.drop(['chrom', 'start','end'], axis=1, inplace=True)
hbcc_wg1kb_bed.set_index(['position',], inplace=True)

hbcc_wg1kb_bed.rename(columns=lambda x: x.replace('avgMod_', '') if isinstance(x, str) and x.startswith('avgMod_') else x, inplace=True)




""" NABEC whole genome 1kb windows Age Regression with covariates: Gender, PMI 
        age ~ meth + sex + PMI
"""


# drop sex chromosomes
nabec_wg1kb_bed_f = nabec_wg1kb_bed[~nabec_wg1kb_bed.index.str.contains('X|Y')]


nabec_wg1kb_age = nabec_meta_num.join(nabec_wg1kb_bed_f.T, how = 'inner')

# Perform linear regression for each column separately 
""" regression with Age PMI Sex """

nabec_wg1kb_cov_lr_df = regression_covs(nabec_wg1kb_age, "nabec_wg1kb", "11252024")

# Apply Benjamini-Hochberg correction
alpha = 0.05  # significance level
reject, p_values_corrected, _, fdr = multipletests(nabec_wg1kb_cov_lr_df['p'], alpha=alpha, method='fdr_bh')

nabec_wg1kb_cov_lr_df['bh_corrected_p'] = p_values_corrected
nabec_wg1kb_cov_lr_df['neg_log10_bh_corrected_p'] = -np.log10(p_values_corrected)

nabec_wg1kb_cov_lr_df.loc[nabec_wg1kb_cov_lr_df['bh_corrected_p']<alpha].to_csv('nabec_wg1kb_olsRcov_Age_PMI_Gender_sig_11252024.csv',index=True,header=True,sep=",")

""" Reformat the chr position to bed format """
nabec_wg1kb_lr_df_slopeI_noi = nabec_wg1kb_cov_lr_df.reset_index()
df_split = nabec_wg1kb_lr_df_slopeI_noi['position'].str.split('_', expand=True)
df_split.columns = ['chrom', 'start','end']

nabec_wg1kb_lr_df_slopeI_noi['chrom'] = df_split['chrom']
nabec_wg1kb_lr_df_slopeI_noi['start'] = df_split['start']
nabec_wg1kb_lr_df_slopeI_noi['end'] = df_split['end']



# select b-h sig 
""" Write out all the B-H Significant hits; and then those with slope (-.05 > S > 0.05) """
nabec_wg1kb_lr_df_slopeI_noi_bhSig = nabec_wg1kb_lr_df_slopeI_noi.loc[nabec_wg1kb_lr_df_slopeI_noi['bh_corrected_p'] < 0.05]

nabec_wg1kb_lr_df_slopeI_noi_bhSig[['chrom', 'start','end','p']].to_csv('nabec_wg1kb_lr_Age_PMI_Gender_bh_sig.bed',index=False,header=False,sep="\t")
# nabec_wg1kb_lr_df_slopeI_noi_bhSig.loc[(nabec_wg1kb_lr_df_slopeI_noi_bhSig['slope'] > 0.1) | (nabec_wg1kb_lr_df_slopeI_noi_bhSig['slope'] < -0.05)][['chrom', 'start','end','slope']].to_csv('nabec_wg1kb_lr_df_bh_sig_slope1_05.bed',index=False,header=False,sep="\t")

nabec_wg1kb_lr_df_slopeI_noi_bhSig.loc[(nabec_wg1kb_lr_df_slopeI_noi_bhSig['slope'] > 0.1) & (nabec_wg1kb_lr_df_slopeI_noi_bhSig['r2'] > 0.4)][['chrom', 'start','end','slope']].to_csv('nabec_wg1kb_lr_Age_PMI_Gender_bh_sig_slope1_r2.4.bed',index=False,header=False,sep="\t")

""" WG1kb Linear Regression Volcano Plot""" 
# correlations_hbcc.sort_values(ascending=False) cpgI_hbcc_trendline_corr_df

fig, axs = plt.subplots(1, 1, figsize=(8, 10)) 

sns.scatterplot(x='slope', y='neg_log10_bh_corrected_p', data=nabec_wg1kb_lr_df_slopeI, hue='neg_log10_bh_corrected_p', 
                palette='coolwarm', edgecolor='k', s=28, ax=axs)
bh_num_sig = nabec_wg1kb_lr_df_slopeI.loc[(nabec_wg1kb_lr_df_slopeI['bh_corrected_p']<0.05)].shape[0]
axs.axhline(y=-np.log10(alpha), color='red', label=f'FDR=0.05: {bh_num_sig} sig hits')

axs.legend()
axs.set_title('1Kb Whole Genome (50 CpG) Avg Methylation Linear Regression x Age')
axs.set_ylabel("B-H Corrected - log10(p)")
axs.set_xlabel("Slope of Avg Methylation x Age Linear Regression")
axs.grid(True)

# t = plt.xticks(rotation=90)

plt.savefig('WG1Kb_Age_r2_x_LinearRegressionSlopeVolcano_Benjamini-Hochberg_NABEC.png',dpi=300, facecolor='white', transparent=False)





"""  HBCC whole genome 1kb windows Age Regression 
        age ~ meth + sex + PMI

"""

# drop sex chromosomes
hbcc_wg1kb_bed_f = hbcc_wg1kb_bed[~hbcc_wg1kb_bed.index.str.contains('X|Y')]
# remove any non prefrontal cortex samples 
hbcc_wg1kb_bed_f.drop(not_FTX_samples, axis=1, inplace=True)


# merge methylation and metadata 
hbcc_wg1kb_age = hbcc_meta_num.join(hbcc_wg1kb_bed_f.T, how = 'inner')

# Perform linear regression 
""" regression with Age PMI Sex """
hbccc_wg1kb_cov_lr_df = regression_covs(hbcc_wg1kb_age, "hbccc_wg1kb", "11272024")

# Apply Benjamini-Hochberg correction
alpha = 0.05  # significance level
reject, p_values_corrected, _, fdr = multipletests(hbccc_wg1kb_cov_lr_df['P_Value_Age'], alpha=alpha, method='fdr_bh')

hbccc_wg1kb_cov_lr_df['bh_corrected_p'] = p_values_corrected
hbccc_wg1kb_cov_lr_df['neg_log10_bh_corrected_p'] = -np.log10(p_values_corrected)

hbccc_wg1kb_cov_lr_df.loc[hbccc_wg1kb_cov_lr_df['bh_corrected_p']<alpha].to_csv('hbcc_wg1kb_olsRcov_Age_PMI_Gender_11272024_df_bh_sig.csv',index=True,header=True,sep=",")

""" Reformat the chr position to bed format """
hbcc_wg1kb_lr_df_slopeI_noi = hbccc_wg1kb_cov_lr_df.reset_index()
df_split = hbcc_wg1kb_lr_df_slopeI_noi['index'].str.split('_', expand=True)
df_split.drop([3,4], axis=1, inplace=True)
df_split.columns = ['chrom', 'start','end']

hbcc_wg1kb_lr_df_slopeI_noi['chrom'] = df_split['chrom']
hbcc_wg1kb_lr_df_slopeI_noi['start'] = df_split['start']
hbcc_wg1kb_lr_df_slopeI_noi['end'] = df_split['end']

# select b-h sig 
""" Write out all the B-H Significant hits; and then those with slope (-.05 > S > 0.05) """
hbcc_wg1kb_lr_df_slopeI_noi_bhSig = hbcc_wg1kb_lr_df_slopeI_noi.loc[hbcc_wg1kb_lr_df_slopeI_noi['bh_corrected_p'] < 0.05]

hbcc_wg1kb_lr_df_slopeI_noi_bhSig[['chrom', 'start','end','slope']].to_csv('hbcc_wg1kb_olsRcov_Age_PMI_Gender_11272024_df_bh_sig.bed',index=False,header=False,sep="\t")


""" WG1kb Linear Regression """ 

fig, axs = plt.subplots(1, 1, figsize=(8, 10)) 

sns.scatterplot(x='slope', y='neg_log10_bh_corrected_p', data=hbccc_wg1kb_cov_lr_df, hue='neg_log10_bh_corrected_p', 
                palette='coolwarm', edgecolor='k', s=28, ax=axs)
bh_num_sig = hbccc_wg1kb_cov_lr_df.loc[(hbccc_wg1kb_cov_lr_df['bh_corrected_p']<0.05)].shape[0]
axs.axhline(y=-np.log10(alpha), color='red', label=f'FDR=0.05: {bh_num_sig} sig hits')

axs.legend()
axs.set_title('HBCC 1Kb Whole Genome (50 CpG) Avg Methylation Linear Regression x Age')
axs.set_ylabel("B-H Corrected - log10(p)")
axs.set_xlabel("Slope of Avg Methylation x Age Linear Regression")
axs.grid(True)

# t = plt.xticks(rotation=90)

plt.savefig('WG1Kb_Age_r2_x_LinearRegressionSlopeVolcano_Benjamini-Hochberg_HBCC.png',dpi=300, facecolor='white', transparent=False)


### Totals and Intersections: 
nabec_cov_ols = pd.read_csv('nabec_wg1kb_olsRcov_Age_PMI_Gender_sig_11252024.csv', sep='\t')
nabec_cov_ols.rename(columns={"Unnamed: 0":"position"},inplace=True)
# nabec_cov_ols.set_index("position", inplace=True)

hbcc_cov_ols = pd.read_csv('hbcc_wg1kb_olsRcov_Age_PMI_Gender_sig_11272024.csv', sep='\t')
hbcc_cov_ols.drop("Unnamed: 0.1", axis=1, inplace=True)
hbcc_cov_ols.rename(columns={"Unnamed: 0":"position"}, inplace=True)

# How many windows are significantly associated with age in both cohorts? 
cohort_merged = pd.merge(hbcc_cov_ols, nabec_cov_ols, on="position", how="inner", suffixes=("_h","_n"))
# 26,707 positions

# How many shared positions are increasing with age? 
print('increasing',cohort_merged.loc[(cohort_merged['Age_Coeff_h']>0) & (cohort_merged['Age_Coeff_n']>0)].shape)
# 3,556

# How many shared positions are decreasing with age? 
print('decreasing',cohort_merged.loc[(cohort_merged['Age_Coeff_h']<0) & (cohort_merged['Age_Coeff_n']<0)].shape)
# 23,135

# How many merged positions have a slope >1
print('slope>0.1',cohort_merged.loc[(cohort_merged['Age_Coeff_h']>0.1) & (cohort_merged['Age_Coeff_n']>0.1)].shape)
# 159


######### CGI ###############
""" NABEC CGI Regression """

# load in NABEC CGI bed
dataCols = ['chrom', 'start','end','name','length','cpgNum','gcNum','perCpg','perGc','obsExp']
cpgI_nabec_bed = pd.read_csv('nabec_cpg_islands_bmap_avgmod.bed',  delimiter='\t', na_values=['.'])
cpgI_nabec_metadata = cpgI_nabec_bed[dataCols]
cpgI_nabec_bed['position'] = cpgI_nabec_bed['chrom'] + '_' + cpgI_nabec_bed['start'].astype(str) + '_' + cpgI_nabec_bed['end'].astype(str)
cpgI_nabec_bed.drop(['chrom', 'start','end','name','length','cpgNum','gcNum','perCpg','perGc','obsExp'], axis=1, inplace=True)
cpgI_nabec_bed.set_index(['position',], inplace=True)

# drop sex chromosomes
cpgI_nabec_bed_f = cpgI_nabec_bed[~cpgI_nabec_bed.index.str.contains('X|Y')]

# merge methylation and metadata 
nabec_cgi_age = nabec_meta.join(cpgI_nabec_bed_f.T, how = 'inner')

# Perform linear regression for each column separately after removing NaN values
nabec_cgi_cov_lr_df = regression_covs(nabec_cgi_age, "nabec_cgi", "12062024")

# Apply Benjamini-Hochberg correction
alpha = 0.05  # significance level

reject, p_age_values_corrected, _, fdr = multipletests(nabec_cgi_cov_lr_df['P_Value_Age'], alpha=alpha, method='fdr_bh')

nabec_cgi_cov_lr_df['bh_corrected_age_p'] = p_age_values_corrected
nabec_cgi_cov_lr_df['neg_log10_bh_corrected_age_p'] = -np.log10(p_age_values_corrected)

nabec_cgi_cov_lr_df.loc[nabec_cgi_cov_lr_df['bh_corrected_age_p']<alpha].to_csv('nabec_cgi_olsRcov_Age_PMI_Gender_sig_12062024.csv',index=True,header=True,sep="\t")

bhsig = nabec_cgi_cov_lr_df.loc[nabec_cgi_cov_lr_df['bh_corrected_age_p']<alpha]
bhsig.reset_index(inplace=True)

split_coordinates = bhsig['index'].str.split("_", expand=True)
# split_coordinates
bhsig[['chr','start','end']] = split_coordinates
bhsig.drop('index', axis=1, inplace=True)
bhsig.head()
bhsig[['chr','start','end','r2','Age_Coeff','intercept','bh_corrected_age_p']].to_csv('nabec_cgi_olsRcov_Age_PMI_Gender_sig_12062024.bed',
                                                                                            index=False,header=False,sep="\t") 
# Plot the results: 
fig, axs = plt.subplots( 2,3, figsize=(15,10) ) 

axs[0,0].hist(bhsig['bh_corrected_age_p'])
axs[0,0].set_title("NABEC Age BH Pvalue")

axs[0,1].hist(bhsig['Age_Coeff'], log=True)
axs[0,1].set_title("Age coeff")
axs[0,2].hist(bhsig['PMI_Coeff'], log=True)
axs[0,2].set_title("PMI coeff")
axs[1,0].hist(bhsig['Gender_coeff'], log=True)
axs[1,0].set_title("Gender coeff")
axs[1,1].hist(bhsig['r2'], log=True)
axs[1,1].set_title("R2")

axs[1,2].scatter( bhsig['Age_Coeff'], 
                 -np.log10(bhsig['bh_corrected_age_p']),
                  s = 2, alpha=0.6)
axs[1,2].set_title("Age_Coefbhsigf x -log10(p) NABEC")

plt.savefig("NABEC_CGI_volcano_age_cov_12062024.png",dpi=300)




""" HBCC CGI Regression """

"""  HBCC CGI group BED """
# load in HBCC CGI bed
hbcc_cgi_bed = pd.read_csv('hbcc_cpg_islands_bmap_avgmod.bed',  delimiter='\t', na_values=['.'])

dataCols = ['chrom', 'start','end']
cgi_hbcc_metadata = hbcc_cgi_bed[dataCols]
hbcc_cgi_bed['position'] = hbcc_cgi_bed['chrom'] + '_' + hbcc_cgi_bed['start'].astype(str) + '_' + hbcc_cgi_bed['end'].astype(str)
hbcc_cgi_bed.drop(['chrom', 'start','end'], axis=1, inplace=True)
hbcc_cgi_bed.set_index(['position',], inplace=True)

hbcc_cgi_bed.rename(columns=lambda x: x.replace('avgMod_', '') if isinstance(x, str) and x.startswith('avgMod_') else x, inplace=True)

# drop sex chromosomes
cpgI_hbcc_bed_f = cpgI_hbcc_bed[~cpgI_hbcc_bed.index.str.contains('X|Y')]

# drop samples that are not FTX 
hbcc_meta_FTX = hbcc_meta_num.loc[hbcc_meta_num['Region']=="PFC"]
not_FTX_samples = hbcc_meta_num.loc[hbcc_meta_num['Region']!="PFC"].index.tolist()


# # merge methylation and metadata 
hbcc_cgi_age = hbcc_meta.join(cpgI_hbcc_bed_f.T, how = 'inner')

# Perform linear regression for each column separately after removing NaN values
hbcc_cgi_cov_lr_df = regression_covs(hbcc_cgi_age, "hbcc_cgi", "12062024")


# Apply Benjamini-Hochberg correction
alpha = 0.05  # significance level

reject, p_age_values_corrected, _, fdr = multipletests(hbcc_cgi_cov_lr_df['P_Value_Age'], alpha=alpha, method='fdr_bh')

hbcc_cgi_cov_lr_df['bh_corrected_age_p'] = p_age_values_corrected
hbcc_cgi_cov_lr_df['neg_log10_bh_corrected_age_p'] = -np.log10(p_age_values_corrected)

hbcc_cgi_cov_lr_df.loc[hbcc_cgi_cov_lr_df['bh_corrected_age_p']<alpha].to_csv('hbcc_cgi_olsRcov_Age_PMI_Gender_sig_12062024.csv',index=True,header=True,sep="\t")

hbcc_bhsig = hbcc_cgi_cov_lr_df.loc[hbcc_cgi_cov_lr_df['bh_corrected_age_p']<alpha]
hbcc_bhsig.reset_index(inplace=True)

split_coordinates = hbcc_bhsig['index'].str.split("_", expand=True)
# split_coordinates
hbcc_bhsig[['chr','start','end']] = split_coordinates
hbcc_bhsig.drop('index', axis=1, inplace=True)
hbcc_bhsig[['chr','start','end','r2','Age_Coeff','intercept','bh_corrected_age_p']].to_csv('hbcc_cgi_olsRcov_Age_PMI_Gender_sig_12062024.bed',
                                                                                            index=False,header=False,sep="\t")

# Plot the results: 
fig, axs = plt.subplots( 2,3, figsize=(15,10) ) 

axs[0,0].hist(hbcc_bhsig['bh_corrected_age_p'])
axs[0,0].set_title("NABEC Age BH Pvalue")

axs[0,1].hist(hbcc_bhsig['Age_Coeff'], log=True)
axs[0,1].set_title("Age coeff")
axs[0,2].hist(hbcc_bhsig['PMI_Coeff'], log=True)
axs[0,2].set_title("PMI coeff")
axs[1,0].hist(hbcc_bhsig['Gender_coeff'], log=True)
axs[1,0].set_title("Gender coeff")
axs[1,1].hist(hbcc_bhsig['r2'], log=True)
axs[1,1].set_title("R2")

axs[1,2].scatter( hbcc_bhsig['Age_Coeff'], 
                 -np.log10(hbcc_bhsig['bh_corrected_age_p']),
                  s = 2, alpha=0.6)
axs[1,2].set_title("Age_Coeff x -log10(p) HBCC")

plt.savefig("HBCC_CGI_volcano_age_cov_12062024.png",dpi=300)


## Totals and Intersections: 

# CGI NABEC + HBCC
nabec_sig_cgi = pd.read_csv("nabec_cgi_olsRcov_Age_PMI_Gender_sig_12062024.csv", sep="\t")
nabec_sig_cgi.rename(columns={"Unnamed: 0":'position'},inplace=True)
hbcc_sig_cgi = pd.read_csv("hbcc_cgi_olsRcov_Age_PMI_Gender_sig_12062024.csv", sep="\t")
hbcc_sig_cgi.rename(columns={"Unnamed: 0":'position'},inplace=True)

merged_cgi_cohorts = pd.merge(nabec_sig_cgi, hbcc_sig_cgi, on=['position'], suffixes=("_h","_n") )

nabec_sig_cgi.shape
# 3890
nabec_sig_cgi.loc[nabec_sig_cgi['Age_Coeff']>0].shape
# 2214
nabec_sig_cgi.loc[nabec_sig_cgi['Age_Coeff']<0].shape
#1676

hbcc_sig_cgi.shape
# 5778
hbcc_sig_cgi.loc[hbcc_sig_cgi['Age_Coeff']>0].shape
# 5302
hbcc_sig_cgi.loc[hbcc_sig_cgi['Age_Coeff']<0].shape
# 476

# how many shared significant CGIs? 
print(merged_cgi_cohorts.shape[0]/27949)
# 1,706: 6%

# How many shared positions are decreasing with age? 
print('decreasing',merged_cgi_cohorts.loc[(merged_cgi_cohorts['Age_Coeff_h']<0) & (merged_cgi_cohorts['Age_Coeff_n']<0)].shape)
# 205

# How many merged positions have a slope >1
print('slope>0.1',merged_cgi_cohorts.loc[(merged_cgi_cohorts['Age_Coeff_h']>0.1) & (merged_cgi_cohorts['Age_Coeff_n']>0.1)].shape)
# 170

# how many shared CGIs in the same direction
print( (205+170)/27949 )
# 1% of CGIs 

# Plot metylation data for regressions 

NPos_wg1kb = ['chr5_141408236_141409617','chr5_141409618_141410619']
genes = NPos_wg1kb

# shape subplots 
rows = math.ceil(len(genes)/6)
if rows<2:
    rows=2
cols = 6 
if len(genes)<cols:
    cols = len(genes) 
if cols<2:
    cols=2

# plot 
fig, axs = plt.subplots(rows, cols, figsize=(10, 7))  

xaxisField = "Age"


hbcc = hbcc_meta.join(hbcc_wg1kb_bed.loc[genes].T, how = 'inner')
nabec = nabec_meta.join(nabec_wg1kb_bed.loc[genes].T, how = 'inner')


# Color map for Gender
color_map = {'male': 'blue', 'female': 'red'}
colors = nabec['Gender'].map(color_map)


color_map2 = {'Male': 0, 'Female': 1}
hbcc_colors = hbcc['Sex_at_Birth'].map(color_map2)

counter = 0 
for i in np.arange(0,rows):
    if counter > len(genes)-1:
            break
    for j in np.arange(0,cols):

        print(i,j,counter, genes[counter])        
        # promoter vs Age scatter plot
        axs[i,j].scatter(hbcc['AgeDeath'], hbcc[genes[counter]], c=hbcc_colors, alpha=0.4, s=8, cmap='PiYG')
        axs[i,j].scatter(nabec[xaxisField], nabec[genes[counter]], c=colors, alpha=0.75, s=8)
        axs[i,j].set_title(genes[counter]+' vs Age')
        axs[i,j].set_xlabel(xaxisField, fontsize=16)
        axs[i,j].grid(True)
    
        # Add trend line
        x = hbcc['AgeDeath'].to_numpy()
        y = hbcc[genes[counter]].to_numpy()
        p = Polynomial.fit(x, y, 1)  
        hbcc_coefficients = p.convert().coef
        hcorr = np.corrcoef(x, y)[0, 1]

        x_fit = np.linspace(min(x), 100, 100) #max(x), 100)
        y_fit = p(x_fit)
    
        # # plot fitted polynomial
        axs[i,j].plot(x_fit, y_fit, color='#ff7f0e', label=f'HBCC: {hbcc_coefficients[1]:.4f}; {hcorr:.2f}')
    
        # Add trend line
        x = nabec['Age']
        y = nabec[genes[counter]]
        p = Polynomial.fit(x, y, 1)  
        nabec_coefficients = p.convert().coef
        ncorr = np.corrcoef(x, y)[0, 1]
        
        # plot the line
        x_fit = np.linspace(min(x), 100, 100) #max(x), 100)
        y_fit = p(x_fit)
        axs[i,j].plot(x_fit, y_fit, color='#1f77b4', label=f'NABEC: {nabec_coefficients[1]:.4f}; {ncorr:.2f}')
        axs[i,j].tick_params(axis='both', which='major', labelsize=18)
        axs[i,j].legend()  #loc='lower right')
        axs[i,j].set_ylim(20,60)
        counter+=1
        if counter > len(genes)-1:
            break
        
        
axs[0,1].set_yticklabels([])


plt.tight_layout()

plt.savefig('WG1kg_Methylation_HighAgeCorr_and_LR_0.1Slope_NABEC_PCDHA_B-HCorr_&_HBCC..png',dpi=300, facecolor='white', transparent=False)

