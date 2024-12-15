# functions for merging unphased bedmethyls from bedtools map
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


def clean_df_write_out_files(cohort_dfm, cohort, num_samples, regional_bed_columns, 
                             directory_path, extrastr="", write_files=True):
    """ Removes duplicated input columns in group dataframe. 
        Writes out a group bed file with desired regiona_bed_columns and avgMod_, 
            and avgMod_hap1 and avgMod_hap2 for each individual. 
        Writes out the group dataframe as a tsv. 
        
        Returns: A reduced bedmethyl group dataframe for individuals in the terra table
    """
    print('cohort DataFrame shape', cohort_dfm.shape, cohort_dfm.shape[1]/num_samples, 'columns per sample' )
    
    # remove the repeated bed columns
    cohort_dfm_no_dups = remove_duplicate_columns(cohort_dfm)

    print('after removing duplicate columns:\n', cohort_dfm_no_dups.shape, 
          (cohort_dfm_no_dups.shape[1])/num_samples, 'columns per sample' )
    
    if write_files:
        # write out the whole thing as a tsv
        write_out_dataframe(cohort_dfm_no_dups, directory_path, cohort+extrastr, suffix="_calcMeth.tsv", 
                            includeIndex = False)
    
    # isolate avgMod columns
    cohort_dfm_avgmod= isolate_avgMod_cols(cohort_dfm_no_dups, cohort,regional_bed_columns )

    print('just avgmod per sample:\n', cohort_dfm_avgmod.shape, 
          (cohort_dfm_avgmod.shape[1])/num_samples, 'columns per sample')
    if write_files:
        write_out_dataframe(cohort_dfm_avgmod, directory_path, cohort+extrastr)

### Reigonal Methylation calculation: 

def load_in_regionalMeth_promoters(filePath, sample):
    """ Human promoter fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'agg_pct_mean_'+sample, 'avgMod_'+sample, 'region_std_'+sample, 'totalCov_'+sample,
           'num_cpgs_'+sample, 'pass_min_cpgs_'+sample, 'numModReads_'+sample, 'numCanonReads_'+sample, 'sample', 'extra'] )

    return df.iloc[:,:-2]


def load_in_bedtoolsMap_promoters(filePath, sample):
    """ Human promoter fields and modbamtools calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine='python',
    names=['chrom', 'start', 'end', 'promoter_name','number', 'strand', 'ratio_'+sample, 'avgMod_'+sample] )

    return df

def load_in_bedtoolsMap_islands(filePath, sample):
    """ CpG Islands fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'name','length','cpgNum','gcNum','perCpg','perGc','obsExp',
               'ratio_'+sample, 'avgMod_'+sample] )

    return df

def load_in_bedtoolsMap_geneBodies(filePath, sample):
    """ GeneBody fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'name','dot','strand',
               'ratio_'+sample, 'avgMod_'+sample] )

    return df

def load_in_bedtoolsMap_cagePeaksTSS(filePath, sample):
    """ TSS cagePeaks fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'geneInfo','num','strand','start2','end2',
               'ratio_'+sample, 'avgMod_'+sample],
    usecols=['chrom', 'start', 'end', 'geneInfo','num','strand','ratio_'+sample, 'avgMod_'+sample])

    return df

def load_in_bedtoolsMap_VISTAenhancers(filePath, sample):
    """ VISTA Enhancer fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'name','num','strand','start2','end2','color',
               'ratio_'+sample, 'avgMod_'+sample],
    usecols=['chrom', 'start', 'end', 'name','num','strand','ratio_'+sample, 'avgMod_'+sample])

    return df

def load_in_bedtoolsMap_ENCODE_cCRES(filePath, sample):
    """ human_ENCODE_cCRES fields and bedtoolsMap wdl calc meth fields"""
    df = pd.read_csv(filePath,
    sep="\t", header=0, engine="c",
    names=['chrom', 'start', 'end', 'name','num','strand','start2','end2','color',
               'ratio_'+sample, 'avgMod_'+sample],
    usecols=['chrom', 'start', 'end', 'name','num','strand','ratio_'+sample, 'avgMod_'+sample])

    return df




def plot_heatmap(d, cohort, directory_path, cov, region="", 
                 measurement="Average Mod", col_prefix="avgMod_", xls="chrom",
                 savePlots=True):
    """
        plots a heatmap using matplotlib imshow, including a color bar.         
    """


    # Select only columns that match {col_prefix}+{COHORT}
    cols = [f for f in d.columns if f[:7+len(cohort)] == col_prefix+cohort.upper()]

    fig, ax = plt.subplots(figsize=(12,12))

    im = ax.imshow(d.loc[:,tuple(cols)].T , aspect='auto')

    
    xlabels = [d.index.get_level_values(xls)[int(t)] for t in ax.get_xticks()[1:-1]  ]
    
    ax.figure.colorbar(im, ax=ax, label=str(directory_path)+" "+str(measurement), orientation="vertical")

    ax.set_yticks(np.arange(len(cols)), labels=[("_").join(c.split("_")[1:]) for c in cols], fontsize=5)
    ax.set_xticklabels(["0"]+xlabels+["00"])
    
    title_string = str(directory_path)+" "+str(measurement)+" "+cohort.upper()+" "+region
    ax.set_title(title_string, fontsize=16)
    
    if savePlots:
        if len(region)>0:
            fig.savefig(directory_path+"/"+cohort.upper()+"_"+directory_path+"_"+region+"_heatmap.png", dpi=300)
            fig.savefig(directory_path+"/"+cohort.upper()+"_"+directory_path+"_"+regio+"_heatmap.svg", dpi=300)
        else:
            fig.savefig(directory_path+"/"+cohort.upper()+"_"+directory_path+"_heatmap.png", dpi=300)
            fig.savefig(directory_path+"/"+cohort.upper()+"_"+directory_path+"_heatmap.svg", dpi=300)

def get_bed_cols(data, cohort, prefix="avg"):
    """ returns the columns not created by modbamtools calcmeth """
    return [f for f in data.columns if not f.startswith(prefix) and len(f.split(cohort.upper()))==1]



def calc_var_std(df_vs, cohort, directory_path, field , high_var_cutoff=0.1):   
    
    # New index of gene 'ids' to add
    new_index = ['gene_'+str(i) for i in np.arange(df_vs.shape[0])]
    df_vs['gene'] = new_index

    df_vs = df_vs.set_index('gene', append=True)
    field='gene'
    
    transposed_df1 = df_vs.T

    # variance and standard deviation for each column
    variance_per_column = transposed_df1.var()
    std_deviation_per_column = transposed_df1.std()
    # Median line
    median_var = variance_per_column.median()
    median_std = std_deviation_per_column.median()
    v_std_df = pd.DataFrame({'variance':variance_per_column, 
                             'std_deviation':std_deviation_per_column},
                           index=variance_per_column.index)

    # Choose columns with the top 10% highest variance
    top_10_percent_columns = variance_per_column.nlargest(int(len(variance_per_column) * high_var_cutoff))
    print(variance_per_column.shape, int(len(variance_per_column) * high_var_cutoff))
    # print("Columns with the Top 10% Highest Variance:",top_10_percent_columns.index)

    chosen_columns = df_vs.loc[top_10_percent_columns.index]

    print(chosen_columns.var())
    # List of 'chrom' values to color
    chrom_colors = {
        'chr1': 'red',
        'chr2': 'green',
        'chr3': 'blue',
        'chr4': 'orange',
        'chr5': 'purple',
        'chr6': 'brown',
        'chr7': 'pink',
        'chr8': 'gray',
        'chr9': 'cyan',
        'chr10': 'yellow',
        'chr11': 'magenta',
        'chr12': 'olive',
        'chr13': 'lime',
        'chr14': 'teal',
        'chr15': 'blueviolet',
        'chr16': 'gold',
        'chr17': 'navy',
        'chr18': 'indigo',
        'chr19': 'darkred',
        'chr20': 'darkgreen',
        'chr21': 'darkblue',
        'chr22': 'darkorange',
        'chrX': 'darkcyan',
        'chrY': 'orangered',
        'chrM': 'black'
    }
    
    v_std_df['numeric_chrom'] = v_std_df.index.get_level_values('chrom').str.strip("chr")
    v_std_df['numeric_chrom'] = pd.to_numeric(v_std_df['numeric_chrom'], errors='coerce')
    v_std_df_sorted = v_std_df.sort_values(by='numeric_chrom')
    
    # Scatter plot with lines for median and mean
    fig, (ax_combined, ax_variance, ax_std_deviation) = plt.subplots(3, 1, figsize=(10, 12), sharex=False)

    # Scatter plot with points colored by 'chrom'
    for chrom, group in v_std_df_sorted.groupby(level='chrom'):
        ax_combined.scatter(group['variance'], group['std_deviation'], color=chrom_colors[chrom], alpha=0.3, label=f'{chrom})')
        
    ax_combined.axvline(median_var, color='red', linestyle='--', label=f'Median Variance: {median_var:.2f}')
    ax_combined.axhline(median_std, color='green', linestyle='--', label=f'Median Std Deviation: {median_std:.2f}')
    ax_combined.set_ylabel('Standard Deviation')
    ax_combined.legend(prop ={'size': 6})

    # Scatter plot for variance
    var_xtix = []
    var_xtix_labels = []
    last_group_end = 0
    for chrom, group in v_std_df_sorted.groupby(level='chrom'):
        ax_variance.scatter(group.index.get_level_values(field), group['variance'], color=chrom_colors[chrom], alpha=0.7, s=2)
        var_xtix.append( (last_group_end+ax_variance.get_xticks()[-2])/2)
        var_xtix_labels.append( chrom )
        last_group_end = ax_variance.get_xticks()[-1]
        
    ax_variance.axhline(median_var, color='red', linestyle='--', label=f'Median Variance: {median_var:.2f}')
    ax_variance.set_ylabel('Variance')
    ax_variance.set_xticks(var_xtix)
    ax_variance.set_xticklabels(var_xtix_labels,rotation=90)
    ax_variance.legend(prop ={'size': 6})

    # Scatter plot for standard deviation
    std_xtix = []
    std_xtix_labels = []
    last_group_end = 0
    for chrom, group in v_std_df_sorted.groupby(level='chrom'):
        ax_std_deviation.scatter(group.index.get_level_values(field), group['std_deviation'], color=chrom_colors[chrom], alpha=0.7, s=2)
        std_xtix.append( (last_group_end+ax_std_deviation.get_xticks()[-2])/2)
        std_xtix_labels.append(chrom)
        last_group_end = ax_std_deviation.get_xticks()[-1]
        
    ax_std_deviation.axhline(median_std, color='green', linestyle='--', label=f'Median Std Deviation: {median_std:.2f}')
    ax_std_deviation.set_xlabel('Genes')
    ax_std_deviation.set_ylabel('Standard Deviation')
    ax_std_deviation.set_xticks(std_xtix)
    ax_std_deviation.set_xticklabels(std_xtix_labels,rotation=90,)
    ax_std_deviation.legend()

    # Set title for the combined plot
    ax_combined.set_title('Combined Scatter Plot with Median and Mean Lines')

    plt.tight_layout()
    plt.show() 
    plt.savefig(directory_path+"/"+cohort.upper()+"_"+directory_path+"_variance_std.png",dpi=300)



def make_regional_bed_directories(directory_path, cohort):
    """ check for existance of regionally named directory """
    # verify existance 
    if not (os.path.exists(directory_path) and os.path.isdir(directory_path)):
        print('making', directory_path, 'directory for cohorts', cohort)
        !mkdir {directory_path}
        !mkdir {directory_path}/{cohort}
    else:
        # just print the contents of directory 
        print(directory_path, 'directory contents for cohort', cohort,"\n\t", ('\n\t').join(os.listdir(directory_path)) )
        print(cohort+"_"+directory_path+"_calcMeth.tsv")
        if os.path.exists(directory_path+"/"+cohort+"_"+directory_path+"_calcMeth.tsv"):
            return True
    return False

def write_out_dataframe(df, directory_path, cohort, delimiter="\t", suffix="_avgmod.bed", 
                        includeHeader = True, includeIndex = True):
    """ Make tsv's or beds from bedmethyls """
    output_path = directory_path+"/"+cohort+"_"+directory_path+suffix
    print('writing to :', output_path)
    df.to_csv(output_path, 
                           sep=delimiter,
                           header=includeHeader, 
                           index=includeIndex)

def write_coverage_cutoff_dataframe(cov_dfm, directory_path, cohort, bed_columns, cov_suffix="cov5"):
    """ takes the coverage dataframe and outputs 
    """

    write_out_dataframe(isolate_avgMod_cols(cov_dfm, cohort, bed_columns ), 
                        directory_path, cohort, suffix="_"+cov_suffix+".bed")
          
def remove_duplicate_columns(df):
    """
    Remove duplicate columns from the DataFrame.
    """
    transposed_df = df.transpose()

    # Drop duplicate rows 
    unique_transposed_df = transposed_df.drop_duplicates()
    
    # back to the original shape  
    result_df = unique_transposed_df.transpose()

    return result_df
        
def isolate_avgMod_cols(dfm1, cohort, regional_bed_columns, field="avgMod_" ):
    # isolate avgMod columns
    prefix = field + cohort.upper()
    cols = regional_bed_columns+[f for f in dfm1.columns 
                                 if f[:len(prefix)] == prefix ]
    
    cohort_dfm_avgmod = dfm1.loc[:,tuple(cols)]
    
    # move bed table information into the dataframe index
    if len(regional_bed_columns)>0:       
        # move bed columns to index
        print('move cols to index',regional_bed_columns)
        cohort_dfm_avgmod.set_index(regional_bed_columns, inplace=True)
        
        # change the datatype of all calls back to floats 
        all_cols = (f for f in cohort_dfm_avgmod.columns)
        cohort_dfm_avgmod.loc[:,all_cols]=cohort_dfm_avgmod.loc[:,all_cols].astype({element: float for element in all_cols})
    
    return cohort_dfm_avgmod

def clean_df_write_out_files(cohort_dfm, cohort, num_samples, regional_bed_columns, 
                             directory_path, extrastr="", write_files=True):
    """ Removes duplicated input columns in group dataframe. 
        Writes out a group bed file with desired regiona_bed_columns and avgMod_, 
            and avgMod_hap1 and avgMod_hap2 for each individual. 
        Writes out the group dataframe as a tsv. 
        
        Returns: A reduced bedmethyl group dataframe for individuals in the terra table
    """
    print('cohort DataFrame shape', cohort_dfm.shape, cohort_dfm.shape[1]/num_samples, 'columns per sample' )
    
    # remove the repeated bed columns
    cohort_dfm_no_dups = remove_duplicate_columns(cohort_dfm)

    print('after removing duplicate columns:\n', cohort_dfm_no_dups.shape, 
          (cohort_dfm_no_dups.shape[1])/num_samples, 'columns per sample' )
    
    if write_files:
        # write out the whole thing as a tsv
        write_out_dataframe(cohort_dfm_no_dups, directory_path, cohort+extrastr, suffix="_calcMeth.tsv", 
                            includeIndex = False)
        
    return cohort_dfm_no_dups
        
def localize(ttable_df, field, directory_path, cohort, read_func):
    """ Bring files to this workspace and load into a dataframe 
        Returns: A dictionary of bedmethyl dataframes for individuals in the terra table
    """
    
    # Create a dictionary to hold individual sample DataFrames 
    gBEDdict = {}
    # set the location to store files 
    directory = directory_path+'/'+cohort+'/'
    # loop through all the files for the field of interest
    for i, row in ttable_df.iterrows():
        # confirm that the file exists in this cell of the table
        b = row[field]
        if type(b)!=float and len(b)>0: 
            # isolate the file name from the gs link
            fname = b.split("/")[-1]

            # if the file alread exists store it in a sample dataframe
            if os.path.isfile(directory+fname):
                gBEDdict[str(i)] = read_func(directory+fname, i)
            else:
                print('localizing:', fname)
                # quietly copy and get read & asm data from shasta log, then delete
                !gsutil -q cp {b} {directory_path}/{cohort}/{fname}
                gBEDdict[str(i)] = read_func(directory+fname, i)

    print(cohort,'samples:',len(gBEDdict.keys()) )
    # wise to remove the cohort directories here
#     !rm -rf {directory_path}/{cohort}
    
    return gBEDdict

def coverage_and_avgMeth_distributions(dfm, cohort, directory_path):
    # how many genes are covered by at least x reads
    cov_shapes=[]

    for i in np.arange(-1,40):
        # set Coverage minmum 
        min_value = i
        field = 'totalCov_'+cohort.upper()

        # make query string
        field_filter = ['`'+str(f)+'`'+'>'+str(min_value) for f in dfm.columns if f[:len(field)] == field]
        filter_string = ' & '.join(field_filter)

        dfm_tcov = dfm.query(filter_string)
        cov_shapes.append(dfm_tcov.shape)
    
    # average methylation for each bed entry
    avgmod_shapes = []

    for i in np.arange(-1,100):
        # set Methylation minimum 
        min_value = i
        nfield = 'avgMod_'+cohort.upper()

        # make query string
        field_filter = ['`'+str(f)+'`'+'>='+str(min_value) for f in dfm.columns if f[:len(nfield)] == nfield]
        filter_string = ' & '.join(field_filter)
        dfm_tmod = dfm.query(filter_string)
        avgmod_shapes.append(dfm_tmod.shape)
    
    # maybe write this data out to a file also
    plot_cov_mod(cov_shapes, avgmod_shapes, directory_path, cohort)
    
    return cov_shapes, avgmod_shapes

def ratio_and_avgMeth_distributions(dfm, cohort, directory_path):
    # how many genes are covered by at least x reads
    cov_shapes=[]

    for i in np.arange(-1,1.25,0.25):
        # set Coverage minmum 
        min_value = i
        field = 'ratio_'+cohort.upper()

        # make query string
        field_filter = ['`'+str(f)+'`'+'<='+str(min_value) for f in dfm.columns if f[:len(field)] == field]
        filter_string = ' & '.join(field_filter)
        if len(filter_string) > 0:
            dfm_tcov = dfm.query(filter_string,engine='python')
            cov_shapes.append(dfm_tcov.shape)
        else:
            cov_shapes.append([0,0])
    
    # average methylation for each bed entry
    avgmod_shapes = []

    for i in np.arange(-1,100):
        # set Methylation minimum 
        min_value = i
        nfield = 'avgMod_'+cohort.upper()

        # make query string
        field_filter = ['`'+str(f)+'`'+'>='+str(min_value) for f in dfm.columns if f[:len(nfield)] == nfield]
        filter_string = ' & '.join(field_filter)
        if len(filter_string) > 0:
            dfm_tmod = dfm.query(filter_string,engine='python')
            avgmod_shapes.append(dfm_tmod.shape)
        else:
            avgmod_shapes.append(0)
        
    
    # maybe write this data out to a file also
    plot_cov_mod(cov_shapes, avgmod_shapes, directory_path, cohort)
    
    return cov_shapes, avgmod_shapes

def plot_cov_mod(cov_shapes, mod_shapes, directory_path, cohort):
    fig, ax = plt.subplots(1,2,figsize=(12,6))

    # Check coverage cutoff ratios for this region
    ax[0].plot(np.arange(-1,1.25,0.25),[s[0] for s in cov_shapes])
    ax[0].set_ylabel("Number "+directory_path)
    ax[0].set_xlabel("Ratio's for Region")
    ax[0].set_title("Number of "+directory_path+" by Filter Ratios "+cohort.upper() )
    ax[0].grid()
    ax[0].legend([cohort.upper()], prop ={'size': 10})

    # modbatools avg meth for promoters.
    ax[1].plot(np.arange(-1,(len(mod_shapes)-1)),[s[0] for s in mod_shapes])
    ax[1].set_ylabel("Number of "+directory_path)
    ax[1].set_xlabel("Average Mod for Region")
    ax[1].set_title("Number of "+directory_path+" by Average Modification "+cohort.upper())
#     ax[1].set_yscale("log")
    ax[1].grid()
    ax[1].legend([cohort.upper()], prop ={'size': 10})
    plt.tight_layout()

    fig.savefig(directory_path+"/"+cohort.upper()+"_cov_avgMod.png",dpi=250)
    plt.close(fig)
    
def coverage_kneePoint_derivitive_analysis(coverage_data, cohort, directory_path):
    # irst numbers from each tuple
    first_numbers = [item[0] for item in coverage_data]

    # derivative of the first numbers
    derivative = np.gradient(first_numbers)

    # index corresponding to the maximum absolute derivative
    knee_point_index = np.argmax(np.abs(derivative))
    
    # Plot 
    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.plot(first_numbers, label='Data')
    plt.scatter(knee_point_index, first_numbers[knee_point_index], color='red', label='Knee Point')
    plt.title(f"The knee point is at coverage {knee_point_index}, corresponding to the value {first_numbers[knee_point_index]}.")
    plt.xlabel('Index')
    plt.ylabel('First Number')
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(derivative, label='Derivative', color='green')
    plt.scatter(knee_point_index, derivative[knee_point_index], color='red', label='Knee Point')
    plt.title('Derivative of the Data and Knee Point')
    plt.xlabel('Coverage')
    plt.ylabel('Derivative')
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)  # Horizontal line at y=0
    plt.legend()

    plt.tight_layout()
    plt.savefig(cohort.upper()+"_"+directory_path+"_coverageCutoff_calculation.png", dpi=300)
    plt.close()
    
    return knee_point_index
    
def coverage_cutoffs(dfm, cohort, cov_min = 5, prefix="totalCov_"):
    """ limits to coverage of 5. 
        Input DataFrame must include totalCov field
        
    """
    field = prefix+cohort.upper()
    print('field:',field)
    # make query string for coverage cutoffs
    field_filter = ['`'+str(f)+'`'+'>'+str(cov_min) for f in dfm.columns if f[:len(field)] == field]
    filter_string = ' & '.join(field_filter)
    print('field_filter',field_filter)

    return dfm.query(filter_string)

def ratio_cutoffs(dfm, cohort, ratio_min = 0.9, prefix="ratio_"):
    """ limits to promoters that had coverage of 5 and a minimum of 10 cpgs in the region. 
        Input DataFrame must include ratio_ field
        
    """
    field = prefix+cohort.upper()
    print('field:',field)
    # make query string for coverage cutoffs
    field_filter = ['`'+str(f)+'`'+'>='+str(ratio_min) for f in dfm.columns if f[:len(field)] == field]

    filter_string = ' & '.join(field_filter)
    print(filter_string)


    return dfm.query(filter_string,engine='python')

def regional_sample_coverage(dfm, filter_threshold = 0.75):
    """ 
        Filter out regions that aren't covered by a minimum number of samples
        default 75% of samples in cohort. 
        Input dataframe must be only numeric values. 
    """
    
    
    # Check if the DataFrame contains numerical data
    if pd.api.types.is_numeric_dtype(dfm.dtypes):
        raise ValueError("Input DataFrame must contain numerical data.")
        
    # ensure that all datatypes are floats befor converting -1's to nan's
    all_cols = (f for f in dfm.columns)
    dfm.loc[:,all_cols] = dfm.loc[:,all_cols].astype({element: float for element in all_cols})

    # Set a threshold of the number of samples that need to cover a region to be includede
    threshhold = filter_threshold * dfm.shape[1]
    # Replace -1's with np.nan floats
    dfm.replace(-1, np.nan, inplace=True)
    # Filter rows based on threshold ratio
    filtered_dfm = dfm.dropna(thresh=threshhold)
    
    print('removed',dfm.shape[0]-filtered_dfm.shape[0], 'regions')
    
    return filtered_dfm
    
def plot_entries_per_chr(bed_df, directory_path, cohort):
    
    # Create a mapping between chromosome strings and integer values
    chromosome_mapping = {'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8,
                              'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
                              'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22,
                              'chrX': 23, 'chrY': 24, 'chrM':25}

    # Create a reverse mapping for labeling the y-axis
    reverse_chromosome_mapping = {v: k for k, v in chromosome_mapping.items()}
    
    # Calculate the number of entries for each chromosome
    entries_per_chromosome = bed_df.loc[bed_df['chrom'].str.len() <= 6]['chrom'].value_counts()

    # Plot a histogram using matplotlib
    plt.figure(figsize=(8, 3))

    plt.bar([chromosome_mapping[chrm] for chrm in entries_per_chromosome.index if chrm in chromosome_mapping.keys()], entries_per_chromosome, color='skyblue', edgecolor='black')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Entries')
    plt.title('Histogram of '+directory_path+' per Chromosome '+cohort)
    plt.xticks(list(reverse_chromosome_mapping.keys()), list(reverse_chromosome_mapping.values()))
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    plt.show()
    plt.savefig(directory_path+"/"+cohort.upper()+"_chr_entryCounts.png",dpi=250)


def analyze_regional_calcmeth_beds(cohort, ttable_df, directory_path, field, 
                                   read_func, regional_bed_columns, default_ratio_value = 0.9,
                                   autosomeHeatmaps=False, write_files=True, do_cov_work=True,
                                   variance_plots=True):
    """
    Remove duplicate columns from a DataFrame.
    """

    # Set up directories for localizing files
    merged_df_exists = make_regional_bed_directories(directory_path, cohort)
    print("tsv already exists?", merged_df_exists,'\n',directory_path, cohort)
    if not merged_df_exists:
        samples_dict = localize(ttable_df, field, directory_path, cohort, read_func)

        cohort_dfm = pd.concat(samples_dict.values(), axis=1)


        num_samples = ttable_df.shape[0]

        dfm = clean_df_write_out_files(cohort_dfm, 
                                                      cohort, 
                                                      num_samples, 
                                                      regional_bed_columns, 
                                                      directory_path, 
                                                      write_files=True)
    else:
        print('reading tsv')
        dfm = pd.read_csv(directory_path+"/"+cohort+"_"+directory_path+"_calcMeth.tsv", 
                          sep="\t", dtype="object")
        
    # ensure the ratios and avgMeths are floats    
    firstCol = len(regional_bed_columns)-1
    for col in dfm.columns[firstCol:]:
        dfm[col] = pd.to_numeric(dfm[col], errors='coerce')
        

    print('plot coverage and methyl dist')
    cov_shapes, avgmod_shapes = ratio_and_avgMeth_distributions(dfm, 
                                                                cohort,
                                                                directory_path)

    dfm_min_cov = ratio_cutoffs(dfm, cohort, 
                                ratio_min = default_ratio_value, prefix="ratio_")
    print('ratio cutoff',default_ratio_value,dfm_min_cov.shape )
    

    print('isolate avgMod cols and write df')
    cohort_dfm_min_cov_avgmod= isolate_avgMod_cols(dfm_min_cov, 
                                cohort, regional_bed_columns )
    write_out_dataframe(cohort_dfm_min_cov_avgmod, directory_path, cohort)
    
    

    plot_entries_per_chr(dfm_min_cov, directory_path, cohort)
    
    # plot heatmap of avgs after coverage coutoff
    plot_heatmap(regional_sample_coverage(cohort_dfm_min_cov_avgmod), cohort, 
                 directory_path, default_ratio_value)
    
    if do_cov_work:

        # Write out group beds after coverage cutoffs 
        print('write cov cutoff',default_min_cov)
        write_coverage_cutoff_dataframe(dfm_min_cov, directory_path, cohort, 
                                        regional_bed_columns, cov_suffix="cov5")

        # calculate second coverage cutoff
        second_cov_cutoff = coverage_kneePoint_derivitive_analysis(cov_shapes, cohort, directory_path)
        print('cov cutoff',second_cov_cutoff)
        dfm_min_cov = coverage_cutoffs(dfm, cohort, cov_min = second_cov_cutoff, prefix="totalCov_")
        cohort_dfm_min_cov_avgmod2= isolate_avgMod_cols(dfm_min_cov, cohort, get_bed_cols(dfm_min_cov, cohort) )
        # plot counts of bed entries per chr
        plot_entries_per_chr(dfm_min_cov, directory_path, cohort+"_"+str(second_cov_cutoff))
        
        # plot heatmap of avgs after coverage coutoff
        plot_heatmap(regional_sample_coverage(cohort_dfm_min_cov_avgmod2), cohort, 
                     directory_path, second_cov_cutoff)

        # Write out group beds after coverage cutoffs 
        print('write cov cutoff',second_cov_cutoff)
        write_coverage_cutoff_dataframe(dfm_min_cov, directory_path, cohort, 
                                        regional_bed_columns, cov_suffix="cov"+str(second_cov_cutoff))
    
    
        
    # pca, plot variance of each region across samples, correlation, methylation variation per gene

    if autosomeHeatmaps:
        #autosomes heatmaps
        chromosome_names = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']
        for cn in chromosome_names:
            if cn in cohort_dfm_min_cov_avgmod.index:
                plot_heatmap(cohort_dfm_min_cov_avgmod.loc[cn], 
                             cohort, 
                             directory_path,
                             default_min_cov,
                             region=cn,
                             xls="start")
    
    if variance_plots:
        print('regional_bed_columns',regional_bed_columns)
        # it was regional_bed_columns[3] to plot variance of each promoter, now I'm using the 'start col' [1]
        calc_var_std(cohort_dfm_min_cov_avgmod, cohort, directory_path, 
                     regional_bed_columns[1], high_var_cutoff=0.1)

    
    
    
#     return result_df




