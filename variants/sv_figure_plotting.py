%%capture
import terra_notebook_utils as tnu
import terra_pandas as tp
import os
import io
# import gzip
# import math
import pandas as pd
import numpy as np
# import tarfile
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pysam
import bgzip
import datetime
import seaborn as sns
from PIL import Image
from adjustText import adjust_text
import matplotlib.lines as mlines
from collections import defaultdict


Set up workspace
# Get the Google billing project name and workspace name
PROJECT = os.environ['WORKSPACE_NAMESPACE']
WORKSPACE =os.path.basename(os.path.dirname(os.getcwd()))
bucket = os.environ['WORKSPACE_BUCKET'] + "/"

# Verify that we've captured the environment variables
print("Billing project: " + PROJECT)
print("Workspace: " + WORKSPACE)
print("Workspace storage bucket: " + bucket)
Billing project: AUTH_ANVIL_CARD_LR_WGS_ANALYSI
Workspace: AnVIL_NIA_CARD_LR_WGS_Analysis
Workspace storage bucket: gs://fc-secure-578c07ca-f188-44a1-a254-b524a5e73ecf/
hapdiff_color = (0/256.0, 170/256.0, 231/256.0)
sniffles_color = (95/256.0, 51/256.0, 139/256.0)
NABEC VCFs and Truvari
Access the HBCC data table and move all the data to this workspace
# get the NABEC Table of links
nabec_cohort = tp.table_to_dataframe("NABEC_cohort", workspace=WORKSPACE, workspace_namespace=PROJECT)
hbcc_cohort = tp.table_to_dataframe("HBCC_Cohort", workspace=WORKSPACE, workspace_namespace=PROJECT)
# set up the bucket loction to store files 
base_gs_url="gs://fc-secure-578c07ca-f188-44a1-a254-b524a5e73ecf/NABEC_plots"
# make the samples a column in the df instead of just the index
nabec_cohort['sample_ID'] = nabec_cohort.index
hbcc_cohort['sample_ID'] = hbcc_cohort.index

#---------------------------
#SVS per individual per cohort
#---------------------------

# !rm nabec_harmonized_vcfs/nabec_hapdiff_harmony_smallerThan50bps_numRecords.txt
!rm nabec_harmonized_vcfs/nabec_hapdiff_harmony_biggerThan50bps_numRecords.txt
# how many variants per individual are there? 
!for f in nabec_harmonized_vcfs/NABEC*hapdiff_harmonized.vcf.gz.50.vcf.gz; do \
    bcftools stats -i 'INFO/SVLEN <= 50 && INFO/SVLEN >= -50' $f \
    | grep -w -e "^ID" -e "number of records:" - | awk '{print $NF}' - >> nabec_harmonized_vcfs/nabec_hapdiff_harmony_smallerThan50bps_numRecords.txt; done

    
!for f in nabec_harmonized_vcfs/NABEC*hapdiff_harmonized.vcf.gz.50.vcf.gz; do \
    bcftools stats -e 'INFO/SVLEN >= 50 && INFO/SVLEN <= -50' $f \
    | grep -w -e "^ID" -e "number of records:" - | awk '{print $NF}' - >> nabec_harmonized_vcfs/nabec_hapdiff_harmony_biggerThan50bps_numRecords.txt; done


!rm hbcc_harmonized_svs/hbcc_hapdiff_harmony_smallerThan50bps_numRecords.txt
!rm hbcc_harmonized_svs/hbcc_hapdiff_harmony_biggerThan50bps_numRecords.txt
# # # how many variants per individual are there? 

!for f in hbcc_harmonized_svs/HBCC*hapdiff_harmonized.vcf.gz.50.vcf.gz; do \
    bcftools stats -i 'INFO/SVLEN <= 50 && INFO/SVLEN >= -50' $f \
    | grep -w -e "^ID" -e "number of records:" - | awk '{print $NF}' - >> hbcc_harmonized_svs/hbcc_hapdiff_harmony_smallerThan50bps_numRecords.txt ; done


!for f in hbcc_harmonized_svs/HBCC*hapdiff_harmonized.vcf.gz.50.vcf.gz; do \
    bcftools stats -i 'INFO/SVLEN >= 50 | INFO/SVLEN <= -50' $f \
    | grep -w -e "^ID" -e "number of records:" - | awk '{print $NF}' - >> hbcc_harmonized_svs/hbcc_hapdiff_harmony_biggerThan50bps_numRecords.txt ; done


svs_per_indiv_NABEC_SNF('nabec_snf', 'NABEC')

svs_per_indiv_HBCC_SNF('hbcc_snf', 'HBCC')


!echo bigger
!head -3 sv_merges_all/nabec_harmonized_vcfs/nabec_hapdiff_harmony_biggerThan50bps_numRecords.txt
# !head -3 hbcc_harmonized_svs/hbcc_hapdiff_harmony_biggerThan50bps_numRecords.txt

!echo snf
!head -3 nabec_snf/sv_num_records_file_NABEC.txt 


# bigger
# nabec_harmonized_vcfs/NABEC_KEN-1066_FTX_hapdiff_harmonized.vcf.gz.50.vcf.gz
# 22148
# nabec_harmonized_vcfs/NABEC_KEN-1069_FTX_hapdiff_harmonized.vcf.gz.50.vcf.gz
# snf
# nabec_snf/NABEC_KEN-1066_FTX.sniffles.vcf.gz
# 23664
# nabec_snf/NABEC_KEN-1069_FTX.sniffles.vcf.gz    

nabec_color = '#1f77b4'
hbcc_color = '#ff7f0e'

# Read sniffles num svs file
file_path = 'nabec_snf/sv_num_records_file_NABEC.txt'

linenum = 0
sample = []
snf_num_svs = []
with open(file_path, 'r') as file:
    for lines in file:
#     lines = file.readline() #s()
        if linenum%2==0:
            sample.append( ('_').join(lines.strip().split('/')[1].split('_')[0:2]) )

        else:
            snf_num_svs.append( int(lines.strip()))
        linenum+=1

# Truvari Harmonized SVs > 50bps
file_path = 'sv_merges_all/nabec_harmonized_vcfs/nabec_hapdiff_harmony_biggerThan50bps_numRecords.txt'
linenum = 0
large_svs = []
big_samples = []
with open(file_path, 'r') as file:
    for lines in file:
#     lines = file.readline() #s()
        if lines.startswith('nabec'):  #linenum%2==0:
            sname = ('_').join(lines.strip().split('/')[1].split('_')[0:2])
            
            if sname in big_samples:
                break
            else:
                big_samples.append( ('_').join(lines.strip().split('/')[1].split('_')[0:2]) )
        else:
            large_svs.append( int(lines.strip()))
        linenum+=1




# Create a DataFrame of NABEC hapdiff 
data = {'sample': big_samples, 'hapdiff_num_svs':large_svs}
hapdup_nabec_df = pd.DataFrame(data)

# remove These samples: 
# 143  NABEC_UMARY-1859        0   # there are no hapdiff SVs.. ?
# 176  NABEC_UMARY-4915    16349   # This sample has a low quality assembly 
hapdup_nabec_df.drop(index=143, inplace=True)
hapdup_nabec_df.drop(index=176, inplace=True)

# Display the DataFrame
print(hapdup_nabec_df.sort_values(by="hapdiff_num_svs").tail(15))

# Create a DataFrame of NABEC sniffles 
snf_data = {'sample':sample, 'snf_num_svs':snf_num_svs}
nabec_snf_df = pd.DataFrame(snf_data)

print(nabec_snf_df.sort_values(by="snf_num_svs").tail(5))


"""
               sample  hapdiff_num_svs
117  NABEC_UMARY-1543            22673
193   NABEC_UMARY-544            22689
20     NABEC_KEN-5015            22689
110  NABEC_UMARY-1486            22697
171  NABEC_UMARY-4786            22698
95   NABEC_UMARY-1326            22701
158  NABEC_UMARY-4549            22720
65     NABEC_SH-97-53            22741
173  NABEC_UMARY-4842            22751
30     NABEC_SH-02-08            22760
170  NABEC_UMARY-4782            22773
103  NABEC_UMARY-1442            22796
162  NABEC_UMARY-4640            22819
106  NABEC_UMARY-1461            22921
126  NABEC_UMARY-1672            23410
               sample  snf_num_svs
47     NABEC_SH-08-04        24324
126  NABEC_UMARY-1672        24411
41     NABEC_SH-06-25        24435
188  NABEC_UMARY-5116        24450
51     NABEC_SH-94-35        27417
"""

#------------------------
# HBCC individual SVs 
# Read sniffles num svs file
file_path = 'hbcc_snf/sv_num_records_file_HBCC.txt'

linenum = 0
sample = []
snf_num_svs = []
with open(file_path, 'r') as file:
    for lines in file:
#     lines = file.readline() #s()
        if linenum%2==0:
            sample.append( ('_').join(lines.strip().split('/')[1].split('_')[0:2]) )

        else:
            snf_num_svs.append( int(lines.strip()))
        linenum+=1

file_path = 'hbcc_harmonized_svs/hbcc_hapdiff_harmony_biggerThan50bps_numRecords.txt'
hb_large_svs = []
hb_big_samples = []
with open(file_path, 'r') as file:
    for lines in file:
#     lines = file.readline() #s()
        if lines.startswith('hbcc'):  #linenum%2==0:
            sname = ('_').join(lines.strip().split('/')[1].split('_')[0:2])
            
            if sname in hb_big_samples:
                break
            else:
                hb_big_samples.append( ('_').join(lines.strip().split('/')[1].split('_')[0:2]) )
        else:
            hb_large_svs.append( int(lines.strip()))
        linenum+=1

# Create a DataFrame
data = {'sample': hb_big_samples, 'num_svs':hb_large_svs}
hapdup_hbcc_df = pd.DataFrame(data)
# hdf.drop(index=143, inplace=True)

# Display the DataFrame
print(hapdup_hbcc_df.sort_values(by="num_svs").tail(5))

snf_data = {'sample':sample, 'snf_num_svs':snf_num_svs}
hbcc_snf_df = pd.DataFrame(snf_data)

print(hbcc_snf_df.sort_values(by="snf_num_svs").tail(5))

#          sample  num_svs
# 91   HBCC_82021    26380
# 52   HBCC_81981    26392
# 42   HBCC_81971    26408
# 118  HBCC_82048    26448
# 84   HBCC_82014    26494
#          sample  snf_num_svs
# 7    HBCC_81933        28974
# 144  HBCC_82075        28996
# 39   HBCC_81968        29106
# 58   HBCC_81987        29262
# 68   HBCC_81997        29284





#------------------------
# Svs per indiv

nabec_color = '#1f77b4'
hbcc_color = '#ff7f0e'

# Violin Swarm plot 
def violin_swarm(x,y,data,ax,cohort_color,swarm_pt_size = 3):
    sns.violinplot(x=x, y=y, data=data, cut=0.25,
               inner=None, color=cohort_color, alpha = 0.4, ax=ax, edgecolor='black')
    sns.swarmplot(x=x, y=y, data=data, color=cohort_color, 
              s=swarm_pt_size, alpha=1, ax=ax)

def violin_median(x,y,data,ax,cohort_color):
    sns.violinplot(x=x, y=y, data=data, color=cohort_color, cut=0.25, ax=ax)
    medianValue = np.median(data[y])
    ax.scatter(x[0], medianValue, color='red', zorder=3) #, label='Median')

# Create the violin plot
fig = plt.figure(figsize=(12, 12))
gs = fig.add_gridspec(1,1, height_ratios=[1])



ax1 = fig.add_subplot(gs[0, 0])

""" Hapdiff SVs """
violin_swarm([1]*hapdup_nabec_df.shape[0], 
             'hapdiff_num_svs', hapdup_nabec_df, ax1, nabec_color)

violin_swarm([2]*hapdup_hbcc_df.shape[0], 
             'num_svs', hapdup_hbcc_df, ax1, hbcc_color)

""" Sniffles2.2 SVs """
violin_swarm([3]*nabec_snf_df.shape[0], 
             'snf_num_svs', nabec_snf_df, ax1, nabec_color)

violin_swarm([4]*hbcc_snf_df.shape[0], 
             'snf_num_svs', hbcc_snf_df, ax1, hbcc_color)


# Add title and labels
ax1.set_title('SVs Per Individual (>50 bps)')
ax1.set_ylabel('Number of SVs')
ax1.set_xticklabels(['Hapdiff', 'Hapdiff', 'Sniffles','Sniffles'])

# # Create custom lines for custom legend
custom_line1 = mlines.Line2D([], [], marker='o', color=nabec_color, label='NABEC')
custom_line2 = mlines.Line2D([], [], marker='o', color=hbcc_color, label='HBCC')
ax1.legend(handles=[custom_line1, custom_line2], loc='upper left')




#------------------------
# SVs per cohort Merge
#------------------------
""" so Hapdiff has more SV calls in the merged file across cohorts than Sniffles2.3
Maybe this is because there are less similar SVs to merge in Hapdiff. 
"""

""" Combined cohort VCFs """
vcf_snf_nabec = 'sv_merges_all/nabec_directory/nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz'
vcf_hapdiff_nabec = 'sv_merges_all/nabec_directory/nabec_harmonized_truvari_merge_genotyped_pct75_07022024.sorted.norm.ordered.vcf.gz'

""" Load in the VCF """
print('hapdiff')
# # Open the VCF file using pysam
hapdiff_combined_vcf_in = pysam.VariantFile(vcf_hapdiff_nabec) 

# Extract SV lengths
nabec_hapdiff_sv_lengths = []
nabec_hapdiff_small_sv_lengths = []

nabec_hapdiff_numCollapsed = []
nabec_hapdiff_numConsolidated = []

nabec_hapdiff_numRecords =0
nabec_hapdiff_svTypes = {}

for record in hapdiff_combined_vcf_in:
    nabec_hapdiff_numRecords+=1
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
#         for key in record.info:
#             print(key,record.info[key])
        svlen = record.info['SVLEN']
        # SVLEN can be a list, so take the first element if it's a list
        if isinstance(svlen, list):
            svlen = int(svlen[0])
        if svlen<50 and svlen>-50:
            nabec_hapdiff_small_sv_lengths.append((svlen))
        else:
            nabec_hapdiff_sv_lengths.append((svlen))
        
            svtype = record.info['SVTYPE']

            if svtype not in nabec_hapdiff_svTypes.keys():
                print(svtype)
                nabec_hapdiff_svTypes[svtype]=1
            else:
                nabec_hapdiff_svTypes[svtype]+=1
    
    if ('NumCollapsed' in record.info) and (len(record.chrom)<6):
        nabec_hapdiff_numCollapsed.append(record.info['NumCollapsed'])
    if ('NumConsolidated' in record.info) and (len(record.chrom)<6):
        nabec_hapdiff_numConsolidated.append(record.info['NumConsolidated'])



# DataFrame of SV lengths
# hapdiff_combined_df = pd.DataFrame([small_sv_lengths, sv_lengths], columns=['smallSVs','SVLength'])

""" sniffles """
print('sniffles')
snf_vcf_in = pysam.VariantFile(vcf_snf_nabec)

nabec_sniffles_sv_lengths = []
nabec_sniffles_small_sv_lengths = []

nabec_sniffles_numRecords =0
nabec_sniffles_svTypes = {}
for record in snf_vcf_in:
    nabec_sniffles_numRecords+=1
    if ('SVLEN' in record.info) and (len(record.chrom)<6):

        svlen = record.info['SVLEN']
        # SVLEN can be a list, so take the first element if it's a list
        if isinstance(svlen, list):
            svlen = int(svlen[0])
        if svlen<50 and svlen>-50:
            nabec_sniffles_small_sv_lengths.append((svlen))
        else:
            nabec_sniffles_sv_lengths.append((svlen))
        
            svtype = record.info['SVTYPE']

            if svtype not in nabec_sniffles_svTypes.keys():
                print(svtype)
                nabec_sniffles_svTypes[svtype]=1
            else:
                nabec_sniffles_svTypes[svtype]+=1




#------------------------
# plot SV lengths Hapdiff
#------------------------
""" Individual Sample Hapdiff VCFs"""
# vcf_file = "nabec_harmonized_vcfs/NABEC_KEN-1070_FTX_hapdiff_harmonized.vcf.gz.50.vcf.gz"
# vcf_file = "hbcc_harmonized_svs/HBCC_81930_FTX_hapdiff_harmonized.vcf.gz.50.vcf.gz"
# vcf_file = "nabec_r9_r10_hbcc_named_truvari_genotyped4_merge_pct95_05242024.vcf.gz"

""" Hapdiff Cohort VCFs"""
# vcf_file = 'nabec_truvari_merge_pct95_03212024_biggerthan50bps.sorted.vcf.gz'
# vcf_file = 'nabec_truvari_merge_pct95_03212024_biggerthan50bps.sorted_noBlacklist_06012024.vcf.gz'
# black_hbcc_vcf_file = "hbcc_harmonized_truvari_merge_pct95_04012024_biggerthan50bps.sorted_noBlacklist_04022024.vcf.gz"

# nabec_hapdiff_black_vcf_file = "nabec_truvari_genotype_merge_pct95_03212024_biggerthan50bps.vcf_noBlacklist_03282024.vcf.gz"
# hbcc_hapdiff_black_vcf_file = 'hbcc_truvari_genotype_merge_pct95_01022024_biggerthan50bps_03212024.vcf_noBlacklist_03292024.vcf.gz'

vcf_hapdiff_nabec = 'sv_merges_all/nabec_directory/nabec_harmonized_truvari_merge_genotyped_pct75_07022024.sorted.norm.ordered.vcf.gz'
vcf_hapdiff_hbcc = 'sv_merges_all/hbcc_directory/hbcc_harmonized_truvari_merge_genotyped_pct75_07012024.sorted.sampleFiltered.vcf.gz'

""" Sniffles2.3 cohort VCFs """
vcf_snf_nabec = 'sv_merges_all/nabec_directory/nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz'
vcf_snf_hbcc  = 'sv_merges_all/hbcc_directory/hbcc_snf2.3_multisample_07012024.sampleFiltered.vcf.gz'

# vcf_file = 'NABEC_snifles2_3_multisample_biggerthan50bps.sorted_noBlacklist_03252024.vcf.gz'
# vcf_file = 'HBCC_snifles2_3_multisample_biggerthan50bps.sorted.vcf.gz'


""" Load in the VCFs """
# NABEC VCF files 
nabec_hapdiff_vcf_in = pysam.VariantFile(vcf_hapdiff_nabec)


# Extract SV lengths
nabec_hapdiff_sv_lengths = []
for record in nabec_hapdiff_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        # SVLEN can be a list, so take the first element if it's a list
        if isinstance(svlen, list):
            svlen = svlen[0]
        nabec_hapdiff_sv_lengths.append((svlen))

# HBCC Hapdiff VCF file 
hbcc_hapdiff_vcf_in = pysam.VariantFile(vcf_hapdiff_hbcc) 


# Extract SV lengths
hbcc_hapdiff_sv_lengths = []
for record in hbcc_hapdiff_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        # SVLEN can be a list, so take the first element if it's a list
        if isinstance(svlen, list):
            svlen = svlen[0]
        hbcc_hapdiff_sv_lengths.append((svlen))
        

# DataFrame of SV lengths
nabec_hapdiff_sv_lengths_df = pd.DataFrame(nabec_hapdiff_sv_lengths, columns=['SV Length'])
hbcc_hapdiff_sv_lengths_df = pd.DataFrame(hbcc_hapdiff_sv_lengths, columns=['SV Length'])


# plot the SV lengths 
fig, axs = plt.subplots(2, 1, figsize=(6,4))

nabec_hpd_b,nabec_hpd_l,nabec_hpd_h = axs[0].hist(
                nabec_hapdiff_sv_lengths_df.loc[(nabec_hapdiff_sv_lengths_df['SV Length']<7000) 
                & (nabec_hapdiff_sv_lengths_df['SV Length']>-7000)]['SV Length'], bins=14000, color=nabec_color)
hbcc_hpd_b,hbcc_hpd_l,hbcc_hpd_h = axs[1].hist(
                hbcc_hapdiff_sv_lengths_df.loc[(hbcc_hapdiff_sv_lengths_df['SV Length']<7000) 
                & (hbcc_hapdiff_sv_lengths_df['SV Length']>-7000)]['SV Length'], bins=14000, color=hbcc_color)


### 
axs[0].set_title('Histogram of NABEC Hapdiff Merged (75%) SV Lengths')
axs[1].set_title('Histogram of HBCC Hapdiff Merged (75%) SV Lengths')
for ax in axs:
    ax.set_xlabel('Length (bp)')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
plt.tight_layout()
plt.show()


""" plot the lengths of all SVs in the merged VCFs 75% """
plt.figure(figsize=(12, 6))
# Deletions

plt.scatter(nabec_hpd_l[:6950],nabec_hpd_b[:6950],s=2, label = 'NABEC', color=nabec_color)
plt.scatter(hbcc_hpd_l[:6950],hbcc_hpd_b[:6950],s=2, label = 'HBCC', alpha=0.45, color=hbcc_color)

# Insertions
plt.scatter(nabec_hpd_l[7050:-1],nabec_hpd_b[7050:],s=2, color=nabec_color)
plt.scatter(hbcc_hpd_l[7050:-1],hbcc_hpd_b[7050:],s=2, alpha=0.45, color=hbcc_color)

positions=[-700,-300,300,700] # -50,50,,6000
# Add vertical lines and text labels at specified positions
x1,x2,y1,y2 = plt.axis()  
plt.axis((-1000,1000,y1,y2))

ymax = plt.axis()[-1]  # Get the max y-value of the current plot
for pos in positions:
    plt.axvline(x=pos, color='grey', linestyle='--', linewidth=1.5, alpha=0.4)
    if abs(pos)!=50:
        if abs(pos)==300:
            plt.text(pos, ymax * 0.65, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)

        else:
            plt.text(pos, ymax * 0.15, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)

# plt.axes()

# plt.ylim(1,500)
# plt.yscale('log')
# plt.xscale('log')
plt.legend()
# Plot legend.
# lgnd = plt.legend(loc="upper left", numpoints=1, fontsize=10)
lgnd = plt.legend(loc="upper left", scatterpoints=1, fontsize=10)

for handle in lgnd.legendHandles:
    handle.set_sizes([26.0])
    
plt.title('Histogram of Hapdiff Merged Structural Variant Lengths (75%)')
plt.xlabel('SV Length (bp)')
plt.ylabel('Frequency')





#------------------------
# plot SV lengths Hapdiff 2 panel size
#------------------------
"""
2 panel of SV lengths 1Kb and 7KB 
With vertical lines at 300(SINE), 700(?), 2000 (SVA), 6000(LINE) bps

Plots both NABEC and HBCC HapDiff Truvari SV Length dfs

06262024

"""
def plot_hapdiff_sv_size_2panels(axs):
    sizes = [1000,7000]
    size_positions = [[-300,300],
                 [-6000,-2000,-300,300,2000,6000]]
    for i,ax in enumerate(axs):

        # Deletions 
        ax.scatter(nabec_hpd_l[:6950],nabec_hpd_b[:6950],s=2, label = 'NABEC', color=nabec_color)  
        ax.scatter(hbcc_hpd_l[:6950],hbcc_hpd_b[:6950],s=2, label = 'HBCC', alpha=0.45, color=hbcc_color)

        # Insertions 
        ax.scatter(nabec_hpd_l[7050:-1],nabec_hpd_b[7050:],s=2, color=nabec_color)
        ax.scatter(hbcc_hpd_l[7050:-1],hbcc_hpd_b[7050:],s=2, alpha=0.45, color=hbcc_color)

        if i>-1:
            ax.set_yscale('log')
        positions= size_positions[i]
        # Add vertical lines and text labels at specified positions
        x1,x2,y1,y2 = ax.axis() 
        x1 = (-1*sizes[i])
        x2 = sizes[i]
        ax.axis((x1,x2,y1,y2))

        ymax = ax.axis()[-1]  # Get the max y-value of the current plot
        for pos in positions:
            ax.axvline(x=pos, color='grey', linestyle='--', linewidth=1.5, alpha=0.15)
            if abs(pos)!=50:
                if abs(pos)==300:
                    ax.text(pos, ymax * 0.35, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)

                else:
                    ax.text(pos, ymax * 0.05, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)
        
        if i==0:
            ax.legend()
            # Plot legend.
            lgnd = ax.legend(loc="upper left", scatterpoints=1, fontsize=10)

            for handle in lgnd.legendHandles:
                handle.set_sizes([26.0])

        ax.set_title('Histogram of HapDiff+Truvari Merged (75%) Structural Variant Lengths')
        ax.set_xlabel('SV Length (bp)')
        ax.set_ylabel('Frequency')
        
fig, axes = plt.subplots(2, 1, figsize=(16, 14))
plot_hapdiff_sv_size_2panels(axes)



#------------------------
#New SNF Plot
#------------------------

""" Sniffles2.3 cohort VCFs """
# snf_vcf_file_nabec = 'NABEC_snifles2_3_multisample_biggerthan50bps.sorted_noBlacklist_03252024.vcf.gz'
# snf_vcf_file_hbcc = 'HBCC/HBCC_phase22_snifles2_3_multisample_04042024_biggerthan50bps.sorted_noBlacklist_04042024.vcf.gz'
snf_vcf_file_nabec = 'sv_merges_all/nabec_directory/nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz'
snf_vcf_file_hbcc  = 'sv_merges_all/hbcc_directory/hbcc_snf2.3_multisample_07012024.sampleFiltered.vcf.gz'


""" Combined cohort VCFs """
# vcf_snf_combined = 'NABEC_HBCCphase22_snifles2_3_multisample_04082024_biggerthan50bps.sorted_noBlacklist_04092024.vcf.gz'
# vcf_hapdiff_combined = 'nabec_hbccHarmony_truvari_genotype_merge_pct95_04222024_biggerthan50bps.sorted_noBlacklist_04252024.vcf.gz'

""" Load in the VCFs """
# NABEC Sniffles2.2 VCF file 
nabec_snf_black_vcf_in = pysam.VariantFile(snf_vcf_file_nabec) 

# Extract SV lengths
nabec_snf_black_sv_lengths = []
for record in nabec_snf_black_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        if isinstance(svlen, list):
            svlen = svlen[0]
        nabec_snf_black_sv_lengths.append((svlen))

# HBCC Sniffles2.2 VCF file 
hbcc_snf_black_vcf_in = pysam.VariantFile(snf_vcf_file_hbcc) 

# Extract SV lengths
hbcc_snf_black_sv_lengths = []
for record in hbcc_snf_black_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        if isinstance(svlen, list):
            svlen = svlen[0]
        hbcc_snf_black_sv_lengths.append((svlen))
        

# DataFrame of SV lengths
nabec_snf_black_sv_lengths_df = pd.DataFrame(nabec_snf_black_sv_lengths, columns=['SV Length'])
hbcc_snf_black_sv_lengths_df = pd.DataFrame(hbcc_snf_black_sv_lengths, columns=['SV Length'])



# plot the SV lengths 
""" Sniffles 2.2 / 2.3"""
fig, axs = plt.subplots(2, 1, figsize=(6,4))

nabec_snf_b,nabec_snf_l,nabec_snf_h = axs[0].hist(
                nabec_snf_black_sv_lengths_df.loc[(nabec_snf_black_sv_lengths_df['SV Length']<7000) 
                & (nabec_snf_black_sv_lengths_df['SV Length']>-7000)]['SV Length'], bins=14000, color=nabec_color)
hbcc_snf_b,hbcc_snf_l,hbcc_snf_h = axs[1].hist(
                hbcc_snf_black_sv_lengths_df.loc[(hbcc_snf_black_sv_lengths_df['SV Length']<7000) 
                & (hbcc_snf_black_sv_lengths_df['SV Length']>-7000)]['SV Length'], bins=14000, color=hbcc_color)


axs[0].set_title('Histogram of NABEC Sniffles Merged Structural Variant Lengths')
axs[1].set_title('Histogram of HBCC Sniffles Merged Structural Variant Lengths')
for ax in axs:
    ax.set_xlabel('Length (bp)')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
plt.tight_layout()
plt.show()


"""
2 panel of SV lengths 1Kb and 7KB 
With vertical lines at 300(SINE), 700(?), 2000 (SVA), 6000(LINE) bps

Plots both NABEC and HBCC Sniffles SV Length dfs

062620204

"""

def plot_snf_sv_size_2panels(axs):
    sizes = [1000,7000]
    size_positions = [[-300,300],
                 [-6000,-2000,-300,300,2000,6000]]
    for i,ax in enumerate(axs):

        # Deletions
        ax.scatter(nabec_snf_l[:6950],nabec_snf_b[:6950],s=2, label = 'NABEC', color=nabec_color)
        ax.scatter(hbcc_snf_l[:6950],hbcc_snf_b[:6950],s=2, label = 'HBCC', alpha=0.45, color=hbcc_color)


        # Insertions nabec_hpd_l[7050:-1]
        ax.scatter(nabec_snf_l[7050:-1],nabec_snf_b[7050:],s=2, color=nabec_color)
        ax.scatter(hbcc_snf_l[7050:-1],hbcc_snf_b[7050:],s=2, alpha=0.45, color=hbcc_color)

        if i>-1:
            ax.set_yscale('log')
        positions= size_positions[i]
        # Add vertical lines and text labels at specified positions
        x1,x2,y1,y2 = ax.axis() 
        x1 = (-1*sizes[i])
        x2 = sizes[i]
        ax.axis((x1,x2,y1,y2))

        ymax = ax.axis()[-1]  # Get the max y-value of the current plot
        for pos in positions:
            ax.axvline(x=pos, color='grey', linestyle='--', linewidth=1.5, alpha=0.15)
            if abs(pos)!=50:
                if abs(pos)==300:
                    ax.text(pos, ymax * 0.35, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)

                else:
                    ax.text(pos, ymax * 0.05, f'{pos}', color='black', ha='left', va='bottom', fontsize=10, rotation=90)

        if i==0:
            ax.legend()
            lgnd = ax.legend(loc="upper left", scatterpoints=1, fontsize=10)

            for handle in lgnd.legendHandles:
                handle.set_sizes([26.0])

        ax.set_title('Histogram of Sniffles2.3 Merged Structural Variant Lengths')
        ax.set_xlabel('SV Length (bp)')
        ax.set_ylabel('Frequency')
        
# from matplotlib.lines import Line2D
fig, axes = plt.subplots(2, 1, figsize=(16, 14))
plot_snf_sv_size_2panels(axes)




""" Combined SV caller cohort VCFs """
# vcf_snf_combined = 'NABEC_HBCCphase22_snifles2_3_multisample_04082024_biggerthan50bps.sorted_noBlacklist_04092024.vcf.gz'
# vcf_hapdiff_combined = 'nabec_hbccHarmony_truvari_genotype_merge_pct95_04222024_biggerthan50bps.sorted_noBlacklist_04252024.vcf.gz'
nabec_combined = "sv_merges_all/nabec_directory/nabec_harmonized_pct75_snf2.3_merge_pct75_07182024.symbolicQualResolved.keepMaxQual.sorted.vcf.gz"
hbcc_combined = "sv_merges_all/hbccc_directory/hbcc_harmonized_pct75_snf2.3_merge_pct75_07182024.symbolicQualResolved.keepMaxQual.vcf"

################## todo ####################
""" Load in the VCFs """
# NABEC Sniffles2.2 VCF file 
nabec_snf_black_vcf_in = pysam.VariantFile(snf_vcf_file_nabec) 

# Extract SV lengths
nabec_snf_black_sv_lengths = []
for record in nabec_snf_black_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        if isinstance(svlen, list):
            svlen = svlen[0]
        nabec_snf_black_sv_lengths.append((svlen))

# HBCC Sniffles2.2 VCF file 
hbcc_snf_black_vcf_in = pysam.VariantFile(snf_vcf_file_hbcc) 

# Extract SV lengths
hbcc_snf_black_sv_lengths = []
for record in hbcc_snf_black_vcf_in:
    if ('SVLEN' in record.info) and (len(record.chrom)<6):
        svlen = record.info['SVLEN']
        if isinstance(svlen, list):
            svlen = svlen[0]
        hbcc_snf_black_sv_lengths.append((svlen))
        

# DataFrame of SV lengths
nabec_snf_black_sv_lengths_df = pd.DataFrame(nabec_snf_black_sv_lengths, columns=['SV Length'])
hbcc_snf_black_sv_lengths_df = pd.DataFrame(hbcc_snf_black_sv_lengths, columns=['SV Length'])











