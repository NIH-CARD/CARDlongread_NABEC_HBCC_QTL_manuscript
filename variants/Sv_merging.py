

# !python3 -m pip install truvari==4.2.1

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pysam
from pysam import VariantFile
import bgzip
from adjustText import adjust_text
import matplotlib.lines as mlines
from collections import defaultdict

"""
Localize NABEC and HBCC Files
Localize Hapdiff Harmonized SV Files ( >= 50bps )
NABEC Harmoinzed Hapdiff VCFs
"""

# NABEC copy all the vcf files to this working directory, rename the files with their sample id
# store all the file names for bcftools merge
datatable_field = "harmonized_GRCh38_vcf"
nabec_directory = "sv_merges_all/nabec_harmonized_vcfs"

bad_qc_samples = []

nabec_vcfFiles = []
df = nabec_cohort.loc[:,[datatable_field,'sample_ID']]

for index, row in df.iterrows():
    
    if index in bad_qc_samples:
        print(index, 'skipping')
        continue
    
    # filepath and file renamed with sample id
    svVcf = row[datatable_field]
    outname = str(index)+"_hapdiff_harmonized.vcf.gz"
    
    # if the file for this sample is already in the directory 
    if os.path.isfile(nabec_directory+"/"+outname+".50.vcf.gz"):
        nabec_vcfFiles.append(nabec_directory+"/"+outname+".50.vcf.gz")
    elif type(svVcf)!=float and len(svVcf)>0:
        print(index)
        !gsutil -q cp {svVcf} {nabec_directory}/{outname}
        !bcftools view -i 'INFO/SVLEN >= 50 | INFO/SVLEN <= -50' {nabec_directory}/{outname} | bgzip  > {nabec_directory}/{outname}.50.vcf.gz
        !tabix {nabec_directory}/{outname}.50.vcf.gz
        nabec_vcfFiles.append(nabec_directory+"/"+outname+".50.vcf.gz")
        !rm {nabec_directory}/{outname}
    else:
        print('missing',index)

print(len(nabec_vcfFiles), 'samples')

# HBCC copy all the vcf files to this working directory, rename the files with their sample id
# store all the file names for bcftools merge
datatable_field = "harmonized_GRCh38_vcf"
hbcc_directory = "sv_merges_all/hbcc_harmonized_vcfs"
hbcc_vcfFiles = []
df = hbcc_cohort.loc[:,[datatable_field,'sample_ID']]

for index, row in df.iterrows():
    
    # filepath and file renamed with sample id
    svVcf = row[datatable_field]
    outname = str(index)+"_hapdiff_harmonized.vcf.gz"
    
    # if the file for this sample is already in the directory 
    if os.path.isfile(hbcc_directory+"/"+outname+".50.vcf.gz"):
        hbcc_vcfFiles.append(hbcc_directory+"/"+outname+".50.vcf.gz")
        if not os.path.isfile(hbcc_directory+"/"+outname+".50.vcf.gz.tbi"):
            print('tabix', hbcc_directory+"/"+outname+".50.vcf.gz")
            !tabix {hbcc_directory}/{outname}.50.vcf.gz
        
            
    elif type(svVcf)!=float and len(svVcf)>0:
        print(index)
        # localize harmonized SNV and SV VCF and select 
        !gsutil -q cp {svVcf} {hbcc_directory}/{outname}      
        !bcftools view -i 'INFO/SVLEN >= 50 | INFO/SVLEN <= -50' {hbcc_directory}/{outname} | bgzip  > {hbcc_directory}/{outname}.50.vcf.gz
        !tabix {hbcc_directory}/{outname}.50.vcf.gz
        hbcc_vcfFiles.append(hbcc_directory+"/"+outname+".50.vcf.gz")
        !rm {hbcc_directory}/{outname}
    else:
        print('missing',index)

print(len(hbcc_vcfFiles), 'samples')


"""
Localize Sniffles2.2 Phase SV Files
NABEC Sniffles VCFs and SNF Files
"""

# NABEC copy all the vcf files to this working directory, rename the files with their sample id
# store all the file names for bcftools merge
nabec_snf_vcfFiles = []
nabec_snfFiles = []
snf_field = 'snifflesSnf_22_GRCh38'
vcf_field = 'snifflesVcf_22_GRCh38'
nabec_snf = "nabec_snf"
df = nabec_cohort.loc[:,[snf_field, vcf_field,'sample_ID']]

bad_qc_sample = ['NABEC_UMARY-4915_FTX']

for index, row in df.iterrows():
    svVcf = row[vcf_field]
    svSnf = row[snf_field]
    vcfOutname = str(index)+".sniffles.vcf"
    snfOutname = str(index)+".snf"
    
    if os.path.isfile(nabec_snf+"/"+snfOutname):
        nabec_snfFiles.append(nabec_snf+"/"+snfOutname)
    else:
        print('not found?',nabec_snf+"/"+snfOutname)
        !gsutil -q cp {svSnf} nabec_snf/{snfOutname}
        nabec_snfFiles.append("nabec_snf/"+snfOutname)

print(len(nabec_snfFiles), 'samples')


"""
HBCC Sniffles VCFs and SNF Files
"""
# HBCC copy all the vcf files to this working directory, rename the files with their sample id
# store all the file names for bcftools merge
hbcc_vcfFiles = []
hbcc_snfFiles = []
snf_field = 'snifflesSnf_GRCh38'
vcf_field = 'snifflesVcf_GRCh38'
hbcc_snf = "hbcc_snf"
df = hbcc_cohort.loc[:,[snf_field, vcf_field,'sample_ID']]

for index, row in df.iterrows():
    svVcf = row[vcf_field]
    svSnf = row[snf_field]
    vcfOutname = str(index)+".sniffles.vcf"
    snfOutname = str(index)+".snf"
    
    if os.path.isfile(hbcc_snf+"/"+snfOutname):
        hbcc_snfFiles.append(hbcc_snf+"/"+snfOutname)
    else:
        print('localizing:',hbcc_snf+"/"+snfOutname)
        !gsutil -q cp {svSnf} hbcc_snf/{snfOutname}
        hbcc_snfFiles.append("hbcc_plots/"+snfOutname)

print(len(hbcc_snfFiles), 'samples')

# Reference for Merging
# copy the reference here 
ref_path = "gs://fc-secure-578c07ca-f188-44a1-a254-b524a5e73ecf/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ref_indx = "gs://fc-secure-578c07ca-f188-44a1-a254-b524a5e73ecf/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"

!gsutil -q cp {ref_path} grch38_ref.fna
!gsutil -q cp {ref_indx} grch38_ref.fna.fai


"""
Merge Samples
Merge NABEC Hapdiff + Truvari VCFs
"""

""" Truvari instructions:
https://github.com/acenglish/truvari/wiki/collapse"""
# make a space separated list of file names for bcftools merge 
vcf_files = " ".join(nabec_vcfFiles)
!bcftools merge --force-samples -m none {vcf_files} | bgzip > nabec_harmonized_merge.vcf.gz
!tabix nabec_harmonized_merge.vcf.gz


# run truvari collapse with relaxed --pctsize 0.75 --pctseq 0.75
!truvari collapse \
    --pctsize 0.75 \
    --pctseq 0.75 \
    --passonly \
    -i nabec_harmonized_merge.vcf.gz \
    -o nabec_harmonized_truvari_merge_pct75_07022024.vcf \
    -c nabec_harmonized_truvari_collapsed_pct75_07022024.vcf \
    -f grch38_ref.fna


"""
Merge HBCC Hapdiff + Truvari VCFs
"""
# make a space separated list of file names for bcftools merge 
vcf_files = " ".join(hbcc_vcfFiles)
!bcftools merge --force-samples -m none {vcf_files} | bgzip > hbcc_merge.vcf.gz
!tabix hbcc_merge.vcf.gz

# run truvari collapse with relaxed --pctsize 0.75 --pctseq 0.75
!truvari collapse \
    --pctsize 0.75 \
    --pctseq 0.75 \
    --passonly \
    -i hbcc_merge.vcf.gz \
    -o hbcc_harmonized_truvari_merge_pct75_09262024.vcf \
    -c hbcc_harmonized_truvari_collapsed_pct75_09262024.vcf \
    -f grch38_ref.fna


"""
Fill in genotypes For Truvari
"""
datatable_field = 'structuralVariantsAsm_confBed'
conf_directory = 'conf_beds'

df = nabec_cohort.loc[:,[datatable_field,'sample_ID']]
for index, row in df.iterrows():
    
    # filepath and file renamed with sample id
    bed = row[datatable_field]
    outname = str(index)+"_confident_regions.bed"
    
    # if the file for this sample is already in the directory 
    if os.path.isfile(conf_directory+"/"+outname):
        continue
    elif type(bed)!=float and len(bed)>0:
        print(index)
        !gsutil -q cp {bed} {conf_directory}/{outname}

        
df = hbcc_cohort.loc[:,[datatable_field,'sample_ID']]

for index, row in df.iterrows():
    
    # filepath and file renamed with sample id
    bed = row[datatable_field]
    outname = str(index)+"_confident_regions.bed"
    
    # if the file for this sample is already in the directory 
    if os.path.isfile(conf_directory+"/"+outname):
        continue
    elif type(bed)!=float and len(bed)>0:
        print(index)
        !gsutil -q cp {bed} {conf_directory}/{outname}

%%bash

# IN_VCF="NABEC_truvari_merge_biggerthan50bps_02082024.vcf"
IN_VCF="nabec_harmonized_truvari_merge_pct75_07022024.vcf"
OUT_VCF="nabec_harmonized_truvari_merge_genotyped_pct75_09262024.vcf"
IT=0


for b in conf_beds/NABEC*confident_regions.bed
do
echo $b
python genotype_truvari.py $b $IN_VCF $OUT_VCF
IN_VCF=$OUT_VCF
OUT_VCF=$IN_VCF
echo #in vcf
echo $IT #$IN_VCF

(( IT++ ))
done

%%bash


# IN_VCF="hbcc_harmonized_truvari_merge_pct75_07012024.vcf"
IN_VCF="hbcc_harmonized_truvari_merge_pct75_09262024.vcf"
OUT_VCF="hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.vcf"
IT=0


for b in conf_beds/HBCC*confident_regions.bed
do
echo $b
python genotype_truvari.py $b $IN_VCF $OUT_VCF
IN_VCF=$OUT_VCF
OUT_VCF=$IN_VCF
echo #in vcf
echo $IT #$IN_VCF

(( IT++ ))
done

""" Merge NABEC SNFs """
# install newest Sniffles2.3.3 version

!pip install --upgrade sniffles

# make tsv of inputs for Sniffles 
"""
Sniffles merge
sniffles --input snf_files_list.tsv --vcf multisample.vcf

"""
# NABEC
# write out snf file paths one per line for both cohorts
with open('snf_files_nabec.tsv','w') as nfile:
    for snf_nfile in nabec_snfFiles:
        nfile.write(snf_nfile+'\t'+snf_nfile.split('/')[1].split('.')[0]+'\n')
        
        
# HBCC write file paths to a file
with open('snf_files_hbcc.tsv','w') as hfile:
    for snf_hfile in hbcc_snfFiles:
        hfile.write(snf_hfile+'\n')

# NABEC 
print('nabec merging', datetime.datetime.now())

!sniffles --threads 64 --input snf_files_nabec.tsv --vcf nabec_snf2.3_multisample_07012024.vcf 

print('\nnabec done', datetime.datetime.now())


# HBCC 
print('HBCC merging', datetime.datetime.now())

!sniffles --threads 64 --input snf_files_hbcc.tsv --vcf hbcc_snf2.3_multisample_07012024.vcf

print('\nHBCC done', datetime.datetime.now())

"""
Merge SV Callers per Cohort
NABEC Merge
Merge: sort, normalize variants, and order the samples to match both SV Caller VCFs
"""

%%bash
#  concatenating the 2 vcf's and then truvari collapse to merge the 2 SV callers

bcftools sort nabec_harmonized_truvari_merge_genotyped_pct75_09262024.vcf | bgzip > nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz
tabix nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz

bcftools sort hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.vcf | bgzip > hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz
tabix hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz

%%bash

# Normalize the VCF files
bcftools norm -m -both -Oz -o nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.norm.vcf.gz nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz
# bcftools norm -m -both -Oz -o nabec_snf2.3_multisample_07012024.norm.vcf.gz nabec_snf2.3_multisample_07012024.vcf

bcftools index nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.norm.vcf.gz
# bcftools index nabec_snf2.3_multisample_07012024.norm.vcf.gz

# Merge the normalized VCF files
# !rm nabec_harmonized_pct75_snf2.3_07022024.vcf.gz
# !bcftools merge --force-samples -Oz -o nabec_harmonized_pct75_snf2.3_07022024.vcf.gz nabec_harmonized_truvari_merge_genotyped_pct75_07022024.sorted.norm.vcf.gz nabec_snf2.3_multisample_07012024.norm.vcf.gz

# reorder the samples in truvari 
!bcftools view --samples {samples}    -o nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.norm.ordered.vcf.gz nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.norm.vcf.gz

!bcftools index nabec_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.norm.ordered.vcf.gz
!tabix nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz

# remove poor QC sample
!bcftools view --samples {qc_sample}  -o nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz nabec_snf2.3_multisample_07012024.norm.vcf.gz
!bcftools index nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz

!tabix nabec_snf2.3_multisample_07012024.sampleFiltered.norm.vcf.gz


# !rm nabec_harmonized_pct75_snf2.3_07032024.vcf.gz
!bcftools concat --allow-overlaps -o nabec_harmonized_pct75_snf2.3_09262024.vcf.gz nabec_harmonized_truvari_merge_genotyped_pct75_09262024

# Change Merged Hapdiff None Quality Values to 0 & Change Symbolic SVs to Actual Sequence
# https://github.com/ACEnglish/truvari/discussions/216
!python truvari_resolve_symbolic_svs.py nabec_harmonized_pct75_snf2.3_09262024.vcf.gz grch38_ref.fna | \
    bgzip  > nabec_harmonized_pct75_snf2.3_09262024.symbolicQualResolved.vcf.gz

!bgzip nabec_harmonized_pct75_snf2.3_09262024.symbolicQualResolved_biggerthan50bps.sorted.vcf
!tabix nabec_harmonized_pct75_snf2.3_09262024.symbolicQualResolved_biggerthan50bps.sorted.vcf.gz

""" Merge the non-symbolic Snf and Qual 0 Hapdiff SVs"""

# try truvari merge on the combined VCF 
# run truvari collapse with relaxed --pctsize 0.75 --pctseq 0.75
!truvari collapse \
    --pctsize 0.75 \
    --pctseq 0.75 \
    --passonly --keep maxqual --chain \
    -i nabec_harmonized_pct75_snf2.3_09262024.symbolicQualResolved_biggerthan50bps.sorted.vcf.gz \
    -o nabec_harmonized_pct75_snf2.3_merge_pct75_09262024.symbolicQualResolved.keepMaxQual.vcf \
    -c nabec_harmonized_pct75_snf2.3_collapsed_pct75_09262024.symbolicQualResolved.keepMaxQual.vcf \
    -f grch38_ref.fna

!bcftools sort nabec_harmonized_pct75_snf2.3_merge_pct75_07032024.vcf.gz |bgzip > nabec_harmonized_pct75_snf2.3_merge_pct75_07032024.sorted.vcf.gz
!tabix nabec_harmonized_pct75_snf2.3_merge_pct75_07032024.sorted.vcf.gz


"""
HBCC merge
Merge: sort, normalize variants, and order the samples to match both SV Caller VCFs
"""

%%bash

# # Normalize the VCF files
# bcftools norm -m -both -Oz -o hbcc_harmonized_truvari_merge_genotyped_pct75_07012024.sorted.norm.vcf.gz hbcc_harmonized_truvari_merge_genotyped_pct75_07012024.sorted.vcf.gz
# bcftools norm -m -both -Oz -o hbcc_snf2.3_multisample_07012024.norm.vcf.gz hbcc_snf2.3_multisample_07012024.vcf

# rm nabec_harmonized_truvari_merge_genotyped_pct75_07022024.sorted.norm.vcf.gz
# rm hbcc_harmonized_truvari_merge_genotyped_pct75_07012024.sorted.norm.vcf.gz

# # Merge the normalized VCF files
# bcftools merge -Oz -o $merged_vcf $normalized_vcf1 $normalized_vcf2


# remove outlier samples from hapdiff merged VCF
!bcftools view --samples {qc_samples}   -Oz -o hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.sampleFiltered.vcf.gz hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.vcf.gz 

!tabix hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.sampleFiltered.vcf.gz

# remove outlier samples from sniffles meged VCF
!bcftools view --samples {qc_samples}   -Oz -o hbcc_snf2.3_multisample_07012024.sampleFiltered.vcf.gz hbcc_snf2.3_multisample_07012024.vcf

!tabix hbcc_snf2.3_multisample_07012024.sampleFiltered.vcf.gz

# Combine the 2 SV caller VCFs

!bcftools concat --allow-overlaps -o hbcc_harmonized_pct75_snf2.3_09262024.vcf.gz hbcc_harmonized_truvari_merge_genotyped_pct75_09262024.sorted.sampleFiltered.vcf.gz hbcc_snf2.3_multisample_07012024.sampleFiltered.vcf.gz        

!tabix hbcc_harmonized_pct75_snf2.3_09262024.vcf.gz

!python truvari_resolve_symbolic_svs.py hbcc_harmonized_pct75_snf2.3_09262024.vcf.gz grch38_ref.fna | \
    bgzip  > hbcc_harmonized_pct75_snf2.3_09262024.symbolicQualResolved.vcf.gz

# Merge the non-symbolic Snf and Qual 0 Hapdiff SVs
!truvari collapse \
    --pctsize 0.75 \
    --pctseq 0.75 \
    --passonly --keep maxqual --chain \
    -i hbcc_harmonized_pct75_snf2.3_09262024.symbolicQualResolved_biggerthan50bps.sorted.vcf.gz \
    -o hbcc_harmonized_pct75_snf2.3_merge_pct75_09262024.symbolicQualResolved.keepMaxQual.vcf \
    -c hbcc_harmonized_pct75_snf2.3_collapsed_pct75_09262024.symbolicQualResolved.keepMaxQual.vcf \
    -f grch38_ref.fna

# # combine VCFs fo each cohort 
!bcftools merge --force-samples -m none \
    nabec_harmonized_pct75_snf2.3_merge_pct75_09262024.symbolicQualResolved.keepMaxQual.sorted.vcf.gz \
    hbcc_harmonized_pct75_snf2.3_merge_pct75_09262024.symbolicQualResolved.keepMaxQual.sorted.vcf.gz | \
    bgzip > nabec_hbcc_allSV_merge_09262024.keepMaxQual.vcf.gz

!tabix nabec_hbcc_allSV_merge_09262024.keepMaxQual.vcf.gz

!truvari collapse \
    --pctsize 0.75 \
    --pctseq 0.75 \
    --passonly --keep maxqual --chain --sizemax 60000 \
    -i nabec_hbcc_allSV_merge_09262024.keepMaxQual.vcf.gz \
    -o nabec_hbcc_allSV_merge_pct75_09262024.keepMaxQual.vcf \
    -c nabec_hbcc_allSV_collapsed_pct75_09262024.keepMaxQual.vcf \
    -f grch38_ref.fna

# Extract multi-allelic SVs and Merge pct=0

def read_bed(bed_file):
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start'])
    return bed_df

# check if a variant is in a BED region
def is_variant_in_bed(chrom, pos, bed_df):
    region = bed_df[(bed_df['chrom'] == chrom) & (bed_df['start'] == pos) ]
    return not region.empty


# split VCF by bed 
def filter_vcf_by_bed(vcf_file, bed_file, output_vcf_file, kept_vcf_file):
    bed_df = read_bed(bed_file)
    
    vcf_in = pysam.VariantFile(vcf_file)
    vcf_out = pysam.VariantFile(output_vcf_file, 'w', header=vcf_in.header)
    vcf_keep_out = pysam.VariantFile(kept_vcf_file, 'w', header=vcf_in.header)
    
    for record in vcf_in:
        chrom = record.chrom
        pos = record.pos
        
        if is_variant_in_bed(chrom, pos, bed_df):
            vcf_out.write(record)
        else: 
            vcf_keep_out.write(record)
    
    vcf_in.close()
    vcf_out.close()
    vcf_keep_out.close()
    print(f'Filtered VCF written to {output_vcf_file}')
    print(f'Kept SV VCF written to {kept_vcf_file}')


bed_file="grep_svcount_chr.txt"
vcf_file="nabec_hbcc_allSV_merge_pct75_09302024.keepMaxQual.60kmax_nbl.vcf"

output_vcf_file = 'nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.vcf' 
kept_vcf_file = 'nabec_hbcc_allSV_merge_pct75_09302024.kept_output.vcf' 


!bcftools sort nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.vcf | bgzip > nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.sorted.vcf.gz
!tabix nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.sorted.vcf.gz

# merge all multi-allelic
!truvari collapse \
    --pctsize 0 \
    --pctseq 0 \
    --passonly --keep maxqual --chain --sizemax 60000 \
    -i nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.sorted.vcf.gz \
    -o nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.60kmax_nbl.vcf \
    -c nabec_hbcc_allSV_collapsed_pct75_09302024.filtered_output.keepMaxQual.60kmax_nbl.vcf \
    -f grch38_ref.fna


!bcftools sort nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.60kmax_nbl.vcf | bgzip > nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.60kmax_nbl.sorted.vcf.gz
!tabix nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.60kmax_nbl.sorted.vcf.gz

!truvari collapse \
    --pctsize 0 \
    --pctseq 0 --typeignore \
    --keep maxqual --chain --sizemax 180000 --refdist 2000\
    -i nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.100kmax_typeIgnore.sorted.vcf.gz \
    -o nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.vcf \
    -c nabec_hbcc_allSV_collapsed_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.vcf \
    -f grch38_ref.fna

!bgzip nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.vcf
!tabix nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.vcf.gz

# Combine the multi-allelic Merged VCF with the non-multi-allelic VCF
# This is the multi-allelic fixed VCF, lets limit the SVs to < 200kb

!bcftools view -i 'SVLEN >= -200000 & SVLEN <= 200000' nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.vcf.gz | \
   bcftools sort > nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.200kbmax.sorted.vcf


!bcftools concat --allow-overlaps -o nabec_hbcc_allSV_merge_pct75_09302024.uniqLoci.vcf \
    nabec_hbcc_allSV_merge_pct75_09302024.kept_output.sorted.vcf.gz \
    nabec_hbcc_allSV_merge_pct75_09302024.filtered_output.keepMaxQual.180kmax_typeIgnore.refDist_2kb.200kbmax.sorted.vcf.gz

!bcftools sort nabec_hbcc_allSV_merge_pct75_09302024.uniqLoci.vcf | bgzip > nabec_hbcc_allSV_merge_pct75_09302024.uniqLoci.sorted.vcf.gz 
!tabix nabec_hbcc_allSV_merge_pct75_09302024.uniqLoci.sorted.vcf.gz
