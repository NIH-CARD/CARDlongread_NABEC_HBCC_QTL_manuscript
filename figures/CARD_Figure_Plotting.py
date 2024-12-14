

nabec_color = '#1f77b4'
hbcc_color = '#ff7f0e'

def violin_swarm(x,y,data,ax,cohort_color,swarm_pt_size = 2.5):
    sns.violinplot(x=x, y=y, data=data, cut=0.25,
               inner=None, color=cohort_color, alpha = 0.4, ax=ax, edgecolor='black')
    sns.swarmplot(x=x, y=y, data=data, color=cohort_color, 
              s=swarm_pt_size, alpha=.75, ax=ax)
    
def violinBox_swarm(x,y,data,ax,cohort_color,swarm_pt_size = 2.5, innerplot="box"):
    
    sns.swarmplot(x=x, y=y, data=data, color=cohort_color, 
              s=swarm_pt_size, alpha=1, ax=ax)
    sns.violinplot(x=x, y=y, data=data, cut=0.25,
               inner=innerplot, color=cohort_color, alpha = 0.4, ax=ax, edgecolor='black',
                  inner_kws=dict(box_width=30, whis_width=2, color=".8"))


def lighten_color(color, amount=0.5):
    try:
        c = mcolors.cnames[color]
    except:
        c = color
    c = mcolors.to_rgba(c)
    return [1 - (1 - x) * amount for x in c]

def get_name_of_unit(n):
    if n == 1_000_000_000:
        return "Gbp"
    if n == 1_000_000:
        return "Mbp"
    if n == 1_000:
        return "Kbp"
    if n == 0:
        return "bp"

def plot_ngx(fig, axes, name, title, color, lengths, genome_size):

    total = float(sum(lengths))
    genome_size = float(genome_size)

    if total > genome_size:
        sys.stderr.write("WARNING: observed total sequence length is greater than genome size\n")

    unit = float(10)**float(max(0,round((np.log10(genome_size) - 3) / 3)*3))

    ng50 = 0.0
    n50 = 0.0

    x = list()
    y = list()

    x_prev = 0.0
    s_prev = 0.0
    for i,s in enumerate(sorted(lengths, reverse=True)):
        s = float(s)/float(unit)

        x0 = x_prev
        x1 = x_prev + s

        if x0*unit <= float(genome_size)/2.0 < x1*unit:
            ng50 = s*unit
#             print("ng50\t%.0f" % ng50)

        if x0*unit <= float(total)/2.0 < x1*unit:
            n50 = s*unit
#             print("n50\t%.0f" % n50)

        if i > 0:
            # pyplot.plot([x0,x0],[s_prev,s], color=color)
            x.extend([x0,x0])
            y.extend([s_prev,s])

        x.extend([x0,x1])
        y.extend([s,s])

        x_prev = x1
        s_prev = s

    # Bring down to 0 for easier comparison
    x.extend([x_prev,x_prev])
    y.extend([s_prev,0])

    axes.plot(x,y,color=color,label=name, alpha=0.15)

    axes.axvline(genome_size/2/unit, linestyle='--', linewidth=0.6, color='gray')
    axes.ticklabel_format(style='plain')


    axes.set_xlim([0,float(genome_size)/float(unit)])

    axes.set_xlabel("Cumulative coverage (%s)" % get_name_of_unit(unit)) #, fontsize = 12)
    axes.set_ylabel("Length (%s)" % get_name_of_unit(unit)) #, fontsize = 12)
    axes.set_title(title) #, fontsize = titlefontSize)

nabec_dual_lengths = ['nabec_dual_ngx/'+a for a in os.listdir('nabec_dual_ngx') if a[:5] == "NABEC"]
hbcc_dual_lengths = ["hbcc_dual_ngx/"+a for a in os.listdir('hbcc_dual_ngx') if a[:4] == "HBCC"]

nabec_phased_lengths = ['all_phased_ngx/'+a for a in os.listdir('all_phased_ngx') if a[:5] == "NABEC"]
hbcc_phased_lengths = ['all_phased_ngx/'+a for a in os.listdir('all_phased_ngx') if a[:4] == "HBCC"]


# Original colors
nabec_color = '#1f77b4'  # blue
hbcc_color = '#ff7f0e'  # orange

# Lightened colors
light_nabec_color = lighten_color(nabec_color, 0.7)  # Lighten by 70%
light_hbcc_color = lighten_color(hbcc_color, 0.7)


# go through the length defined wambam dfs and get the mean and median for each indiv
nabec_asm_1Mb_summary2_dfs = []

for name, wb_asum in nabec_cohort['asmDual1_wambam_alignedSummary_2'].dropna().items():
#     print(name)
    tmpdf = pd.read_csv(wb_asum, sep="\t")
    tdf = pd.DataFrame( {'id': [name],
                          'mean': [tmpdf.loc[tmpdf['inferred_len']>=1000000]['identity'].mean()],
                          'median': [tmpdf.loc[tmpdf['inferred_len']>=1000000]['identity'].median()]
                        })
    tdf.set_index('id', inplace=True)
    nabec_asm_1Mb_summary2_dfs.append(tdf)

nabec_hp1_1Mb_dual_id_2 = pd.concat(nabec_asm_1Mb_summary2_dfs)

print(nabec_hp1_1Mb_dual_id_2['median'].mean(), nabec_hp1_1Mb_dual_id_2['mean'].mean())

nabec_hp1_1Mb_dual_id_2.to_csv("NABEC_Grch38_Dual1_1Mb_2mismatchIndel_Identity.csv", header=True, index=True)

# go through the length defined wambam dfs and get the mean and median for each indiv
hbcc_asm_1Mb_summary2_dfs = []

for name, wb_asum in hbcc_cohort['asmDual1_wambam_alignedSummary_2'].dropna().items():
#     print(name)
    tmpdf = pd.read_csv(wb_asum, sep="\t")
    tdf = pd.DataFrame( {'id': [name],
                          'mean': [tmpdf.loc[tmpdf['inferred_len']>=1000000]['identity'].mean()],
                          'median': [tmpdf.loc[tmpdf['inferred_len']>=1000000]['identity'].median()]
                        })
    tdf.set_index('id', inplace=True)
    hbcc_asm_1Mb_summary2_dfs.append(tdf)

hbcc_hp1_1Mb_dual_id_2 = pd.concat(hbcc_asm_1Mb_summary2_dfs)

hbcc_hp1_1Mb_dual_id_2.to_csv("HBCC_Grch38_Dual1_1Mb_2mismatchIndel_Identity.csv", header=True, index=True)


# load in data
hbcc_all = pd.read_csv("HBCC_Grch38_readStats.csv")
hbcc_all = hbcc_all.set_index('sample')

nabec_all = pd.read_csv("NABEC_Grch38_readStats.csv")
nabec_all = nabec_all.set_index('sample')

hbcc_all_dual1_1mb = pd.read_csv("HBCC_Grch38_Dual1_1Mb_2mismatchIndel_Identity.csv")
hbcc_all_dual1_1mb = hbcc_all_dual1_1mb.set_index('id')

hbcc_all_dual1 = pd.read_csv("HBCC_Grch38_Dual1_Stats.csv")
hbcc_all_dual1 = hbcc_all_dual1.set_index('sample')                             


nabec_all_dual1_1mb = pd.read_csv("NABEC_Grch38_Dual1_1Mb_2mismatchIndel_Identity.csv")
nabec_all_dual1_1mb = nabec_all_dual1_1mb.set_index('id')

nabec_all_dual1 = pd.read_csv("NABEC_Grch38_Dual1_Stats.csv")
nabec_all_dual1 = nabec_all_dual1.set_index('sample')

gencode_pc_genes_asmbld = pd.read_csv("hg38.gencode.v44._proteinCoding_numberAssembled_100pct.csv")
gencode_pc_genes_asmbld = gencode_pc_genes_asmbld.set_index('ID')

nabec_all_filt = nabec_all.drop(bad_nabec_samples)
hbcc_all_filt = hbcc_all.drop(bad_hbcc_samples)

nabec_all_dual1_filt = nabec_all_dual1.drop(bad_nabec_samples)
hbcc_all_dual1_filt = hbcc_all_dual1.drop(bad_hbcc_samples)

nabec_all_dual1_1mb = nabec_all_dual1_1mb.drop(bad_nabec_samples)
hbcc_all_dual1_1mb = hbcc_all_dual1_1mb.drop(bad_hbcc_samples)

gencode_pc_genes_asmbld_filt = gencode_pc_genes_asmbld


nabec_all_filt['n50_kb'] = nabec_all_filt['n50']/1000
hbcc_all_filt['n50_kb'] = hbcc_all_filt['n50']/1000

nabec_all_dual1_1mb['pct_med.identity'] = nabec_all_dual1_1mb['median']*100
hbcc_all_dual1_1mb['pct_med.identity'] = hbcc_all_dual1_1mb['median']*100

gencode_pc_genes_asmbld_filt['percent_PC_Genes'] = gencode_pc_genes_asmbld_filt['pct_PC_Genes']*100

# Figure 1 Ages and Reads
# Create the figure and GridSpec
fig = plt.figure(figsize=(12,9))
gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.25], height_ratios=[1, 1, 1])

nabec_color = '#1f77b4'  # blue
hbcc_color = '#ff7f0e'  # orange
titlefontsize = 18



""" col1 all rows """
ax1 = fig.add_subplot(gs[0:3, 0])
violin_swarm(['NABEC']*nabec_cohort.shape[0], 
             'Age', nabec_cohort, ax1, nabec_color) 

violin_swarm(['HBCC']*hbcc_cohort.shape[0], 
             'AgeDeath', hbcc_cohort, ax1, hbcc_color)
ax1.set_title('Cohort Ages \n(Age at Death)', fontsize = titlefontsize)
ax1.set_ylabel('')
ax1.tick_params(axis='y', labelsize=16)

""" Read Stats: N50, Read ID, Coverage"""
ax2 = fig.add_subplot(gs[0, 1])
violin_swarm('n50_kb', ['NABEC']*nabec_all_filt.shape[0], 
              nabec_all_filt, ax2, nabec_color, swarm_pt_size = 1.5)
violin_swarm('n50_kb', ['HBCC']*hbcc_all_filt.shape[0], 
              hbcc_all_filt, ax2, hbcc_color, swarm_pt_size = 1.5)
nMed = nabec_all_filt['n50_kb'].median()
hMed = hbcc_all_filt['n50_kb'].median()
for i, m in enumerate([nMed, hMed]):
    ax2.plot([m,m], [i-0.25,i+0.25], color='firebrick', linestyle='-', alpha=0.75)
    if i==0:
        ax2.scatter(y=i, x=m, marker='s', s=10 , color='firebrick', label='Median')
    else:
        ax2.scatter(y=i, x=m, marker='s', s=10 , color='firebrick')

ax2.legend(loc='center left', fontsize='small') #, bbox_to_anchor=(-0.1, 0.5))
ax2.set_title('Read N50s (Kbp)', fontsize = titlefontsize)
ax2.set_xlabel('')


""" prop.length.geq.95 """
ax3 = fig.add_subplot(gs[1, 1])

violin_swarm('prop.length.geq.95',['NABEC']*nabec_all_filt.shape[0], 
              nabec_all_filt, ax3, nabec_color, swarm_pt_size = 1.5)
violin_swarm('prop.length.geq.95',['HBCC']*hbcc_all_filt.shape[0], 
              hbcc_all_filt, ax3, hbcc_color, swarm_pt_size = 1.5)
nMed = nabec_all_filt['prop.length.geq.95'].median()
hMed = hbcc_all_filt['prop.length.geq.95'].median()
for i, m in enumerate([nMed, hMed]):
    ax3.plot([m,m], [i-0.25,i+0.25], color='firebrick', linestyle='-', alpha=0.75)
    ax3.scatter(y=i, x=m, marker='s', s=10 , color='firebrick')
ax3.set_xlabel('')
ax3.set_title('Proportion of Reads >10Kb', fontsize = titlefontsize)


""" genome Cov """
ax4 = fig.add_subplot(gs[2, 1])
violin_swarm('gCov', ['NABEC']*nabec_all_filt.shape[0], 
             nabec_all_filt, ax4, nabec_color, swarm_pt_size = 1.5)
violin_swarm('gCov', ['HBCC']*hbcc_all_filt.shape[0], 
              hbcc_all_filt, ax4, hbcc_color, swarm_pt_size = 1.5)
nMed = nabec_all_filt['gCov'].median()
hMed = hbcc_all_filt['gCov'].median()
for i, m in enumerate([nMed, hMed]):
    ax4.plot([m,m], [i-0.25,i+0.25], color='firebrick', linestyle='-', alpha=0.75)
    ax4.scatter(y=i, x=m, marker='s', s=10 , color='firebrick')
ax4.set_xlabel('')
ax4.set_title('Read GRCh38 Coverage', fontsize = titlefontsize)
# ax6.set_ylabel('Read Coverage')

for ax in [ax2, ax3, ax4]:
    
    ax.yaxis.tick_right()          # Move ticks to the right
    ax.yaxis.set_label_position("right")  # Move labels to the right
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='x', labelsize=16)

    
plt.savefig('figures/figure1_ages_readStats_12122024.png',dpi=300, transparent=True)

# Figure 2 Reads and Assembly
# Create the figure and GridSpec
fig = plt.figure(figsize=(14, 14))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1])

nabec_color = '#1f77b4'  # blue
hbcc_color = '#ff7f0e'  # orange
titlefontsize = 18
labelfontsize = 18

""" Row 3 """
# NGX
ax1 = fig.add_subplot(gs[0, 0:2]) #0:2])
# fig.add_subplot(gs[1, 0:2])
genome_size = 3100000000

for n in nabec_dual_lengths:
    lengths = []
    with open(n,'r') as file:
        for line in file:
            lengths.append(float(line.strip()))
    
    plot_ngx(
            name='NABEC', title="",
            color=nabec_color,
            fig=fig,
            axes=ax1,
            lengths=lengths,
            genome_size=genome_size
        )
    
for n in hbcc_dual_lengths:
    lengths = []
    with open(n,'r') as file:
        for line in file:
            lengths.append(float(line.strip()))
    
    plot_ngx(
            name='HBCC', title="",
            color=hbcc_color,
            fig=fig,
            axes=ax1,
            lengths=lengths,
            genome_size=genome_size
        )

    
custom_line1 = mlines.Line2D([], [], color=nabec_color, label='NABEC')
custom_line2 = mlines.Line2D([], [], color=hbcc_color, label='HBCC')
ax1.legend(handles=[custom_line1, custom_line2], loc='upper right', fontsize='large',
          prop={'size': 14})
ax1.tick_params(axis='both', labelsize=18)
ax1.set_title("Hapdup Dual Assembly NGx", fontsize=titlefontsize)
ax1.set_xlabel("Cumulative coverage (Mbp)" , fontsize = labelfontsize)
ax1.set_ylabel("Length (Mbp)" , fontsize = labelfontsize)

# # Dual Assembly n50
ax2 = fig.add_subplot(gs[0, 2])
phased_merged_df = pd.read_csv('all_phased_ngx/ng50.txt', header=None, names = ['sampleId','ngx'])

n_phased_merged_df = phased_merged_df.loc[phased_merged_df['sampleId'].str.startswith('N')]
h_phased_merged_df = phased_merged_df.loc[phased_merged_df['sampleId'].str.startswith('H')]

violin_swarm(['NABEC']*n_phased_merged_df.shape[0], 
             'ngx', n_phased_merged_df, ax2, nabec_color)

violin_swarm(['HBCC']*h_phased_merged_df.shape[0], 
             'ngx', h_phased_merged_df, ax2, hbcc_color)

nMed = n_phased_merged_df['ngx'].median()
hMed = h_phased_merged_df['ngx'].median()
for i, m in enumerate([nMed, hMed]):
    ax2.plot([i-0.25,i+0.25], [m,m], color='firebrick', linestyle='-', alpha=0.75)
    ax2.scatter(x=i, y=m, marker='s', s=20 , color='firebrick')
    
ax2.set_title('Assembly Phase Block \nNG50s (Mbp)', fontsize=titlefontsize)
ax2.set_ylabel('Mbp', fontsize = labelfontsize)
ax2.tick_params(axis='y', labelsize=18)
ax2.set_ylabel('')
# ax8.set_facecolor('#f0f0f0')

""" total.gbp """
ax3 = fig.add_subplot(gs[1, 0])
violin_swarm(['NABEC']*nabec_all_dual1_filt.shape[0], 
             'total.gbp', nabec_all_dual1_filt, ax3, nabec_color)
violin_swarm(['HBCC']*hbcc_all_dual1_filt.shape[0], 
             'total.gbp', hbcc_all_dual1_filt, ax3, hbcc_color)
nMed = nabec_all_dual1_filt['total.gbp'].median()
hMed = hbcc_all_dual1_filt['total.gbp'].median()
for i, m in enumerate([nMed, hMed]):
    ax3.plot([i-0.25,i+0.25], [m,m], color='firebrick', linestyle='-', alpha=0.75)
    ax3.scatter(x=i, y=m, marker='s', s=20 , color='firebrick')
    
ax3.tick_params(axis='y', labelsize=18)
ax3.set_title('Assembled Base Pairs (Gbp)', fontsize=titlefontsize)
ax3.set_ylabel('')
# ax9.set_facecolor('#f0f0f0')

""" median.identity """
ax4 = fig.add_subplot(gs[1, 1])
violin_swarm(['NABEC']*nabec_all_dual1_1mb.shape[0], 
             'pct_med.identity', nabec_all_dual1_1mb, ax4, nabec_color)
violin_swarm(['HBCC']*hbcc_all_dual1_1mb.shape[0], 
             'pct_med.identity', hbcc_all_dual1_1mb, ax4, hbcc_color)
nMed = nabec_all_dual1_1mb['pct_med.identity'].median()
hMed = hbcc_all_dual1_1mb['pct_med.identity'].median()
for i, m in enumerate([nMed, hMed]):
    ax4.plot([i-0.25,i+0.25], [m,m], color='firebrick', linestyle='-', alpha=0.75)
    ax4.scatter(x=i, y=m, marker='s', s=20 , color='firebrick')
    
ax4.tick_params(axis='y', labelsize=18)
ax4.set_title('Assembly Identity CHM13', fontsize=titlefontsize)
ax4.set_ylabel('')

""" Protein Coding Genes Assembled """
ax5 = fig.add_subplot(gs[1, 2])
violin_swarm(['NABEC']*gencode_pc_genes_asmbld_filt.iloc[146::,:].shape[0], 
             'percent_PC_Genes', gencode_pc_genes_asmbld_filt.iloc[146::,:], ax5, nabec_color)
violin_swarm(['HBCC']*gencode_pc_genes_asmbld_filt.iloc[0:146,:].shape[0], 
             'percent_PC_Genes', gencode_pc_genes_asmbld_filt.iloc[0:146,:], ax5, hbcc_color)
nMed = gencode_pc_genes_asmbld_filt.iloc[146::,:]['percent_PC_Genes'].median()
hMed = gencode_pc_genes_asmbld_filt.iloc[0:146,:]['percent_PC_Genes'].median()
for i, m in enumerate([nMed, hMed]):
    ax5.plot([i-0.25,i+0.25], [m,m], color='firebrick', linestyle='-', alpha=0.75)
    ax5.scatter(x=i, y=m, marker='s', s=20 , color='firebrick')
    
ax5.set_title('Percent of Genbank Protein Coding \nGenes Fully Assembled', fontsize=titlefontsize)
ax5.tick_params(axis='y', labelsize=18)
ax5.set_ylabel('')
# ax10.set_facecolor('#f0f0f0')

# # Add letters to the panels
axes = [ax1, ax2, ax3, ax4, ax5]
labels = ['a', 'b', 'c', 'd', 'e']

for ax, label in zip(axes, labels):
    ax.text(-0.1, 1.25, label, transform=ax.transAxes, 
            fontsize=26, va='top', ha='right')
    
plt.subplots_adjust(wspace=0.6, hspace=0.6)
# plt.tight_layout()

plt.savefig('figures/figure2_assemblyStats_filtered_12122024.png',dpi=300)


