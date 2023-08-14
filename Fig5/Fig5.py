import pandas as pd
import pyranges as pr
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import plot, from_contents
import matplotlib.cm as cm

# v1 cCRE peak annot, make the annotation consistent with Bin Ren and Wang's paper. delete genebody and TTS peak and
# add distance information to nearest gene.
# we have the following rules:
# 1. The goal of peak annotation is to map peak to gene symbols, which is the union of all transcripts of a given gene.
# 2. A peak can be mapped to multiple genes.
# 3. A peak can only be one type of peaks for a given gene, which means a peak cannot be annotated as both a promoter peak and a distal peak of the same gene.
# 4. Only protein coding genes are included for annotation.
# 6. If a peak overlaps with promoter region (-200 bp, +200 bp) of any transcription start site (TSS), it is annotated as a promoter peak of the gene.
# 7. If a peak overlaps with proximal region (-2000 ~ -200 bp, +200 ~ + 2000 bp) of any transcription start site (TSS), it is annotated as a promoter peak of the gene.
# 8. If a peak is within +-200 kb to +- 2kb of the closest TSS, and if it is not a promoter or proximal peak of the gene of the closest TSS, it will be annotated as a distal peak of that gene.
# 9. If a peak overlaps the body of a transcript, and it is not a promoter nor a proximal nor a distal peak of the gene, it will be annotated as a distal peak of that gene with distance set as zero.
# 10. If a peak has not been mapped to any gene at the step, it will be annotated as an intergenic peak without a gene symbol assigned.
EnsDb_Hsapiens_v86_seqlen = {
    "chr1":248956422, "chr2":242193529, "chr3":198295559, "chr4":190214555, "chr5":181538259,
    "chr6":170805979, "chr7":159345973, "chr8":145138636, "chr9":138394717, "chr10":133797422,
    "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX":156040895, 'chrY':57227415
}
class cCREAnnot(object):
    def __init__(self, df, gene_annot_df, genomic_region_ranges=None, genomic_regions_idx=None):
        """
        :param df: dataframe with first three columns: chr/start/end
        :param gene_annot_df: dataframe for gene reference
               gene_annot_df.columns = ['Chromosome', 'Start', 'End', 'Strand', 'GeneName']
        :param genomic_region_ranges:
        :param genomic_regions_idx:
        """
        tmp_df = df.iloc[:,:3]
        tmp_df.columns = ['Chromosome','Start','End']
        self.df = tmp_df
        self.gene_annot_df = gene_annot_df
        if genomic_region_ranges is None:
            self.genomic_regions_ranges = {
                # prom:-1000 ~ 1000 of TTS
                'Promoter': [-200, -200, 200, 200],
                'Genebody': [200, 0],
                'Proximal': [-2000, 200, -200, 2000],
                'Distal': [-200000, 200, -200, 200000]
            }
        else:
            self.genomic_regions_ranges = genomic_region_ranges

        if genomic_regions_idx is None:
            # genomic location: promoter > proximal> diatal > genebody > intergenic
            self.genomic_regions_idx = {'Promoter':2**4,'Proximal':2**3,'Distal': 2**2,'Genebody':2**1}
        else:
            self.genomic_regions_idx = genomic_regions_idx

    def _GenerateGenomicRegion(self):
        genomic_regions = {}
        for genoloc in self.genomic_regions_ranges:
            gene_annot_cp = self.gene_annot_df.copy()
            current_params = self.genomic_regions_ranges[genoloc]
            if genoloc not in ['Genebody','Promoter']:
                gene_annot_cp_up = gene_annot_cp.copy()
                gene_annot_cp_dw = gene_annot_cp.copy()
                # upstream direction
                gene_annot_cp_up[f'{genoloc}_start'] = np.where(gene_annot_cp_up['Strand'] == '+',
                                                                gene_annot_cp_up['Start'] + current_params[0],
                                                                gene_annot_cp_up['End'] + current_params[1])
                gene_annot_cp_up[f'{genoloc}_end'] = np.where(gene_annot_cp_up['Strand'] == '+',
                                                              gene_annot_cp_up['Start'] + current_params[2],
                                                              gene_annot_cp_up['End'] + current_params[3])
                # downstream direction
                gene_annot_cp_dw[f'{genoloc}_start'] = np.where(gene_annot_cp_dw['Strand'] == '+',
                                                                gene_annot_cp_dw['Start'] - current_params[2],
                                                                gene_annot_cp_dw['End'] - current_params[3])
                gene_annot_cp_dw[f'{genoloc}_end'] = np.where(gene_annot_cp_dw['Strand'] == '+',
                                                              gene_annot_cp_dw['Start'] - current_params[0],
                                                              gene_annot_cp_dw['End'] - current_params[1])
                gene_annot_cp = pd.concat([gene_annot_cp_up,gene_annot_cp_dw]).reset_index(drop=True)
            elif genoloc == 'Promoter':
                gene_annot_cp[f'{genoloc}_start'] = np.where(gene_annot_cp['Strand'] == '+',
                                                             gene_annot_cp['Start'] + current_params[0],
                                                             gene_annot_cp['End'] + current_params[1])
                gene_annot_cp[f'{genoloc}_end'] = np.where(gene_annot_cp['Strand'] == '+',
                                                           gene_annot_cp['Start'] + current_params[2],
                                                           gene_annot_cp['End'] + current_params[3])
            elif genoloc == 'Genebody':
                gene_annot_cp[f'{genoloc}_start'] = gene_annot_cp['Start'] + current_params[0]
                gene_annot_cp[f'{genoloc}_end'] = gene_annot_cp['End'] + current_params[1]
            else:
                print('Not possible condition')
            # restrict within chromosome length
            gene_annot_cp.loc[gene_annot_cp[f'{genoloc}_start'] < 0, f'{genoloc}_start'] = 1
            gene_annot_cp.loc[gene_annot_cp[f'{genoloc}_end'] < 0, f'{genoloc}_end'] = 1
            gene_annot_cp['maxsize'] = gene_annot_cp['Chromosome'].map(EnsDb_Hsapiens_v86_seqlen)
            gene_annot_cp.loc[gene_annot_cp[f'{genoloc}_start'] > gene_annot_cp['maxsize'], f'{genoloc}_start'] = gene_annot_cp['maxsize']
            gene_annot_cp.loc[gene_annot_cp[f'{genoloc}_end'] > gene_annot_cp['maxsize'], f'{genoloc}_end'] = gene_annot_cp['maxsize']

            gene_annot_cp = gene_annot_cp[['Chromosome', f'{genoloc}_start', f'{genoloc}_end', 'GeneName']]
            gene_annot_cp.columns = ['Chromosome', 'Start', 'End', 'GeneName']
            genomic_regions[genoloc] = gene_annot_cp
        return genomic_regions

    def _genoloc_annot_dict(self, regions_pr, gene_region_pr):
        regions_annot = regions_pr.join(gene_region_pr).as_df()
        # in case no overlaps found
        if len(regions_annot) == 0:
            return {}
        else:
            regions_set = set(
                regions_annot[['Chromosome', 'Start', 'End', 'GeneName']].itertuples(index=False, name=None))
            regions_annot_dict = defaultdict(list)
            for x in regions_set:
                regions_annot_dict['_'.join([x[0], str(x[1]), str(x[2])])].append(x[3])
            regions_annot_dict = {key: ','.join(value) for key, value in regions_annot_dict.items()}
        return regions_annot_dict

    def annot(self):
        df_annot = self.df.copy()
        df_annot['key'] = df_annot['Chromosome'] + '_' + df_annot['Start'].astype(str) + '_' + df_annot['End'].astype(str)
        print("Generating genomic region range from reference")
        genomic_regions = self._GenerateGenomicRegion()
        print("Annotating")
        df_annot_pr = pr.PyRanges(df_annot)
        for genoloc in genomic_regions:
            genomic_regions_pr = pr.PyRanges(genomic_regions[genoloc])
            genolocdict = self._genoloc_annot_dict(df_annot_pr, genomic_regions_pr)
            df_annot[f'{genoloc}'] = df_annot['key'].map(genolocdict)
            df_annot[f'{genoloc}_idx'] = df_annot[f'{genoloc}'].notnull().astype('int') * self.genomic_regions_idx[genoloc]

        idx_cols = [f'{genoloc}_idx' for genoloc in genomic_regions]
        df_annot['Intergenic'] = 'None'
        df_annot['sum'] = df_annot[idx_cols].sum(axis=1)
        # initial final_annot, order is important here
        df_annot['genomeLoc_annot'] = 'Intergenic'
        # see rule 9
        df_annot.loc[df_annot['sum'] >= 2, 'genomeLoc_annot'] = 'Distal'
        df_annot.loc[df_annot['sum'] >= 8, 'genomeLoc_annot'] = 'Proximal'
        df_annot.loc[df_annot['sum'] >= 16, 'genomeLoc_annot'] = 'Promoter'

        # Apply the function to create the new column 'f'
        df_annot['genes'] = df_annot.apply(lambda x:x[x['genomeLoc_annot']], axis=1)

        gene_df_annot_final = df_annot[['Chromosome', 'Start', 'End', 'genomeLoc_annot','genes']]
        return gene_df_annot_final

# load annotation file
gene_annot_df = pd.read_csv('/Users/kfang/Documents/lab/Jin_lab/refseq/hg38/gene_hg38_ranges.csv')
gene_annot_df = gene_annot_df[['seqnames','start','end','strand','symbol']]
gene_annot_df.columns = ['Chromosome', 'Start', 'End', 'Strand', 'GeneName']
peaks_df = pd.read_csv('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scATAC/results/TTs_peaks_filt_df.csv')
outdir = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig4/v2'
# annot the peaks
ccre_annot = cCREAnnot(peaks_df,gene_annot_df)
peaks_annot_df = ccre_annot.annot()

peaks_annot_df = pd.merge(peaks_annot_df,peaks_df[['seqnames','start','end','peak_called_in']],left_on=['Chromosome','Start','End'],
                          right_on=['seqnames','start','end'])
peaks_annot_df = peaks_annot_df[['Chromosome','Start','End','genomeLoc_annot','genes','peak_called_in']]
peaks_annot_df['cell_states'] = peaks_annot_df['peak_called_in'].str.split(',')
peaks_annot_df['region'] = peaks_annot_df['Chromosome'] + '-' + peaks_annot_df['Start'].astype(str) + '-' + peaks_annot_df['End'].astype(str)
peaks_annot_expand_df = peaks_annot_df.explode('cell_states')
peaks_annot_df.drop(columns=['cell_states'],inplace=True)
peaks_annot_df.to_csv('/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scATAC/results/TTs_peaks_filt_annot_df.csv')

peaks_count_mat = pd.crosstab(columns=peaks_annot_expand_df['genomeLoc_annot'],index=peaks_annot_expand_df['cell_states'])
peaks_count_mat = peaks_count_mat[['Distal','Proximal','Promoter','Intergenic']]

fig, ax = plt.subplots(1,1)
# create stacked bar
peaks_count_mat.plot(kind='bar', stacked=True, color=['chocolate', 'darkorange', 'orange', 'grey'], ax=ax).legend(
    bbox_to_anchor=(1.0, 1.0)
)
# labels for x & y axis
ax.set_ylabel('count')
ax.set_yscale('log', basey=2)
ax.spines[['right', 'top']].set_visible(False)
plt.tight_layout()
# title of plot
plt.savefig(f'{outdir}/peaks_compoistion.png',dpi=300)
plt.close()

# load differential peaks for rt_cs1-3
cs_list = ['RT_CS1','RT_CS2','RT_CS3']
genomeLocs = ['Distal','Proximal','Promoter','Intergenic']
dpeaks_f = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig4/v2/CS_dpeaks_macs2.xlsx'
# scATAC differential peaks
cs_dp = {}
cs_dp_expand = {}
for cs in cs_list:
    print(f"Loading scATAC {cs}")
    rt_cs_dp = pd.read_excel(dpeaks_f,sheet_name=cs)
    rt_cs_dp = pd.concat([rt_cs_dp['region'].str.split('-',expand=True),rt_cs_dp],axis=1)
    rt_cs_dp.columns = ['Chromosome','Start','End'] + list(rt_cs_dp.columns[3:])
    rt_cs_dp[['Start','End']] = rt_cs_dp[['Start','End']].astype(int)
    # Annot gene
    rt_cs_dp = pd.merge(rt_cs_dp, peaks_annot_df[['Chromosome','Start','End','genomeLoc_annot','genes']],
                        on=['Chromosome','Start','End'])
    rt_cs_dp['type'] = np.where(rt_cs_dp['avg_log2FC']>0,'UpA','DownA')
    cs_dp[cs] = rt_cs_dp

    # Split the rows based on commas and explode the resulting values
    rt_cs_dp['geneslist'] = rt_cs_dp['genes'].str.split(',')
    rt_cs_dp_expand = rt_cs_dp.explode('geneslist')
    # Reset the index of the dataframe
    rt_cs_dp_expand = rt_cs_dp_expand.reset_index(drop=True)
    cs_dp_expand[cs] = rt_cs_dp_expand

# save sync.all.genes
all_dp_writer = pd.ExcelWriter(f'{outdir}/Dp_all.xlsx')
for cs in cs_list:
    cs_dp[cs].to_excel(all_dp_writer,cs)
all_dp_writer.close()
# Profiling annotation
fig,axs = plt.subplots(3, 1, figsize=(3,6), sharex=True)
for idx, cs in enumerate(cs_list):
    sns.countplot(data=cs_dp[cs],x='genomeLoc_annot',hue='type',hue_order=['UpA','DownA'],order=genomeLocs,ax=axs[idx],
                  )
    axs[idx].legend().set_visible(False)
    axs[idx].spines[['right', 'top']].set_visible(False)
    axs[idx].set_ylabel(cs)
    axs[idx].set_xlabel("")
lines_labels = [axs[0].get_legend_handles_labels()]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines, labels,loc='right',bbox_to_anchor=(1.0, 0.9))
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f'{outdir}/dp_genomeLoc_count.png',dpi=300)
plt.close()

# RNA
deg_f = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig3/RT_RPT_CS_marker_all.xlsx'
deg_module_f = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/results/Fig3/RT_RPT_CS_module_all.xlsx'
cs_de_all = {}
cs_de_module = {}
for cs in cs_list:
    print(f"Loading scRNA {cs}")
    rt_cs_de_all = pd.read_excel(deg_f,sheet_name=f'{cs}_all')
    rt_cs_de_module = pd.read_excel(deg_module_f,sheet_name=f'{cs}_module')
    cs_de_all[cs] = rt_cs_de_all
    cs_de_module[cs] = rt_cs_de_module

# find uniq and common genes between dps and de_all
# disgard intergenic region
cs_dp_de_both = {}
for cs in cs_list:
    print(f"Processing scATAC scRNA {cs}")
    rt_cs_dp_expand = cs_dp_expand[cs]
    rt_cs_dp_expand_upg = rt_cs_dp_expand[rt_cs_dp_expand['avg_log2FC'] > 0]
    rt_cs_dp_expand_downg = rt_cs_dp_expand[rt_cs_dp_expand['avg_log2FC'] < 0]

    rt_cs_de_all = cs_de_all[cs]
    rt_cs_de_all_upg = rt_cs_de_all[rt_cs_de_all['avg_log2FC']>0]
    rt_cs_de_all_downg = rt_cs_de_all[rt_cs_de_all['avg_log2FC']<0]

    rt_cs_dp_de_both_up = list(set(rt_cs_dp_expand_upg['genes']) & set(rt_cs_de_all_upg['gene']))
    rt_cs_dp_de_both_down = list(set(rt_cs_dp_expand_downg['genes']) & set(rt_cs_de_all_downg['gene']))
    rt_cs_dp_de_both_upg = rt_cs_dp_expand_upg[rt_cs_dp_expand_upg['genes'].isin(rt_cs_dp_de_both_up)]
    rt_cs_dp_de_both_downg = rt_cs_dp_expand_downg[rt_cs_dp_expand_downg['genes'].isin(rt_cs_dp_de_both_down)]
    rt_cs_dp_de_both = pd.concat([rt_cs_dp_de_both_upg,rt_cs_dp_de_both_downg])
    de_avgfc = dict(zip(rt_cs_de_all['gene'],rt_cs_de_all['avg_log2FC']))
    rt_cs_dp_de_both['avg_log2FC_scRNA'] = rt_cs_dp_de_both['genes'].map(de_avgfc)
    rt_cs_dp_de_both['sync_score'] = np.sign(rt_cs_dp_de_both['avg_log2FC_scRNA'] * rt_cs_dp_de_both['avg_log2FC']) * \
                                     np.sqrt(np.abs(rt_cs_dp_de_both['avg_log2FC_scRNA']) *
                                             np.abs(rt_cs_dp_de_both['avg_log2FC']))
    cs_dp_de_both[cs] = rt_cs_dp_de_both.sort_values('sync_score', ascending=False)

# save sync.all.genes
all_dp_de_writer = pd.ExcelWriter(f'{outdir}/Dp_De_all.xlsx')
for cs in cs_list:
    cs_dp_de_both[cs].to_excel(all_dp_de_writer,cs)
all_dp_de_writer.close()
# sync scattterplot
sync_scatter_list = []
for cs in cs_list:
    rt_cs_dp_expand = cs_dp_expand[cs].copy()
    rt_cs_de_all = cs_de_all[cs]
    de_all_avgfc = dict(zip(rt_cs_de_all['gene'],rt_cs_de_all['avg_log2FC']))
    rt_cs_dp_expand['avg_log2FC_scRNA'] = rt_cs_dp_expand['genes'].map(de_all_avgfc)
    rt_cs_dp_expand.dropna(inplace=True)
    rt_cs_dp_expand['cell_states'] = cs
    rt_cs_out = rt_cs_dp_expand[['avg_log2FC','avg_log2FC_scRNA','cell_states','genomeLoc_annot']]
    rt_cs_out.columns = ['avg_log2FC_scATAC','avg_log2FC_scRNA','cell_states','region']
    sync_scatter_list.append(rt_cs_out)

# sync scatterplot
sync_scatter_df = pd.concat(sync_scatter_list)
sync_scatter_df['sync_score'] = np.sign(sync_scatter_df['avg_log2FC_scRNA'] * sync_scatter_df['avg_log2FC_scATAC']) * \
                                     np.sqrt(np.abs(sync_scatter_df['avg_log2FC_scRNA']) *
                                             np.abs(sync_scatter_df['avg_log2FC_scATAC']))
sync_avg = sync_scatter_df['sync_score'].mean()
fig, ax = plt.subplots(1,1)
g = sns.scatterplot(data=sync_scatter_df,x='avg_log2FC_scRNA',y='avg_log2FC_scATAC',hue='cell_states',s=6,style='region',
                palette={'RT_CS1':'orangered','RT_CS2':'lime','RT_CS3':'cyan'},  linewidth = 0 ,ax=ax)
ax.spines[['right', 'top']].set_visible(False)
g.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1)
plt.tight_layout()
plt.savefig(f'{outdir}/sync.scatter.png',dpi=300)
plt.close()


# find uniq and common genes between dps and de_module
module_summary_dict = {'counts':[], 'cell_states':[], 'type':[],'region':[]}
cs_dp_de_module = {}
for cs in cs_list:
    print(f"Processing scATAC scRNA {cs}")
    rt_cs_dp_expand = cs_dp_expand[cs]
    rt_cs_dp_expand_upg = rt_cs_dp_expand[rt_cs_dp_expand['avg_log2FC'] > 0]
    rt_cs_de_module = cs_de_module[cs]
    rt_cs_dp_de_module_upg = rt_cs_dp_expand_upg[rt_cs_dp_expand_upg['genes'].isin(rt_cs_de_module['gene'])]
    rt_cs_dp_de_module_rest = rt_cs_dp_expand[rt_cs_dp_expand['genes'].isin(list(set(rt_cs_de_module['gene']) -
                                                                                set(rt_cs_dp_de_module_upg['genes'])))]
    rt_cs_dp_de_module = pd.concat([rt_cs_dp_de_module_upg,rt_cs_dp_de_module_rest])
    de_avgfc = dict(zip(rt_cs_de_module['gene'], rt_cs_de_module['avg_log2FC']))
    rt_cs_dp_de_module['avg_log2FC_scRNA'] = rt_cs_dp_de_module['genes'].map(de_avgfc)
    rt_cs_dp_de_module['sync_score'] = np.sign(rt_cs_dp_de_module['avg_log2FC_scRNA'] * rt_cs_dp_de_module['avg_log2FC']) * \
                                     np.sqrt(np.abs(rt_cs_dp_de_module['avg_log2FC_scRNA']) *
                                             np.abs(rt_cs_dp_de_module['avg_log2FC']))
    cs_dp_de_module[cs] = rt_cs_dp_de_module
    # for region in genomeLocs[:3]:
    #     rt_cs_dp_de_module_upg_region = rt_cs_dp_de_module_upg[rt_cs_dp_de_module_upg['genomeLoc_annot']==region]
    #     rt_cs_dp_de_module_rest_region = rt_cs_dp_de_module_rest[rt_cs_dp_de_module_rest['genomeLoc_annot']==region]
    num_sync_gene = len(rt_cs_dp_de_module_upg['genes'].unique())
    num_nosync_gene =  len(rt_cs_dp_de_module_rest['genes'].unique())
    module_summary_dict['counts'].append(num_sync_gene)
    module_summary_dict['cell_states'].append(cs)
    module_summary_dict['type'].append('Sync')
    module_summary_dict['region'].append(region)

    module_summary_dict['counts'].append(num_nosync_gene)
    module_summary_dict['cell_states'].append(cs)
    module_summary_dict['type'].append('Not.Sync')
    module_summary_dict['region'].append(region)

    module_summary_dict['counts'].append(len(rt_cs_de_module)-num_sync_gene-num_nosync_gene)
    module_summary_dict['cell_states'].append(cs)
    module_summary_dict['type'].append('No.Dpeaks')
    module_summary_dict['region'].append(region)

module_summary_df = pd.DataFrame(module_summary_dict)
# construct data for violin
violindata_list = []
for cs in cs_list:
    current_df = cs_dp_de_module[cs][['type','sync_score']]
    current_df['cell_states'] = cs
    violindata_list.append(current_df)
violindata_df = pd.concat(violindata_list)
fig,axs = plt.subplots(2,1,sharex=True)
sns.barplot(data=module_summary_df,x='cell_states',y='counts',hue='type',ax=axs[0])
sns.move_legend(axs[0], "lower center", bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False)
axs[0].spines[['right', 'top']].set_visible(False)
sns.violinplot(data=violindata_df,x='cell_states',y='sync_score',palette={'RT_CS1':'orangered',
                                                                          'RT_CS2':'lime',
                                                                          'RT_CS3':'cyan'},ax=axs[1])
axs[1].axhline(sync_avg,c='r',ls='--')
axs[1].spines[['right', 'top']].set_visible(False)
plt.tight_layout()
plt.savefig(f'{outdir}/sync_score_rt_cs.png',dpi=300)
plt.close()

# save module
module_dp_de_writer = pd.ExcelWriter(f'{outdir}/Dp_De_module.xlsx')
for cs in cs_list:
    cs_dp_de_module[cs].to_excel(module_dp_de_writer,cs,index=False)
module_dp_de_writer.close()

def plot_clustered_stacked(data, x, y, value, hue, hue_order=None, axe=None, labels, title="",  H="/", **kwargs):

    if hue_order is None:
        hue_order = data[hue].unique()
    else:
        hue_order = hue_order
    n_df = len(hue_order)
    n_col = len(data[y].unique())
    ind = data[x].unique()
    n_ind = len(ind)
    if axe is None:
        axe = plt.subplot(111)
    else:
        axe = axe

    for h in hue_order : # for each data frame
        df = data[data[hue]==h].pivot_table(columns=y, index=x, values=value).reset_index()
        print(df)
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1])
    axe.add_artist(l1)
    return axe


fig, ax = plt.subplots(1,1)
ax = plot_clustered_stacked(data=module_summary_df,x='cell_states',y='type',value='counts',hue='region',
                            hue_order=['Distal','Proximal','Promoter'],axe=ax)
plt.tight_layout()
plt.savefig(f'{outdir}')



