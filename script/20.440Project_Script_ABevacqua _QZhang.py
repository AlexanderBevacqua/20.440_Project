#!/usr/bin/env python
# coding: utf-8

# In[ ]:


### 20.440 Analysis of Biological Networks:
### Project - Script

# by Alexander Bevacqua and Qianhe Zhang
# May 9, 2022


# In[9]:


### Tumor vs normal tissue comparison ###


# In[4]:


# Import modules
# Make sure to have statsmodels, bioinfokit and pyvenn installed
import csv
import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
import scipy.stats as sp
from statsmodels.stats.multitest import fdrcorrection
from bioinfokit import visuz
import matplotlib.pyplot as plt
from venn import venn
from collections import Counter


# In[5]:


# Make nested dictionaries of selected tumor and tissue types

cancer_DIR = '../data/TCGA_eRNAexp/'
tissue_DIR = '../data/m1_eRNAcsv/'

cancer_types = ['BLCA', 'BRCA', 'KICH', 'KIRC', 'KIRP', 
                'LIHC', 'LUAD', 'LUSC', 'SKCM', 'STAD']
tissue_types = ['Bladder', 'Breast', 'Kidney', 'Liver', 
                'Lung', 'Skin_no_sun', 'Stomach']

cancer_dictionary = {}
tissue_dictionary = {}

for cancer in cancer_types:
    cancer_csv = csv.reader(open(cancer_DIR + cancer + 
                                 '.eRNA_exp.csv', "r"))
    cancer_sampleID = next(cancer_csv) 
    cancer_type_dict = {}
    for line in cancer_csv:
        cancer_type_eRNA = {line[0]: line[1:]}
        cancer_type_dict.update(cancer_type_eRNA)    
    if cancer not in cancer_dictionary:
        cancer_dictionary[cancer] = cancer_type_dict

for tissue in tissue_types:
    tissue_csv = csv.reader(open(tissue_DIR + tissue + '.csv', "r"))
    tissue_sampleID = next(tissue_csv)
    tissue_type_dict = {}
    for line in tissue_csv:
        tissue_type_eRNA = {line[0]: line[1:]}
        tissue_type_dict.update(tissue_type_eRNA)  
    if tissue not in tissue_dictionary:
        tissue_dictionary[tissue] = tissue_type_dict


# In[6]:


# Generate dataframe needed
# Input: dictionary of eRNA exp of a cancer type/tissue type
# Output: dataframe to generate volcano plot (columns: eRNA, log2FC, 
#         p-value, q-valuw, -log10qval, rejected or not)
# Workflow: identify matched eRNA in cancer and normal tissue data;
#           normalize data using Min-max normalization;
#           perform mannwhitneyu to generate p value;
#           calculate log2 fold change of mean expression;
#           perform FDR correction to generate corrected p value
#           generate -log10qval
def volcanodf(cancer_data, tissue_data):
    volcano_dict = {}
    volcano_dict['eRNA'] = []
    volcano_dict['log2FC'] = []
    p_value = []
    for i in cancer_data:
        if i in tissue_data.keys():
            volcano_dict['eRNA'].append(i)
            cancer_exp = np.array(cancer_data[i], dtype='float')
            cancer_exp_norm = minmax_scale(cancer_exp, feature_range=(0,1))
            cancer_mean = np.mean(cancer_exp_norm)
            tissue_exp = np.array(tissue_data[i], dtype='float')
            tissue_exp_norm = minmax_scale(tissue_exp, feature_range=(0,1))
            tissue_mean = np.mean(tissue_exp_norm)
            sig, pvals = sp.mannwhitneyu(cancer_exp_norm, tissue_exp_norm)
            p_value.append(pvals)
            volcano_dict['log2FC'].append(np.log2(cancer_mean/tissue_mean))
    
    rejected, qvals = fdrcorrection(p_value)
    volcano_dict['p_value'] = p_value
    volcano_dict['q_value'] = qvals
    volcano_dict['-log10qval'] = -np.log10(qvals)
    volcano_dict['rejected'] = rejected
    volcano_df = pd.DataFrame.from_dict(volcano_dict)
    volcano_df = volcano_df.sort_values(by='-log10qval', ascending=False)
    volcano_df = volcano_df.reset_index(drop=True)
    
    return volcano_df


# In[7]:


# Write a function to make volcano plot for tumor types that all eRNAs are downregulated
def volscat(df):
    df_false = df.drop(df[df.q_value < 0.05].index)
    df_true = df.drop(df[df.q_value >= 0.05].index)
    df_up = df_true.drop(df_true[df_true.log2FC < 1.5].index)
    df_down = df_true.drop(df_true[df_true.log2FC > -1.5].index)
    df_nochange = df_true.drop(df_true[df_true.log2FC < -1.5].index)
    df_nochange = df_nochange.drop(df_nochange[df_nochange.log2FC > 1.5].index)
    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax1 = fig.add_subplot(111)
    ax1.scatter(df_up.log2FC, df_up['-log10qval'], s=8, 
                color='yellow', alpha=0.75, label='higher')
    ax1.scatter(df_false.log2FC, df_false['-log10qval'], s=8, 
                color='grey', alpha=0.75, label='insignificant')
    ax1.scatter(df_nochange.log2FC, df_nochange['-log10qval'], 
                s=8, alpha=0.75, color='grey')
    ax1.scatter(df_down.log2FC, df_down['-log10qval'], s=8, 
                color='purple', alpha=0.75, label='lower')
    plt.axvline(x=-1.5, ls='--', color='grey')
    plt.axhline(y=(-np.log10(0.05)), ls='--', color='grey')
    plt.xlabel('log2(fold change)')
    plt.ylabel('-log10(corrected p-value)')
    plt.legend(bbox_to_anchor=(1.46, 0.7), markerscale=2)
    return ax1


# In[ ]:


# volcano plot for BLCA
'''
BLCA_df = volcanodf(cancer_dictionary['BLCA'], tissue_dictionary['Bladder'])
display(BLCA_df)
BLCA_df.to_csv('..table/BLCA_df.csv')
visuz.GeneExpression.volcano(df=BLCA_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_BLCA")
'''


# In[ ]:


# volcano plot for BRCA
'''
BRCA_df = volcanodf(cancer_dictionary['BRCA'], tissue_dictionary['Breast'])
display(BRCA_df)
BRCA_df.to_csv('..table/BRCA_df.csv')
visuz.GeneExpression.volcano(df=BRCA_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_BRCA")
'''


# In[10]:


# volcano plot for KICH
KICH_df = volcanodf(cancer_dictionary['KICH'], tissue_dictionary['Kidney'])
display(KICH_df)
KICH_df.to_csv('../table/KICH_df.csv')
'''
visuz.GeneExpression.volcano(df=KICH_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_KICH")
'''

volscat(KICH_df)
plt.savefig('../figure/volcano_KICH.png', bbox_inches='tight')


# In[11]:


# volcano plot for KIRC
KIRC_df = volcanodf(cancer_dictionary['KIRC'], tissue_dictionary['Kidney'])
display(KIRC_df)
KIRC_df.to_csv('../table/KIRC_df.csv')
'''
visuz.GeneExpression.volcano(df=KIRC_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_KIRC")
'''
volscat(KIRC_df)
plt.savefig('../figure/volcano_KIRC.png', bbox_inches='tight')


# In[12]:


# volcano plot for KIPR
KIRP_df = volcanodf(cancer_dictionary['KIRP'], tissue_dictionary['Kidney'])
display(KIRP_df)
KIRP_df.to_csv('../table/KIRP_df.csv')
visuz.GeneExpression.volcano(df=KIRP_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_KIRP")


# In[13]:


# volcano plot for LIHC
LIHC_df = volcanodf(cancer_dictionary['LIHC'], tissue_dictionary['Liver'])
display(LIHC_df)
LIHC_df.to_csv('../table/LIHC_df.csv')
'''
visuz.GeneExpression.volcano(df=LIHC_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_LIHC")
'''
volscat(LIHC_df)
plt.savefig('../figure/volcano_LIHC.png', bbox_inches='tight')


# In[14]:


# volcano plot for LUAD
LUAD_df = volcanodf(cancer_dictionary['LUAD'], tissue_dictionary['Lung'])
display(LUAD_df)
LUAD_df.to_csv('../table/LUAD_df.csv')
visuz.GeneExpression.volcano(df=LUAD_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_LUAD")


# In[15]:


# volcano plot for LUSC
LUSC_df = volcanodf(cancer_dictionary['LUSC'], tissue_dictionary['Lung'])
display(LUSC_df)
LUSC_df.to_csv('../table/LUSC_df.csv')
visuz.GeneExpression.volcano(df=LUSC_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_LUSC")

# geneid='eRNA', genenames=("ENSR00000324472", "ENSR00000224069"), gfont=8, gstyle=2,


# In[16]:


# volcano plot for SKCM
SKCM_df = volcanodf(cancer_dictionary['SKCM'], tissue_dictionary['Skin_no_sun'])
display(SKCM_df)
SKCM_df.to_csv('../table/SKCM_df.csv')
'''
visuz.GeneExpression.volcano(df=SKCM_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_SKCM")
'''

volscat(SKCM_df)
plt.savefig('../figure/volcano_SKCM.png', bbox_inches='tight')


# In[ ]:


# volcano plot for STAD
'''
STAD_df = volcanodf(cancer_dictionary['STAD'], tissue_dictionary['Stomach'])
display(STAD_df)
STAD_df.to_csv('..table/STAD_df.csv')
visuz.GeneExpression.volcano(df=STAD_df, lfc='log2FC', pv='q_value', 
                             lfc_thr=(1.5, 1.5),
                             color=("yellow","grey","purple"), valpha=0.75, 
                             ar=0, sign_line=True,
                             axtickfontsize=10, axlabelfontsize=10,
                             axtickfontname='DejaVu Sans', axlabelfontname='DejaVu Sans',
                             legendlabels=['higher', 'insignificant', 'lower'],
                             axxlabel='log2(fold change)',
                             axylabel='-log10(corrected p-value)',
                             plotlegend=True, 
                             legendanchor=(1.46,0.7), figname="../figure/volcano_STAD")
'''


# In[17]:


# Generate venn plot for all eRNAs

# generate set for all matched eRNAs between tumor and tissue
KICH_set=set(KICH_df['eRNA'].values.tolist())
KIRP_set=set(KIRP_df['eRNA'].values.tolist())
KIRC_set=set(KIRC_df['eRNA'].values.tolist())
LIHC_set=set(LIHC_df['eRNA'].values.tolist())
LUSC_set=set(LUSC_df['eRNA'].values.tolist())
LUAD_set=set(LUAD_df['eRNA'].values.tolist())
SKCM_set=set(SKCM_df['eRNA'].values.tolist())

# venn plot for KICH, KIRC, LIHC, LUAD, SKCM (all eRNAs)
venn({'KICH': KICH_set, 'KIRC': KIRC_set, 'LIHC': LIHC_set, 
      'LUAD': LUAD_set, 'SKCM': SKCM_set})
plt.savefig('../figure/venn_5tumor_alleRNA.png', bbox_inches='tight')


# In[ ]:


# venn plot for KICH, KIRC, KIRP, LUAD, LUSC (all eRNAs)
'''
venn({'KICH': KICH_set, 'KIRC': KIRC_set, 'KIRP': KIRP_set, 
      'LUAD': LUAD_set, 'LUSC': LUSC_set})
plt.savefig('../figure/venn_5tumor2_alleRNA.png', bbox_inches='tight')
'''


# In[18]:


# Generate venn plot for dysregulated eRNA hits

# function to generate set for dysregulated eRNA hits
def hitset(df):
    df_true = df.drop(df[df.q_value >= 0.05].index)
    df_up = df_true.drop(df_true[df_true.log2FC < 1.5].index)
    df_down = df_true.drop(df_true[df_true.log2FC > -1.5].index)
    df_change = df_up.append(df_down, ignore_index=True)
    hitsetdf = set(df_change['eRNA'].values.tolist())
    return hitsetdf

# set of all dysregulated eRNA hits for all cancer types
KICH_hitset=hitset(KICH_df)
KIRP_hitset=hitset(KIRP_df)
KIRC_hitset=hitset(KIRC_df)
LIHC_hitset=hitset(LIHC_df)
LUSC_hitset=hitset(LUSC_df)
LUAD_hitset=hitset(LUAD_df)
SKCM_hitset=hitset(SKCM_df)

# venn plot for KICH, KIRC, LIHC, LUAD, SKCM (dysregulated eRNA hits)
venn({'KICH': KICH_hitset, 'KIRC': KIRC_hitset, 'LIHC': LIHC_hitset, 
      'LUAD': LUAD_hitset, 'SKCM': SKCM_hitset})
plt.savefig('../figure/venn_5tumor_hiteRNA.png', bbox_inches='tight')


# In[ ]:


# venn plot for BLCA, LUAD, LUSC, SKCM, STAD (dysregulated eRNA hits)
'''
venn({'BLCA': BLCA_hitset, 'LUAD': LUAD_hitset, 'LUSC': LUSC_hitset, 
      'SKCM': SKCM_hitset, 'STAD': STAD_hitset})
plt.savefig('../figure/venn_5tumor2_hiteRNA.png', bbox_inches='tight')
'''


# In[ ]:


# venn plot for LUAD, LUSC (dysregulated eRNA hits)
'''
venn({'LUAD': LUAD_hitset, 'LUSC': LUSC_hitset})
plt.savefig('../figure/venn_2tumor_hiteRNA.png', bbox_inches='tight')
'''


# In[ ]:


# venn plot for KICH, KIRC, KIPR (dysregulated eRNA hits)
'''
venn({'KICH': KICH_hitset, 'KIRC': hitset(KIRC_df), 'KIPR': KIRP_hitset})
plt.savefig('../figure/venn_3tumor_hiteRNA.png', bbox_inches='tight')
'''


# In[23]:


# function generate dataframe for dysregulated eRNA hits
def hit(df):
    df_true = df.drop(df[df.q_value >= 0.05].index)
    df_up = df_true.drop(df_true[df_true.log2FC < 1.5].index)
    df_down = df_true.drop(df_true[df_true.log2FC > -1.5].index)
    df_change = df_up.append(df_down, ignore_index=True)
    return df_change

# dataframe of all dysregulated eRNA hits for all cancer types
KICH_hits=hit(KICH_df)
KIRC_hits=hit(KIRC_df)
KIRP_hits=hit(KIRP_df)
LIHC_hits=hit(LIHC_df)
LUAD_hits=hit(LUAD_df)
LUSC_hits=hit(LUSC_df)
SKCM_hits=hit(SKCM_df)

# generate a full list of all the hits
merged_hits = pd.concat([KICH_hits, KIRC_hits, KIRP_hits, LIHC_hits, 
                         LUAD_hits, LUSC_hits, SKCM_hits])
Tissue=(['Kidney']*KICH_hits.shape[0]+['Kidney']*KIRC_hits.shape[0]
        +['Kidney']*KIRP_hits.shape[0]+['Liver']*LIHC_hits.shape[0]
        +['Lung']*LUAD_hits.shape[0]+['Lung']*LUSC_hits.shape[0]
        +['Skin']*SKCM_hits.shape[0])
merged_hits['Tissue']=Tissue
merged_hits.to_csv('../table/merged_hits_in_tumors.csv')
eRNA_hits = merged_hits['eRNA'].values.tolist()
# Count max apperance of all eRNAs
counter = Counter(eRNA_hits).most_common()


# In[36]:


# Count apperance of all eRNAs
a1=0
b2=0
c3=0
d4=0
e5=0
counter_dict = {}
for i in range(len(counter)):
    counter_dict[counter[i][0]] = counter[i][1]
    if counter[i][1] == 1:
        a1=a1+1
    elif counter[i][1] == 2:
        b2=b2+1
    elif counter[i][1] == 3:
        c3=c3+1
    elif counter[i][1] == 4:
        d4=d4+1
    elif counter[i][1] == 5:
        e5=e5+1     
print(a1, b2, c3, d4, e5)
counter_df=pd.DataFrame.from_dict(counter_dict, orient='index',columns=['Counts'])
counter_df.to_csv('../table/counts of eRNA hits in tumors.csv')

# Plot distribution of eRNAs
x=[5,4,3,2,1]
size=[200,160,120,80,40]
y=[e5,d4,c3,b2,a1]
fig, ax = plt.subplots()
bars = ax.barh(x,y,height=0.5,alpha=0.5)
ax.bar_label(bars)
plt.xlabel('Number of eRNAs')
plt.ylabel('Frequency')
plt.title('Distribution of eRNA hits in seven tumor types')
plt.savefig('../figure/Distribution of eRNA hits in tumors')


# In[ ]:


### Human - Mouse comparison ###


# In[ ]:


### Human - Mouse comparison ###
### The tissue types as vectors ###


human_tissue_types = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland',
                      'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 
                      'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 
                      'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 
                      'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
                      'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)',
                      'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
                      'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes',
                      'Cells - Leukemia cell line (CML)', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix',
                      'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
                      'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage',
                      'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland',
                      'Minor Salivary Gland', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate',
                      'Skin - Not Sun Exposed', 'Skin - Sun Exposed', 'Small Intestine - Terminal Ileum',
                      'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']

expression_file_names = ['Adipose_sub', 'Adipose_vis', 'Adrenal', 'Artery_aor', 'Artery_cor', 'Artery_tib', 
                         'Bladder', 'Brain_Amy', 'Brain_Ant', 'Brain_Cau', 'Brain_Cer', 'Brain_Cerebellum', 
                         'Brain_Cortex', 'Brain_Fro', 'Brain_Hippo', 'Brain_Hypo', 'Brain_Nuc', 'Brain_Put', 
                         'Brain_Spi', 'Brain_Sub', 'Breast', 'Cells_ebv', 'Cells_fibroblast', 'Cells_leukemia', 
                         'Cervix_ecto', 'Cervix_endo', 'Colon_sigmoid', 'Colon_transverse', 'Esophagus_gast', 
                         'Esophagus_muc', 'Esophagus_mus', 'Fallopian_tube', 'Heart_Atr', 'Heart_Left', 'Kidney', 
                         'Liver', 'Lung', 'Minor_salivary_gland', 'Muscle_skeletal', 'Nerve', 'Ovary', 'Pancreas', 
                         'Pituitary', 'Prostate', 'Skin_no_sun', 'Skin_sun', 'Small_intestine', 'Spleen', 
                         'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_blood']

mouse_tissue_types = ['Abdominal Aorta', 'Abdominal Muscle', 'Atrium', 
                'Blood', 'Brain', 'Diaphragm', 'Duodenum', 'Epididymal White Adipose', 
                'Extensor Digitorum Longus Muscle', 'Extraocular Eye Muscle', 
                'Flexor Digitorum Brevis Muscle', 'Forebrain', 'Gastrocnemius Muscle',
                'Gonad Adipose', 'Heart', 'Hemisphaere', 'Hindbrain', 'Ileum', 'Kidney',
                'Left Ventricle', 'Liver', 'Lung', 'Masseter Muscle', 'Ovary', 'Plantaris Muscle',
                'Quadricep Muscle', 'Right Ventricle', 'Skin', 'Soleus Muscle', 'Spleen',
                'Testis', 'Thoracic Aorta', 'Tibialis Anterior Muscle', 'Tongue', 'Total Aorta']


# In[25]:


human_files = ['Kidney', 'Liver', 'Lung', 'Skin_no_sun']
mouse = ['Kidney', 'Liver', 'Lung', 'Skin']


# In[26]:


### Creation of a Human eRNA expression dictionary ###

# Runtime = About 6 seconds

import os
import csv
import scipy.stats as sp
DIR = '../data/m1_eRNAcsv/'


# Tissue Types that we are interested in
#human_files = [expression_file_names[3], expression_file_names[-1], expression_file_names[11], expression_file_names[34], 
#               expression_file_names[33], expression_file_names[35], expression_file_names[44], expression_file_names[47]]



# Human: Empty dictionary of dictionaries
human_dictionary = {}

# Going through each CSV file
for tissue_type in human_files:
    
    # Opening the CSV file out of m1_eRNArds/
    csv_file = csv.reader(open(DIR + tissue_type + '.csv', "r"))
    
    # A blank subdictionary
    a_new_dictionary = {}
    
    # Going through each line, which corresponds to an individual eRNA
    for line in csv_file:
        
        # If the line corresponds to an eRNA
        if 'ENSR' in line[0]:
            
            # A new key value pair is added to the subdictionary
            new_human_eRNA = {line[0]: line[1:]}
            a_new_dictionary.update(new_human_eRNA)
    
    # Once the subdictionary is complete, it is added to the grand dictionary
    
    if tissue_type not in human_dictionary:
        human_dictionary[tissue_type] = a_new_dictionary
        
# human_dictionary is a dictionary of dictionaries
    

# To retrieve an entry, 
# a = human_dictionary.get('A_Tissue_Type')
# print(a.get('An_eRNA_ID'))

# Examples:
# a = human_dictionary.get('Spleen')
# print(a.get('ENSR00000000587'))
#
# a = human_dictionary.get('Skin_no_sun')
# print(a.get('ENSR00000000321'))
# You get a list of all the recorded expression data


# In[27]:


### Creation of a Mouse eRNA expression dictionary ###

# Runtime = About 6.5 seconds

# Tissue Types that we are interested in
# mouse = [mouse_tissue_types[0], mouse_tissue_types[3], mouse_tissue_types[4], mouse_tissue_types[18], 
#         mouse_tissue_types[19], mouse_tissue_types[20], mouse_tissue_types[26], mouse_tissue_types[28]]

# Opening the singular mouse file
csv_file2 = csv.reader(open('../data/mouse.eRNA.csv', "r"))

# Initial conditions
mouse_dictionary = {}
a_new_dictionary = {}
previous_tissue_type = 'Abdominal Aorta'

# Going through each line, which corresponds to an individual eRNA in a certain tissue
for line in csv_file2:
    
    # The tissue type is here
    tissue_type = line[3]
    
    # Ignore header line and skip line if the tissue type is not in mouse
    if tissue_type == 'tissue' or tissue_type not in mouse:
        continue
    
    # If the previous tissue type is not the current tissue type, then 
    # the current dictionary is assigned to the grand mouse dictionary
    if previous_tissue_type != tissue_type:
        mouse_dictionary[previous_tissue_type] = a_new_dictionary
        a_new_dictionary = {}
    
    # Creating a subdictionary
    if tissue_type in mouse and 'chr' in line[4]:
        new_mouse_eRNA = {line[4]: line[5:8]}
        a_new_dictionary.update(new_mouse_eRNA)
        
    # Updating the previous tissue type to the current tissue type
    previous_tissue_type = tissue_type


mouse_dictionary[previous_tissue_type] = a_new_dictionary

# To retrieve an entry, 
# a = human_dictionary.get('A_Tissue_Type')
# print(a.get('An_eRNA_ID'))

# Examples:
# a = mouse_dictionary.get('Liver')
# print(a.get('chr7:105876398-105882398'))
#
# You get a list of all the recorded expression data


# In[28]:


### Creation of a Dataframe containing human and non-human eRNA hits based on Sequence Similarity ###

# Runtime varies depending on cutoff selected

import pandas as pd

# Opening and reading the M8 file

# (1,2) identifiers for query and target sequences/profiles, 
# (3) sequence identity, (4) alignment length, 
# (5) number of mismatches, (6) number of gap openings, 
# (7-8, 9-10) domain start and end-position in query and in target, 
# (11) E-value, and (12) bit score.

columns = ['Human', 'Mouse', 'Sequence Identity', 'Alignment Length']
df = pd.DataFrame(columns=columns)


# Crucial number!!! The sequence identity is used for sequence alignment cutoff
cutoff_sequence_identity = 0.975

f = open("../data/RESULTS_v1.m8", "r")
# 16388620 entries to go through

for line in f:
    
    entries = line.split()
    
    if float(entries[2]) >= cutoff_sequence_identity:
        df.loc[len(df.index),:] = entries[0:4]

f.close()
df = df[df.Mouse != 'chr16:9-95006001']
df = df.drop_duplicates(subset = ['Human'], ignore_index = True)
display(df)


# In[29]:


### Bringing all of the pieces together ###
### Accessing sequence alignment hits, collecting expression data necessary for statistical testing ###
### Performing a Mann–Whitney U test ###

# Runtime varies depending on the number of hits

with open('../table/Results_V3_0.975.csv', 'w') as csvfile:
    fieldnames = ['Human eRNA', 'Mouse eRNA', 'Tissue', 'p-value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writerow({'Human eRNA': 'Human eRNA', 'Mouse eRNA': 'Mouse eRNA', 
                     'Tissue': 'Tissue', 'p-value': 'p-value'})

    for i in range(len(df)):

        entry = df.loc[i, :].values.tolist()

        # retrieve human eRNA first from the dataframe
        if 'ENSR' in entry[0]:
            human_eRNA = entry[0].partition('_')[0]
            mouse_eRNA = entry[1]

            # Extract human eRNA expression data from dictionary given an eRNA ID
            for k in range(len(human_files)):
                a = human_dictionary.get(human_files[k])

                ### Human expression data for statistical test ###
                human_expression_data = a.get(human_eRNA)

                if human_expression_data == None:
                    continue

                else: # Continue on the mouse path
                    mouse_enhancer = entry[1].partition(':')
                    chromosome = mouse_enhancer[0]
                    mouse_enhancer_reg = mouse_enhancer[2].partition('-')
                    #print([chromosome, mouse_enhancer_reg[0], mouse_enhancer_reg[2]])
                    mouse_range_enh = range(int(mouse_enhancer_reg[0]), int(mouse_enhancer_reg[2]))

                    b = mouse_dictionary.get(mouse[k])
                    mouse_expression_data = []

                    for key, values in b.items():

                        if chromosome == values[0]:

                            middle = int(values[1])

                            if middle in mouse_range_enh:
                                # print(human_files[k], values)
                                mouse_expression_data.append(values[-1])

                    if not mouse_expression_data:
                        continue
                    
                    import numpy as np
                    human_expression_data = np.array(human_expression_data, dtype=float)
                    mouse_expression_data = np.array(mouse_expression_data, dtype=float)
                    
                    ### Mouse expression data for statistical test ###
                    # print([chromosome, mouse_enhancer_reg[0], mouse_enhancer_reg[2]], mouse_expression_data)

                    U, p_value = sp.mannwhitneyu(human_expression_data,mouse_expression_data)
                    writer.writerow({'Human eRNA': human_eRNA, 'Mouse eRNA': mouse_eRNA, 
                                     'Tissue': human_files[k], 'p-value': p_value})

                    print(f'Human: {human_eRNA}, Mouse: {mouse_eRNA}, Tissue: {human_files[k]}') 
                    print(f'Significance: {sp.mannwhitneyu(human_expression_data,mouse_expression_data)}')
                    print('-----------------------------------------------------------------------------')
                


# In[30]:


### Copy of the Above Code Block, slightly modified, ###
### for Opening the animalDB authors' alignment data ###

# Runtime varies depending on cutoff selected

columns = ['Mouse', 'Human', 'Sequence Identity', 'Alignment Length']
df_authors = pd.DataFrame(columns=columns)


# Crucial number!!! The sequence identity is used for sequence alignment cutoff
cutoff_sequence_identity = 97.5

f = open("../data/mouse2human.txt", "r")
i = 0
for line in f:
    
    entries = line.split()
    
    if float(entries[2]) >= cutoff_sequence_identity:
        df_authors.loc[len(df_authors.index),:] = entries[0:4]
        
    #if (i % 100000 == 0):
        # print(i)

    i += 1
        
f.close()
print(i)
df_authors = df_authors.drop_duplicates(subset = ['Human'], ignore_index = True)
display(df_authors)


# In[31]:


### Bringing all of the pieces together ###
### Accessing sequence alignment hits, collecting expression data necessary for statistical testing ###
### Performing a Mann–Whitney U test ###

# Runtime varies depending on the number of hits

with open('../table/Results_V4_0.975Authors.csv', 'w') as csvfile2:
    
    fieldnames = ['Human eRNA', 'Mouse eRNA', 'Tissue', 'p-value']
    writer = csv.DictWriter(csvfile2, fieldnames=fieldnames)
    writer.writerow({'Human eRNA': 'Human eRNA', 'Mouse eRNA': 'Mouse eRNA', 
                     'Tissue': 'Tissue', 'p-value': 'p-value'})

    for i in range(len(df_authors)):

        entry = df_authors.loc[i, :].values.tolist()

        # retrieve human eRNA first from the dataframe
        if 'ENSR' in entry[1]:
            human_eRNA = entry[1].partition('_')[0]
            mouse_eRNA = entry[0]

            # Extract human eRNA expression data from dictionary given an eRNA ID
            for k in range(len(human_files)):
                a = human_dictionary.get(human_files[k])

                ### Human expression data for statistical test ###
                human_expression_data = a.get(human_eRNA)

                if human_expression_data == None:
                    continue

                else: # Continue on the mouse path
                    mouse_enhancer = entry[0].partition(':')
                    chromosome = mouse_enhancer[0]
                    mouse_enhancer_reg = mouse_enhancer[2].partition('-')
                    #print([chromosome, mouse_enhancer_reg[0], mouse_enhancer_reg[2]])
                    mouse_range_enh = range(int(mouse_enhancer_reg[0]), int(mouse_enhancer_reg[2]))

                    b = mouse_dictionary.get(mouse[k])
                    mouse_expression_data = []

                    for key, values in b.items():

                        if chromosome == values[0]:

                            middle = int(values[1])

                            if middle in mouse_range_enh:
                                # print(human_files[k], values)
                                mouse_expression_data.append(values[-1])

                    if not mouse_expression_data:
                        continue

                    human_expression_data = np.array(human_expression_data, dtype=float)
                    mouse_expression_data = np.array(mouse_expression_data, dtype=float)
                    
                    ### Mouse expression data for statistical test ###
                    # print([chromosome, mouse_enhancer_reg[0], mouse_enhancer_reg[2]], mouse_expression_data)

                    U, p_value = sp.mannwhitneyu(human_expression_data,mouse_expression_data)
                    writer.writerow({'Human eRNA': human_eRNA, 'Mouse eRNA': mouse_eRNA, 
                                     'Tissue': human_files[k], 'p-value': p_value})

                    print(f'Human: {human_eRNA}, Mouse: {mouse_eRNA}, Tissue: {human_files[k]}') 
                    print(f'Significance: {sp.mannwhitneyu(human_expression_data,mouse_expression_data)}')
                    print('-----------------------------------------------------------------------------')


# In[40]:


### Merging the benign-tumor and cross-species hits together ###

# from the Benjamini-Hochberg Procedure
df_benj_mouse_human = pd.read_csv('../table/Results_V4_0.975Authors.csv')
new_alpha=0.01768969
df_benj_mouse_human = df_benj_mouse_human.drop(df_benj_mouse_human[df_benj_mouse_human['p-value']<new_alpha].index)
benj_mouse_human = df_benj_mouse_human['Human eRNA'].tolist()
tissue_benj_mouse_human = df_benj_mouse_human['Tissue'].tolist()

new_list = [a + '_' + b for a, b in zip(benj_mouse_human, tissue_benj_mouse_human)]


# from the eRic hits
df_cancer = pd.read_csv('../table/merged_hits_in_tumors.csv')
cancer = df_cancer['eRNA'].tolist()
tissue_cancer = df_cancer['Tissue'].tolist()

new_list2 = [a + '_' + b for a, b in zip(cancer, tissue_cancer)]

# The final hits at the intersection of these two branches, tissue-specific as well
final_eRNA_hits = [x for x in new_list if x in new_list2]
print(final_eRNA_hits)
print(len(final_eRNA_hits))

# Counting the number of final hits present in each tissue type

# ['Kidney', 'Liver', 'Lung', 'Skin']
count = [0, 0, 0, 0]
for entry in final_eRNA_hits:
    
    if 'Kidney' in entry:
        count[0] += 1
    
    if 'Liver' in entry:
        count[1] += 1
        
    if 'Lung' in entry:
        count[2] += 1
    
    if 'Skin' in entry:
        count[3] += 1
    
print(count)


# In[ ]:




