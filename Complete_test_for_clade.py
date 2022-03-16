#!/usr/bin/env python
import sys, getopt, os


def usage():
    print('''
    use to find the sepecific orthologs in each clade
    -t adjusted tree file with node names
    -T targeted clades to find specific orthologs, clades should be separated by comma, targets containing more than one clade will be merged by plus
    -m matrix data with completeness
    -s species_id.txt
    -o output
    ''')



import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from scipy import stats
from scipy import mean
import scipy.stats as ss
import statsmodels.api as sa
import scikit_posthocs as sp





opts,args=getopt.getopt(sys.argv[1:],'-h-T:-t:-s:-m:-o:',['help'])
for k, v in opts:
    if k in ("-h", "--help"):
        usage()
        sys.exit()
    elif k in ("-t"):
        tree = v
    elif k in ("-s"):
        spemap = v
    elif k in ("-m"):
        matrix = v
    elif k in ("-T"):
        target = v
    elif k in ("-o"):
        output = v
    else:
        usage()
        sys.exit()

if len(opts)>0:
    pass
else:
    usage()
    sys.exit()

from ete3 import EvolTree

def get_leaves(target=target,tree=tree):
    
    t = EvolTree(tree,format=1)
    node = t&target
    TarLeafList = node.get_leaf_names() 
    return(TarLeafList)

def get_interlist(TarLeafList, mapper):
    #seqlist = speRep.keys()
    if not(set(mapper) & set(TarLeafList)):
        return(False)
    else:
        return(list(set(mapper) & set(TarLeafList)))


def id_mapper(SpeMap):
    mapper = {}
    name = {}
    IN=open(SpeMap)
    for lines in IN.readlines():
        ID = lines.split()[0].replace(":","")
        SeqID = lines.split()[1].replace(".faa",".kegg.m8")
        mapper[SeqID] = "T_"+ID
        TID = "T_" + ID
        name[TID] = SeqID
    return(mapper, name)

IDmap, nameMap = id_mapper(spemap)



def get_targets(string=target):
    target.strip()
    taxa_dict = {}
    tar = string.split(",")
    cladeName = {}
    for clade in tar:
        
        lists = []
        if re.search("\+",clade):
            for item in clade.split("+"):
                Leaves = get_leaves(item)
                lists.extend(Leaves) # = list(set(get_leaves(item)) | set(lists))
                for tips in Leaves:
                    cladeName[tips] = clade
        else:
            lists = get_leaves(clade)  
            for tips in lists:
                cladeName[tips] = clade

        taxa_dict[clade] = get_interlist(lists, IDmap.values())
    return(taxa_dict, cladeName)


def get_diffList(TarLeafList, speRep):
    seqlist = speRep.keys()
    if not(set(seqlist) & set(TarLeafList)):
        return(False)
    else:
        return(list(set(seqlist) - set(TarLeafList)))

def data_extract(row, Lists):
    sel = list(row[lists])
    return(sel)


df = pd.read_table(matrix,index_col=0,sep="\t")

df.rename(columns=IDmap, inplace=True)


taxonList, cladeName = get_targets(target)


df =df[cladeName.keys()]

import statsmodels.stats.multitest as smt
Name_final = {}
for item in cladeName:
    printName = nameMap[item] + "_" + cladeName[item]
    Name_final[item] = printName

def anova_test(row,taxonList=taxonList):
    data = {}
    for tar in taxonList:
        data[tar] = list(row[taxonList[tar]])
    try:
        p = ss.kruskal(*data.values())[1]
    except ValueError:
        p = 1
    #comparison = MultiComparison(df_melt['value'], df_melt['treatment'])
    #tukey = comparison.tukeyhsd(0.05)
    return(p)

def posthoc(row,taxonList=taxonList):
    data = {}
    for tar in taxonList:
        data[tar] = row[taxonList[tar]]
    df_anova = pd.DataFrame.from_dict(data, orient='index')
    df_anova = df_anova.transpose()
    df_melt = pd.melt(df_anova).dropna()
    df_out = sp.posthoc_conover(df_melt, val_col='value', group_col='variable', p_adjust = 'holm')
    return(df_out)


df["anova_p"] = df.apply(anova_test,axis=1)
df["anova_p"][pd.isnull(df["anova_p"])] = 1

df["anova_padj"] = smt.multipletests(df["anova_p"],method="fdr_bh",alpha=0.1)[1]


df_bool = df > 0.75


df_sel = df.loc[(df["anova_padj"]<=0.05) & (df_bool.sum(axis=1) > 0)]
df_sel.rename(columns=Name_final,inplace=True)
df_sel.to_csv(f"{output}.tsv", sep="\t")

for index,row in  df_sel.iterrows(): 
    row = df.loc[index]
    row = row[0:-2]
    df_posthoc = posthoc(row)
    index = index[0:6]
    df_posthoc.to_csv(f"{index}.tsv", sep="\t")


heatDF = df_sel.iloc[:,:-2]
import seaborn as sns

fontsize = 9
plt.rcParams['ytick.labelsize'] = fontsize

fontsize_pt = plt.rcParams['ytick.labelsize']

dpi = 72.27
matrix_width_pt = fontsize_pt * (heatDF.shape[1]+80)
matrix_width_in = float(matrix_width_pt) / dpi
left_margin = 0.05
right_margin = 0.15
figure_width = matrix_width_in / (1 - left_margin - right_margin)
matrix_height_pt = fontsize_pt * heatDF.shape[0]
matrix_height_in = float(matrix_height_pt) / dpi
top_margin = 0.1  # in percentage of the figure height
bottom_margin = 0.15 # in percentage of the figure height
figure_height = matrix_height_in / (1 - top_margin - bottom_margin)


maps = sns.clustermap(data=heatDF,
               row_cluster=True, #行方向不聚类
               col_cluster=True, #列方向聚类
               cmap = "vlag",

               xticklabels = True,
               yticklabels = True,
               figsize = (figure_width,(figure_height)),
              )

maps.ax_heatmap.set_xticklabels(maps.ax_heatmap.get_xmajorticklabels(), fontsize = fontsize)
maps.ax_heatmap.set_yticklabels(maps.ax_heatmap.get_ymajorticklabels(), fontsize = fontsize)
maps.savefig(f"{output}.pdf")



