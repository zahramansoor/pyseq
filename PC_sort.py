# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 19:11:39 2021

@author: zahhr

Script to sort top genes expressed in P21 in Purkinje cells
"""

import pandas as pd, os

#path to file
pth = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Purkinje_cell_gene_list\Purkinje_RNA-seq_combined_ver6.xlsx"
dst = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Purkinje_cell_gene_list"
df = pd.read_excel(pth)

#sort by p28 (might just be mislabelled p21)
df_sort = df.sort_values("P28_ave", ascending=False)
#only get genes with nonzero RNA counts
df_sort = df_sort[(df_sort["P28_ave"]>0) & (df_sort["P14_ave"]>0) 
                  & (df_sort["P8_ave"]>0) & (df_sort["PO_ave"]>0)
                  & (df_sort["P4_ave"]>0)]
#get top 500 genes
df_top500 = df_sort[:500]
#names of top 500 genes
top500 = df_top500["Gene_symbol"].values
#export for blasting to mgi database
df_top500["Gene_symbol"].to_csv(os.path.join(dst, r"top500_PC_P28_genes_Cliffordetal_2019.csv"), index=None)

#get top 500 for P21 PCs from https://www.omicsdi.org/dataset/geo/GSE86824 (Clifford et al., 2019, Nature)
#dowloaded raw files from https://www.omicsdi.org/dataset/geo/GSE86824
#in triplicate
fls = [os.path.join(dst, r"GSE86824_RAW\GSM2309055_P21_R1.txt"),
       os.path.join(dst, r"GSE86824_RAW\GSM2309056_P21_R2.txt"),
       os.path.join(dst, r"GSE86824_RAW\GSM2309057_P21_R3.txt")]
dfs = []
#concatenate all replicates and find ones only that exist in all 3?
for fl in fls:
    df = pd.read_csv(fl, sep="	", header=None)
    df.columns = ["geneID", "RNAcounts"]
    dfs.append(df)
bigdf = pd.concat(dfs)
#get all unique geneIDs
geneIDs = bigdf.geneID.unique()
#collect geneIDs common to only all the replicates
replicate_geneIDs = [gID for gID in geneIDs if len(bigdf[bigdf.geneID == gID]) == 3]
#now only get counts for the genes that are in all 3 replicates in the bigdf
replicate_bigdf = bigdf[bigdf.geneID.isin(replicate_geneIDs)]
rep_bigdf_mean = replicate_bigdf.groupby("geneID", as_index=False).mean().round(2)
#now sort by highest expressing genes
rep_bigdf_mean_sort = rep_bigdf_mean.sort_values("RNAcounts", ascending=False)
#remove first few 'no_feature', 'ambigous' rows
rep_bigdf_mean_sort=rep_bigdf_mean_sort[2:]
top500 = rep_bigdf_mean_sort.geneID[:500]
#export ensembl ids
rep_bigdf_mean_sort[rep_bigdf_mean_sort["RNAcounts"]>1].to_csv(os.path.join(dst, r"RNAcounts1plus_PC_P21_genes_Cliffordetal_2019.csv"), index=None)
rep_bigdf_mean_sort[:500].to_csv(os.path.join(dst, r"top500_PC_P21_genes_Cliffordetal_2019.csv"), index=None)
#converted to gene symbol IDs and added column for gene symbol
#%%
#find genes with GO pathways related to tubulin, mitochondrial transport, and axonal transport, etc. 
#run list of genes from Corbo's dataframe in MGI browser to get GO annotations
pth = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\MGI_analysis\MGIBatchReport_20210118_124556_Purkinje_RNA-seq_combined_ver6.xlsx"
dst = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021"
df = pd.read_excel(pth)
#convert df to strings
df = df.applymap(str)
#first, combine all the GO term columns
df_t = df.groupby(["Input", "Input Type", "MGI Gene/Marker ID", "Symbol", "Name",
       "Feature Type"])["Term"].apply(",".join).reset_index()
df_go = df.groupby(["Input", "Input Type", "MGI Gene/Marker ID", "Symbol", "Name",
       "Feature Type"])["GO ID"].apply(",".join).reset_index()
#add go id column to dataframe
df_t["GO ID"] = df_go["GO ID"]

#search terms
terms = ["tubulin","posttranslational","motility","axonal transport","transport along",
         "microtubule","kinesin","dynein","cili","mitochondrion transport"]
dfs_s = []
for term in terms:
    df_s = df_t[["Term"]].apply(lambda x: x.str.contains(term))
    #make boolean mask of original dataframe
    m = df.astype(bool)
    m["Term"] = df_s
    #get converted df
    df_s = df_t[m]
    df_s = df_s.dropna(0)
    #rank based on search terms
    if term == "tubulin" or term == "microtubule" or term == "mitochondrion transport":
        df_s["Rank"] = [3]*len(df_s)
    elif term == "axonal transport" or term == "transport along":
        df_s["Rank"] = [2]*len(df_s)
    else:
        df_s["Rank"] = [1]*len(df_s)
    dfs_s.append(df_s)
#concat
bigdf_s = pd.concat(dfs_s)

#add RNA expression data from original dataframe as reference
src = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Purkinje_cell_gene_list\Purkinje_RNA-seq_combined_ver6.xlsx"
orgdf = pd.read_excel(src)
for gene in bigdf_s.Input.values:
    if gene in orgdf.Gene_symbol.values:
        #add exp data for all time points
        bigdf_s.loc[bigdf_s.Input == gene, "RNA_exp_ave_P21"] = orgdf.loc[orgdf.Gene_symbol == gene, "P28_ave"].values[0]
        bigdf_s.loc[bigdf_s.Input == gene, "RNA_exp_ave_P0"] = orgdf.loc[orgdf.Gene_symbol == gene, "PO_ave"].values[0]
        bigdf_s.loc[bigdf_s.Input == gene, "RNA_exp_ave_P4"] = orgdf.loc[orgdf.Gene_symbol == gene, "P4_ave"].values[0]
        bigdf_s.loc[bigdf_s.Input == gene, "RNA_exp_ave_P8"] = orgdf.loc[orgdf.Gene_symbol == gene, "P8_ave"].values[0]
        bigdf_s.loc[bigdf_s.Input == gene, "RNA_exp_ave_P14"] = orgdf.loc[orgdf.Gene_symbol == gene, "P14_ave"].values[0]

#sort by 1) RNA exp at p21, 2) rank
sortdf = bigdf_s.sort_values(by=["RNA_exp_ave_P21"], ascending=False)
sortdf = sortdf[(sortdf.Rank > 1) | ((sortdf.Rank == 1) & (sortdf.RNA_exp_ave_P21 > 200))]
#drop duplicate gene symbols
sortdf = sortdf.drop_duplicates(subset=["Symbol"])
#save
sortdf.to_excel(os.path.join(dst, r"Purkinje_cell_gene_list\Purkinje_RNA-seq_combined_ver6_matching_GO_search_terms.xlsx"),index=None)

#%%
