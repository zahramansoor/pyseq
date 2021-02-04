# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 23:30:12 2021

@author: zahhr
"""

import pandas as pd,os

src = r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Purkinje_cell_gene_list"
oligo_lib = pd.read_excel(os.path.join(src,"STable 22 Brie.xlsx"))
gene_list = pd.read_excel(os.path.join(src,"Purkinje_RNA-seq_combined_ver6_matching_GO_search_terms.xlsx"))

oligos = oligo_lib[oligo_lib["Target Gene Symbol"].isin(gene_list.Symbol.values)]
repeating_oligos = list(oligos["Target Gene Symbol"].values)
oligos["Target Gene Symbol"] = list(map(lambda x: x[1] + "-" + str(repeating_oligos[:x[0]].count(x[1]) + 1) if repeating_oligos.count(x[1]) > 1 else x[1], enumerate(repeating_oligos)))
#read in flanks df
flanks = pd.read_excel(os.path.join(src,"Otx2GFP_screen.xlsx"), sheet_name="HiFidelity clone")
oligos["HiFidelity 5' homology"] = [flanks["5' homology"].unique()[0]]*len(oligos)
oligos["HiFidelity 3' homology"] = [flanks["3 ' homology"].unique()[0]]*len(oligos)
#get oligo sequence and base pair #
for gene in oligos["Target Gene Symbol"].values:
    oligos.loc[oligos["Target Gene Symbol"]==gene, "oligo"] = oligos.loc[oligos["Target Gene Symbol"]==gene, 
                    "HiFidelity 5' homology"]+oligos.loc[oligos["Target Gene Symbol"]==gene,
                    "sgRNA Target Sequence"]+oligos.loc[oligos["Target Gene Symbol"]==gene,"HiFidelity 3' homology"]
    oligos.loc[oligos["Target Gene Symbol"]==gene, "total_bp"] = len(oligos.loc[oligos["Target Gene Symbol"]==gene, "oligo"].values[0])
    
oligos.to_csv(os.path.join(src,"Purkinje_RNA-seq_combined_ver6_oligo_library.csv"),index=None)

#find ones that did not match gene lsit
genes = gene_list.Symbol.values
unmatched = [xx for xx in genes if xx not in oligo_lib["Target Gene Symbol"].values]
pd.DataFrame(unmatched).to_csv(os.path.join(src,"Purkinje_RNA-seq_combined_ver6_unmatched.csv"),index=None)
