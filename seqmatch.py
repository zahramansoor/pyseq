# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 12:38:54 2021

@author: zahhr
"""

import pandas as pd

#extract sequences 
f = open(r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Otx2-GFP_gene_list\Corb0_DMZD_1-26-21_SIC_index_0889_GACCAGAGGA_S144_L001_R1_001_seq.txt")

seq = [l[25:45] for l in f]

#original oligos w reverse complement
df = pd.read_csv(r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Otx2-GFP_gene_list\OligosTest2-rc.csv")

#get matches and # of matches
matches = [xx for xx in df.revComp if xx in seq]
nonmatches = [xx for xx in seq if xx not in df.revComp.values]

for xx in df.revComp.values:
    df.loc[df.revComp==xx, "seq coverage"] = list(seq).count(xx)
#export   
df.to_csv(r"C:\Users\zahhr\OneDrive - Washington University in St. Louis\corbo_lab_spring2021\Otx2-GFP_gene_list\OligosTest2-rc_with_seq_coverage.csv",
          index=None)
