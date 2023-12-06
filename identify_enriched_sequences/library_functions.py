#! /usr/bin/env python
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from math import log
import pandas as pd
import numpy as np
import subprocess
import seaborn
import glob
import re

def get_seq_param(df, seq):
    df['Net charge'] = [count_net_charge(seq) for seq in df['rep_CDS']]
    df['Isoelectric point'] = [IsoelectricPoint(seq).pi() for seq in df['rep_CDS']]
    df['length'] = [len(seq) for seq in df['rep_CDS']]
    #print ([seq.replace('X','').replace('*','') for seq in df['rep_CDS']])
    df['GRAVY'] = [ProteinAnalysis(seq.replace('X','').replace('*','')).gravy() for seq in df['rep_CDS']]
    df['% His'] = [float(seq.count('H'))/len(seq)*100.0 for seq in df['rep_CDS']]
    df['% Thr'] = [float(seq.count('T'))/len(seq)*100.0 for seq in df['rep_CDS']]
    df['% Arg'] = [float(seq.count('R'))/len(seq)*100.0 for seq in df['rep_CDS']]
    df['% Asp'] = [float(seq.count('D'))/len(seq)*100.0 for seq in df['rep_CDS']]
    df['% Cys'] = [float(seq.count('C'))/len(seq)*100.0 for seq in df['rep_CDS']]
    return df

def count_net_charge(seq, pH=8):
    return ProteinAnalysis(seq).charge_at_pH(pH)

def name_to_scaffold(name):
    # print (scaffold)
    scaffold = re.sub(r'^(.*?)_[^_]+$', r'\1', name.split('_surf')[0])
    scaffold = scaffold.replace('_Val10', '')
    scaffold = scaffold.replace('_Scaffold','')
    scaffold = scaffold.split('_seq')[0]
    # print (scaffold)
    return scaffold

def name_to_surface(name):
    surface = re.sub(r'^.*?_([^_]+)$', r'\1', name.split('_surf')[0])
    surface = surface.replace('Al50', 'Val10Al50')
    surface = surface.replace('0001', 'scaffold')
    return surface 

def shuffle_dataframe_for_plotting(df, legend='subpool', select=False):
    legend_sorter = []

    if select:
        elements = sorted(set(df[legend]))
        skip = [ele for ele in elements if ele != select]
        shift = len(df[legend])
        selections = []
        
        if type(select) == str:
            for i, l in enumerate(df[legend]):
                if l == select:
                    legend_sorter.append(i+shift)
                    selections.append(select)
                else:
                    legend_sorter.append(i) 
                    selections.append('not selected')
        elif type(select) == list:
            for i, l in enumerate(df[legend]):
                if l in select:
                    legend_sorter.append(i+shift)
                    selections.append(l)
                else:
                    legend_sorter.append(i) 
                    selections.append('not selected')            
            
        df['select'] = selections
        df['legend_sort'] = legend_sorter
        df.sort_values(by='legend_sort', inplace=True)
        return df
        
    elif legend == 'subpool':
        order = ['SP1', 'SP2', 'SP3', 'SP4', 'SP5', 'SP6', 'SP7', 'SP8', 'SP9', 'SP10', 'SP11', 'SP12', 'SP13', 'SP14', 'SP15', 'SP16']
        skip = []
    else:
        order = sorted(set(df[legend]))
        skip = []
    
    df = df.sample(frac=1)
    seen = []
    for i, l in enumerate(df[legend]):
        if l in seen:
            legend_sorter.append(i+len(order))
        else:
            seen.append(l)
            legend_sorter.append(order.index(l))
            
    df['legend_sort'] = legend_sorter
    df.sort_values(by='legend_sort', inplace=True)
    return df
