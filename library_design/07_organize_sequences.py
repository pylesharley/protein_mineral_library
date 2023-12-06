#! /usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import subprocess
import argparse
import random
import sys
import os
import re

comp_category_list = ['group1', 'group2', 'group3', 'group4', 'group5', 'group6', 'group7']

comp_category_dict = {  'group1':['NegGre45',    'greNegDEY',  'DEAroMidGr', 'lsNeGreDiv'],
                        'group3':['H40E20D20',  'Val10_Al50', 'His20Cys20', 'H40LIM40' ],
                        'group5':['PosGre45',   'PosMidGr', 'PosMidGrDv', 'posGreDiv'],
                        'group2':['Neg60',  'Neg90',    'NegDEqQN', 'PolyAs'],
                        'group4':['Q65N25',   'S35T35', 'Thr90', 'Thr60'],
                        'group6':['Pos60',  'Pos90',    'PolyLys'],
                        'group7':['Scaffold' ]}

def gs_pad(string, target_len):
    start_len = len(string)
    while len(string) < target_len:
        string += 'GS'
        
    string = string[:target_len]
    
    # if start_len < target_len and string.endswith('S'):
    #     string = string[:-1]
    #     string += 'G'
    
    assert len(string) == target_len, string +' '+ str(target_len)

    return string


def main():
# if 1:
    parser = argparse.ArgumentParser(description='Sorts sequences based on HARDCODED dictionary')
    # parser.add_argument('-pdb', type=str, help='input pdbs', nargs='+', required=True)
    parser.add_argument('-tsv', type=str, help='table merged_output table from 03_extract_output.py OR 04_remove_redundant_sequences.py', required=True)
    parser.add_argument('-barcodes', type=str, help='text file with barcodes to add to N term of sorted designs', default=False)
    # parser.add_argument('-output', '-out', '-o', type=str, help='output name stem (used for .fas and .tsv output)', default='organized')
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep='\t')
    scaffolds = sorted(list(set(df['scaffold'])))
    print (df)
    ALL_REDUNDANTS_TO_REMOVE = []

    lengths = []

    for i, row in df.iterrows():
        length = len(row['sequence'])
        if length > 154:
            row['sequence'] = row['sequence'][:154]
            length = len(row['sequence'])
        lengths.append(length)
        # print (length, row['scaffold'])

    df['seq_len'] = lengths
    df = df.sort_values('seq_len')

    print ('df sorted', df)
    len_set = list(set(df['seq_len']))
    print ('length set', sorted(len_set))

    fig,ax = plt.subplots()
    df['seq_len'].hist(bins=64, color='black')
    plt.axvline(x=93, c='#e6194b')
    plt.axvline(x=130, c='#3cb44b')
    plt.axvline(x=154, c='#4363d8')
    plt.xlabel('Sequence length')
    plt.ylabel('Number of designs')
    plt.title('Library design lengths')
    plt.autoscale()
    plt.tight_layout()
    fig.savefig('Library_design_lengths.png', dpi=450)
    plt.close('all')

    padded_sequences = []
    padded_lengths = []
    for i, row in df.iterrows():
        if row['seq_len'] < 93:
            target = 92
        elif row['seq_len'] < 130:
            target = 129
        else:
            target = 154

        padded_seq = gs_pad(row['sequence'], target)
        padded_sequences.append(padded_seq)
        padded_lengths.append(len(padded_seq))

    df['padded_sequence'] = padded_sequences
    df['padded_len'] = padded_lengths

    fig,ax = plt.subplots()
    df['padded_len'].hist(bins=58, color='black')
    plt.axvline(x=93, c='#e6194b')
    plt.axvline(x=130, c='#3cb44b')
    plt.axvline(x=154, c='#4363d8')
    plt.xlabel('Sequence length')
    plt.ylabel('Number of designs')
    plt.title('Adjusted library sequences lengths')
    plt.autoscale()
    plt.tight_layout()
    fig.savefig('Library_binned_lengths.png', dpi=450)
    print (df.columns)

    size_categories = df['padded_len'].unique()

    all_aa_comps = []
    [all_aa_comps.extend(comp_category_dict[group]) for group in comp_category_list ]
    print (all_aa_comps, 'all_aa_comps')

    for comp in df['aa_comp']:
        assert comp in all_aa_comps, '{0} not found in hardcoded aa_comp groups dictionary'.format(comp)
    
    final_group_info_table = ['chip_group\tlength\taa_comps']
    final_group_ticker = 0 


    number_of_seqs = len(df.index)

    if args.barcodes:
        with open(args.barcodes, 'r') as barcode_file:
            barcode_lines = barcode_file.readlines()
        barcodes = [line.strip() for line in barcode_lines if len(line.strip())]
        assert number_of_seqs < len(barcodes), 'Not enough barcodes for the number of sequences! \n\t# barcodes = {0} \t# sequences = {1}'.format(len(barcodes), number_of_seqs)
        # randomized_barcodes = random.sample(barcodes, number_of_seqs)
        randomized_barcodes = random.sample(barcodes, len(barcodes))

    for size in size_categories:
        df_size = df[df['padded_len'] == size]

        for group in comp_category_list:
            final_group_ticker += 1
            aa_comp_list = comp_category_dict[group]    
            df_size_plus_aa_comp_group = df_size[df_size['aa_comp'].isin(aa_comp_list)]
            print (len(df_size_plus_aa_comp_group), final_group_ticker)
            
            chip_group = 'chip_group_{0}'.format(str(final_group_ticker).rjust(2, '0'))
            if not os.path.isdir(chip_group):
                os.mkdir(chip_group)

            with open("{0}.fas".format(chip_group), 'w') as output_fasta:
                for name, sequence in zip(df_size_plus_aa_comp_group['description'], df_size_plus_aa_comp_group['padded_sequence']):
                    print('>'+name, file=output_fasta)
                    if args.barcodes: barcode = randomized_barcodes.pop()
                    else: barcode = ''
                    print(barcode+sequence, file=output_fasta)            
                    subprocess.check_output(['cp', 'non_redundant/'+name+'.pdb', chip_group+'/'])

            final_group_info_table.append( '\t'.join([chip_group, str(size)]+aa_comp_list) )

    with open("chip_group_info.tsv", 'w') as output_table_file:
        print ('\n'.join(final_group_info_table), file=output_table_file )
    

    if args.barcodes: 
        barcode_output = re.sub(r'(.*)\.(.*)', r'\1_unused.\2', args.barcodes)
        assert barcode_output != args.barcodes
        with open(barcode_output, 'w') as unused_barcode_file:
            print ('\n'.join(randomized_barcodes), file=unused_barcode_file )


if __name__ == "__main__":
    sys.exit(main())
