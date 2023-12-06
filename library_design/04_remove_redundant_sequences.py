#! /usr/bin/env python

import pandas as pd
import argparse
import subprocess
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def main():
# if 1:
    parser = argparse.ArgumentParser(description='run like:\n\t04_remove_redundant_sequences.py -tsv merged_output.tsv -threshold 5 -output output_fasta.fas')
    # parser.add_argument('-pdb', type=str, help='input pdbs', nargs='+', required=True)
    parser.add_argument('-tsv', type=str, help='table merged_output table from 03_extract_output.py', default='merged_output.tsv')
    parser.add_argument('-threshold', type=int, help='threshold number of differences', default=5)
    parser.add_argument('-output', '-out', '-o', type=str, help='output name stem (used for .fas and .tsv output)', default='non_redundant_output')
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep='\t')
    scaffolds = sorted(list(set(df['scaffold'])))

    print (df)
    ALL_REDUNDANTS_TO_REMOVE = []

    for scaffold in scaffolds:
        df_subset = df.loc[df['scaffold'] == scaffold]
        # df_subset = df_subset.reset_index(drop=True)
        print (df_subset)
        # for i, row1 in df_subset.iterrows():
        #     print (index)
        all_v_all_differences = []
        redundants_to_remove = []

        number_of_sequences = len(df_subset['sequence'])
        for i, row1 in df_subset.iterrows():
            for j, row2 in df_subset.iterrows():
                if i == j:
                    continue
                seq1 = df_subset['sequence'][i]
                seq2 = df_subset['sequence'][j]
                count = sum(1 for a, b in zip(seq1, seq2) if a != b) + abs(len(seq1) - len(seq2))
                all_v_all_differences.append(count)

                if count < args.threshold:
                    if not i in redundants_to_remove:
                        if not j in redundants_to_remove:
                            redundants_to_remove.append(j)
        if len(all_v_all_differences):
            all_v_all_differences.sort()
            print (all_v_all_differences[:20])
            print (all_v_all_differences[-20:])
            print (len(all_v_all_differences))
            print (number_of_sequences)
            print (len(redundants_to_remove))
            print ()
            fig,ax = plt.subplots()
            n, bins, patches = plt.hist(all_v_all_differences, max(all_v_all_differences), density=False)

            plt.xlabel('Number of sequence differences')
            plt.ylabel('Number of sequence pairs')
            plt.title('{0} based design sequence diversity'.format(scaffold))
            plt.axvline(x=args.threshold, c='r')

            plt.grid(True)
            fig.savefig('{0}_sequence_differences.png'.format(scaffold), dpi=450)
            plt.close('all')

            ALL_REDUNDANTS_TO_REMOVE.extend(redundants_to_remove)
            # break 
        
    print (len(ALL_REDUNDANTS_TO_REMOVE))
    print (len(list(set(ALL_REDUNDANTS_TO_REMOVE))))
    print (ALL_REDUNDANTS_TO_REMOVE[:100])

    non_redundant_df = df.drop(df.index[ALL_REDUNDANTS_TO_REMOVE])
    non_redundant_df.to_csv("{0}.tsv".format(args.output), sep="\t", index=False)

    with open("{0}.fas".format(args.output), 'w') as output_fasta:
        for name, sequence in zip(non_redundant_df['description'], non_redundant_df['sequence']):
            print('>'+name, file=output_fasta)
            print(sequence, file=output_fasta)

    if not os.path.isdir('non_redundant'):
        os.mkdir('non_redundant')

    for name in non_redundant_df['description']:
        try:
            subprocess.check_output(['cp', 'merged_output/'+name+'.pdb', 'non_redundant/'])
        except:
            print ('{0} is missing, probably due to bad silent file'.format(name))

if __name__ == "__main__":
    sys.exit(main())
