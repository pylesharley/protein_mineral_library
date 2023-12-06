#! /usr/bin/env python

import glob
import re
import sys
import itertools
import numpy as np
import math

import pyrosetta
pyrosetta.init( extra_options = "-mute all -corrections::beta_nov16" )

import pandas as pd
import matplotlib.pyplot as plt
import seaborn
import argparse


def main():
    parser = argparse.ArgumentParser(description='')
    # parser.add_argument('-pdb', type=str, help='input pdbs', nargs='+', required=True)
    parser.add_argument('-tsv', type=str, help='non_redundant_output.tsv', default='non_redundant_output.tsv')
    parser.add_argument('-selections', type=str, help='surface_cooked_selections.tsv', default='./scaffolds/surface_cooked_selections.tsv')
    # parser.add_argument('-output', '-out', '-o', type=str, help='output name stem (used for .fas and .tsv output)', default='non_redundant_output')
    args = parser.parse_args()

    with open(args.selections, 'r') as cooked_selections:
        selection_lines = cooked_selections.readlines()

    dhr_surfaces = {}
    for i, line in enumerate(selection_lines):
        try: pdb, repeat, rep_n, select_a, select_b = line.split()
        except ValueError: print ('ERROR: could not parse line #{0}: {1}'.format(i+1,line)); continue
        repeat = int(repeat)
        select_a = [int(a) for a in select_a.split('+')]
        select_b = [int(b) for b in select_b.split('+')]

        dhr_surfaces[re.sub(r'(.*)\.pdb',r'\1',pdb)] = (select_a, select_b)


    dhr_surfaces['RiAFP_4DT5'] = ([41,62,82,129,109,43,64,84,131,111,45,66,86,133,113,47,68,88,135,115],[41,62,82,129,109,43,64,84,131,111,17,22,45,66,86,133,113,15,24,47,68,88,135,115,13,26])
    dhr_surfaces['FI998143'] = ([2,4,6,8,30,32,34,36,60,62,64,66],[14,15,18,22,25,43,44,47,51,54])
    dhr_surfaces['FI998252'] = ([2,4,6,8,32,34,36,38,64,66,68,70],[13,17,20,24,27,43,47,50,54,57])
    dhr_surfaces['FI2000161'] = ([3,5,7,9,11,36,38,40,42,44,70,72,74,76,78],[17,18,21,24,28,31,32,50,51,54,57,61,64,65])

    # results_path = '/home/pylesh/designs/repeats/library/july8th_results'
    df = pd.read_csv(args.tsv, sep='\t')

    amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    # COUNT FRACTIONS OF RESIDUES IN SURFACES
    # import random
    all_aa_proportions = []
    for AA in amino_acids:
        proportions = []
        for scaffold, name, seq in zip(df['scaffold'], df['description'], df['sequence']):
            if 'surfA' in name:
                subset = dhr_surfaces[scaffold][0]
            elif 'surfB' in name:
                subset = dhr_surfaces[scaffold][1]
            else:
                print ('WARNING! Could not parse surface from pdb name: {0}'.format(name))
                subset = None
            if subset: seq = [res for i, res in enumerate(seq) if i+1 in subset]
            proportions.append(float(seq.count(AA))/float(len(seq)))
    #         if not random.randint(0,50):
    #             print (AA, name, float(seq.count(AA))/float(len(seq)), '+'.join([str(x) for x in subset]))
    #             break
        all_aa_proportions.append(proportions)

    # ADD TO DATAFRAME
    for AA, proportion in zip(amino_acids, all_aa_proportions):
        try:
            df[AA]
        except KeyError:
            df[AA] = proportion
            print('added surface fractions of '+AA)

    aa_compositions = sorted(list(set(df['aa_comp'])))
    # aa_compositions

    table = []
    row_names = []


    for aa_comp in aa_compositions:
        average_aa_fractions = []

        aa_comp_df = df[df['aa_comp'] == aa_comp]
        for aa in amino_acids:
            average_aa_fractions.append(np.mean(aa_comp_df[aa]))

        row_names.append(aa_comp)
        table.append(average_aa_fractions)

    aa_fraction_df = pd.DataFrame(table, row_names, amino_acids)


    plt.rcParams["figure.figsize"] = (8,8)
    seaborn.color_palette("mako", as_cmap=True)

    ax = seaborn.heatmap(aa_fraction_df, annot=False, cmap="mako")
    # ax = seaborn.heatmap(aa_fraction_df, annot=aa_fraction_df.values.round(2), cmap="rocket")
    plt.tight_layout()


    figure = ax.get_figure()
    figure.savefig('AA_composition_report2.png', dpi=300)


if __name__ == "__main__":
    sys.exit(main())
