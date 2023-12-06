#! /usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import seaborn
import math
import sys
import os

sashacolors = ['#000000', '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#800000', '#fffac8', '#42aaf5', '#ce42f5', '#eba534']
seaborn.set_palette(sashacolors)
seaborn.palplot(seaborn.color_palette())
print (len(sashacolors))
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams["figure.figsize"] = (15,10)

def divide_categories(categories):
    number_of_plots = math.ceil(len(categories) / 21 )
    cats_per_plot = math.ceil(len(categories) / number_of_plots )
    for i in range(0, len(categories), cats_per_plot): 
        yield categories[i:i + cats_per_plot]

def plot_metrics(df, x_values, y_values, color_coding=None, filename=False, prefix='', suffix='', close=1, title=False, size_factor=1):
    if color_coding:
        if len(set(df[color_coding])) == 2:
            red_df = df.loc[df[color_coding] == color_coding]
            grey_df = df.loc[df[color_coding] != color_coding]
            fig,ax = plt.subplots()
            ax.scatter(grey_df[x_values], grey_df[y_values], s=(6*size_factor), c='black', alpha=0.4)
            ax.scatter(red_df[x_values], red_df[y_values], s=(15*size_factor), c='red', alpha=1)

        elif len(set(df[color_coding])) > 24:
            categories = sorted(list(set(df[color_coding])))
            count = 1
            for category_set in divide_categories(categories):
                df_subset = df.loc[df[color_coding].isin(category_set)]
                df_subset = df_subset.reset_index(drop=True)
#                 print ('df_subset',df_subset)
                plot_metrics(df_subset, x_values, y_values, color_coding=color_coding, filename=filename, prefix=prefix, suffix=suffix+'__plot{0}'.format(count))
                count += 1
            return
#             cats_per_plot = math.ceil(len(categories) / number_of_plots )
#             print (number_of_plots, len(categories), cats_per_plot)
        else:
            if type(df[color_coding][0]) == str:
                #https://sashamaps.net/docs/resources/20-colors/
                fig,ax = plt.subplots()
                seaborn.scatterplot(data=df, hue=color_coding, x=x_values, y=y_values, s=(12*size_factor))
                handles, labels = ax.get_legend_handles_labels()
                # sort both labels and handles by labels
                labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
                ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #             plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            else:
                df.plot.scatter(x=x_values, y=y_values, c=color_coding, s=(12*size_factor), sharex=False);
    else:
        df.plot.scatter(x=x_values, y=y_values);
    plt.ylabel(y_values)
    plt.xlabel(x_values)
    plt.autoscale()
    plt.tight_layout()
    
    if title:
        plt.title(title)
    if not type(filename) == str:
        if color_coding: filename = '{0}{1}_vs_{2}__{3}_colored{4}.png'.format(prefix, x_values, y_values, color_coding, suffix)
        else: filename = '{0}{1}_vs_{2}{3}.png'.format(prefix, x_values, y_values, suffix)
    if filename:
        fig.savefig(filename, dpi=450)
    if close:
        plt.close('all')


def main():
# if 1:
    parser = argparse.ArgumentParser(description='run like:\n06_plot_metrics.py -tsv merged_output.tsv -threshold 5 -output output_fasta.fas')
    # parser.add_argument('-pdb', type=str, help='input pdbs', nargs='+', required=True)
    parser.add_argument('-tsv', type=str, help='table merged_output table from 03_extract_output.py OR 04_remove_redundant_sequences.py', default='non_redundant_output.tsv')
    parser.add_argument('-x_values', '-x', type=str, help='x values', nargs='+', default=['charge_pH_7'])
    parser.add_argument('-y_values', '-y', type=str, help='y values', nargs='+', default=['score_resnorm','SAP_resnorm'])
    # parser.add_argument('-threshold', type=int, help='threshold number of differences', default=5)
    # parser.add_argument('-output', '-out', '-o', type=str, help='output name stem (used for .fas and .tsv output)', default='non_redundant_output')
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep='\t')
    
    print ('Columns that could be plotted:')
    print(df.columns.tolist())

    x_values = ['charge_pH_7']#,'mean_seq_diff']
    y_values = ['score_resnorm','SAP_resnorm','charge_pH_7']

    for x in args.x_values:
        for y in args.y_values:
            if x == y: continue
            plot_metrics(df, x, y, color_coding='aa_comp', filename=1, close=False)
            plot_metrics(df, x, y, color_coding='scaffold', filename=1, close=False)



    all_aa_comps = list(set(df['aa_comp']))
    # all_aa_comps

    average_dict = {'aa_comp':[]}
    for x in args.x_values:
        average_dict[x] = []
    for y in args.y_values:
        average_dict[y] = []

    for aa_comp in all_aa_comps:
        df_subset = df[df['aa_comp'] == aa_comp]
        average_dict['aa_comp'].append(aa_comp)
        for x in args.x_values:
            average_dict[x].append(np.mean(df_subset[x]))
        for y in args.y_values:
            average_dict[y].append(np.mean(df_subset[y]))
            # df_subset[x]
    print (average_dict)
    average_df = pd.DataFrame(average_dict)
    for x in args.x_values:
        for y in args.y_values:
            plot_metrics(average_df, x, y, color_coding='aa_comp', filename=1, prefix='averaged_values__', close=False, size_factor=24)




    for aa_comp in all_aa_comps:
        ingroup_outgroup_list = []
        for value in df['aa_comp']:
            if value == aa_comp:
                ingroup_outgroup_list.append(aa_comp)
            else:
                ingroup_outgroup_list.append('other')
        df[aa_comp] = ingroup_outgroup_list

    plt.rcParams["figure.figsize"] = (10,8)

    x_values = ['charge_pH_7']#,'mean_seq_diff']
    y_values = ['score_resnorm','SAP_resnorm','charge_pH_7']

    for aa_comp in all_aa_comps:
        if not os.path.isdir(aa_comp):
            os.mkdir(aa_comp)
        for x in args.x_values:
            for y in args.y_values:
                if x == y: continue
                plot_metrics(df, x, y, color_coding=aa_comp, filename=1, prefix='{0}/{0}'.format(aa_comp), close=1, title=aa_comp)



if __name__ == "__main__":
    sys.exit(main())