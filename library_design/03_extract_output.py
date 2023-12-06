#! /usr/bin/env python

import matplotlib.pyplot as plt
from PIL import Image 
import matplotlib
import numpy as np
import subprocess
import argparse
import glob
import time
import sys
import os
import re                         

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import copy

def compare_one_to_list(seq1, seq_list):
    differences = [ count_differences(seq1, seq2) for seq2 in seq_list ]
    differences.remove(0) 
    return np.mean(differences)

def count_differences(seq1, seq2):
    #print ([seq1, seq2])
    if seq1 is None:
        return float('NaN')
    if seq2 is None:
        return float('NaN')
    try:
        assert len(seq1) == len(seq2)
    except AssertionError:
        print (seq1, seq2)
        assert len(seq1) == len(seq2)
    number = 0
    for r1, r2 in zip(seq1, seq2):
        if r1 != r2:
            number += 1
    return number 

def add_to_dict(dictionary, key, value):
    if type(value) == list:
        try: dictionary[key].extend(value)
        except KeyError: dictionary[key] = value
    elif type(value) == str:
        try: dictionary[key].append(value)
        except KeyError: dictionary[key] = [value]
    else:
        assert type(value) == str or type(value) == list

# def infer_silent_name(pdb):
#     pdb = pdb.replace('_RELAX','')
#     silent = re.sub(r'(.*)_\d\d\d\d', r'\1.silent', pdb)
#     return silent

def main():
# if 1:
    parser = argparse.ArgumentParser(description='')
    # parser.add_argument('-pdb', type=str, help='input pdbs', nargs='+', required=True)
    parser.add_argument('-glob','-g', type=str, help='path to silent and score files', default='*/*/*')
    parser.add_argument('-sortby', '-sort', type=str, help='column to sort by', default='total_score' )
    parser.add_argument('-score', type=str, help='column to report uncst scores', default='BetaNov16' )
    # parser.add_argument('-title', type=str, help='title table file(tsv)', default=1 )
    # parser.add_argument('-off', type=str, help="query string to identify 'off' surface state, aka APO models", default=False )
    parser.add_argument('-extract', type=int, help='boolean for whether to extract pdbs or not (turn off to replot)', default=1)
    parser.add_argument('-silent_extracter', type=str, help='boolean for whether to extract pdbs or not (turn off to replot)', default='/software/rosetta/latest/main/source/bin/extract_pdbs.hdf5.linuxgccrelease')
    parser.add_argument('-database', type=str, help='boolean for whether to extract pdbs or not (turn off to replot)', default='/software/rosetta/latest/main/database')
    # parser.add_argument('-range', type=int, help='score_min score_max', nargs='+', default=False ) 
    # parser.add_argument('-combine', type=int, help='concatenate output plots', default=1 ) 
    parser.add_argument('-plot', type=int, help='plot or not?', default=1 ) 
    parser.add_argument('-autorange', type=int, help='automatically calculate energy range for plotting?', default=1 ) 
    parser.add_argument('-mute', type=int, help='print internal variables for debugging?', default=1 ) 
    # parser.add_argument('-combine', type=int, help='concatenate output plots', default=1 ) 
    # parser.add_argument('-col', type=int, help='columns in combined plots', default=4 ) 
    args = parser.parse_args()

    print ('{0}_score.sc'.format(args.glob))
    score_files = glob.glob( '{0}_score.sc'.format(args.glob) )
    if not len(score_files):
        print ('{0}.sc'.format(args.glob))
        score_files = glob.glob( '{0}.sc'.format(args.glob) )

    # silent_files = glob.glob('{0}.silent'.format(args.glob) )
    os.system('cat {0}.silent > merged_output.silent'.format(args.glob))
    # subprocess.command(['cat', '{0}.silent'.format(args.glob), '>', 'merged_output.silent'])
    
    design_keyed_scorelines = {}
    designs = []

    for score_file in score_files:
        # print ('score_file', score_file)
        with open(score_file, 'r') as open_score_file:
            score_lines = open_score_file.readlines()
        
        # ADD SILENT FILE TO SCORE LINE
        header = score_lines[1].split()
        silent = re.sub(r'(.*)_score.sc', r'\1.silent', score_file)
        score_lines = [line.split()+[silent] for line in score_lines[2:]]

        design_name = re.sub(r'(.*)_thread\d*_score.sc', r'\1', score_file)
        designs.append(design_name)

        add_to_dict(design_keyed_scorelines, design_name, score_lines)
        # try: design_keyed_scorelines[design_name].extend(score_lines)
        # except KeyError: design_keyed_scorelines[design_name] = [header] + score_lines

    unique_designs = list(set(designs))
    ALL_PDBS_TO_EXTRACT = []
    ALL_EXTRACTED_SCORE_LINES = []
    ALL_EXTRACTED_SILENTS = []
    ALL_EXTRACTED_SEQUENCES = []
    # print ('unique_designs', unique_designs)

    for design_name in unique_designs:
        if not args.mute: print ('design_name:',design_name)
        # if design_name == 'FI_AP_6MRR_DEAroMidGr_surfA_antirep_thread1_0001':
        #     continue
        # header = design_keyed_scorelines[design_name][0]
        score_lines = design_keyed_scorelines[design_name][1:]
        # SORT SCORE LINES HERE
        sort_index = header.index(args.sortby)
        score_index = header.index(args.score)
        score_lines.sort(key=lambda line: float(line[sort_index]))

        # print(score_lines[:10])
        all_structures = [ line[-2] for line in score_lines ]
        #print (all_structures[:10], all_structures[:10])
        all_scores = [ float(line[score_index]) for line in score_lines ]
        all_silents = [ line[-1] for line in score_lines ]

        quarter_percental_index = max(1, int(len(all_scores)*0.25))

        sequences = []
        silents = []
        for structure, silent_file in zip(all_structures[:quarter_percental_index], all_silents[:quarter_percental_index]):
            # print (' '.join(['grep', 'ANNOTATED_SEQUENCE.*{0}'.format(structure), silent_file]))
            silents.append(silent_file)
            try:
                sequences.append(subprocess.check_output(['grep', 'ANNOTATED_SEQUENCE.*{0}'.format(structure), silent_file]))
            except:
                sequences.append(None)
        sequences = [str(sequence).replace('[HIS_D]','') for sequence in sequences if not sequence is None]
        sequences = [re.sub(r"b'ANNOTATED_SEQUENCE: ([A-Z])\[.*\]([A-Z]+)\[.*", r'\1\2', sequence) for sequence in sequences]        
        sequences = [re.sub(r"b'ANNOTATED_SEQUENCE: ([A-Z])\[.*\]([A-Z]+)\[.*", r'\1\2', sequence) for sequence in sequences]        
        sequences = [re.sub(r"b'ANNOTATED_SEQUENCE: ([A-Z])\[.*\]([A-Z]+).*", r'\1\2', sequence) for sequence in sequences]        

        ALL_PDBS_TO_EXTRACT.append(all_structures[0]) # sorted already
        ALL_EXTRACTED_SCORE_LINES.append(score_lines[0])
        ALL_EXTRACTED_SEQUENCES.append(sequences[0])
        ALL_EXTRACTED_SILENTS.append(silents[0])

        diff_to_low_nrg = []
        for sequence in sequences:
            diff_to_low_nrg.append(count_differences(sequences[0], sequence))
        if max(diff_to_low_nrg) > 3: # If any structures in quarter lowest nrg have different sequences than lowest energy extract that too
            most_diff = diff_to_low_nrg.index(max(diff_to_low_nrg))
            ALL_PDBS_TO_EXTRACT.append(all_structures[most_diff])
            ALL_EXTRACTED_SCORE_LINES.append(score_lines[most_diff])# + [scaffold, aa_comp, seq])
            ALL_EXTRACTED_SEQUENCES.append(sequences[most_diff])
            ALL_EXTRACTED_SILENTS.append(silents[most_diff])


    print ('ALL_PDBS_TO_EXTRACT',ALL_PDBS_TO_EXTRACT[:10])
    print ('ALL_EXTRACTED_SILENTS',ALL_EXTRACTED_SILENTS[:10])

    if not os.path.isdir('merged_output'): os.mkdir('merged_output')
    failures = []
    if args.extract:
        print ('Extracting', len(ALL_PDBS_TO_EXTRACT), 'pdbs')    
        for silent, tag in zip(ALL_EXTRACTED_SILENTS, ALL_PDBS_TO_EXTRACT):
            extract_command = f'{args.silent_extracter} -database {args.database} -in:file:silent {0} -out:prefix merged_output/ -in:file:tags {1}'.format(silent, tag)
            if not args.mute:
                print(extract_command)
            try:
                subprocess.check_output(extract_command.split())
            except:
                failures.append((silent, tag))
                print ('fail:', silent, tag)

    print ('Failed to extract {0} designs'.format(len(failures)))

    average_difference = {}
    ALL_EXTRACTED_SCAFFOLDS = [re.sub(r'([^/]+)/.*(surf[ABCDE]).*', r'\1_\2', LINE[-1]) for LINE in ALL_EXTRACTED_SCORE_LINES]
    for SCAFFOLD, SEQUENCE in zip(ALL_EXTRACTED_SCAFFOLDS, ALL_EXTRACTED_SEQUENCES):
        other_scaffold_sequences = [seq for scaff, seq in zip(ALL_EXTRACTED_SCAFFOLDS, ALL_EXTRACTED_SEQUENCES) if scaff == SCAFFOLD]
        average_difference[SEQUENCE] = compare_one_to_list(SEQUENCE, other_scaffold_sequences)

    HEADER = header[:-1] + ['charge_pH_7', 'isoelectric_point', 'GRAVY_value', 'mean_seq_diff', 'scaffold', 'aa_comp', 'description', 'sequence']
    df = pd.DataFrame(columns=HEADER)
        
    for LINE, SEQUENCE in zip(copy.deepcopy(ALL_EXTRACTED_SCORE_LINES), ALL_EXTRACTED_SEQUENCES):
        silent = LINE.pop()
        description = LINE.pop()
        scaffold = re.sub(r'([^/]+)/.*(_surf[ABCDE]).*', r'\1', silent)
        aa_comp = re.sub(r'[^/]+/([^/]+)_surf.*', r'\1', silent)
        
        bio_prot = ProteinAnalysis(SEQUENCE)
        ph_value = 7
        charge_at_ph = bio_prot.charge_at_pH(ph_value)
        # Protein GRAVY returns the GRAVY (grand average of hydropathy) value for the protein sequences you enter. 
        # Algorithm and values from Bjellqvist and Tabb [https://biopython.org/docs/1.75/api/Bio.SeqUtils.IsoelectricPoint.html]
        isoelectric_point = bio_prot.isoelectric_point()
        # The GRAVY value is calculated by adding the hydropathy value for each residue and dividing by the length of the sequence (Kyte and Doolittle,1982)
        gravy_value = bio_prot.gravy()
        
        df.loc[len(df)] = LINE + [charge_at_ph, isoelectric_point, gravy_value, average_difference[SEQUENCE], scaffold, aa_comp, description, SEQUENCE]
    
    df = df.round(decimals=4)
    df.to_csv("merged_output.tsv", sep="\t", index=False)

if __name__ == "__main__":
    sys.exit(main())
