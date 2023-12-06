#! /usr/bin/env python
import Bio
import Bio.Seq
import Bio.Align
from Bio.Align import AlignInfo

import pandas as pd
import subprocess
import argparse
import random
import math
import time
import glob
import sys
import os
import re


def gap_free_alignment(design, read):
    best_alignment = None
    best_score = 0
    
    for shift in range(0,len(design)):
        alignment = '-'*shift + read
        
        score = 0
        for d, r in zip(design, alignment):
            if r == d:
                score += 1
        if score > best_score:
            best_alignment = alignment  
            best_score = score        
    return best_alignment

def translate(align):
   trim = math.ceil(align.count('-')/3)*3
   CDS = align[trim:]
   CDS = CDS[:math.floor(len(CDS)/3)*3]
   protein = str(Bio.Seq.Seq(CDS).translate(to_stop=False))
   protein = '-'*int(trim/3)+protein
   return protein

def load_fasta(fas_filename):
    fas_list = []
    with open(fas_filename, 'r') as fas_file:
        fas_lines = fas_file.readlines()
    for i, line in enumerate( fas_lines ):
        if line.startswith('>'):
            header = re.sub(r'>(.+)\n$', r'\1', line)
            seq = re.sub(r'(.+)\n$', r'\1', fas_lines[i+1])
            fas_list.append( (header, seq) )    
    return fas_list

def main():
   arg_parser = argparse.ArgumentParser(description="  ")
   arg_parser.add_argument('-fasta', '-fas', '-f', type=str, nargs='+', help="fasta files of reads matched to a given category", required=True )   
   arg_parser.add_argument('-design_fasta', '-design', '-d', type=str, help="fasta file with (complete/assemblied) design sequence DNA", required=True )    
   arg_parser.add_argument('-nterm_anchor', '-nterm', '-anchor', '-n', type=str, help="N-terminal anchor", default='TCAGGAAGTTCTGGCGGC' )
   arg_parser.add_argument('-print_alignment', '-print', '-p', type=int, help="print DNA and protein alignments along with consensus sequence (which will always be print)?", default=1 )
   arg_parser.add_argument('-max_reads', '-max', '-m', type=int, help="max number of reads to look at for any design category", default=1000 )
   args = arg_parser.parse_args()

   for Suffix in ['_CONSENSUS.fas','_PROTEIN.fas','_ALIGN.fas']:
      args.fasta = [file for file in args.fasta if not file.endswith(Suffix)]

   for read_file in args.fasta:
      consensus_file = read_file.replace('.fas', '_CONSENSUS.fas')
      assert consensus_file != read_file
      protein_file = read_file.replace('.fas', '_PROTEIN.fas')
      assert protein_file != read_file
      alignment_file = read_file.replace('.fas', '_ALIGN.fas')
      assert alignment_file != read_file 

      design_name = read_file.split('.')[0].split('/')[-1]
      design_name = design_name.replace('_reads', '')
      
      should_translate = False
      if 'unassigned' in design_name:
         consensus = '-'  
      elif 'chimera' in design_name:   
         chimera_names = design_name.split('_chimera_')
         scaffold1 = re.sub(r'^(.*)_\w+_surf.*$', r'\1', chimera_names[0])
         scaffold2 = re.sub(r'^(.*)_\w+_surf.*$', r'\1', chimera_names[1])

         if scaffold1 == scaffold2:
            name_to_grep = chimera_names[0]
            should_translate = True
         else:
            consensus = '-'   

      else: #'if "design"':
         name_to_grep = design_name
         should_translate = True

      if should_translate:
         #print (f'TRANSLATING\t{design_name}'*8)
         #print (f'^>{name_to_grep}\n')
         #print (args.design_fasta)
         grepped_design_sequences = subprocess.check_output(['grep', '-A', '1', f'^>{name_to_grep}$', args.design_fasta]).decode('UTF-8')   
         # print (grepped_design_sequences,'grepped_design_sequences')
         design_DNA = grepped_design_sequences.split('\n')[1]

         if args.print_alignment:            
            seq_to_keep = grepped_design_sequences.split('\n')[0:10]
            # assert grepped_design_sequences.count('>') < 12
            # print (grepped_design_sequences,'grepped_design_sequences')

            with open(alignment_file, 'w') as output:
               print ('\n'.join(seq_to_keep), file=output)
            with open(protein_file, 'w') as output:
               pass
         
         Reads = [read for name, read in load_fasta(read_file)]
         Read1s = Reads[0::2]
         Read2s = Reads[1::2]
         Align1s = []
         Align2s = []
         Protein1s = []
         Protein2s = []

         for read1, read2 in zip(Read1s[:args.max_reads], Read2s[:args.max_reads]):
            if args.nterm_anchor not in read1: continue
            # trim read 1
            read1tr = args.nterm_anchor + read1.split(args.nterm_anchor)[-1]
            # rev comp read 2
            read2rc = str(Bio.Seq.Seq(read2).reverse_complement())
            alignment1 = gap_free_alignment(design_DNA, read1tr)
            alignment2 = gap_free_alignment(design_DNA, read2rc)
            Align1s.append( alignment1 )
            Align2s.append( alignment2 )
            Protein1s.append( translate(alignment1) )
            Protein2s.append( translate(alignment2) )

         print (f'\t{len(Protein1s)}/{len(Read1s)} translated')

         if not len(Protein1s):
            consensus = '-'
         else:
            max_len = max([len(p) for p in Protein1s] + [len(p) for p in Protein2s])
            Protein1s = [p.ljust(max_len,'-') for p in Protein1s]
            Protein2s = [p.ljust(max_len,'-') for p in Protein2s]

            if args.print_alignment:
               with open(alignment_file,'a') as output:
                  for align1, align2 in zip(Align1s, Align2s):
                     print (f'>{design_name}_r1\n{align1}', file=output)
                     print (f'>{design_name}_r2\n{align2}', file=output)
               with open(protein_file, 'a') as output:
                  for protein1, protein2 in zip(Protein1s, Protein2s):
                     print (f'>{design_name}_r1\n{protein1}', file=output)
                     print (f'>{design_name}_r2\n{protein2}', file=output)

            # Bio.Align.MultipleSeqAlignment import MultipleSeqAlignment
            protein_align = Bio.Align.MultipleSeqAlignment([])
            for protein1, protein2, n in zip(Protein1s, Protein2s, range(1,len(Protein1s)+1,1)):
               protein_align.add_sequence(f"{n}", protein1)
               protein_align.add_sequence(f"{n}", protein2)

            summary_align = AlignInfo.SummaryInfo(protein_align)
            consensus = str(summary_align.dumb_consensus(0.45))
            consensus = consensus.split('*')[0]
            # print (f'{design_name}_consenCDS')
            cMYC = 'EQKLISEEDL'
            # print (consensus)
            '''Remove flaking features from display construct'''
            consensus = re.sub(r'^[GS]+(.*)$',        r'\1',   consensus)
            # print (consensus)
            consensus = re.sub(r'^(.*)EQKL.{1,40}$',  r'\1',   consensus)
            # print (consensus)
            consensus = re.sub(r'^(.*)LEGGG.{1,20}$', r'\1',   consensus)
            # print (consensus)
            consensus = re.sub(r'^(.*)LEG$',          r'\1',   consensus)
            # print (consensus)            
            consensus = re.sub(r'^(.*?)[GS]+$',        r'\1',   consensus)

            if consensus.count('X') > 10:
               consensus = '-'
      else:
         print (f'skipping\t\t{design_name}')

      print (consensus)
      with open(consensus_file, 'w') as output:
         print (f'>{design_name}_consenCDS\n{consensus}', file=output) 



if __name__ == "__main__":
   sys.exit(main())