#! /usr/bin/env python
from datetime import datetime
from Bio.Seq import Seq
import pandas as pd
import subprocess
import argparse
import random
import time
import glob
import sys
import os
import re

def blast(seq, full=0, short=0, return_seq=0, max_hits=10, subject='./chip_ordered_sequences/chipmunk_sequences.fas'):
   #print ('blasting seq:\n',seq)
   p1 = subprocess.Popen(['echo', '-e', '>query\n'+seq], stdout=subprocess.PIPE)
   if short:
      p2 = subprocess.Popen(['blastn', '-subject', subject, '-task', 'blastn-short', '-evalue', '0.00001', '-perc_identity', '100' ], stdin=p1.stdout, stdout=subprocess.PIPE)
   else:
      p2 = subprocess.Popen(['blastn', '-subject', subject], stdin=p1.stdout, stdout=subprocess.PIPE)

   blast_output = p2.communicate()[0].decode('UTF-8')
   blast_lines = blast_output.split('\n')
   # print ('blast_lines', '\n'.join(blast_lines))
   if full:
      print (blast_output)
   num_hits = 0
   for i, line in enumerate(blast_lines):
      if line.startswith('>'):
         name = line.split()[-1]
         identity_line = blast_lines[i+4]
         #print ('identity_line', identity_line)
         identity_ratio = re.sub(r'^ Identities = (\d+\/\d+).*\n?', r'\1', identity_line)
         #print ('identity_ratio', identity_ratio)
         matches = int(identity_ratio.split('/')[0])
         identity = float(identity_ratio.split('/')[0]) / float(identity_ratio.split('/')[1])
         length = int(identity_ratio.split('/')[1])
         if num_hits < max_hits:
            num_hits += 1
            if not return_seq:
               yield (matches, length, identity, name)
            else:
               query = blast_lines[i+7].split()[2]
               yield (matches, length, identity, name, query.upper() )

def parse_fastq(fastq):
   with open(fastq, 'r') as fastq_file:
      fastq_lines = fastq_file.readlines()
   sequences = fastq_lines[1:len(fastq_lines):4]
   sequences = [seq.strip() for seq in sequences]
   return(sequences)

def make_seq_count_dict(fasta):
   design_counts = {}
   with open(fasta, 'r') as fasta_file:
      fasta_lines = fasta_file.readlines()
   for line in fasta_lines:
      if line.startswith('>'):
         name = line.strip()[1:]
         design_counts[name] = []
   return design_counts

def count_fraction_of_reads(library_sequences, r1_fastq, r2_fastq, fraction, min_identity = 0.90, min_length = 150, debug=0, outdir='.', ):
   r1_seqs = parse_fastq(r1_fastq)
   r2_seqs = parse_fastq(r2_fastq)
   assert len(r1_seqs) == len(r2_seqs), f'There should be the same number of sequences in read1 and read2 files\nread1:{r1_fastq}\nread2:{r2_fastq}'
   # if debug: print (f'loaded {len(r1_seqs)} reads...')
   design_count_dict = make_seq_count_dict(library_sequences)
   design_count_dict['unassigned'] = []
   #print ('len(design_count_dict):', len(design_count_dict))
   blasted_counter = 0
   total_sequences_to_blast = len(r1_seqs)
   print ('total_sequences', total_sequences_to_blast)
   paired_reads_subset = []
   for seq1, seq2 in zip(r1_seqs, r2_seqs): 
      if random.random() < fraction:
         paired_reads_subset.append((seq1, seq2))

   # print ('len(paired_reads_subset)', len(paired_reads_subset))
   # print (paired_reads_subset[:10])
   for seq1, seq2 in paired_reads_subset: 
   # for seq1, seq2 in zip(r1_seqs, r2_seqs): 
      blasted_counter += 1
      hits1 = [hit for hit in blast(seq1, subject=library_sequences)]
      hits2 = [hit for hit in blast(seq2, subject=library_sequences)]

      # sort hits by sequence identity
      hits1.sort()
      hits1.reverse()
      hits2.sort()
      hits2.reverse()
      
      # if len(hits1) or len(hits2):            
      best_hit1 = 'unassigned'
      best_non_chimera1 = 'unassigned'
      for m1, l1, i1, name1 in hits1:
         if i1 > min_identity and l1 > min_length:
            best_hit1 = name1
            if not 'CHIMERA' in best_hit1:
               best_non_chimera1 = best_hit1
               break
         if i1 > min_identity and l1 > min_length and not 'CHIMERA' in name1:
            best_non_chimera1 = name1
         if best_hit1 and best_non_chimera1:
            break

      best_hit2 = 'unassigned'
      best_non_chimera2 = 'unassigned'
      for m2, l2, i2, name2 in hits2:
         if i2 > min_identity and l2 > min_length:
            best_hit2 = name2
            if not 'CHIMERA' in best_hit2:
               best_non_chimera2 = best_hit2
               break
         if i2 > min_identity and l2 > min_length and not 'CHIMERA' in name2:
            best_non_chimera2 = name2
         if best_hit2 and best_non_chimera2:
            break

      if best_hit1 == best_hit2:
         complete_name = best_hit1
      # If conflicted assignment, add to chimera 
      else:
         ## XX No chimeras of CHIMERAS allowed XX
         complete_name = best_non_chimera1+'_chimera_'+best_non_chimera2

      try:
         design_count_dict[complete_name].append(seq1)
         design_count_dict[complete_name].append(seq2)
      except KeyError: 
         design_count_dict[complete_name] = [seq1, seq2]
         
      #    # else:
      #    #    #debug = 3
      #    #    design_count_dict['unassigned'].append(seq1)
      #    #    design_count_dict['unassigned'].append(seq2)
      # else:
      #    #debug = 4
      #    design_count_dict['unassigned'].append(seq1)
      #    design_count_dict['unassigned'].append(seq2)

      # FOR DEV / DEBUGGING ONLY      
      # if debug:
      #    print ('debug:',debug)
      #    #FOR DEBUGGING
      #    for hit in hits1:
      #       print ('h1', hit)
      #    print ()
      #    for hit in hits2:
      #       print ('h2', hit)
      #    print ('\n'*3)
      
      if not blasted_counter % 1000:
         # print (f'{blasted_counter} out of {total_sequences_to_blast} sequences blasted', datetime.now())
         # print ('printing to fastas')
         for design in design_count_dict:
            if len(design_count_dict[design]):
               with open(f'{outdir}/{design}_reads.fas', 'w') as outfile:
                  for i, seq in enumerate(design_count_dict[design]):
                     print (f'>{design}_read{str(i).rjust(6,"0")}\n{seq}', file=outfile)
   #print ('len(design_count_dict):', len(design_count_dict))
   # print (f'{blasted_counter} out of {total_sequences_to_blast} sequences blasted', datetime.now())
   # print ('printing to fastas')
   for design in design_count_dict:
      if len(design_count_dict[design]):
         with open(f'{outdir}/{design}_reads.fas', 'w') as outfile:
            for i, seq in enumerate(design_count_dict[design]):
               print (f'>{design}_read{str(i).rjust(6,"0")}\n{seq}', file=outfile)

   designs = []
   counts = []
   for design in design_count_dict:
      designs.append(design)
      counts.append( len(design_count_dict[design])/2 ) # (Two reads added for each match)
   
   design_count_df = pd.DataFrame()
   design_count_df['design'] = designs
   design_count_df['count'] = counts
   design_count_df.sort_values('count', inplace=True, ascending=False)
   return design_count_df

def main():
   arg_parser = argparse.ArgumentParser(description="  ")
   arg_parser.add_argument('-fastq_stem', '-fastq', '-f', type=str, help="fastq stem for MiSeq data", required=True )    
   arg_parser.add_argument('-fraction', type=float, help="what fraction of sequence reads to blast?", default=0.002 )    
   arg_parser.add_argument('-tag', '-t', type=str, help="tag for experimental condition of fastq files (i.e. expressed/sorted)", required=True )    
   arg_parser.add_argument('-design_fasta', '-design', '-d', type=str, help="fasta file with (complete/assemblied) design sequence DNA", required=True )    
   arg_parser.add_argument('-subpool', '-s', type=str, help="what subpool (subdir) should read fastas go into", required=True )    
   args = arg_parser.parse_args()

   if not os.path.exists( 'matched_read_subset' ): os.mkdir('matched_read_subset')
   full_path = f'matched_read_subset/{args.subpool}_{args.tag}'
   if not os.path.exists(full_path): os.mkdir(full_path)

   reads1 = args.fastq_stem + '_R1_001.fastq'
   reads2 = args.fastq_stem + '_R2_001.fastq' 
   # print ('counting_reads...')
   count_df = count_fraction_of_reads(args.design_fasta, reads1, reads2, args.fraction, debug=0, outdir=full_path)
   print (count_df,'count_df')
   design_pool_name_stem = args.design_fasta.split('.')[0]
   design_pool_name_stem = design_pool_name_stem.split('/')[-1]
   tsv_name = f'matched_read_subset/{design_pool_name_stem}_{args.tag}_counts.tsv'

   count_df.to_csv(tsv_name, sep='\t', index=False)
   

if __name__ == "__main__":
   sys.exit(main())