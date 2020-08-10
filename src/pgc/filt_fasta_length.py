#!/usr/bin/env python

from tinyfasta import FastaParser
import argparse

#usage: python ~/metagenome_data/filt_fasta_length.py --input /home/wzk/metagenome_data/assembly/HMP_MOCK_SRR2726667_8.25M/HMP_MOCK_SRR2726667_8.25M.contigs.fa --length 500 --output /home/wzk/metagenome_data/assembly/HMP_MOCK_SRR2726667_8.25M/HMP_MOCK_SRR2726667_8.25M.contigs_cut500.fa

#Author: Zhikun Wu
#Date: 2017.11.21
#Function: Filt the fasta file with threshold of sequence length

def cut_fasta(fasta_file, cut_len, fasta_out):
    cut_len = int(cut_len)
    out_h = open(fasta_out, 'w')
    for record in FastaParser(fasta_file):
        des = record.description
        seq = record.sequence
        length = len(seq)
        if length >= cut_len:
            out_h.write('%s\n%s\n' % (des, seq))
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description='Output the fasta file with length above the threshold.')
    parser.add_argument('-i', '--input', help='The input sequence with fasta format.')
    parser.add_argument('-l', '--length', help='The threshold for filted sequence length.')
    parser.add_argument('-o', '--output', help='The output file.')
    args = parser.parse_args()
    cut_fasta(args.input, args.length, args.output)

if __name__ == '__main__':
   main()
