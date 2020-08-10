import argparse
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.22'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='anchor NRS from raw bed file')
parser.add_argument('-bed', '--bed_raw', type=str, help='path to the raw bed of coordinate')
parser.add_argument('-out', '--bed_new', type=str, help='path to the new bed of coordinate')
parser.add_argument('-fasta', '--NRS_fasta_nc', type=str, help='non-contaminant NRS')


def main():
    args = parser.parse_args()
    bed_raw = args.bed_raw
    bed_new = args.bed_new
    nrs_fasta_nc = args.NRS_fasta_nc

    raw_record = {}
    with open(bed_raw, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            raw_record[line[3]] = line

    f_shuchu = open(bed_new, 'w')
    for seq in SeqIO.parse(nrs_fasta_nc, 'fasta'):
        name = '_'.join(seq.id.split('_')[0:3])
        if name in raw_record:
            newline = raw_record[name]
            newline[3] = seq.id
            f_shuchu.write('\t'.join(newline) + '\n')
    f_shuchu.close()


if __name__ == '__main__':
    main()

