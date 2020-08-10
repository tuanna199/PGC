import re
import argparse
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.01.19'


# python lastz_32_filter.py -m {maf} -i {unmap_contigs} -o {filter_result} -l 1000 -c 0.1
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='filter unmap contigs region')
parser.add_argument('-m', '--maf_file', type=str, help='path to the maf')
parser.add_argument('-i', '--unmap_path', type=str, help='path to the unmap contigs')
parser.add_argument('-o', '--filter_result', type=str, help='path to the filter result')
parser.add_argument('-l', '--min_length', type=int, help='minimum length(bp) reserved')
parser.add_argument('-c', '--min_coverage', type=float, help='minimum coverage reserved')


def main():
    args = parser.parse_args()
    maf_file = args.maf_file
    unmap_path = args.unmap_path
    filter_result = args.filter_result
    min_length = args.min_length
    min_coverage = args.min_coverage
    # maf_file = 'hg38_vs_m416-cns-unmap.maf'
    # unmap_path = 'm416_cns_unmap.fa'
    # filter_result = 'm416_cns_unmap_lastz_0.1.fa'
    # min_length = 1000
    # min_coverage = 0.1
    read_unmap2filter(maf_file, unmap_path, filter_result, min_length, min_coverage)


def read_unmap2filter(maf_file, unmap_path, filter_result, min_length, min_coverage):
    unmap_contigs = SeqIO.to_dict(SeqIO.parse(unmap_path, 'fasta'))
    idy_list = []
    pattern_str = re.compile(r'\S+')
    f_lastz_32 = open(maf_file, 'r')
    for line in f_lastz_32:
        if not line.startswith('#'):
            if 'ctg' in line:
                line = line.strip('\n')
                temp = re.findall(pattern_str, line)
                idy = int(temp[3]) / int(temp[5])
                if idy > min_coverage:
                    idy_list.append(temp[1])
    f_lastz_32.close()
    # print(len(set(idy_list)))
    partly_align = 0
    partly_align_base = 0
    non_align = 0
    non_align_base = 0
    shuchu = []
    for seq_name, seq_unmap in unmap_contigs.items():
        if seq_name not in idy_list:
            if len(seq_unmap) >= min_length:
                shuchu.append(seq_unmap)
                if 'none_none' in seq_name:
                    non_align += 1
                    non_align_base += len(seq_unmap)
                else:
                    partly_align += 1
                    partly_align_base += len(seq_unmap)
    print('Partly alignment: %s, %s bp' % (partly_align, partly_align_base))
    print('No alignment at all: %s, %s bp' % (non_align, non_align_base))
    SeqIO.write(shuchu, filter_result, 'fasta')


if __name__ == '__main__':
    main()
