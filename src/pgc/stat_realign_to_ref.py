import argparse
from collections import defaultdict
from Bio import SeqIO
import os


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.21'


# need bedtools in path
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='statistic of realign to ref genome using minimap2')
parser.add_argument('-NRS_one', '--NRS_one', type=str, help='path to NRS one')
parser.add_argument('-in_paf', '--input_paf_part', type=str, help='non reference realign to the ref')
parser.add_argument('-out_fasta', '--out_fasta', type=str, help='the fasta output')


def main():
    args = parser.parse_args()
    nrs_one = args.NRS_one
    input_paf_part = args.input_paf_part
    out_fasta = args.out_fasta
    exist = paf2parse(input_paf_part)
    shuchu_fa(nrs_one, out_fasta, exist)


def paf2parse(input_paf_part):
    length_dict = {}
    total_dict = defaultdict(list)
    with open(input_paf_part, 'r') as f:
        for raw_line in f:
            line = raw_line.strip('\n').split('\t')
            info = float(line[15].split(':')[-1])
            if info <= 0.1:
                length = int(line[3]) - int(line[2]) + 1
                if length >= 100:
                    length_dict[line[0]] = int(line[1])
                    total_dict[line[0]].append((int(line[2]), int(line[3])))

    sample_name = input_paf_part.split('/')[-1].split('_')[0]
    temp = open('%s_temp.bed' % sample_name, 'w')
    name_order = total_dict.keys()
    name_order = sorted(name_order)
    for i in name_order:
        new_per = sorted(total_dict[i], key=lambda x: (x[0], x[1]))
        for per in new_per:
            temp.write('\t'.join([i, str(per[0]), str(per[1])]) + '\n')
    temp.close()
    os.system('bedtools merge -i %s_temp.bed > %s_temp_new.bed' % (sample_name, sample_name))
    total_merge = defaultdict(list)
    with open('%s_temp_new.bed' % sample_name, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            total_merge[line[0]].append((int(line[1]), int(line[2])))
    os.system('rm %s_temp.bed %s_temp_new.bed' % (sample_name, sample_name))
    
    shuchu = set()
    for per in total_merge:
        length = sum(int(i[1]) - int(i[0]) + 1 for i in total_merge[per])
        cov = length / length_dict[per]
        if cov >= 0.8:
            shuchu.add(per)
            # print(per, cov)
    print(len(shuchu))
    return shuchu


def shuchu_fa(nrs_one, out_fasta, exist):
    shuchu = []
    exist_length = 0
    for i in SeqIO.parse(nrs_one, 'fasta'):
        if len(i.seq) >= 1000:
            if i.id not in exist:
                shuchu.append(i)
            else:
                exist_length += len(i.seq)
    print('exist number: %s' % len(exist))
    print('exist length: %s' % exist_length)
    SeqIO.write(shuchu, out_fasta, 'fasta')


if __name__ == '__main__':
    main()

