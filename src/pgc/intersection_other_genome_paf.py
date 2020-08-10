import argparse
from collections import defaultdict
from Bio import SeqIO
import os


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.24'


# need bedtools in path
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='intersection with other genome using minimap2')
parser.add_argument('-cpg', '--chinese_pan', type=str, help='chinese pan genome')
parser.add_argument('-in_paf', '--input_paf', type=str, help='using other genome as reference')
parser.add_argument('-pav_out', '--pav_record_out', type=str, help='the cpg pav record output')
parser.add_argument('-paf_part', '--paf_part', type=str, help='the part of paf file')


def main():
    args = parser.parse_args()
    chinese_pan = args.chinese_pan
    input_paf = args.input_paf
    pav_record_out = args.pav_record_out
    paf_part = args.paf_part
    exist = paf2parse(input_paf, paf_part)
    shuchu_pav(chinese_pan, pav_record_out, exist)


def paf2parse(input_paf, paf_part):
    f_paf_part = open(paf_part, 'w')
    length_dict = {}
    total_dict = defaultdict(list)
    with open(input_paf, 'r') as f:
        for raw_line in f:
            line = raw_line.strip('\n').split('\t')
            info = float(line[15].split(':')[-1])
            if info <= 0.1:
                length = int(line[3]) - int(line[2])
                if length >= 100:
                    f_paf_part.write(raw_line)
                    length_dict[line[0]] = int(line[1])
                    total_dict[line[0]].append((int(line[2]), int(line[3])))
    f_paf_part.close()

    temp = open('temp.bed', 'w')
    name_order = total_dict.keys()
    name_order = sorted(name_order)
    for i in name_order:
        new_per = sorted(total_dict[i], key=lambda x: (x[0], x[1]))
        for per in new_per:
            temp.write('\t'.join([i, str(per[0]), str(per[1])]) + '\n')
    temp.close()
    os.system('bedtools merge -i temp.bed > temp_new.bed')
    total_merge = defaultdict(list)
    with open('temp_new.bed', 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            total_merge[line[0]].append((int(line[1]), int(line[2])))
    os.system('rm temp.bed temp_new.bed')

    shuchu = set()
    for per in total_merge:
        length = sum(int(i[1]) - int(i[0]) + 1 for i in total_merge[per])
        cov = length / length_dict[per]
        if cov >= 0.8:
            shuchu.add(per)
            print(per, cov)
    print(len(shuchu))
    return shuchu


def shuchu_pav(chinese_pan, pav_record_out, exist):
    id_order = []
    pan_length = {}
    for seq in SeqIO.parse(chinese_pan, 'fasta'):
        id_order.append(seq.id)
        pan_length[seq.id] = len(seq.seq)

    exist_length = 0
    pav_order = []
    for i in id_order:
        if i in exist:
            pav_order.append(1)
            exist_length += pan_length[i]
        else:
            pav_order.append(0)
    f_pav = open(pav_record_out, 'w')
    for name, idx in zip(id_order, pav_order):
        f_pav.write('\t'.join([name, str(idx)]) + '\n')
    f_pav.close()
    print('exist number: %s' % len(exist))
    print('exist length: %s' % exist_length)


if __name__ == '__main__':
    main()

