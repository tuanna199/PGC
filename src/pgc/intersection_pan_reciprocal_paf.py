import argparse
from collections import defaultdict
import networkx as nx
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.24'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='intersection with pan genome using reciprocal strategy')
parser.add_argument('-cpg', '--chinese_pan', type=str, help='chinese pan genome')
parser.add_argument('-cpg_ref', '--cpg_reference', type=str, help='using cpg as reference')
parser.add_argument('-other_ref', '--other_reference', type=str, help='using other pan genome as reference')
parser.add_argument('-pav_out', '--pav_record_out', type=str, help='the cpg pav record output')
parser.add_argument('-paf_part_cpg', '--paf_part_cpg_ref', type=str, help='the part of paf file (raw paf too large)')
parser.add_argument('-paf_part_other', '--paf_part_other_ref', type=str, help='the part of paf file (other ref)')


def main():
    args = parser.parse_args()
    chinese_pan = args.chinese_pan
    cpg_reference = args.cpg_reference
    other_reference = args.other_reference
    pav_record_out = args.pav_record_out
    paf_part_cpg_ref = args.paf_part_cpg_ref
    paf_part_other_ref = args.paf_part_other_ref
    exist1 = other2cpg(cpg_reference, paf_part_cpg_ref)
    exist2 = cpg2other(other_reference, paf_part_other_ref)
    shuchu_pav(chinese_pan, pav_record_out, exist1, exist2)


def other2cpg(cpg_reference, paf_part_cpg_ref):
    f_paf_part_cpg_ref = open(paf_part_cpg_ref, 'w')
    length_dict = {}
    total_dict = defaultdict(list)
    with open(cpg_reference, 'r') as f:
        for raw_line in f:
            line = raw_line.strip('\n').split('\t')
            info = float(line[15].split(':')[-1])
            if info <= 0.1:
                length = int(line[3]) - int(line[2])
                if length >= 100:
                    f_paf_part_cpg_ref.write(raw_line)
                    length_dict[line[0]] = int(line[1])
                    total_dict[(line[0], line[5])].append((int(line[2]), int(line[3])))
    f_paf_part_cpg_ref.close()

    total_merge = {}
    for per_name, per_query in total_dict.items():
        new_per = []
        g = nx.Graph()
        for per1 in per_query:
            for per2 in per_query:
                if overlap_fun(per1[0], per1[1], per2[0], per2[1]):
                    g.add_edge(per1, per2)
        for cluster in nx.algorithms.connected_components(g):
            if len(cluster) == 1:
                new_per.append(list(cluster)[0])
            else:
                cluster_s = sorted(list(cluster), key=lambda x: x[0])
                cluster_e = sorted(list(cluster), key=lambda x: x[1])
                new_per.append((cluster_s[0][0], cluster_e[-1][1]))
        new_per = sorted(new_per, key=lambda x: (x[0], x[1]))
        total_merge[per_name] = new_per

    shuchu = set()
    for per in total_merge:
        length = sum(int(i[1]) - int(i[0]) + 1 for i in total_merge[per])
        cov = length / length_dict[per[0]]
        if cov >= 0.8:
            shuchu.add(per[1])
            print(per, cov)
    print('other2cpg', len(shuchu))
    return shuchu


def cpg2other(other_reference, paf_part_other_ref):
    f_paf_part_other_ref = open(paf_part_other_ref, 'w')
    length_dict = {}
    total_dict = defaultdict(list)
    with open(other_reference, 'r') as f:
        for raw_line in f:
            line = raw_line.strip('\n').split('\t')
            info = float(line[15].split(':')[-1])
            if info <= 0.1:
                length = int(line[3]) - int(line[2])
                if length >= 100:
                    f_paf_part_other_ref.write(raw_line)
                    length_dict[line[0]] = int(line[1])
                    total_dict[(line[0], line[5])].append((int(line[2]), int(line[3])))
    f_paf_part_other_ref.close()

    total_merge = {}
    for per_name, per_query in total_dict.items():
        new_per = []
        g = nx.Graph()
        for per1 in per_query:
            for per2 in per_query:
                if overlap_fun(per1[0], per1[1], per2[0], per2[1]):
                    g.add_edge(per1, per2)
        for cluster in nx.algorithms.connected_components(g):
            if len(cluster) == 1:
                new_per.append(list(cluster)[0])
            else:
                cluster_s = sorted(list(cluster), key=lambda x: x[0])
                cluster_e = sorted(list(cluster), key=lambda x: x[1])
                new_per.append((cluster_s[0][0], cluster_e[-1][1]))
        new_per = sorted(new_per, key=lambda x: (x[0], x[1]))
        total_merge[per_name] = new_per

    shuchu = set()
    for per in total_merge:
        length = sum(int(i[1]) - int(i[0]) + 1 for i in total_merge[per])
        cov = length / length_dict[per[0]]
        if cov >= 0.8:
            shuchu.add(per[0])
            print(per, cov)
    print('cpg2other', len(shuchu))
    return shuchu


def shuchu_pav(chinese_pan, pav_record_out, exist1, exist2):
    id_order = []
    pan_length = {}
    for seq in SeqIO.parse(chinese_pan, 'fasta'):
        id_order.append(seq.id)
        pan_length[seq.id] = len(seq.seq)

    exist_length = 0
    pav_order = []
    for i in id_order:
        if i in exist1 or i in exist2:
            pav_order.append(1)
            exist_length += pan_length[i]
        else:
            pav_order.append(0)
    f_pav = open(pav_record_out, 'w')
    for name, idx in zip(id_order, pav_order):
        f_pav.write('\t'.join([name, str(idx)]) + '\n')
    f_pav.close()
    print('exist number: %s' % len(exist1 | exist2))
    print('exist length: %s' % exist_length)


def overlap_fun(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)


if __name__ == '__main__':
    main()

