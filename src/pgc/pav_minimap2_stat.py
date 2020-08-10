from collections import defaultdict
from Bio import SeqIO
import networkx as nx
import argparse


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.24'


# python pav_minimap2_stat.py -paf {paf} -pan {CPG} -mincov 0.8 -mindiv 0.1 -pav {result}
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='presence and absence analyze')
parser.add_argument('-paf', '--paf_file', type=str, help='path to the paf')
parser.add_argument('-pan', '--pan_genome', type=str, help='path to the pan_genome')
parser.add_argument('-two_end', '--two_end_placed_fasta', type=str, help='path to the two end placed fasta')
parser.add_argument('-mincov', '--min_coverage', type=float, help='minimum coverage reserved')
parser.add_argument('-mindiv', '--min_divergence', type=float, help='minimum divergence reserved')
parser.add_argument('-pav', '--pav_result', type=str, help='path to the pav result')
parser.add_argument('-one_clstr', '--one_clstr', type=str, help='path to the one end clstr')
parser.add_argument('-two_clstr', '--two_clstr', type=str, help='path to the two end clstr')


def main():
    args = parser.parse_args()
    one_clstr = args.one_clstr
    two_clstr = args.two_clstr
    p_paf = args.paf_file
    sample = p_paf.split('/')[-1].split('_')[0]
    p_seq_pan = args.pan_genome
    two_end_placed_fasta = args.two_end_placed_fasta
    p_pav = args.pav_result
    min_coverage = args.min_coverage
    min_divergence = args.min_divergence
    pan_align = paf2region(p_paf, min_divergence)
    total_merge = merge2region(pan_align)
    exist1 = clstr2parse(one_clstr, sample)
    exist2 = clstr2parse(two_clstr, sample)
    pav_shuchu(total_merge, p_seq_pan, two_end_placed_fasta, min_coverage, sample, p_pav, exist1, exist2)


def paf2region(p_paf, min_divergence):
    pan_align = defaultdict(list)
    with open(p_paf, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            info = float(line[15].split(':')[-1])
            if info <= min_divergence:
                length = int(line[3]) - int(line[2])
                if length >= 100:
                    pan_align[(line[5], line[0])].append((int(line[7]), int(line[8])))
    return pan_align


def overlap_fun(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)


def merge2region(pan_align):
    total_merge = {}
    for per_name, per_query in pan_align.items():
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
    return total_merge


def clstr2parse(clstr, sample):
    cluster_num = 0
    represent = {}
    cluster_dict = defaultdict(list)
    with open(clstr, 'r') as f:
        for line in f:
            if 'Cluster' in line:
                cluster_num = line.strip('\n').split(' ')[1]
            if '*' in line:
                seq_name = line.split('...')[0].split('>')[-1]
                represent[cluster_num] = seq_name
            if 'at' in line:
                seq_sample = line.split('...')[0].split('>')[-1].split('_')[0]
                cluster_dict[cluster_num].append(seq_sample)
    exist = []
    for idx, per_list in cluster_dict.items():
        if sample in per_list:
            exist.append(represent[idx])
    return exist


def pav_shuchu(total_merge, p_seq_pan, two_end_placed_fasta, min_coverage, sample, p_pav, exist1, exist2):
    seq_dict = SeqIO.to_dict(SeqIO.parse(two_end_placed_fasta, 'fasta'))
    id_order = []
    pan_length = {}
    shuchu = set()
    for seq in SeqIO.parse(p_seq_pan, 'fasta'):
        id_order.append(seq.id)
        pan_length[seq.id] = len(seq.seq)
        if sample in seq.id:
            shuchu.add(seq.id)
    for per in total_merge:
        length = sum(int(i[1]) - int(i[0]) + 1 for i in total_merge[per])
        cov = length / pan_length[per[0]]
        if cov >= min_coverage:
            if per[0] not in seq_dict:
                shuchu.add(per[0])
                # print(per, total_merge[per], length, pan_length[per[0]])
    # print(shuchu)
    pav_order = []
    for i in id_order:
        if i in shuchu or i in exist1 or i in exist2:
            pav_order.append(1)
        else:
            pav_order.append(0)
    f_pav = open(p_pav, 'w')
    for name, idx in zip(id_order, pav_order):
        f_pav.write('\t'.join([name, str(idx)]) + '\n')
    f_pav.close()


if __name__ == '__main__':
    main()

