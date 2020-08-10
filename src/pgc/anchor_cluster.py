import argparse
from collections import defaultdict
import networkx as nx
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.23'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='cluster by anchor')
parser.add_argument('-b', '--bed_merge', type=str, default='', help='merged of bed')
parser.add_argument('-f', '--fa_merge', type=str, default='', help='merge of fa')
parser.add_argument('-two_end_clstr', '--two_end_clstr', type=str, help='clstr of two end anchor')
parser.add_argument('-one_end_clstr', '--one_end_clstr', type=str, help='clstr of one end anchor')
parser.add_argument('-two_fa', '--two_end_merge_fa', type=str, help='fasta of two end anchor cluster')
parser.add_argument('-one_fa', '--one_end_merge_fa', type=str, help='fasta of one end anchor cluster')


def main():
    args = parser.parse_args()
    bed_merge = args.bed_merge
    fa_merge = args.fa_merge
    two_end_clstr = args.two_end_clstr
    one_end_clstr = args.one_end_clstr
    two_end_merge_fa = args.two_end_merge_fa
    one_end_merge_fa = args.one_end_merge_fa
    raw_record, two_e_a, one_e_a, raw_fa, raw_fa_len = read2bed2fa(bed_merge, fa_merge)
    clique_two = bed2clique(two_e_a)
    clique_one = bed2clique(one_e_a)
    clique2pav(clique_two, raw_fa, raw_fa_len, two_end_clstr, two_end_merge_fa)
    clique2pav(clique_one, raw_fa, raw_fa_len, one_end_clstr, one_end_merge_fa)


def read2bed2fa(bed_merge, fa_merge):
    raw_fa = {}
    raw_fa_len = {}
    for seq in SeqIO.parse(fa_merge, 'fasta'):
        raw_fa[seq.id] = seq
        raw_fa_len[seq.id] = len(seq.seq)
    raw_record = {}
    two_e_a = defaultdict(list)
    one_e_a = defaultdict(list)
    with open(bed_merge, 'r') as f:
        for raw_line in f:
            line = raw_line.strip('\n').split('\t')
            raw_record[line[3]] = raw_line
            if line[1] != '.' and line[2] != '.':
                if int(line[1]) < int(line[2]):
                    mid = int(line[1]) + (int(line[2]) - int(line[1])) / 2
                    two_e_a[line[0]].append((line[3], mid - 250, mid + 250, line[5]))
                elif int(line[1]) > int(line[2]):
                    mid = int(line[2]) + (int(line[1]) - int(line[2])) / 2
                    two_e_a[line[0]].append((line[3], mid - 250, mid + 250, line[5]))
                else:
                    mid = int(line[1])
                    two_e_a[line[0]].append((line[3], mid - 250, mid + 250, line[5]))
            if line[1] == '.':
                one_e_a[line[0]].append((line[3], int(line[2]) - 250, int(line[2]) + 250, line[5]))
            if line[2] == '.':
                one_e_a[line[0]].append((line[3], int(line[1]) - 250, int(line[1]) + 250, line[5]))
    # print(len(raw_record))
    return raw_record, two_e_a, one_e_a, raw_fa, raw_fa_len


def bed2clique(two_e_a):
    clique_two = {}
    for chrom in two_e_a:
        per_chrom = sorted(two_e_a[chrom], key=lambda x: (x[1], x[2]))
        # print(per_chrom)
        g = nx.Graph()
        for per1 in per_chrom:
            # g.add_node(per1)
            for per2 in per_chrom:
                if overlap_fun(per1[1], per1[2], per2[1], per2[2]):
                    g.add_edge(per1, per2)
        max_clique_list = []
        temp_nodes = set(g.nodes)
        # print(temp_nodes)
        while temp_nodes:
            temp_graph = nx.Graph.subgraph(g, list(temp_nodes))
            clique_list = [(i, len(i)) for i in nx.algorithms.clique.find_cliques(temp_graph)]
            clique_list = sorted(clique_list, key=lambda x: x[1], reverse=True)
            clique_max = clique_list[0][0]
            max_clique_list.append(tuple(clique_max))
            temp_nodes -= set(clique_max)
            # print(clique_max)
        clique_two[chrom] = max_clique_list
    return clique_two


def clique2pav(clique_two, raw_fa, raw_fa_len, two_end_clstr, two_end_merge_fa):
    f_shuchu = open(two_end_clstr, 'w')
    n = 0
    shuchu_seq = []
    for chrom in clique_two:
        for clique in clique_two[chrom]:
            clique = [(i[0], raw_fa_len[i[0]]) for i in clique]
            clique = sorted(clique, key=lambda x: x[1], reverse=True)
            shuchu_seq.append(raw_fa[clique[0][0]])
            f_shuchu.write('>Cluster %s\n' % n)
            n += 1
            m = 0
            f_shuchu.write('%s\t%snt, >%s... *\n' % (m, raw_fa_len[clique[0][0]], clique[0][0]))
            for idx in range(1, len(clique)):
                m += 1
                f_shuchu.write('%s\t%snt, >%s... at ./%.2f%s\n' %
                               (m, raw_fa_len[clique[idx][0]], clique[idx][0], 0.00, '%'))
    print(len(shuchu_seq))
    SeqIO.write(shuchu_seq, two_end_merge_fa, 'fasta')
    f_shuchu.close()


def overlap_fun(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)


if __name__ == '__main__':
    main()

