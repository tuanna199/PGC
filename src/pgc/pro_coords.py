import re
import argparse
import networkx as nx
from collections import defaultdict, Counter
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.01.15'


# python pro_coords.py -c {coords} -i {contigs} -o {unmap} -log {prefix} -l 1000
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='extract unmap contigs region')
parser.add_argument('-c', '--coords_file', type=str, help='path to the coords')
parser.add_argument('-i', '--contigs_path', type=str, help='path to the contigs')
parser.add_argument('-o', '--unmap_out_path', type=str, help='path to the unmap region')
parser.add_argument('-log', '--log_name', type=str, help='prefix of log file')
parser.add_argument('-l', '--min_length', type=int, help='minimum length reserved')


def main():
    args = parser.parse_args()
    coords_file = args.coords_file
    min_length = args.min_length
    query_dict, query_length_dict, query_most_freq_chr = coords2parse(coords_file)
    total_merge = sort2merge(query_dict)
    total_unmap = merge2unmap(total_merge, query_length_dict, min_length)

    log_name = args.log_name
    shuchu_log(log_name, query_dict, query_length_dict, query_most_freq_chr, total_merge, total_unmap)

    contigs_path = args.contigs_path
    unmap_out_path = args.unmap_out_path
    read_all_contig(unmap_out_path, contigs_path, total_unmap, query_most_freq_chr, query_length_dict, min_length)


def coords2parse(coords_file):
    query_most_freq_chr = defaultdict(list)
    query_dict = defaultdict(list)
    query_length_dict = {}
    idy_list = []
    f_file = open(coords_file, 'r')
    for line in f_file:
        if '\t' in line:
            line = line.strip('\n').split('|')
            pattern_num = re.compile(r'\d+')
            # Start(End) of the alignment region in the query sequence
            query_se = re.findall(pattern_num, line[1])
            if int(query_se[0]) < int(query_se[1]):
                query_s, query_e = int(query_se[0]), int(query_se[1])
            else:
                query_s, query_e = int(query_se[1]), int(query_se[0])
            # Start(End) of the alignment region in the reference sequence
            ref_se = re.findall(pattern_num, line[0])
            ref_s, ref_e = int(ref_se[0]), int(ref_se[1])
            # Percent identity of the alignment
            idy = line[3].replace(' ', '')
            idy_list.append(float(idy))
            # The query FastA ID
            query_name = line[6].split('\t')[1]
            # The reference FastA ID
            ref_name = line[6].split('\t')[0].replace(' ', '')
            # Length of the query sequence
            length_q = int(re.findall(pattern_num, line[4])[1])

            query_dict[query_name].append((query_s, query_e, ref_name, ref_name, ref_s, ref_e))
            query_most_freq_chr[query_name].append(ref_name)
            query_length_dict[query_name] = length_q

    for _, per_query in query_most_freq_chr.items():
        per_query = Counter(per_query)
        per_query = per_query.most_common(1)[0][0]
        query_most_freq_chr[_] = per_query

    for _, per_query in query_dict.items():
        per_query = sorted(per_query, key=lambda x: (x[0], x[1]))
        query_dict[_] = per_query
    return query_dict, query_length_dict, query_most_freq_chr


def overlap_fun(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)


def sort2merge(query_dict):
    total_merge = {}
    for per_name, per_query in query_dict.items():
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
                cluster_s = sorted(list(cluster), key=lambda x: (x[0], x[2]))
                cluster_e = sorted(list(cluster), key=lambda x: (x[1], x[3]))
                new_per.append((cluster_s[0][0], cluster_e[-1][1], cluster_s[0][2], cluster_e[-1][3],
                                cluster_s[0][4], cluster_e[-1][5]))
        new_per = sorted(new_per, key=lambda x: (x[0], x[1]))
        total_merge[per_name] = new_per
    return total_merge


def merge2unmap(total_merge, query_length_dict, min_length):
    total_unmap = {}
    for per_name, per_query in total_merge.items():
        unmap_region = []

        if per_query[0][0] >= min_length:
            unmap_region.append((1, per_query[0][0], 'none', per_query[0][2], 'none', per_query[0][4]))

        for idx in range(len(per_query) - 1):
            region1 = per_query[idx]
            region2 = per_query[idx + 1]
            if region2[0] - region1[1] >= min_length:
                unmap_region.append((region1[1], region2[0], region1[3], region2[2], region1[5], region2[4]))

        if query_length_dict[per_name] - per_query[-1][1] >= min_length:
            unmap_region.append((per_query[-1][1], query_length_dict[per_name],
                                 per_query[-1][3], 'none', per_query[-1][5], 'none'))
        total_unmap[per_name] = unmap_region
    return total_unmap


def shuchu_log(log_name, query_dict, query_length_dict, query_most_freq_chr, total_merge, total_unmap):
    shuchu_contigs = open('%s_contigs.log' % log_name, 'w')
    for n, per in query_dict.items():
        shuchu_contigs.write(n + '\t' + query_most_freq_chr[n] + '\t' + str(query_length_dict[n]) + '\n')
        for m in per:
            shuchu_contigs.write('\t'.join([str(m[0]), str(m[1]), m[2], m[3], str(m[4]), str(m[5])]) + '\n')
    shuchu_contigs.close()

    shuchu_merge = open('%s_merge_region.log' % log_name, 'w')
    for n, per in total_merge.items():
        shuchu_merge.write(n + '\t' + query_most_freq_chr[n] + '\t' + str(query_length_dict[n]) + '\n')
        for m in per:
            shuchu_merge.write('\t'.join([str(m[0]), str(m[1]), m[2], m[3], str(m[4]), str(m[5])]) + '\n')
    shuchu_merge.close()

    shuchu_unmap = open('%s_unmap_region.log' % log_name, 'w')
    for n, per in total_unmap.items():
        shuchu_unmap.write(n + '\t' + query_most_freq_chr[n] + '\t' + str(query_length_dict[n]) + '\n')
        for m in per:
            shuchu_unmap.write('\t'.join([str(m[0]), str(m[1]), m[2], m[3], str(m[4]), str(m[5])]) + '\n')
    shuchu_unmap.close()


def read_all_contig(unmap_out_path, contigs_path, total_unmap, query_most_freq_chr, query_length_dict, min_length):
    whole_align = 0
    partly_align = 0
    partly_align_base = 0
    non_align = 0
    non_align_base = 0
    all_contigs = SeqIO.to_dict(SeqIO.parse(contigs_path, 'fasta'))
    shuchu = []
    for query_name, per_query in total_unmap.items():
        for region in per_query:
            temp = all_contigs[query_name][region[0]-1:region[1]-1]
            temp.id = '%s_%s_%s_%s_%s_%s_%s_%s_%s' \
                      % (query_name, query_most_freq_chr[query_name], str(query_length_dict[query_name]),
                         str(region[0]), str(region[1]), region[2], region[3], str(region[4]), str(region[5]))
            temp.name = ''
            temp.description = ''
            shuchu.append(temp)
            temp_length = region[1] - region[0]
            partly_align_base += temp_length
        if len(per_query) == 0:
            whole_align += 1
            # print('There are no unmaped regions: %s' % query_name)
        else:
            partly_align += 1
    print('Completely alignment: %s' % whole_align)
    print('Partly alignment: %s, %s bp' % (partly_align, partly_align_base))
    for query_name, per_query in all_contigs.items():
        if query_name not in total_unmap:
            if len(per_query) >= min_length:
                per_query.id = '%s_none_%s_none_none_none_none_none_none' % (per_query.id, str(len(per_query)))
                per_query.name = ''
                per_query.description = ''
                shuchu.append(per_query)
                non_align += 1
                non_align_base += len(per_query)
    SeqIO.write(shuchu, unmap_out_path, 'fasta')
    print('No alignment at all: %s, %s bp' % (non_align, non_align_base))


if __name__ == '__main__':
    main()
