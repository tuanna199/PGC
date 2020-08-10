from collections import defaultdict
import os
import argparse
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.23'


# need bedtools in path
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='remove redundancy using minimap2.paf')
parser.add_argument('-paf', '--minimap2_paf', type=str, help='path to the minimap2.paf')
parser.add_argument('-paf_part', '--paf_part', type=str, help='path to the part of paf')
parser.add_argument('-i', '--contigs_path', type=str, help='path to the contigs')
parser.add_argument('-o', '--nonredun_fa', type=str, help='path to the non-redundancy')
parser.add_argument('-cl', '--clstr_out', type=str, help='path to the cluster record')
parser.add_argument('-mincov', '--min_coverage', type=float, help='minimum coverage reserved')
parser.add_argument('-mindiv', '--min_divergence', type=float, help='minimum divergence reserved')


def main():
    args = parser.parse_args()
    cov_threshold = args.min_coverage
    p_paf = args.minimap2_paf
    paf_part = args.paf_part
    out_clstr = args.clstr_out
    input_fa = args.contigs_path
    output_fa = args.nonredun_fa
    min_div = args.min_divergence
    # params
    total_dict, length_dict = read_paf(p_paf, min_div, paf_part)
    merge_dict = merge_record(total_dict)
    merge_dict = merge_region(merge_dict)
    ratio_dict = calcu_ratio(merge_dict, length_dict, cov_threshold)
    like_cdhit(ratio_dict, length_dict, out_clstr, input_fa, output_fa)


def default_fun():
    return [[], []]


def read_paf(p_paf, min_div, paf_part):
    f_shuchu = open(paf_part, 'w')
    total_dict = defaultdict(default_fun)
    length_dict = {}
    f_paf = open(p_paf, 'r')
    for rawline in f_paf:
        line = rawline.strip('\n').split('\t')
        info = float(line[15].split(':')[-1])
        if info <= min_div:
            length = int(line[3]) - int(line[2])
            if length >= 100:
                if line[0] != line[5]:
                    f_shuchu.write(rawline)
                    length_dict[line[0]] = int(line[1])
                    length_dict[line[5]] = int(line[6])
                    total_dict[(line[0], line[5])][0].append((int(line[2]), int(line[3])))
                    total_dict[(line[0], line[5])][1].append((int(line[7]), int(line[8])))
    f_paf.close()
    f_shuchu.close()
    # print(total_dict.keys())
    return total_dict, length_dict


def merge_record(total_dict):
    merge_dict = {}
    all_keys = total_dict.keys()
    for name12, per12 in total_dict.items():
        a_list = []
        b_list = []
        a_list.extend(per12[0])
        b_list.extend(per12[1])
        name21 = (name12[1], name12[0])
        if name21 in all_keys:
            per21 = total_dict[name21]
            a_list.extend(per21[1])
            b_list.extend(per21[0])
        if name12 not in merge_dict:
            if name21 not in merge_dict:
                merge_dict[name12] = [a_list, b_list]
    # print(merge_dict.keys())
    return merge_dict


def merge_region(merge_dict):
    temp1 = open('temp1.bed', 'w')
    temp2 = open('temp2.bed', 'w')
    for perid, per in merge_dict.items():
        per1 = sorted(per[0], key=lambda x: (x[0], x[1]))
        per2 = sorted(per[1], key=lambda x: (x[0], x[1]))
        name = 'fengefu'.join(perid)
        for i in per1:
            temp1.write('\t'.join([name, str(i[0]), str(i[1])]) + '\n')
        for i in per2:
            temp2.write('\t'.join([name, str(i[0]), str(i[1])]) + '\n')
    temp1.close()
    temp2.close()
    os.system('bedtools merge -i temp1.bed > temp_new1.bed')
    os.system('bedtools merge -i temp2.bed > temp_new2.bed')
    total_merge = defaultdict(default_fun)
    with open('temp_new1.bed', 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            name = tuple(line[0].split('fengefu'))
            total_merge[name][0].append((int(line[1]), int(line[2])))
    with open('temp_new2.bed', 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            name = tuple(line[0].split('fengefu'))
            total_merge[name][1].append((int(line[1]), int(line[2])))
    # print(total_merge)
    os.system('rm temp1.bed temp2.bed temp_new1.bed temp_new2.bed')
    return total_merge


def calcu_ratio(merge_dict, length_dict, cov_threshold):
    ratio_dict = defaultdict(list)
    for perid, per in merge_dict.items():
        if perid[0] != perid[1]:
            map_ratio_0 = sum([i[1] - i[0] + 1 for i in per[0]]) / length_dict[perid[0]]
            map_ratio_1 = sum([i[1] - i[0] + 1 for i in per[1]]) / length_dict[perid[1]]
            if map_ratio_0 >= cov_threshold:
                ratio_dict[perid[1]].append((perid[0], map_ratio_0))
            if map_ratio_1 >= cov_threshold:
                ratio_dict[perid[0]].append((perid[1], map_ratio_1))
    # print(ratio_dict)
    return ratio_dict


def like_cdhit(ratio_dict, length_dict, out_clstr, input_fa, output_fa):
    cluster_dict = {}
    wait_remove = set(length_dict.keys())
    while wait_remove:
        max_contig = find_max_contig(wait_remove, length_dict)
        temp_list = [i for i in ratio_dict[max_contig] if i[0] in wait_remove]
        # print(temp_list)
        # Tree structure, one level only
        subtree_list = []
        for per, _ in temp_list:
            for subtree in ratio_dict[per]:
                if subtree[0] in wait_remove:
                    subtree_list.append(subtree)
        temp_list_name = [i[0] for i in temp_list]
        for i in subtree_list:
            if i[0] not in temp_list_name and i[0] != max_contig:
                temp_list.append(i)
                temp_list_name.append(i[0])

        cluster_dict[max_contig] = temp_list
        for name, _ in temp_list:
            wait_remove.remove(name)
        wait_remove.remove(max_contig)
        print(len(wait_remove))
    print(len(cluster_dict))
    f_shuchu = open(out_clstr, 'w')
    n = 0
    for name, per in cluster_dict.items():
        f_shuchu.write('>Cluster %s\n' % n)
        n += 1
        m = 0
        f_shuchu.write('%s\t%snt, >%s... *\n' % (m, length_dict[name], name))
        for ip in per:
            m += 1
            f_shuchu.write('%s\t%snt, >%s... at ./%.2f%s\n' % (m, length_dict[ip[0]], ip[0], ip[1]*100, '%'))

    shuchu_seq = []
    for per_seq in SeqIO.parse(input_fa, 'fasta'):
        if per_seq.id in cluster_dict:
            shuchu_seq.append(per_seq)
    SeqIO.write(shuchu_seq, output_fa, 'fasta')
    f_shuchu.close()


def find_max_contig(wait_remove, length_dict):
    length_list = [(name, length) for name, length in length_dict.items() if name in wait_remove]
    length_list = sorted(length_list, key=lambda x: x[1], reverse=True)
    max_contig = length_list[0][0]
    return max_contig


if __name__ == '__main__':
    main()

