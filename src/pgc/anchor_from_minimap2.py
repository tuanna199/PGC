import argparse
from collections import defaultdict
import re


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.22'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='anchor CPG by minimap2 of flanking sequence')
parser.add_argument('-left', '--left_minimap2', type=str, help='path to the left minimap2')
parser.add_argument('-right', '--right_minimap2', type=str, help='path to the right minimap2')
parser.add_argument('-bed', '--bed_output', type=str, help='path to the bed of coordinate')


def main():
    args = parser.parse_args()
    p_left = args.left_minimap2
    p_right = args.right_minimap2
    p_bed = args.bed_output
    best_match_left = read2minimap2(p_left)
    best_match_right = read2minimap2(p_right)
    anchor, non_anchor = combine(best_match_left, best_match_right)
    bed2output(p_bed, anchor, non_anchor)


def read2minimap2(p_minimap2):
    minimap2_rec_dict = defaultdict(list)
    pattern_dv = re.compile(r'dv:f:(.+)\trl')
    with open(p_minimap2, 'r') as f:
        for line in f:
            div_re = re.findall(pattern_dv, line)
            div_re = float(div_re[0])
            if div_re <= 0.1:
                line = line.strip('\n').split('\t')
                strand = line[4]
                length = abs(int(line[3]) - int(line[2]))
                # fixed parameter
                if length >= 5000:
                    cov = int(length / int(line[1]) * 100)
                    if strand == '+':
                        minimap2_rec_dict[line[0]].append((cov, (1-div_re)*100,
                                                           line[5], int(line[7]), int(line[8]), strand))
                    else:
                        minimap2_rec_dict[line[0]].append((cov, (1 - div_re) * 100,
                                                           line[5], int(line[8]), int(line[7]), strand))
    best_match = {}
    for name, per in minimap2_rec_dict.items():
        per = sorted(per, key=lambda x: (100-x[0], 100-x[1], x[2]))
        name = name.split('_')[:3]
        name = '_'.join(name)
        best_match[name] = per[0]
    return best_match


def combine(best_match_left, best_match_right):
    anchor = {}
    all_name = set(best_match_left.keys()) | set(best_match_right.keys())
    print(len(all_name))
    n = 0
    for i in all_name:
        if i in best_match_left.keys() and i in best_match_right.keys():
            left = best_match_left[i]
            right = best_match_right[i]
            internal = abs(right[3] - left[4])
            # fixed parameter
            if left[5] == right[5] and left[2] == right[2] and internal <= 500:
                anchor[i] = (left[2], left[4], right[3], left[5])
            else:
                # print(i, internal, left, right)
                n += 1
    print(n)
    non_anchor = {}
    for i in all_name:
        if i not in anchor:
            if i in best_match_left.keys() and i in best_match_right.keys():
                temp = [best_match_left[i], best_match_right[i]]
                idx = sorted([0, 1], key=lambda x: (100-temp[x][0], 100-temp[x][1], temp[x][2]))
                if idx[0] == 0:
                    left = temp[0]
                    non_anchor[i] = (left[2], left[4], '.', left[5])
                else:
                    right = temp[1]
                    non_anchor[i] = (right[2], '.', right[3], right[5])
            if i in best_match_left.keys() and i not in best_match_right.keys():
                left = best_match_left[i]
                non_anchor[i] = (left[2], left[4], '.', left[5])
            if i not in best_match_left.keys() and i in best_match_right.keys():
                right = best_match_right[i]
                non_anchor[i] = (right[2], '.', right[3], right[5])
    print(len(anchor))
    # print(non_anchor)
    return anchor, non_anchor


def bed2output(p_bed, anchor, non_anchor):
    f_bed = open(p_bed, 'w')
    for name, per in anchor.items():
        temp = [per[0], str(per[1] - 1), str(per[2] - 1), name, '.', per[3]]
        f_bed.write('\t'.join(temp) + '\n')
    for name, per in non_anchor.items():
        if per[1] == '.':
            temp = [per[0], per[1], str(per[2] - 1), name, '.', per[3]]
            f_bed.write('\t'.join(temp) + '\n')
        if per[2] == '.':
            temp = [per[0], str(per[1] - 1), per[2], name, '.', per[3]]
            f_bed.write('\t'.join(temp) + '\n')
    f_bed.close()


if __name__ == '__main__':
    main()

