import argparse
from Bio import SeqIO
import gzip


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.04.04'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='extract flank sequences')
parser.add_argument('-gff', '--gff_path', type=str, help='path to the gff file')
parser.add_argument('-i', '--unmap_path', type=str, help='path to the unmap fasta')
parser.add_argument('-ctg', '--ctg_path', type=str, help='path to the ctg')
parser.add_argument('-ol', '--left_out_path', type=str, help='path to the left flank out')
parser.add_argument('-or', '--right_out_path', type=str, help='path to the right flank out')


def main():
    args = parser.parse_args()
    p_fa = args.unmap_path
    p_gff = args.gff_path
    p_ctg = args.ctg_path
    left_out_path = args.left_out_path
    right_out_path = args.right_out_path
    flank_dict, partly_align = read2parse(p_gff)
    extract_seq(p_fa, p_ctg, flank_dict, partly_align, left_out_path, right_out_path)


def read2parse(p_gff):
    partly_align = {}
    flank_dict = {}
    with open(p_gff, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            name = line[0]
            ctg_name = line[1]
            ctg_start = int(line[3])
            ctg_end = int(line[4])
            ctg_len = int(line[2])
            partly_align[name] = [ctg_name, ctg_len, ctg_start, ctg_end]

    for per, info in partly_align.items():
        if info[2] - 1 >= 10000:
            left_flank = [info[2] - 10000, info[2] - 1]
        else:
            left_flank = [1, max(info[2] - 1, 1)]
        if info[1] - info[3] >= 10000:
            right_flank = [info[3] + 1, info[3] + 10000]
        else:
            right_flank = [min(info[3] + 1, info[1]), info[1]]
        flank_dict[per] = [left_flank, right_flank]
    return flank_dict, partly_align


def extract_seq(p_fa, p_ctg, flank_dict, partly_align, left_out_path, right_out_path):
    shuchu_left = []
    shuchu_right = []
    unmap_fa = SeqIO.to_dict(SeqIO.parse(p_fa, 'fasta'))
    p_ctg = gzip.open(p_ctg, 'rt')
    ctg_fa_dict = SeqIO.to_dict(SeqIO.parse(p_ctg, 'fasta'))
    p_ctg.close()
    for name in unmap_fa:
        name = name.split('_')
        name = '_'.join(name[0:3])
        if name in partly_align:
            ctg_name = partly_align[name][0]
            left_region = flank_dict[name][0]
            # print(name, ctg_name, left_region, ctg_path)
            if left_region[0] != left_region[1]:
                if abs(left_region[1] - left_region[0] + 1) >= 5000:
                    left_flank_seq = ctg_fa_dict[ctg_name][left_region[0] - 1:left_region[1]]
                    left_flank_seq.id = '%s_flank_left' % name
                    left_flank_seq.name = ''
                    left_flank_seq.description = '%s_%s_%s' % (ctg_name, left_region[0], left_region[1])
                    shuchu_left.append(left_flank_seq)

            right_region = flank_dict[name][1]
            if right_region[0] != right_region[1]:
                if abs(right_region[1] - right_region[0] + 1) >= 5000:
                    right_flank_seq = ctg_fa_dict[ctg_name][right_region[0] - 1:right_region[1]]
                    right_flank_seq.id = '%s_flank_right' % name
                    right_flank_seq.name = ''
                    right_flank_seq.description = '%s_%s_%s' % (ctg_name, right_region[0], right_region[1])
                    shuchu_right.append(right_flank_seq)
    SeqIO.write(shuchu_left, left_out_path, 'fasta')
    SeqIO.write(shuchu_right, right_out_path, 'fasta')


if __name__ == '__main__':
    main()

