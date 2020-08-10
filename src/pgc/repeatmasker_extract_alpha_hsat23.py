import argparse
import re


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.20'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='convert alpha and hsat23 coordinate from repeatmasker')
parser.add_argument('-re_in', '--repeatmasker_in', type=str, help='path to the out file of repeatmasker')
parser.add_argument('-re_bed', '--out_bed', type=str, help='path to the bed file of result')


def main():
    args = parser.parse_args()
    repeatmasker_in = args.repeatmasker_in
    out_bed = args.out_bed
    f_shuchu = open(out_bed, 'w')

    alpha_hsat23_feat = ['HSATII',
                         '(CATTC)n', '(ATTCC)n', '(TTCCA)n', '(TCCAT)n', '(CCATT)n',
                         '(GAATG)n', '(AATGG)n', '(ATGGA)n', '(TGGAA)n', '(GGAAT)n']
    pattern_re = re.compile(r'[^ ]+')
    n = 0
    with open(repeatmasker_in, 'r') as f:
        for line in f:
            n += 1
            if n > 3:
                newline = re.findall(pattern_re, line)
                if newline[9] in alpha_hsat23_feat:
                    f_shuchu.write('\t'.join([newline[4], str(int(newline[5]) - 1), 
                                              str(int(newline[6]) - 1), '1']) + '\n')
                elif newline[9] == 'ALR/Alpha':
                    f_shuchu.write('\t'.join([newline[4], str(int(newline[5]) - 1), 
                                              str(int(newline[6]) - 1), '2']) + '\n')
    f_shuchu.close()


if __name__ == '__main__':
    main()

