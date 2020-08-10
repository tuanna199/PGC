import re
import argparse


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.04.14'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='evaluate protein coding gene in quast')
parser.add_argument('-q', '--quast_file', type=str, help='path to the file of quast')
parser.add_argument('-gff', '--protein_gff', type=str, help='path to the gff of protein')


def main():
    args = parser.parse_args()
    p_quast = args.quast_file
    p_protein = args.protein_gff
    id2complete = {}
    id2ctg = {}
    id2length = {}
    with open(p_quast, 'r') as f:
        for line in f:
            if not line.startswith('ID or #') and not line.startswith('===='):
                line = line.strip('\n').split('\t')
                # empty \t line[4]
                id2complete[line[0]] = line[4]
                id2ctg[line[0]] = line[-1]
                id2length[line[0]] = int(line[3]) - int(line[2]) + 1

    pattern_num = re.compile(r'\d-')
    n = 0
    m = 0
    partial = []
    nottiall = []
    with open(p_protein, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            info = line[8].split(';')
            temp = {}
            for per in info:
                per = per.split('=')
                temp[per[0]] = per[1]
            geneid = temp['ID']
            if geneid in id2complete:
                n += 1
                if id2complete[geneid] == 'complete':
                    m += 1
                if id2complete[geneid] == 'partial':
                    map_len = 0
                    for per in id2ctg[geneid].split(','):
                        per = per.split(':')[1]
                        match = list(re.finditer(pattern_num, per))

                        qian = int(per[0:match[0].end() - 1])
                        hou = int(per[match[0].end():])

                        map_len += abs(qian - hou + 1)
                    ratio = map_len / id2length[geneid]
                    if ratio >= 0.9:
                        m += 1
                    else:
                        partial.append(geneid)
            else:
                nottiall.append(geneid)
    # print(nottiall)
    # print(len(nottiall))
    # print(partial)
    # print(len(partial))
    # print(n, m, m/19982*100)
    print('%s\t%.2f' % (p_quast.split('/')[-1], m/19982*100))


if __name__ == '__main__':
    main()

