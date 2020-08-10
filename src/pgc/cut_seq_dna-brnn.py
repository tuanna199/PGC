import argparse
from collections import defaultdict
from Bio import SeqIO
# import pprint


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.21'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='extract non-alpha and non-hsat2,3')
parser.add_argument('-bed', '--bed_dna_brnn', type=str, help='path to the bed from dna-brnn')
parser.add_argument('-NRS_fasta', '--raw_NRS_fasta', type=str, help='path to the raw NRS')
parser.add_argument('-out_NRS', '--out_NRS', type=str, help='path to the NRS non-alpha and non-hsat2,3')


def main():
    args = parser.parse_args()
    bed_dna_brnn = args.bed_dna_brnn
    raw_nrs_fasta = args.raw_NRS_fasta
    out_nrs = args.out_NRS

    repeat_dict = defaultdict(list)
    sample_name = raw_nrs_fasta.split('/')[-1].split('_')[0]
    with open(bed_dna_brnn, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if sample_name in line[0]:
                repeat_dict[line[0]].append((int(line[1]), int(line[2])))

    seq_dict = SeqIO.to_dict(SeqIO.parse(raw_nrs_fasta, 'fasta'))
    total_unrep = {}
    for per_name, per_query in repeat_dict.items():
        unrep_region = []
        if per_query[0][0] >= 1000:
            unrep_region.append((0, per_query[0][0]))
        for idx in range(len(per_query) - 1):
            region1 = per_query[idx]
            region2 = per_query[idx + 1]
            if region2[0] - region1[1] - 1 >= 1000:
                unrep_region.append((region1[1] + 1, region2[0]))
        if len(seq_dict[per_name].seq) - per_query[-1][1] - 1 >= 1000:
            unrep_region.append((per_query[-1][1] + 1, len(seq_dict[per_name].seq)))
        total_unrep[per_name] = unrep_region
    # pprint.pprint(total_unrep)

    shuchu = []
    for i in seq_dict:
        if i not in total_unrep:
            shuchu.append(seq_dict[i])
        else:
            n = 0
            for per in total_unrep[i]:
                n += 1
                temp = seq_dict[i][per[0]:per[1]]
                temp.id = '%s_%s' % (i, n)
                temp.name = ''
                temp.description = ''
                shuchu.append(temp)
    SeqIO.write(shuchu, out_nrs, 'fasta')


if __name__ == '__main__':
    main()

