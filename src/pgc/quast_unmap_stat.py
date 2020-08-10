import numpy as np
import argparse
from Bio import SeqIO
import gzip


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.29'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='extract unmap using quast')
parser.add_argument('-quast_unalign_info', '--quast_unalign_info', type=str, help='quast unalign info')
parser.add_argument('-assembly_fasta', '--assembly_fasta', type=str, help='assembly genome from wtdbg2')
parser.add_argument('-out_NRS', '--out_NRS', type=str, help='path of NRS fasta')
parser.add_argument('-out_gff_NRS', '--out_gff_NRS', type=str, help='path of NRS gff')
parser.add_argument('-min_contigs_len', '--min_contigs_len', type=str, help='min contigs length in assembly')


def main():
    args = parser.parse_args()
    p_quast_unmap = args.quast_unalign_info
    assembly_fasta = args.assembly_fasta
    out_nrs = args.out_NRS
    out_gff_nrs = args.out_gff_NRS
    min_contigs_len = int(args.min_contigs_len)
    length_list = []
    info_list = []
    with open(p_quast_unmap, 'r') as f:
        for line in f:
            if not line.startswith('Contig'):
                line = line.strip('\n').split('\t')
                info = line[4].split(',')
                for per in info:
                    per = per.split('-')
                    length = int(per[1]) - int(per[0]) + 1
                    if length >= 1000:
                        length_list.append(length)
                        info_list.append((line[0], int(per[0]), int(per[1])))
                        # print(line[0], length, per[0], per[1])
    print(np.array(length_list).mean())
    print(sum(length_list))
    cacul_n50(sum(length_list), length_list)

    n = 0
    shuchu = []
    sample_name = assembly_fasta.split('/')[-1].split('_')[0]
    assembly_fasta = gzip.open(assembly_fasta, 'rt')
    seq_dict = SeqIO.to_dict(SeqIO.parse(assembly_fasta, 'fasta'))
    assembly_fasta.close()
    f_shuchu = open(out_gff_nrs, 'w')
    for i in info_list:
        if len(seq_dict[i[0]].seq) >= min_contigs_len:
            n += 1
            temp = seq_dict[i[0]][i[1]-1:i[2]-1]
            temp.id = '%s_NRS_%s' % (sample_name, n)
            temp.name = ''
            temp.description = ''
            shuchu.append(temp)
            f_shuchu.write('\t'.join([temp.id, i[0], str(len(seq_dict[i[0]].seq)),
                                      str(i[1]), str(i[2]), str(i[2] - i[1] + 1)]) + '\n')
    SeqIO.write(shuchu, out_nrs, 'fasta')
    f_shuchu.close()


def cacul_n50(total_length, length_list):
    n50_pos = total_length / 2
    length_list = sorted(length_list, reverse=True)
    temp = 0
    n50 = 0
    for i in length_list:
        temp += i
        if n50_pos <= temp:
            n50 = i
            break
    # print(length_list)
    print(n50)


if __name__ == '__main__':
    main()

