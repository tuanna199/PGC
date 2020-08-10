from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='separate placed and non-placed')
parser.add_argument('-merge_fa', '--merge_fa', type=str, help='merge NRS fasta')
parser.add_argument('-merge_bed', '--merge_bed', type=str, help='merge NRS anchor bed')
parser.add_argument('-can_anchor', '--can_anchor', type=str, help='the placed fasta')
parser.add_argument('-cant_anchor', '--cant_anchor', type=str, help='the unplaced fasta')


def main():
    args = parser.parse_args()
    p_fasta = args.merge_fa
    p_bed = args.merge_bed
    anchor = set()
    with open(p_bed, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            anchor.add(line[3])
    shuchu = []
    shuchu2 = []
    for seq in SeqIO.parse(p_fasta, 'fasta'):
        if seq.id in anchor:
            shuchu.append(seq)
        else:
            shuchu2.append(seq)
    SeqIO.write(shuchu, args.can_anchor, 'fasta')
    SeqIO.write(shuchu2, args.cant_anchor, 'fasta')


if __name__ == '__main__':
    main()

