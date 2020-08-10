from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='classify non redundant')
parser.add_argument('-all_nrs', '--all_nrs_fasta', type=str, help='path to the total NRS')
parser.add_argument('-all_bed', '--all_anchor_bed', type=str, help='path to the all bed file')
parser.add_argument('-two', '--two_end', type=str, help='path to the two end of NRS')
parser.add_argument('-one', '--one_end', type=str, help='path to the one end of NRS')
parser.add_argument('-unplaced', '--unplaced_fa', type=str, help='path to the unplaced of NRS')


def main():
    args = parser.parse_args()
    p_fasta = args.all_nrs_fasta
    p_bed = args.all_anchor_bed
    two_end = set()
    one_end = set()
    with open(p_bed, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if line[1] != '.' and line[2] != '.':
                two_end.add(line[3])
            else:
                one_end.add(line[3])
    shuchu_two = []
    shuchu_one = []
    shuchu_unp = []
    for seq in SeqIO.parse(p_fasta, 'fasta'):
        if seq.id in two_end:
            shuchu_two.append(seq)
        elif seq.id in one_end:
            shuchu_one.append(seq)
        else:
            shuchu_unp.append(seq)
    SeqIO.write(shuchu_two, args.two_end, 'fasta')
    SeqIO.write(shuchu_one, args.one_end, 'fasta')
    SeqIO.write(shuchu_unp, args.unplaced_fa, 'fasta')


if __name__ == '__main__':
    main()

