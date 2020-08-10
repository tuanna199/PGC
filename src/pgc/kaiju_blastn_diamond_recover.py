import argparse
from collections import defaultdict, Counter
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.07.22'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='filter contaminants base on kaiju result '
                                             'and recover base on blastn or diamond')
parser.add_argument('-dia', '--diamond_result', type=str, help='path to the result of diamond')
parser.add_argument('-bo', '--blastn_result', type=str, help='path to the blastn output')
parser.add_argument('-ko', '--kaiju_result', type=str, help='path to the kaiju output')
parser.add_argument('-taxid', '--taxonomyID', type=str, help='path to the taxonomy ID list')
parser.add_argument('-i', '--input', type=str, help='path to the input fasta')
parser.add_argument('-o', '--output', type=str, help='path to the output fasta')


def main():
    args = parser.parse_args()
    diamond_result = args.diamond_result
    blastn_result = args.blastn_result
    kaiju_result = args.kaiju_result
    taxonomyid = args.taxonomyID
    p_input = args.input
    p_output = args.output

    sample = p_input.split('/')[-1].split('_')[0]
    chordata_taxid = chordata_gain(taxonomyid)

    kaiju_shuchu = read2kaiju(kaiju_result, sample, chordata_taxid)
    blastn_shuchu = read2blastn(blastn_result, sample, chordata_taxid)
    diamond_shuchu = read2diamond(diamond_result, sample, chordata_taxid)

    temp = kaiju_shuchu - blastn_shuchu - diamond_shuchu
    # print(temp)
    # print(len(temp))
    # print(len(kaiju_shuchu))
    # print(len(blastn_shuchu))
    # print(len(diamond_shuchu))
    # print(diamond_shuchu)

    seq_dict = SeqIO.to_dict(SeqIO.parse(p_input, 'fasta'))
    shuchu = []
    for i in seq_dict:
        if i not in temp:
            shuchu.append(seq_dict[i])
    SeqIO.write(shuchu, p_output, 'fasta')


def read2diamond(diamond_result, sample, chordata_taxid):
    contam_dict = defaultdict(list)
    with open(diamond_result, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if sample in line[0]:
                if int(line[4]) > 100 and float(line[5]) > 60:
                    if line[13] in chordata_taxid:
                        contam_dict[line[0]].append((line[12], line[13]))
    comtam_list = []
    for i in contam_dict:
        if len(contam_dict[i]) >= 2:
            comtam_list.append(i)
    return set(comtam_list)


def read2blastn(blastn_result, sample, chordata_taxid):
    taxid_dict = defaultdict(list)
    f_blastn = open(blastn_result, 'r')
    error_taxid = ['111789']
    for line in f_blastn:
        line = line.strip('\n').split('\t')
        name = line[0]
        if sample in name:
            if line[17] not in error_taxid and int(line[4]) > 100:
                # Eukaryotic synthetic construct
                taxid_dict[name].append(line[17])
    f_blastn.close()
    taxid_most = {}
    for name, per in taxid_dict.items():
        temp = Counter(per).most_common(1)
        taxid_most[name] = temp[0][0]
    total_shuchu = []
    for name, taxid in taxid_most.items():
        if taxid in chordata_taxid:
            total_shuchu.append(name)
    return set(total_shuchu)


def read2kaiju(kaiju_result, sample, chordata_taxid):
    taxid_dict = {}
    f_kaiju = open(kaiju_result, 'r')
    error_taxid = ['111789']
    for line in f_kaiju:
        line = line.strip('\n').split('\t')
        if line[0] == 'C':
            name = line[1]
            if sample in name:
                if line[2] not in error_taxid and int(line[3]) >= 100:
                    # Eukaryotic synthetic construct
                    taxid_dict[name] = line[2]
    f_kaiju.close()
    total_shuchu = []
    for name, taxid in taxid_dict.items():
        if taxid not in chordata_taxid:
            total_shuchu.append(name)
    return set(total_shuchu)


def chordata_gain(taxonomyid):
    chordata_taxid = []
    f_chordata = open(taxonomyid, 'r')
    for line in f_chordata:
        line = line.strip('\n')
        chordata_taxid.append(line)
    f_chordata.close()
    return chordata_taxid


if __name__ == '__main__':
    main()

