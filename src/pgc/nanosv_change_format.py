import argparse


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2020.05.01'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='NanoSV vcf re-format')
parser.add_argument('-raw', '--raw_vcf', type=str, default='', help='raw vcf of NanoSV')
parser.add_argument('-new', '--new_vcf', type=str, default='', help='vcf of new format for clique_SV')
parser.add_argument("-s", "--score", default=0, help="The score value threshold for SV.")


def main():
    args = parser.parse_args()
    p_vcf = args.raw_vcf
    new_vcf = args.new_vcf
    score = args.score
    read2vcf(p_vcf, new_vcf, score)


def read2vcf(p_vcf, new_vcf, score):
    score = float(score)

    f_new_vcf = open(new_vcf, 'w')
    with open(p_vcf, 'r') as f:
        for line_raw in f:
            if line_raw.startswith('#'):
                f_new_vcf.write(line_raw)
            else:
                line = line_raw.strip('\n').split('\t')
                scoreValue = line[5]
                if scoreValue == ".":
                    scoreValue = 0
                scoreValue = float(scoreValue)
                ### filt the record using the score threshold
                if scoreValue >= score:
                    
                    info_list = line[7].split(';')

                    if 'IMPRECISE' in info_list:
                        info_list.remove('IMPRECISE')
                    info_dict = {}
                    for i in info_list:
                        ni = i.split('=')
                        info_dict[ni[0]] = ni[1]

                    if '<' in line[4]:
                        chr2 = line[0]
                        strands = '+-'
                        format_genotype = 'GT:DR:DV'
                        genotype = line[-1].split(':')[0]
                        if genotype == './.':
                            genotype = '0/0'
                        dr = line[-1].split(':')[1].split(',')[0]
                        dv = line[-1].split(':')[2].split(',')[0]

                        new_info = ['SVTYPE=%s' % info_dict['SVTYPE'], 'END=%s' % info_dict['END'],
                                    'SVLEN=%s' % info_dict['SVLEN'], 'CHR2=%s' % chr2, 'STRANDS=%s' % strands]
                        new_line = line[0:7]
                        new_line.append(';'.join(new_info))
                        new_line.append(format_genotype)
                        new_line.append('%s:%s:%s' % (genotype, dr, dv))
                        # print(new_line)
                        f_new_vcf.write('\t'.join(new_line) + '\n')
                    else:
                        svtype = 'TRA'
                        svlen = '0'
                        if '[' in line[4]:
                            chr2 = line[4].split('[')[1].split(':')[0]
                            end = line[4].split('[')[1].split(':')[1]
                            if chr2 == line[0]:
                                if line[4].split('[')[-1] in ['A', 'T', 'C', 'G', 'N']:
                                    svtype = 'INV'
                                    svlen = info_dict['SVLEN']
                                if line[4].split('[')[0] in ['A', 'T', 'C', 'G', 'N']:
                                    svtype = 'DEL'
                                    svlen = '-' + info_dict['SVLEN']

                        if ']' in line[4]:
                            chr2 = line[4].split(']')[1].split(':')[0]
                            end = line[4].split(']')[1].split(':')[1]
                            if chr2 == line[0]:
                                if line[4].split(']')[0] in ['A', 'T', 'C', 'G', 'N']:
                                    svtype = 'INV'
                                    svlen = info_dict['SVLEN']
                                if line[4].split(']')[-1] in ['A', 'T', 'C', 'G', 'N']:
                                    svtype = 'DUP'
                                    svlen = info_dict['SVLEN']

                        strands = '+-'
                        format_genotype = 'GT:DR:DV'
                        genotype = line[-1].split(':')[0]
                        if genotype == './.':
                            genotype = '0/0'
                        dr = line[-1].split(':')[1].split(',')[0]
                        dv = line[-1].split(':')[2].split(',')[0]

                        new_info = ['SVTYPE=%s' % svtype, 'END=%s' % end, 'SVLEN=%s' % svlen, 'CHR2=%s' % chr2,
                                    'STRANDS=%s' % strands]

                        new_line = line[0:7]
                        new_line[4] = '<%s>' % svtype
                        new_line.append(';'.join(new_info))
                        new_line.append(format_genotype)
                        new_line.append('%s:%s:%s' % (genotype, dr, dv))
                        # print(new_line)
                        f_new_vcf.write('\t'.join(new_line) + '\n')

    f_new_vcf.close()


if __name__ == '__main__':
    main()
