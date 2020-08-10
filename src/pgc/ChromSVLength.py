#!/usr/bin/python
from BaseFunc import Infor_target_values
import collections
import argparse


#usage: python ~/github/NanoHub/src/NanoHub/ChromSVLength.py --vcf M671-2_read_depth_filt.vcf --out temp.xls


def Chr_SV_length(vcf_file, out_file, cutLength, LenTag):
    TypeList = ["DEL", "INS", "INV", "DUP"]
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

    ChrTypeLength = collections.defaultdict(lambda: collections.Counter())
    cutLength = int(cutLength)

    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start = lines[:2]
            Infor = lines[7]
            Eed = Infor_target_values(Infor, "END")
            SVType = Infor_target_values(Infor, "SVTYPE")

            ### length tag should be "SVLEN" or "AVGLEN" (SURVIVOR merged)
            SVLength = Infor_target_values(Infor, LenTag)
            # SVLength = abs(int(SVLength))
            SVLength = abs(float(SVLength))

            if SVLength > cutLength:
                print(line)
                continue
            else:
                if SVType in TypeList:
                    if Chr in CHRS:
                        ChrTypeLength[Chr][SVType] += SVLength
    vcf_h.close()

    ### output the summary
    out_h = open(out_file, "w")
    out_h.write("Chr\t%s\tTotal\n" % "\t".join(TypeList))

    TotalLenth = collections.defaultdict(list)
    for c in CHRS:
        TypeLEN = 0
        TypeLenList = []
        for sv in TypeList:
            length = ChrTypeLength[c][sv]
            TypeLenList.append(length)
            TypeLEN += length

            TotalLenth[sv].append(length)

        TypeLenList = [str(int(s)) for s in TypeLenList]
        out_h.write("%s\t%s\t%d\n" % (c, "\t".join(TypeLenList), TypeLEN))


    TotalLENs =[]
    for t in TypeList:
        lens = TotalLenth[t]
        totalLen = sum(lens)
        TotalLENs.append(totalLen)

    TotalSUM = sum(TotalLENs)
    TotalLENs = [str(int(t)) for t in TotalLENs]

    sample = vcf_file.split("/")[-1].split("_")[0].split(".")[0]
    out_h.write("%s\t%s\t%d\n" % (sample, "\t".join(TotalLENs), TotalSUM))

    out_h.close()

def main():
    parser = argparse.ArgumentParser(description='Get the summary for SV and length for each chromosome.')
    parser.add_argument('-v', '--vcf', help='The input vcf file.')
    parser.add_argument("-l", "--cutLength", default=100000000, help="Drop the length above this threshould.")
    parser.add_argument("-t", "--LenTag", default="SVLEN", help="The tag name for SV length, should be 'SVLEN' or 'AVGLEN'.")
    parser.add_argument('-o', '--out', help='The output file.')

    args = parser.parse_args()
    Chr_SV_length(args.vcf, args.out, args.cutLength, args.LenTag)



if __name__ == '__main__':
    main()








                