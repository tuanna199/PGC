#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
from BaseFunc import SV_tag
from BaseFunc import Format_tag
import argparse
import sys
import math
import collections
import operator
from tinyfasta import FastaParser

#usage: python ~/github/NanoHub/src/NanoHub/CombineSVSeqID.py  --vcf Sample_common_SV_convert.vcf --sequence /home/wuzhikun/Project/NanoTrio/SVCall/SeqSniffles/M671-2_filt_length.fasta  --out temp.txt > temp

# def Format_tag(Format, geno):
#     """
#     Format: GT:PSV:LN:DR:ST:TY:CO
#     geno: 1/1:NA:137:1,36:+-:INS,INS:1_136673-1_136891,1_136999-1_137209

#     ['NaN-0-NaN', '1_10169-X_449438-0-TRA']
#     """
#     LN = parse_genotype_format(Format, geno, "LN")
#     TY = parse_genotype_format(Format, geno, "TY")
#     CO = parse_genotype_format(Format, geno, "CO")
#     COs = CO.split(",")
#     TYs = TY.split(",")

#     Tags = []
#     COLen = len(CO)

#     ### set length of translocation as 0
#     if TY == "TRA":
#         LN = "0"

#     if COLen == 1:
#         tag = "%s-%s-%s" % (CO, LN, TY)
#         Tags.append(tag)
#     elif COLen > 1:
#         for c, t in zip(COs, TYs):
#             tag = "%s-%s-%s" % (c, LN, t)
#             Tags.append(tag)

#     TagLen = len(Tags)
#     if TagLen == 1:
#         return Tags[0]
#     else:
#         return Tags


def tag_SV_length(GenoTags):
    """
    GenoTags:
    ['1_90320-1_90432-111-INS', '1_90312-1_90395-56-INS', '1_90387-1_90456-119-INS', 'NaN-0-NaN']
    """
    LENS = []
    for t in GenoTags:
        if t == "NaN-0-NaN" or t == "-" or isinstance(t, list):
            l = 0
        else:
            l = int(t.split("-")[2])
        LENS.append(l)
    return LENS


def longest_index(LENS):
    ValueIndex = {}
    for i, v in enumerate(LENS):
        ValueIndex[v] = i
    SortedValues = sorted(ValueIndex.items(), key=operator.itemgetter(0), reverse=True)
    target = SortedValues[0][0]
    targetIndex = SortedValues[0][1]

    return target, targetIndex




def vcf_position_sample(vcf_file):
    TagSampleTag = collections.defaultdict(dict)
    SampleTags = collections.defaultdict(dict)

    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            lines = line.split("\t")
            samples = lines[9:]
            sampleNames = [s.split("/")[-1].split(".")[0] for s in samples]
        else:
            lines = line.split("\t")
            Infor = lines[7]
            Format = lines[8]
            Genos = lines[9:]
            SVType = Infor_target_values(Infor, "SVTYPE")
            GenoTags = [Format_tag(Format, geno) for geno in Genos]
            Tag = SV_tag(lines)
            if SVType == "INS" or SVType == "DEL":
                LENS = tag_SV_length(GenoTags)
                target, targetIndex = longest_index(LENS)
                targetSample = sampleNames[targetIndex]
                targetTag = GenoTags[targetIndex]

                TagSampleTag[Tag][targetSample] = targetTag
                SampleTags[targetSample][targetTag] = 1
    vcf_h.close()
    return TagSampleTag, SampleTags


def tag_target_sample_sequence(vcf_file, seq_files, out_file):
    TagSampleTag, SampleTags = vcf_position_sample(vcf_file)

    SampleTagSeq = collections.defaultdict(dict)

    files = seq_files.split(",")
    files = [f.strip() for f in files]
    fileNames = [f.split("/")[-1].split("_")[0] for f in files]
    for s, f in zip(fileNames, files):
        if s in SampleTags:
            Tags = SampleTags[s]

            getTags = []

            ### parse the fasta file 
            for record in FastaParser(f):
                desc = str(record.description)
                tag = desc.lstrip(">")
                if tag in Tags:
                    seq = str(record.sequence)
                    SampleTagSeq[s][tag] = seq
                    getTags.append(tag)


            ### whether all tags have sequence
            noTags = set(Tags) - set(getTags)
            if len(noTags) != 0:
                print("Please check whether the tags %s can find the sequences for sample %s." % (noTags, s))
        else:
            print("Please check whether the sample %s exist in file %s." % (s, vcf_file))
            sys.exit(1)



    out_h = open(out_file, "w")

    for T in TagSampleTag:
        for s in TagSampleTag[T]:
            t = TagSampleTag[T][s]
            if t in SampleTagSeq[s]:
                seq = SampleTagSeq[s][t]
                out_h.write(">%s:%s:%s\n%s\n" % (T, s, t, seq))
            else:
                print("Please check whether the sequence of tag %s is in file %s for merged Tag %s of file %s." % (t, s, T, vcf_file))
    out_h.close()







def main():
    parser = argparse.ArgumentParser(description="Filt SV file based on the number of support records.")
    parser.add_argument("-v", "--vcf", help="The input sv file with vcf format.")
    parser.add_argument("-s", "--sequence", help="The sequence file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    tag_target_sample_sequence(args.vcf, args.sequence, args.out)



if __name__ == "__main__":
    main()
