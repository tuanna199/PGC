#!/usr/bin/python
import collections
import argparse
import os
import sys


def tag_frequency(geno_freq):
    """
    Tag     Ref_freq        Alt_freq        MAF
    1_10059-1_10060-3242-INS        0.9975  0.0025  0.0025
    1_10380-1_10381-321-INS 0.9988  0.0012  0.0012
    1_10657-1_10760-78-INS  0.9988  0.0012  0.0012
    """

    FreqTags = collections.defaultdict(list)
    TagFreq = {}
    # FreqList = [0.01, 0.05]
    bed_h = open(geno_freq, "r")
    header = bed_h.readline()
    for line in bed_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        freq = lines[2]
        altCount = lines[4]
        altCount = int(altCount)
        freq = float(freq)
        if altCount == 1:
            FreqTags["Singleton"].append(tag)
            TagFreq[tag] = "Singleton"
        if altCount != 1 and freq < 0.01:
            FreqTags["Rare"].append(tag)
            TagFreq[tag] = "Rare"
        elif freq >= 0.01 and freq < 0.05:
            FreqTags["Low"].append(tag)
            TagFreq[tag] = "Low"
        elif freq > 0.05:
            FreqTags["Common"].append(tag)
            TagFreq[tag] = "Common"
    bed_h.close()
    return FreqTags, TagFreq


def split_bed(geno_freq, bed_file, tag_out, out_prefix):
    FreqTags, TagFreq = tag_frequency(geno_freq)

    tag_h = open(tag_out, "w")
    tag_h.write("Category\tTags\n")
    for f in FreqTags:
        tags = FreqTags[f]
        tag_h.write("%s\t%s\n" % (f, ",".join(tags)))
    tag_h.close()

    FreqRecord = collections.defaultdict(list)
    bed_h = open(bed_file, "r")
    for line in bed_h:
        line = line.strip()
        lines = line.split("\t")
        Tag = lines[3]
        Freq = TagFreq[Tag]
        FreqRecord[Freq].append(line)
    bed_h.close()

    if "/" in out_prefix:
        outdir = "/".join(out_prefix.split("/")[:-1])
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    for i in FreqRecord:
        out_h = open("%s_%s.bed" % (out_prefix, i), "w")
        records = FreqRecord[i]
        out_h.write("%s\n" % "\n".join(records))
        out_h.close()



def main():
    parser = argparse.ArgumentParser(description="The the count number of sliding window across the genome.")
    parser.add_argument("-f", "--frequency", help="The input tag and frequency file.")
    parser.add_argument("-b", "--bed", help="The input bed file with tag reocrd.")
    parser.add_argument("-t", "--tag", help="The output file with frequency and tags.")
    parser.add_argument("-o", "--outPrefix", help="The output prefix.")
    args = parser.parse_args()
    split_bed(args.frequency, args.bed, args.tag, args.outPrefix)


if __name__ == "__main__":
    main()