#!/usr/bin/python
from __future__ import division
import collections
import argparse
import sys

#usage: python ~/github/NanoHub/src/NanoHub/SVGenoFreqFilt.py --genotype Sample_SV_genotype.txt --freqThreshold 0.05 --out temp.txt --frequency temp_freq.txt

def genotype_frequency(genotypes, line):
    """
    genotypes:
    ["0/0", "0/1", "1/1"]
    """
    genoLen = len(genotypes)
    refs = 0
    heters = 0
    alts = 0
    altCount = 0
    for ge in genotypes:
        if ge == "0/0":
            refs += 2
        elif ge == "0/1" or ge == "1/0":
            alts += 1
            refs += 1
            altCount += 1
        elif ge == "1/1":
            alts += 2
            altCount += 1
        else:
            print("Please check and make sure that the genotype for SV is '0/0', '0/1', '1/0' or '1/1' for record %s." % line)
            sys.exit(1)
    
    values = refs + alts
    if values != genoLen * 2 :
        print("Please make sure that all genotypes are '0/0', '0/1', '1/0' or '1/1' for record %s." % line)
        sys.exit(1)

    refFeq = refs  / values
    altFeq = alts / values

    return refFeq, altFeq, altCount

     



def SV_genotype_frequency_filt(sv_genotype, out_file, freq_file, freq_threshold):
    """
    sv_genotype:
    Chr1    Pos1    Chr2    Pos2    Type    M509-0  M561-2  M655-1  M668-1  M509-1  M509-2  M561-0  M561-1  M655-0  M655-2  M668-0  M668-2  M425-0  M425-1  M425-2
    1       10022   X       449438  TRA     0/0     0/0     0/0     0/0     0/0     0/0     0/0     1/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       10001   20      64286701        TRA     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0

    freq_file:
    Chr1	Pos1	Chr2	Pos2	Type	Ref_freq	Heter_freq	Alt_freq
    1	10022	X	449438	TRA	0.951	0.000	0.049
    1	10001	20	64286701	TRA	0.992	0.008	0.000
    1	10056	20	64287143	TRA	0.984	0.016	0.000
    """
    freq_threshold = float(freq_threshold)

    sv_h = open(sv_genotype, "r")
    header = sv_h.readline().strip()

    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)

    freq_h = open(freq_file, "w")
    titles = header.split("\t")[:6]
    # freq_h.write("%s\tRef_freq\tAlt_freq\tMAF\tTag\n" % "\t".join(titles))

    freq_h.write("Tag\tRef_freq\tAlt_freq\tMAF\n")

    for line in sv_h:
        line = line.strip()
        lines = line.split("\t")
        ### genotype starts from 6th column
        genotypes = lines[6:]

        ### get the frequency for genotype
        refFeq, altFeq, altCount = genotype_frequency(genotypes, line)


        minFreq = min([refFeq, altFeq])
        if minFreq >= freq_threshold:
            out_h.write("%s\n" % line)
        

        ### output the frequency
        titles = lines[:6]
        Tag = "%s_%s-%s_%s-%s-%s" % tuple(titles)
        # freq_h.write("%s\t%.3f\t%.3f\t%.3f\t%s\n" % ("\t".join(titles), refFeq, altFeq, minFreq, Tag))
        freq_h.write("%s\t%.4f\t%.4f\t%.4f\t%d\n" % (Tag, refFeq, altFeq, minFreq, altCount))


    sv_h.close()
    out_h.close()
    freq_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get and filt the genotype based on frequency.")
    parser.add_argument("-g", "--genotype", help="The input SV genotype file.")
    parser.add_argument("-t", "--freqThreshold", default=0.05, help="The frequency threshold.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-f", "--frequency", help="The output frequency file.")
    args = parser.parse_args()
    SV_genotype_frequency_filt(args.genotype, args.out, args.frequency, args.freqThreshold)

if __name__ == "__main__":
    main()

