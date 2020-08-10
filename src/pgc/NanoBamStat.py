#!/usr/bin/env python
from __future__ import division
import sys
import argparse


def bam_stat_summary(bam_stat, out_file, sample):
    """
    # Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
    SN      raw total sequences:    31865
    SN      filtered sequences:     0
    SN      sequences:      31865
    SN      is sorted:      1
    SN      1st fragments:  31865
    SN      last fragments: 0
    SN      reads mapped:   20412
    SN      reads mapped and paired:        0       # paired-end technology bit set + both mates mapped
    SN      reads unmapped: 11453
    SN      reads properly paired:  0       # proper-pair bit set
    SN      reads paired:   0       # paired-end technology bit set
    SN      reads duplicated:       0       # PCR or optical duplicate bit set
    SN      reads MQ0:      335     # mapped and MQ=0
    SN      reads QC failed:        0
    SN      non-primary alignments: 65295
    SN      total length:   413416059       # ignores clipping
    SN      total first fragment length:    413416059       # ignores clipping
    SN      total last fragment length:     0       # ignores clipping
    SN      bases mapped:   306572781       # ignores clipping
    SN      bases mapped (cigar):   88799666        # more accurate
    SN      bases trimmed:  0
    SN      bases duplicated:       0
    SN      mismatches:     19840449        # from NM fields
    SN      error rate:     2.234293e-01    # mismatches / bases mapped (cigar)
    SN      average length: 12973
    SN      average first fragment length:  12974
    SN      average last fragment length:   0
    SN      maximum length: 100111
    SN      maximum first fragment length:  0
    SN      maximum last fragment length:   0
    SN      average quality:        12.0
    SN      insert size average:    0.0
    SN      insert size standard deviation: 0.0
    SN      inward oriented pairs:  0
    SN      outward oriented pairs: 0
    SN      pairs with other orientation:   0
    SN      pairs on different chromosomes: 0
    SN      percentage of properly paired reads (%):        0.0

    """
    in_h = open(bam_stat, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("SN"):
            lines = line.split("\t")
            tag = lines[1]
            if tag.startswith("raw total sequences"):
                totalRead = lines[2]
            elif tag.startswith("reads mapped:"):
                mappedRead = lines[2]
            elif tag.startswith("total length"):
                totalBase = lines[2]
            elif tag.startswith("bases mapped") and "ignores clipping" in line:
                mappedBase = lines[2]
            elif tag.startswith("bases mapped (cigar)"):
                mappedBaseCigar = lines[2]
            elif tag.startswith("error rate"):
                errorRate = lines[2]
                errorRate = float(errorRate)
    in_h.close()

    try:
        mappedRate = "%.2f"  %  (int(mappedRead) / int(totalRead) * 100)
        mappedBaseRate = "%.2f" % (int(mappedBase) / int(totalBase) * 100)
    except ZeroDivisionError:
        print("Please check 'total sequences' exists in the bam stats file %s." % bam_stat)
        sys.exit(1)
    errorRate = "%.2f" %  (errorRate * 100)

    out_h = open(out_file, "w")
    out_h.write("Sample\tTotal_reads\tMapped_reads\tTotal_bases\tMapped_bases\tMapped_base_rate(%)\tError_rate(%)\tMapped_read_rate(%)\n")
    out_h.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, totalRead, mappedRead, totalBase, mappedBase, mappedBaseRate, errorRate, mappedRate))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the summary of bam mapping statistics.")
    parser.add_argument("-t", "--stat", help="The input stat file derived from 'samtools stats'.")
    parser.add_argument("-o", "--out", help="The ouput file.")
    parser.add_argument("-s", "--sample", default="sample", help="The sample name.")

    args = parser.parse_args()
    bam_stat_summary(args.stat, args.out, args.sample)

if __name__== "__main__":
    main()


