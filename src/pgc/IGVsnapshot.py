#!/usr/bin/env python
import argparse

#usage : python ~/github/TrioWGS/src/TrioWGS/IGVsnapshot.py --bed M625_denovo_SV_proband_test-1.bed --out M625_denovo_SV_proband_test-1_igv.batch --outdir /home/wuzhikun/Project/NanoTrio/IGV --genome /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --bam /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M625-0.bam --sample M625



def IGV_snapshot_bat(bed_file, out_file, genome, bamfile, outdir, sample, height=500):
    height = int(height)

    out_h = open(out_file, "w")
    out_h.write("new\n")
    out_h.write("genome %s\n" % genome)
    out_h.write("snapshotDirectory %s\n" % outdir)
    ### input one or multiple bam files
    bams = bamfile.split(",")
    bams = [b.strip() for b in bams]
    for bam in bams:
        out_h.write("load %s\n" % bam)

    out_h.write("maxPanelHeight %d\n" % height)

    ### read the bed files
    in_h = open(bed_file, "r")
    for line in in_h:
        line = line.strip()
        if line != "":
            lines = line.split("\t")
            Chr, Start, End = lines[:3]
            if len(lines) == 3:
                name = "_".join([sample, Chr, Start, End])
            elif len(lines) >= 4:
                name = "_".join([sample, lines[3]])

            out_h.write("goto %s:%s-%s\n" % (Chr, Start, End))
            out_h.write("snapshot %s\n" % name)

    out_h.write("exit\n")
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Get the script file of IGV based on genome, bam and bed file. ")
    parser.add_argument("-b", "--bed", help="The bed file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-d", "--outdir", help="The out directory used for stroring pictures in batch script.")
    parser.add_argument("-g", "--genome", help="The reference genome for bam file.")
    parser.add_argument("-m", "--bam", help="The input bam files which are separated with ','.")
    parser.add_argument("-s", "--sample", help="The sample name.")
    parser.add_argument("-l", "--height", default=500, help="The height of the picture.")
    args = parser.parse_args()
    IGV_snapshot_bat(args.bed, args.out, args.genome, args.bam, args.outdir, args.sample, height=args.height)



if __name__ == "__main__":
    main()