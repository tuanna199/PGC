#!/usr/bin/python
from __future__ import division
import argparse
import collections
import sys
import math
import os

#usage: python ~/github/NanoHub/src/NanoHub/IGVBatchSplit.py --input M671_igv_Sniffles.batch --outdir temp --number 100

def batch_file_split(batch_file, shotNumber, out_dir):
    """
    batch_file:
    new
    genome /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    snapshotDirectory /home/wuzhikun/Project/NanoTrio/TrioCommon/IGV/M671_igv_Sniffles
    load /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-0.bam
    load /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-1.bam
    load /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-2.bam
    maxPanelHeight 700
    goto 1:820830-821171
    snapshot M671_1_820880_821121_INS
    """
    shotNumber = int(shotNumber)

    fileBase = batch_file.split("/")[-1]

    if not os.path.exists(out_dir):
        # os.mkdir(out_dir)
        cmd = "mkdir -p %s" % out_dir
        os.system(cmd)
    if not out_dir.endswith("/"):
        out_dir += "/"


    batch_h = open(batch_file, "r")
    Common = []
    Record = []

    for line in batch_h:
        line = line.strip()
        if line.startswith("snapshotDirectory"):
            Common.append(line)
            shotDir = line
        elif line.startswith("goto") or line.startswith("snapshot"):
            Record.append(line)
        else:
            Common.append(line)
    batch_h.close()

    RecordNum = len(Record)

    if RecordNum % 2 != 0:
        print("Please check whether the record number of 'goto' and 'snapshot' is even.")
        sys.exit(1)





    eventNum = RecordNum / 2
    if eventNum <= shotNumber:
        outfile = out_dir + fileBase
        out_h = open(outfile, "w")
        out_h.write("%s\n%s\n%s\n" % ("\n".join(Common[:-1]), "\n".join(Record), Common[-1]))
        out_h.close()
    else:
        fileNum = math.ceil(eventNum / shotNumber)
        for i in range(fileNum):
            ### output file
            outfile = out_dir + fileBase + "-" + str(i)
            out_h = open(outfile, "w")

            if i != fileNum:
                events = Record[i* shotNumber* 2 : (i+1) * (shotNumber-1) * 2]
            else:
                events = Record[i* shotNumber* 2 : -1]

            ### change the output file
            newShotDir = shotDir + "-" + str(i)
            NewCommon = []
            for c in Common:
                if c.startswith("snapshotDirectory"):
                    newC = newShotDir
                else:
                    newC = c
                NewCommon.append(newC)

            out_h.write("%s\n%s\n%s\n\n" % ("\n".join(NewCommon[:-1]), "\n".join(events), NewCommon[-1]))
            out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Divide the large batch file to small ones.")
    parser.add_argument("-i", "--input", help="The input batch file.")
    parser.add_argument("-n", "--number", default=100, help="The shot number for each small file.")
    parser.add_argument("-o", "--outdir", help="The output file.")
    args = parser.parse_args()
    batch_file_split(args.input, args.number, args.outdir)


if __name__ == "__main__":
    main()
