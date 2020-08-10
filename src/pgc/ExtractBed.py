#!/usr/bin/python
from BaseFunc import Infor_target_values
import argparse
import sys
import os
import collections



def extract_information(infor_bed, out_prefix):
    """
    infor_bed:
    chr1    786548  786567  chr1    837623  837642  454810  614.95  -       +       DUP     PASS    454810  T       <DUP>
       .       .       .       SVTYPE=DUP;POS=786558;SVLEN=51075;END=837633;STRANDS=-+:31;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0
    ,0;CIEND95=0,0;SU=31;PE=22;SR=9;ALG=PROD;NSAMP=3;MSQ=118.98;AC=3;AN=29246;NFAM=3;AN_NFE=8402;AN_AFR=9872;AN_Other=170;A
    N_AMR=4656;AN_PI=174;AN_FE=2380;AN_EAS=26;AN_SAS=1002;AC_NFE=2;AC_AFR=1;AC_Other=0;AC_AMR=0;AC_PI=0;AC_FE=0;AC_EAS=0;AC
    _SAS=0 
    """


    TypeRecords = collections.defaultdict(list)
    bed_h = open(infor_bed, "r")
    for line in bed_h:
        lines = line.strip().split("\t")
        Chr = lines[0]
        Chr = Chr.lstrip("chr")
        Infor = lines[-2]
        SVType, Start, SVLen, End = Infor_target_values(Infor, "SVTYPE,POS,SVLEN,END")
        SVLen = SVLen.lstrip("-")
        Tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr, End, SVLen, SVType)
        Record = "%s\t%s\t%s\t%s" % (Chr, Start, End, Tag)
        TypeRecords[SVType].append(Record)
    bed_h.close()


    ### output the file
    outdir = "/".join(out_prefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for t in TypeRecords:
        if t == "MEI":
            t = "INS"
        out_h = open("%s_%s.bed" % (out_prefix, t), "w")
        records = TypeRecords[t]
        out_h.write("%s\n" % "\n".join(records))
        out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Extract the SV information from the bed file.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-o", "--outPrefix", help="The output prefix.")
    args = parser.parse_args()
    extract_information(args.bed, args.outPrefix)


if __name__ == "__main__":
    main()


