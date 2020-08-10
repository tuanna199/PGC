#!/usr/bin/python
from tinyfasta import FastaParser
import argparse
import os 

#usage: python ~/github/NanoHub/src/NanoHub/SplitTypeFasta.py --fasta Sample_common_SV_tag_seq.fasta --outPrefix temp/temp

def split_fasta_by_type(fa_file, outPrefix):
    """
    fa_file:
    >1_83968-1_84057-89.0-INS:CN037:1_83968-1_84057-89-INS
    GAAAGAAAGAAAGAGAAAGAGAAAGAAAGAAAGATAGAAAGAAAGAAAAGAAAGAAAGAAAGAAAGAGAAAGAAAGAAAGAAAGATAG
    >1_88684-1_88831-147.0-DEL:CN051:1_88684-1_88831-147-DEL
    CTCGTGTTGGCCGGCAGAGCCGGCCCCCATCTCCTCTGACCTCCTCCCCACCTCTTGCCCTCAGCACCCAGAGTGCTCGTGACGGCCAGCAGAGCCAGCCTCCATCTCCTCTGACCTCCCACCTCTCGCCCTCAGC
    ACCCAGAGTGC
    """
    outdir = "/".join(outPrefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)


    ins_h = open(outPrefix + "_INS.fasta", "w")
    del_h = open(outPrefix + "_DEL.fasta", "w")
    
    for record in FastaParser(fa_file):
        seq = str(record.sequence)
        desc = str(record.description)
        desc = desc.split(":")[0]
        Type = desc.split("-")[-1]
        if Type == "INS":
            ins_h.write("%s\n%s\n" % (desc, seq))
        elif Type == "DEL":
            del_h.write("%s\n%s\n" % (desc, seq))
        else:
            print("Please ckeck whether INS or DEL is in description %s." % desc)
    ins_h.close()
    del_h.close()

def main():
    parser = argparse.ArgumentParser(description="Split the fasta file to DEL and INS files.")
    parser.add_argument("-f", "--fasta", help="The input fasta file.")
    parser.add_argument("-o", "--outPrefix", help="The frefix of output file.")
    args = parser.parse_args()
    split_fasta_by_type(args.fasta, args.outPrefix)



if __name__ == "__main__":
    main()

