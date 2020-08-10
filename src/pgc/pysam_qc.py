#!/usr/bin/python
#-*- coding:utf-8 -*-

import pysam
import re
import os
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/pysam_qc.py --bam M671-2_MHC.bam --out temp.xls

def Bam_stat(bam):

    """
    The output is a list of (operation, length) tuples, such as [(0, 30)]
    
    M   BAM_CMATCH  0
    I   BAM_CINS    1
    D   BAM_CDEL    2
    N   BAM_CREF_SKIP   3
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD    6
    =   BAM_CEQUAL  7
    X   BAM_CDIFF   8
    B   BAM_CBACK   9
    
    """

    bamfile = pysam.AlignmentFile(bam,'rb')

    total_mis = 0
    total_ins = 0
    total_del = 0
    total_S = 0
    total_Lmax = 0
    total_Lseq = 0
    total_fenzi = 0
    total_base = 0

    for line in bamfile:
        if line.is_secondary:
            continue
        else:
            try:
                MD = line.get_tag('MD')
                CIGAR = line.cigar
            except KeyError:
                continue
                
            pat = "[0-9]+[ATGC]+"
            MD_list = re.findall(pat,MD)
            mismatch_MD = 0
            for i in MD_list:
                for j in i:
                    if j == 'A' or j == 'T' or j == 'G' or j == 'C':
                        mismatch_MD += 1

            # match/insertion/deletion/soft clip
            
            match = 0
            ins = 0
            delete = 0
            sclip = 0
            hclip = 0
            other = 0

            for k in CIGAR:
                if k[0] == 0:
                    match += k[1]
                if k[0] == 1:
                    ins += k[1]
                if k[0] == 2:
                    delete += k[1]
                if k[0] == 4:
                    sclip += k[1]
                if k[0] == 5:
                    hclip += k[1]
                else:
                    other += k[1]

            match -= mismatch_MD
            
            sum_base = match + ins + delete + sclip + hclip + other
            fenzi = match + delete + other

            M = match
            S = M + ins + sclip
            L_seq = M + ins
            L_ref = M + delete
            L_max = max(L_seq, L_ref)
            
            total_mis += mismatch_MD
            total_ins += ins
            total_del += delete
            total_S += S
            total_Lmax += L_max
            total_Lseq += L_seq
            total_fenzi += fenzi
            total_base += sum_base
     
    # error rate
    E_ins = int(total_ins)/int(total_Lmax) * 100
    E_del = int(total_del)/int(total_Lmax) * 100
    E_mis = int(total_mis)/int(total_Lmax) * 100

    E_total = E_ins + E_del + E_mis
    E_ins = "%.2f" % E_ins
    E_del = "%.2f" % E_del
    E_mis = "%.2f" % E_mis
    E_total = "%.2f" % E_total
    A = '%.2f' % (int(total_Lseq)/int(total_S) * 100)
        
    BP_rate = int(total_fenzi)/int(total_base) * 100
    BP_align = "%.2f" % BP_rate

    return E_ins,E_del,E_mis,E_total,A,BP_align

def getRes(input_file, output):
    sample = input_file.split("/")[-1].split("_")[0].split(".")[0]
    f = open(output, "w")
    E_ins,E_del,E_mis,E_total,A,BP_align = Bam_stat(input_file)
    head = ["Sample", "E_ins(%)","E_del(%)","E_mis(%)","E_total(%)","Accuracy(%)","BP_align(%)"]
    f.write("%s\n" % "\t".join(head))
    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample,E_ins,E_del,E_mis,E_total,A,BP_align))

def main():
    parser = argparse.ArgumentParser(description="Discard de novo SV for common SV of population.")
    parser.add_argument("-b", "--bam", help="The input bam file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    getRes(args.bam, args.out)


if __name__ == "__main__":
    main()
