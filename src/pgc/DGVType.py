#!/usr/bin/python
import sys
import argparse
import collections
import os

#usage: python ~/github/NanoHub/src/NanoHub/DGVType.py --input GRCh38_hg38_variants_2020-02-25.txt --outPrefix temp/temp

def which_type(Subtype, SUBTYPES):
    targetTypes = []
    for s in SUBTYPES:
        if s in Subtype:
            targetTypes.append(s)
    return targetTypes

def convert_SV_type(DGV_file, outPrefix):
    """
    DGV_file:
    variantaccession        chr     start   end     varianttype     variantsubtype  reference       pubmedid        method  platform        mergedvaria
    dgv1n82 1       10001   22118   CNV     duplication     Sudmant_et_al_2013      23825009        Oligo aCGH,Sequencing                   nsv945697,n
    nsv7879 1       10001   127330  CNV     gain+loss       Perry_et_al_2008        18304495        Oligo aCGH  

    out_file1:
    1   100002154   100003577   esv2567165  deletion
    1   10000942    10005942    nsv4044290  deletion
    1   10001   127330  nsv7879 gain+loss
    1   10001   2368561 nsv482937   loss
    """
    SUBTYPESDICT = {"deletion": "DEL", "insertion": "INS", "duplication": "DUP", "inversion": "INV", "gain": "DUP", "loss": "DEL"}
    SUBTYPES = ["deletion", "insertion", "duplication", "inversion", "gain", "loss"] 

    TypeRecord = collections.defaultdict(list)

    in_h = open(DGV_file, "r")
    header = in_h.readline().strip()

    for line in in_h:
        lines = line.strip().split("\t")
        variantID, Chr, Start, End, Type, Subtype = lines[:6]
        # record = "%s\t%s\t%s\t%s\t%s" % (Chr, Start, End, variantID, Subtype)
        Start = int(Start)
        End = int(End)
        record = [Chr, Start, End, variantID, Subtype]
        ### finall the type key
        targetTypes = which_type(Subtype, SUBTYPES)
        targetLength = len(targetTypes)

        ### it may be not target type or chromosome empty
        if targetLength == 0 or Chr == "" or Start == End:
            continue
        else:
            for s in targetTypes:
                t = SUBTYPESDICT[s]
                TypeRecord[t].append(record)
    in_h.close()


    ### output the file
    outdir = "/".join(outPrefix.split("/")[:-1])
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    targetTypes = sorted(list(TypeRecord.keys()))
    for t in TypeRecord:
        records = sorted(TypeRecord[t])
        out_h = open("%s_%s.bed" % (outPrefix, t), "w")
        for r in records:
            r = [str(n) for n in r]
            record = "\t".join(r)
            out_h.write("%s\n" % record)
        out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Convert the SV types for DGV database.")
    parser.add_argument("-i", "--input", help="The input file.")
    parser.add_argument("-o", "--outPrefix", help="The frefix of output file.")
    args = parser.parse_args()
    convert_SV_type(args.input, args.outPrefix)



if __name__ == "__main__":
    main()
