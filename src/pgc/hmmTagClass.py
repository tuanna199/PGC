import sys
import argparse
import collections

#usage: python ~/github/NanoHub/src/NanoHub/hmmTagClass.py --dfam /home/wuzhikun/database/Dfam/Dfam_type.xls  --input SV_seq_domtblout_1000.txt --out temp.txt

def dfam_type(dfam_file):
    """
    ID      NM      Type    SubType Species
    DF0000001       MIR     SINE    MIR     Mammalia
    DF0000002       AluY    SINE    Alu     Primates
    """
    NameType = {}
    in_h = open(dfam_file, "r")
    header = in_h.readline().strip()
    for line in in_h:
        lines = line.strip().split("\t")
        NM, Type = lines[1:3]

        ### type maybe empty
        if Type == "":
            Type = "Unknown"
            
        NameType[NM] = Type
    in_h.close()

    return NameType





def parse_hmm_table(dfam_file, hmm_file, out_file):
    """
    hmm_file:
    4_24991340-4_24993658-2318-DEL:M509-2:4_24991340-4_24993658-2318-DEL            -           2318 MIR                  DF0000001.4   262   1.2e-80  273.4   0.0   1   5       2.2   2.9e+02    3.2   0.0    89   140   843   891   823   922 0.77 -
    4_24991340-4_24993658-2318-DEL:M509-2:4_24991340-4_24993658-2318-DEL            -           2318 MIR                  DF0000001.4   262   1.2e-80  273.4   0.0   2   5       3.1   4.1e+02    2.7   0.0    21    76  1512  1561  1497  1588 0.70 -
    4_24991340-4_24993658-2318-DEL:M509-2:4_24991340-4_24993658-2318-DEL            -           2318 MIR                  DF0000001.4   262   1.2e-80  273.4   0.0   3   5   3.5e-49   4.6e-47  163.3   0.0     4   262  1603  1860  1602  1863 0.97 -
    4_24991340-4_24993658-2318-DEL:M509-2:4_24991340-4_24993658-2318-DEL            -           2318 MIR                  DF0000001.4   262   1.2e-80  273.4   0.0   4   5   0.00051     0.067   15.1   0.0    65   128  1874  1937  1861  1938 0.85 -
    4_24991340-4_24993658-2318-DEL:M509-2:4_24991340-4_24993658-2318-DEL            -           2318 MIR                  DF0000001.4   262   1.2e-80  273.4   0.0   5   5     7e-24   9.2e-22   80.3   0.0    59   240  1935  2109  1931  2125 0.86 -
    15_83474409-15_83482372-7963-DEL:M537-2:15_83474409-15_83482372-7963-DEL        -           7963 MIR                  DF0000001.4   262   6.8e-70  238.2   4.6   1   3   3.6e-11   4.8e-09   38.6   0.5    60   254  3978  4190  3922  4197 0.74 -

    out_file:
    SVID    SVLength        SVType  DfamName        DfamType
    10_110678753-10_110682142-3393.5-DEL:M518-0:10_110678745-10_110682143-3398-DEL  3393.5  DEL     MIR     SINE
    10_112440059-10_112440147-87.8228-DEL:M636-0:10_112440059-10_112440152-93-DEL   87.8228 DEL     MIR     SINE
    """
    NameType = dfam_type(dfam_file)


    TargetEvalueQue = collections.defaultdict(lambda: collections.defaultdict(set))
    
    in_h = open(hmm_file , "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines = line.split()
            Target = lines[0]
            Query = lines[2]
            Evalue = lines[6]
            Evalue = float(Evalue)
            TargetEvalueQue[Target][Evalue].add(Query)
    in_h.close()



    ### get the query with smallest evalue
    TargetQuery = {}
    for t in TargetEvalueQue:
        evalues = list(TargetEvalueQue[t].keys())
        smallestEvalue = min(evalues)
        querys = TargetEvalueQue[t][smallestEvalue]

        queryLen = len(querys)
        if queryLen == 1:
            query = list(querys)[0]
        else:
            print("Please check whether the target tag %s has two querys %s." % (t, querys))
            query = sorted(list(querys))[0]

        TargetQuery[t] = query

    ### output the file:
    out_h = open(out_file, "w")
    out_h.write("SVID\tSVLength\tSVType\tDfamName\tDfamType\n")
    sortedTargets = sorted(TargetQuery.keys())
    for t in sortedTargets:
        tt = t.split(":")
        Length, SVType = tt[0].split("-")[2:4]

        query = TargetQuery[t]
        if query in NameType:
            DfamType = NameType[query]
        else:
            DfamType = ""
            print("Please check whether the name of Dfam %s has corresponding type in file %s." % (query, dfam_file))
        out_h.write("%s\t%s\t%s\t%s\t%s\n" % (t, Length, SVType, query, DfamType))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the type of repeats for Dfam database.")
    parser.add_argument("-i", "--input", help="The input file with table format of hmm result.")
    parser.add_argument("-d", "--dfam", help="The input file with dfam and type.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    parse_hmm_table(args.dfam, args.input, args.out)


if __name__ == "__main__":
    main()
