#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
import argparse
import re
import collections
from string import digits

#usage: python ~/github/NanoHub/src/NanoHub/AnnoSVClass.py --input Sample_common_SV.tsv --out temp.xls

def Type_freq_number(freq):
    if freq < 0.1: # < 0.1:
        IDPolyType = "Single"
    elif freq >= 0.1 and freq < 0.5:
        IDPolyType = "Poly"
    elif freq >= 0.5 and freq < 1:
        IDPolyType = "Major"
    elif freq == 1:
        IDPolyType = "Common"
    return IDPolyType



def supp_freq(SUPP):
    SuppList = [int(s) for s in SUPP]
    SuppSum = sum(SuppList)
    SuppLength = len(SuppList)
    SuppFreq = SuppSum / SuppLength
    return SuppFreq



def remove_location_digits(string):
    remove_digits = str.maketrans('', '', digits)
    res = string.translate(remove_digits)
    return res


def new_location_code(location):
    locationList = ['exon-exon', 'txStart-exon', 'exon-txEnd', 'txStart-intron', 'intron-txEnd', 'exon-intron', 'intron-exon', 'txStart-txEnd', 'intron-intron']
    locationListReplace = ['Exon-exon', 'Start-exon', 'Exon-end', 'Start-intron', 'Intron-end', 'Exon-intron', 'Intron-exon', 'Start-end', 'Intron-intron']
    if location in locationList:
        locIndex = locationList.index(location)
        newloc = locationListReplace[locIndex]
        return newloc
    else:
        print("Please check whether the location %s in the annotation file." % location)



def assign_SV_type(in_file, out_file):
    """
    {'exon-exon', 'txStart-exon', 'exon-txEnd', 'txStart-intron', 'intron-txEnd', 'exon-intron', 'intron-exon', 'txStart-txEnd', 'intron-intron'}
    """
    PolyList = ["Single", "Poly", "Major", "Common"]


    PolyRegionNumber = collections.defaultdict(lambda:collections.Counter())
    Locations = set()

    IDS = set()

    in_h = open(in_file, "r")
    headers = in_h.readline().strip().split("\t")
    try:
        locationIndex = headers.index("location")
    except ValueError:
        print("Please check whetehr the column name 'location' is in the header of file\n%s." % headers)
        sys.exit(1)

        
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        ID = lines[0]
        Infor = lines[11]

        if ID not in IDS:
            SUPP = Infor_target_values(Infor, "SUPP_VEC")
            SVTYPE = Infor_target_values(Infor, "SVTYPE")
            SuppFreq = supp_freq(SUPP)
            IDPolyType = Type_freq_number(SuppFreq)

            ### get the column of location
            location = lines[locationIndex]
            new_location = remove_location_digits(location)
            newloc = new_location_code(new_location)

            Locations.add(newloc)

            PolyRegionNumber[IDPolyType][newloc] += 1

        IDS.add(ID)
    in_h.close()

    out_h = open(out_file, "w")

    regionAll = []
    for i in PolyRegionNumber:
        regionAll = sorted(list(PolyRegionNumber[i]))
        break

    out_h.write("Content\t%s\n" % "\t".join(regionAll))

    for p in sorted(list(PolyRegionNumber.keys())):
        # regions = sorted(list(PolyRegionNumber[p]))
        number = [PolyRegionNumber[p][r] for r in regionAll]
        number = [str(n)  for n in number]
        out_h.write("%s\t%s\n" % (p, "\t".join(number)))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Convert the BND to SV types for NanoSV results.")
    parser.add_argument("-i", "--input", help="The input annotation file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    assign_SV_type(args.input, args.out)

if __name__ == "__main__":
    main()





