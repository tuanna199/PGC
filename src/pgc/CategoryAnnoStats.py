#!/usr/bin/python
from __future__ import division
import collections
import sys
import os
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/CategoryAnnoStats.py --category /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls --annotation /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap.bed --out temp

def category_tags(category_file):
    """
    category_file:
    Category        Tags
    common  10_101728712-10_101729253-540.3-INS,10_116505797-10_116506120-321.9-DEL,10_124330938-10_124331032-91.7-INS,10_126417903-10_126418217-314.2-
    major   10_100093578-10_100093710-119.6-
    """
    TagCategory = {}
    CategoryCount = collections.Counter()

    in_h = open(category_file, "r")
    header = in_h.readline()
    for line in in_h:
        lines = line.strip().split("\t")
        cate, tag = lines[:2]
        tags = tag.split(",")
        for tag in tags:
            TagCategory[tag] = cate
            CategoryCount[cate] += 1
            CategoryCount["All"] += 1
    in_h.close()
    return TagCategory, CategoryCount
    

def category_annotation_statistics(category_file, anno_file, out_file, stat_file):
    """
    anno_file:
    1       67910   68341   1_67910-1_68341-431.0-DEL       1       68090   69090   OR4F5   ENST00000335137.3       UpStream
    1       88684   88831   1_88684-1_88831-147.0-DEL       1       88294   89294   RP11-34P13.7    ENST00000466430.5       DownStream
    1       88684   88831   1_88684-1_88831-147.0-DEL       1       88550   89550   RP11-34P13.8    ENST00000495576.1       DownStream


    out_file:
    Category        CDS     DownStream      Exon    Intron  UTR3    UTR5    UpStream        Intergenic
    Common  17(2.3) 15(2.1) 6(0.8)  353(48.3)       4(0.5)  7(1.0)  24(3.3) 352(48.2)
    Major   234(1.9)        354(2.9)        170(1.4)        5891(48.0)      81(0.7) 92(0.7) 355(2.9)        6037(49.2)
    Poly    1046(2.1)       1681(3.3)       893(1.8)        24086(47.3)     543(1.1)        543(1.1)        1701(3.3)       25069(49.2)
    Single  1436(2.7)       1996(3.7)       1201(2.2)       26615(49.3)     761(1.4)        777(1.4)        1943(3.6)       25466(47.2)
    """
    CatTagFeatures = collections.defaultdict(lambda: collections.defaultdict(set))
    Features = set()

    TagCategory, CategoryCount = category_tags(category_file)


    in_h = open(anno_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        gene, trans, feature = lines[7:10]

        if "utr" in feature.lower():
            feature = "UTR5/3"

        if "stream" in feature.lower():
            feature = "Up/DownStream"

        if feature == "Exon":
            feature = "NC_Exon"

        ### get all features
        Features.add(feature)

        try:
            cate = TagCategory[tag]
        except KeyError:
            print("Please check whether the tag %s is in file %s." % (tag, category_file))
            sys.exit(1)
        CatTagFeatures[cate][tag].add(feature)
        CatTagFeatures["All"][tag].add(feature)
    in_h.close()


    FeatureList = sorted(list(Features))
    ### output the statistics
    out_h = open(out_file, "w")
    CateAllCounts = {}
    out_h.write("Category\t%s\tOther\tAll\n" % ("\t".join(FeatureList)))

    stat_h = open(stat_file, "w")
    stat_h.write("Category\t%s\tOther\tAll\n" % ("\t".join(FeatureList)))


    Cates = sorted(list(CatTagFeatures.keys()))
    for cate in Cates:
        TagFeatures = CatTagFeatures[cate]

        FeatureCounts = collections.Counter()
        GeneCounts = collections.Counter()

        for tag in TagFeatures:
            GeneCounts[cate] += 1
            features = list(TagFeatures[tag])
            for f in features:
                FeatureCounts[f] += 1


        ### output 
        allTags = CategoryCount[cate]

        geneTags = GeneCounts[cate]
        interTags = allTags - geneTags
        interRatio = interTags / allTags 
        interRatio = "%.3f" % interRatio

        FeatureRecord = []
        FeatureStat = []
        for feature in FeatureList:
            if feature in FeatureCounts:
                count = FeatureCounts[feature]
            else:
                count = 0

            geneRatio = count / allTags
            geneRatio = "%.3f" % geneRatio
            record = "%s(%s)" % (str(count), geneRatio)
            FeatureRecord.append(record)
            FeatureStat.append(str(count))

        out_h.write("%s\t%s\t%d(%s)\t%d\n" % (cate.capitalize(), "\t".join(FeatureRecord), interTags, interRatio, allTags))
        stat_h.write("%s\t%s\t%d\t%d\n" % (cate.capitalize(), "\t".join(FeatureStat), interTags, allTags))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Statistics for category and features of annotation file.")
    parser.add_argument("-a", "--annotation", help="The input annotation file with bed format.")
    parser.add_argument("-c", "--category", help="The input category file of tags.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-s", "--stats", help="The output file without the ratio statistics.")
    args = parser.parse_args()
    category_annotation_statistics(args.category, args.annotation, args.out, args.stats)




if __name__ == "__main__":
    main()





