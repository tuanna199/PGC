#!/usr/bin/python
from __future__ import division
import collections
import argparse
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/hierarchyAnnotation.py --annotation population/Annotation/Sample_common_SV_genepred_overlap.bed --category population/Category/Sample_SV_category_tag.xls --out gene_annotation.txt --stats annotation_stat.xls --ratio annotation_stat_ratio.xls


def hierarchy_feature(features):
    features = [f.lower() for f in features]
    if "cds" in features:
        f = "CDS"
    elif "utr5" in features or "utr3" in features or "utr" in features:
        f = "UTR"
    elif "promoter" in features:
        f = "Promoter"
    elif "intron" in features:
        f = "Intron"
    else:
        print(features)
    return f



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
    


def hierarchy_annotation(anno_file, category_file, out_file, stat_file, ratio_file):
    """
    anno_file:
    1       939439  939492  1_939439-1_939492-52-DEL        1       939274  939460  SAMD11  ENST00000342066.7       CDS
    1       939439  939492  1_939439-1_939492-52-DEL        1       939412  941143  SAMD11  ENST00000617307.4       Intron
    """
    featureList = ["CDS", "Intron", "UTR", "Promoter"]
    featureList = [f.lower() for f in featureList]

    TagGeneFeature = collections.defaultdict(lambda: collections.defaultdict(set))
    anno_h = open(anno_file, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        gene = lines[7]
        feature = lines[-1]
        TagGeneFeature[tag][gene].add(feature)
    anno_h.close()

    out_h = open(out_file, "w")
    TagFeatures = collections.defaultdict(set)
    Tags = sorted(list(TagGeneFeature.keys()))
    for t in Tags:
        genes = sorted(list(TagGeneFeature[t].keys()))
        for g in genes:
            features = TagGeneFeature[t][g]
            f = hierarchy_feature(features)
            out_h.write("%s\t%s\t%s\n" % (t, g, f))
            TagFeatures[t].add(f)
    out_h.close()


    TagCategory, CategoryCount = category_tags(category_file)

    ### statistics for annotation
    CategoryFeatureCount = collections.defaultdict(lambda: collections.Counter())
    AllFeatures = set()
    for t in TagFeatures:
        features = TagFeatures[t]
        f = hierarchy_feature(features)
        AllFeatures.add(f)
        try:
            cat = TagCategory[t]
        except KeyError:
            print("Please check whether %s was in %s." % (t, category_file))
            sys.exit(1)

        CategoryFeatureCount[cat][f] += 1
        CategoryFeatureCount["All"][f] += 1


    sortFeatures = sorted(list(AllFeatures))
    stat_h = open(stat_file, "w")
    stat_h.write("Category\t%s\tIntergenic\tAll\n" % "\t".join(sortFeatures))

    ratio_h = open(ratio_file, "w")
    ratio_h.write("Category\t%s\tIntergenic\tAll\n" % "\t".join(sortFeatures))
    CategoryList = ["All", "Common", "Low", "Rare", "Singleton"]
    for c in CategoryList:
        if c in CategoryFeatureCount:
            features = CategoryFeatureCount[c]
            feaCount = 0
            Counts = []
            for f  in sortFeatures:
                if f in features:
                    count = CategoryFeatureCount[c][f]
                    count = int(count)
                else:
                    count = 0
                feaCount += count
                Counts.append(count)

            AllCount = CategoryCount[c]
            countRatio = ["%.1f" % (c / AllCount * 100) for c in Counts]
            
            ### intergenic region
            intergenic = AllCount - feaCount
            intergenicRatio = "%.1f" % (intergenic / AllCount * 100)

            Counts = [str(c) for c in Counts]

            stats = []
            for n, r in zip(Counts, countRatio):
                stat = "%s (%s)" % (n, r)
                stats.append(stat)

            stat_h.write("%s\t%s\t%d\t%d\n" % (c, "\t".join(Counts), intergenic, AllCount))
            ratio_h.write("%s\t%s\t%d (%s)\t%d\n" % (c, "\t".join(stats), intergenic, intergenicRatio, AllCount))
        else:
            print("Please check whether %s is in file %s." % (c, CategoryFeatureCount))
            sys.exit(1)
    stat_h.close()
    ratio_h.close()




def main():
    parser = argparse.ArgumentParser(description="Hierarchy annotation.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-c", "--category", help="The category file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-s", "--stats", help="Statistics for gene annotation.")
    parser.add_argument("-r", "--ratio", help="The ratio of statistics for annotation")
    args = parser.parse_args()
    hierarchy_annotation(args.annotation, args.category, args.out, args.stats, args.ratio)



if __name__ == "__main__":
    main()


