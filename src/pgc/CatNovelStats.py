#!/usr/bin/python
from __future__ import division
import collections
import argparse
import sys
import os

# usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/CatNovelStats.py --category /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls --overlap /home/wuzhikun/Project/Population/population/bed/SV_overlap_filt_tags_overlap.xls --stats /home/wuzhikun/Project/Population/population/Category/Kown_and_novel_tags_stats.txt --novel /home/wuzhikun/Project/Population/population/Category/Novel_SV_category_tag.xls

def known_tags(overlap_file):
    """
    overlap_file:
    DGV     LRS15   WGS17795        dbVar   gnomAD
    NA      10_100093578-10_100093707-143-INS       NA      NA      NA
    10_100231216-10_100245042-13826-DUP     NA      NA      10_100231216-10_100245042-13826-DUP     NA
    """
    KnownTags = set()
    overlap_h = open(overlap_file, "r")
    header = overlap_h.readline().strip()
    for line in overlap_h:
        line = line.strip()
        lines = line.split("\t")
        tagSet = set(lines) - {"NA"}
        setLen = len(tagSet)
        if setLen == 1:
            tag = list(tagSet)[0]
            KnownTags.add(tag)
        else:
            print("Please check whether the there are multiple tags in the record %s in file %s." % (line, overlap_file))
            sys.exit(1)
    overlap_h.close()
    return KnownTags

def known_and_novel_tags(overlap_file, category_file, stats_file, novel_file):
    """
    stats_file:
    Category        KnownTags       NovelTags       KnownRatio      NovelRatio
    Common  27990   4792    0.854   0.146
    Low     7436    6926    0.518   0.482
    Rare    10517   18434   0.363   0.637
    Singleton       13262   42999   0.236   0.764
    """
    KnownTags = known_tags(overlap_file)

    cat_h = open(category_file, "r")
    header = cat_h.readline().strip()


    KnownCatTags = collections.defaultdict(list)
    NovelCatTags = collections.defaultdict(list)
    AllCatTags = {}

    CATS = set()
    for line in cat_h:
        lines = line.strip().split("\t")
        cat, Tag = lines
        CATS.add(cat)
        tags = Tag.split(",")
        AllCatTags[cat] = len(tags)
        for t in tags:
            if t in KnownTags:
                KnownCatTags[cat].append(t)
            else:
                NovelCatTags[cat].append(t)
    cat_h.close()

    ### output the statistics for known and novel tags
    stats_h = open(stats_file, "w")
    stats_h.write("Category\tKnownTags\tNovelTags\tKnownRatio\tNovelRatio\n")
    cats = sorted(list(CATS))
    for c in cats:
        All = AllCatTags[c]

        if c in KnownCatTags:
            known = len(KnownCatTags[c])
        else:
            known = 0

        if c in NovelCatTags:
            novel = len(NovelCatTags[c])
        else:
            novel = 0

        knownRatio = known / All
        knownRatio = "%.3f" % knownRatio

        novelRatio = novel / All
        novelRatio = "%.3f" % novelRatio
        stats_h.write("%s\t%d\t%d\t%s\t%s\n" % (c, known, novel, knownRatio, novelRatio))
    stats_h.close()

    ### output the novel cattegory tags
    novel_h = open(novel_file, "w")
    novel_h.write("%s\n" % header)
    for c in NovelCatTags:
        tags = NovelCatTags[c]
        novel_h.write("%s\t%s\n" % (c, ",".join(tags)))
    novel_h.close()

def main():
    parser = argparse.ArgumentParser(description="Statistics for known and novel tags of categiries and output novel tags.")
    parser.add_argument("-c", "--category", help="The file with category and tags.")
    parser.add_argument("-e", "--overlap", help="The file with overlapped known tags.")
    parser.add_argument("-n", "--novel", help="Out the novel tags with different categories.")
    parser.add_argument("-s", "--stats", help="The statistics for known and novel tag number for different categories.")
    args = parser.parse_args()
    known_and_novel_tags(args.overlap, args.category, args.stats, args.novel)

if __name__ == "__main__":
    main()




