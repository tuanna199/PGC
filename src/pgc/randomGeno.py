#1/usr/bin/python
import collections
import argparse
import sys
import os
import random
import numpy


#usage: python ~/github/NanoHub/src/NanoHub/randomGeno.py --genotype Sample_SV_genotype.txt --category /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls  --out tempxls --number 6 --repeat 5


def is_geno_diverse(genos):
    genoSet = set(genos)
    setLen = len(genoSet)
    if setLen == 1 and list(genoSet)[0] == "0/0":
        isDiverse = False
    else:
        isDiverse = True
    return isDiverse



def tag_category(cate_file):
    TagCat = {}
    cat_h = open(cate_file, "r")
    header = cat_h.readline().strip()
    for line in cat_h:
        lines = line.strip().split("\t")
        Cat, Tag = lines
        tags = Tag.split(",")
        for tag in tags:
            TagCat[tag] = Cat
    cat_h.close()
    return TagCat



def random_select_genos(geno_file, category_file, number):
    """
    out_file:
    Number  Singleton   Rare    Low Common
    6   735 1236    3119    28457
    """
    number = int(number)

    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    infor = headers[:6]
    samples = headers[6:]
    sampleLen = len(samples)

    if number > sampleLen:
        print("Please make sure that the selected number %s is less or equal to the sample length of file %s." % (number, geno_file))
        sys.exit(1)

    sampleIndex = list(range(sampleLen))
    shuffleIndex = sampleIndex
    random.shuffle(shuffleIndex)

    targetIndex = shuffleIndex[:number]
    targetSample = [samples[i] for i in targetIndex]



    TargetTags = {}
    for line in geno_h:
        lines = line.strip().split("\t")
        infor = lines[:6]
        genos = lines[6:]
        targetGeno = [genos[i] for i in targetIndex]
        isDiverse = is_geno_diverse(targetGeno)
        if isDiverse == True:
            tag = "%s_%s-%s_%s-%s-%s" % tuple(infor)
            TargetTags[tag] = 1
    geno_h.close()


    
    ### category for target tags.
    TagCat = tag_category(category_file)

    CategoryCount = collections.Counter()
    for t in TargetTags:
        try:
            cat = TagCat[t]
            CategoryCount[cat] += 1
        except KeyError:
            print("Please check whether the tag %s is in the tag category file %s." % (t, cat))
            sys.exit(1)
    return CategoryCount




def select_repeat_stats(geno_file, category_file, number, repeat, out_file):
    CateList = ["Singleton", "Rare", "Low", "Common"]
    repeat = int(repeat)

    out_h = open(out_file, "w")
    out_h.write("Number\t%s\tAll\n" % "\t".join(CateList))

    AllCounts = []
    for i in range(repeat):
        CategoryCount = random_select_genos(geno_file, category_file, number)
        Counts = []
        for c in CateList:
            if c in CategoryCount:
                count = CategoryCount[c]
            else:
                count = 0
            Counts.append(count)
        countSum = sum(Counts)
        Counts.append(countSum)
        AllCounts.append(Counts)

    ### statistics for all

    ### transpose the two dim array
    TransConts = [*zip(*AllCounts)]
    CatMean = [int(numpy.mean(s)) for s in TransConts]


    CatMean = [str(c) for c in CatMean]
    out_h.write("%s\t%s\n" % (number, "\t".join(CatMean)))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Shuffle and select target number genotypes.")
    parser.add_argument("-b", "--genotype", help="The input genotype file..")
    parser.add_argument("-c", "--category", help="The category file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-n", "--number", help="The selected number.")
    parser.add_argument("-r", "--repeat", default=10, help="The repeat number.")
    args = parser.parse_args()
    select_repeat_stats(args.genotype, args.category, args.number, args.repeat, args.out)

if __name__ == "__main__":
    main()


