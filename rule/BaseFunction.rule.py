#!/usr/bin/env python
import argparse
import collections
import sys


#########################################
def catagory_SV_distance(cate_file, dist_file, out_file):
    TagCat = {}
    cat_h = open(cate_file, "r")
    header = cat_h.readline()
    for line in cat_h:
        lines = line.strip().split("\t")
        cat, Tag = lines[:2]
        tags = Tag.split(",")
        for t in tags:
            TagCat[t] = cat
    cat_h.close()
    dist_h = open(dist_file, "r")
    header = dist_h.readline()
    out_h = open(out_file, "w")
    out_h.write("Category\tDistance\n")
    for line in dist_h:
        lines = line.strip().split("\t")
        tag, length = lines[:2]
        cat = TagCat[tag]
        out_h.write("%s\t%s\n" % (cat, length))
    dist_h.close()
    out_h.close()




########### length and count of  target tag ###
def target_SV_distribution(freq_file, bed_file, out_file):
    """
    freq_file:
    Tag	Ref_freq	Alt_freq	MAF
    1_10059-1_10060-3242-INS	0.9975	0.0025	0.0025	1
    1_10380-1_10381-321-INS	0.9988	0.0012	0.0012	1

    bed_file:
    1	939439	939444	1_939439-1_939575-122-INS	1	939271	939460	SAMD11	ENST00000622503.4	CDS
    1	939439	939444	1_939439-1_939575-122-INS	1	939274	939460	SAMD11	ENST00000342066.7	CDS
    1	939439	939492	1_939439-1_939492-52-DEL	1	939271	939460	SAMD11	ENST00000622503.4	CDS
    """
    TagCount = {}
    freq_h = open(freq_file, "r")
    header = freq_h.readline()
    for line in freq_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        count = lines[-1]
        TagCount[tag] = count
    freq_h.close()
    bed_h = open(bed_file, "r")
    TagSets = set()
    for line in bed_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        TagSets.add(tag)
    bed_h.close()
    out_h = open(out_file, "w")
    out_h.write("Tag\tLength\tCount\n")
    Tags = sorted(list(TagSets))
    for t in Tags:
        c = TagCount[t]
        length = t.split("-")[-2]
        out_h.write("%s\t%s\t%s\n" % (t, length, c))
    out_h.close()
        




#############  target category ############
def target_category(cate_file, tag_file, out_file):
    import collections
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

    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)


    CatTags = collections.defaultdict(list)
    tag_h = open(tag_file, "r")
    header = tag_h.readline()
    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        if tag in TagCat:
            cat = TagCat[tag]
            CatTags[cat].append(tag)
        else:
            print("Please check whether the tag name %s is in file %s." % (tag, cate_file))
            sys.exit(1)
    tag_h.close()

    for c in CatTags:
        tags = sorted(CatTags[c])
        out_h.write("%s\t%s\n" % (c, ",".join(tags)))
    out_h.close()
    


######################################

def homo_tag_freq(tag_file, out_file, totalNum):
    import collections
    CatTags = collections.defaultdict(list)
    totalNum = int(totalNum)
    tag_h = open(tag_file, "r")
    header = tag_h.readline()
    for line in tag_h:
        lines = line.strip().split("\t")
        tag, Sample = lines
        samples = Sample.split(",")
        sampleNum = len(samples)
        ratio = sampleNum / totalNum
        if sampleNum == 1:
            cat = "Singleton"
        elif sampleNum != 1 and ratio <= 0.01:
            cat = "Rare"
        elif ratio > 0.01 and ratio <= 0.05:
            cat = "Low"
        elif ratio > 0.05:
            cat = "Common"
        CatTags[cat].append(tag)
    tag_h.close()

    out_h = open(out_file, "w")
    out_h.write("Category\tTags\n")
    for c in CatTags:
        tags = sorted(CatTags[c])
        out_h.write("%s\t%s\n" % (c, ",".join(tags)))
    out_h.close()



########### target tag annotation ################
def target_tag_annotation(anno_bed, tag_file, out_file):
    TagSet = {}
    tag_h = open(tag_file, "r")
    header = tag_h.readline().strip()
    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        TagSet[tag] = 1
    tag_h.close()
    
    anno_h = open(anno_bed, "r")
    out_h = open(out_file, "w")
    for line in anno_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[3]
        if tag in TagSet:
            out_h.write("%s\n" % line)
    anno_h.close()
    out_h.close()








############### category genotypes ###############
def genotype_category(category_file, genotype_file, outPrefix):
    import collections
    TagCategory = {}
    cat_h = open(category_file, "r")
    header = cat_h.readline()
    for line in cat_h:
        lines = line.strip().split("\t")
        cat, Tag = lines
        tags = Tag.split(",")
        for t in tags:
            TagCategory[t] = cat
    cat_h.close()

    CategoryRecord = collections.defaultdict(list)
    geno_h = open(genotype_file, "r")
    header = geno_h.readline().strip()
    for line in geno_h:
        line = line.strip()
        lines = line.split("\t")
        Tag = "%s_%s-%s_%s-%s-%s" % tuple(lines[:6])
        if Tag in TagCategory:
            cat = TagCategory[Tag]
            CategoryRecord[cat].append(line)
        else:
            print("Please check whether the tag %s is in category file." % (cat, category_file))
    geno_h.close()

    ### output the file
    outdir = "/".join(outPrefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i in CategoryRecord:
        records = CategoryRecord[i]
        out_h = open("%s_%s.txt" % (outPrefix, i), "w")
        out_h.write("%s\n" % header)
        out_h.write("%s\n" % "\n".join(records))
        out_h.close()




############### PCA add regions for samples ############
def column_index(records, column):
    column = column.lower()
    records = [r.lower() for r in records]
    try:
        columnIndex = records.index(column)
        return columnIndex
    except ValueError:
        print("Please check whether the column name %s is in record %s." % (column, records))
        sys.exit(1)

def PCA_sample_region(meta_file, PCA_file, out_file):
    """
    meta_file:
    Sample	Sex	Age	Province	Region	Group	Type
    M416-1	Male	42	Unknown	Unknown	M416	Parent
    """
    SampleRegion = {}
    meta_h = open(meta_file, "r")
    headers = meta_h.readline().strip().split("\t")
    regionIndex = column_index(headers, "region")
    for line in meta_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        region = lines[regionIndex]
        SampleRegion[sample] = region
    meta_h.close()

    PCA_h = open(PCA_file, "r")
    out_h = open(out_file, "w")
    header = PCA_h.readline().strip()
    out_h.write("%s\tRegion\n" % header)    
    for line in PCA_h:
        line = line.strip()
        lines = line.split("\t")
        sample = lines[0]
        try:
            region = SampleRegion[sample]
            out_h.write("%s\t%s\n" % (line, region))
        except KeyError:
            print("Please check whther the sample %s is in the meta file %s." % (sample, meta_file))
    PCA_h.close()
    out_h.close()




###### target tag frequency  #########
def target_tag_freq(tag_file, freq_file, out_file, method):
    import sys
    Tags = []
    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        if lines[0] != "":
            tag = lines[0]
            Tags.append(tag)
    tag_h.close()

    freq_h = open(freq_file, "r")
    out_h = open(out_file, "w")
    header = freq_h.readline().strip()
    out_h.write("%s\n" % header)

    method = method.lower()
    if method == "include":
        for line in freq_h:
            line = line.strip()
            lines = line.split("\t")
            Tag = lines[0]
            if Tag in Tags:
                out_h.write("%s\n" % line)
    elif method == "exclude":
        for line in freq_h:
            line = line.strip()
            lines = line.split("\t")
            Tag = lines[0]
            if Tag not in Tags:
                out_h.write("%s\n" % line)
    else:
        print("Please make sure that the method parameter is 'include' or 'exclude'.")
        sys.exit(1)

    freq_h.close()
    out_h.close()




########### Merge summary of Multiple files #############
def multiple_file_SV_summary(fileStr, out_file):
    import numpy
    out_h = open(out_file, "w")
    files = fileStr.split(",")
    files = [f.strip() for f in files]
    for i in range(len(files)):
        f = files[i]
        ff = f.split("/")
        fdir = ff[-2]

        in_h = open(f, "r")
        header = in_h.readline().strip()
        if i == 0:
            out_h.write("%s\n" % header)


        Records = []
        for line in in_h:
            lines = line.strip().split("\t")
            record = lines[1:]
            record = [int(r) for r in record]
            Records.append(record)
        in_h.close()
        NewRecord = numpy.transpose(Records)
        print(NewRecord)
        number = [str(sum(i)) for i in NewRecord]
        out_h.write("%s\t%s\n" % (fdir, "\t".join(number)))
    out_h.close()






############# merge record ##############
def merge_record(Sigs, out_file):
    for i in range(len(Sigs)):
        sig = Sigs[i]
        if i == 0:
            cmd = "cat %s > %s" % (sig, out_file)
        else:
            cmd = "sed -n '2p' %s >> %s" % (sig, out_file)
        os.system(cmd)



############## meta target #####################
def meta_target_file(targetList, meta_file, pair_name):
    tumorName = PairInfor(meta_file).pair_tumor(pair_name)
    normalName = PairInfor(meta_file).pair_normal(pair_name)
    for target in targetList:
        if tumorName in target:
            tumorFile = target
        elif normalName in target:
            normalFile = target
    return normalFile, tumorFile



########### merge rcord with one row title #############
def merge_record(Sigs, out_file):
    for i in range(len(Sigs)):
        sig = Sigs[i]
        if i == 0:
            cmd = "cat %s > %s" % (sig, out_file)
        else:
            cmd = "sed -n '2p' %s >> %s" % (sig, out_file)
        os.system(cmd)



############# target file #################
def target_file(samples, files):
    targetFiles = []
    for f in files:
        if isinstance(samples, list):
            for s in samples:
                if s in f:
                    targetFiles.append(f)
        elif isinstance(samples, str):
            if samples in f:
                targetFiles.append(f)
    return targetFiles



################## merge bam summary #####################
def read_file(file):
    in_h = open(file, "r")
    header = in_h.readline().strip()
    records = in_h.readline().strip().split("\t")
    record = "\t".join(records[1:])
    in_h.close()
    return header, record


def merge_bam_stats(bam_files, out_file):
    FileRecord = {}
    files = bam_files.split(",")
    for f in files:
        file_base = f.split("/")[-1].split("_")[0]
        header, record = read_file(f)
        FileRecord[file_base] = record
    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)
    for ff in sorted(list(FileRecord.keys())):
        rec = FileRecord[ff]
        out_h.write("%s\t%s\n" % (ff, rec))
    out_h.close()
    





####################### parse pair for metafile #################
def meta_target_file(targetList, meta_file, pair_name):
    tumorName = PairInfor(meta_file).pair_tumor(pair_name)
    normalName = PairInfor(meta_file).pair_normal(pair_name)
    for target in targetList:
        if tumorName in target:
            tumorFile = target
        elif normalName in target:
            normalFile = target
    return normalFile, tumorFile



######################### remove balsted tag ########################
def remove_short_tag(tag_file, out_file):
    """
    tag_file:
    M416-1:contig1  M426-1:contig1072       7229    19000   0.962
    M416-1:contig3  M426-1:contig394        3000    3000    0.912
    M416-1:contig14 M426-1:contig1684       3000    7680    1.000
    """
    in_h = open(tag_file, "r")
    RemovedTags = set()
    for line in in_h:
        lines = line.strip().split("\t")
        qseqid, sseqid, qlen, slen = lines[:4]
        qlen = int(qlen)
        slen = int(slen)
        if qlen <= slen:
            RemovedTags.add(qseqid)
        else:
            RemovedTags.add(sseqid)
    in_h.close()

    sortedTags = sorted(list(RemovedTags))
    out_h = open(out_file, "w")
    for t in sortedTags:
        out_h.write("%s\n" % t)
    out_h.close()




#######################################################################
def group_samples(meta_file, target):
    import collections
    import sys
    
    meta_h = open(meta_file, "r")
    headers = meta_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]
    try:
        groupIndex = headers.index("group")
        sampleIndex = headers.index("sample")
    except KeyError:
        print("Please check whether the header of file %s contain the column 'group' and 'sample'." % meta_file)
        sys.exit(1)

    GroupSamples = collections.defaultdict(list)

    for line in meta_h:
        lines = line.strip().split("\t")
        group = lines[groupIndex]
        sample = lines[sampleIndex]
        GroupSamples[group].append(sample)
    meta_h.close()

    try:
        samples = GroupSamples[target]
    except KeyError:
        print("Please check whether the target group %s is in the file %s." % (target, meta_file))
        sys.exit(1)
    return samples




##################################################
def change_tag(bed_file, fasta1, fasta2):
    from tinyfasta import FastaParser
    TagInfor = {}
    bed_h = open(bed_file, "r")
    for line in bed_h:
        lines = line.strip().split("\t")
        Chr, Start, End, Tag = lines
        origTag = "%s:%s-%s" % (Chr, Start, End)
        TagInfor[origTag] = Tag

    ### whether file is empty
    if os.stat(fasta1).st_size == 0:
        cmd = "touch %s" % fasta2
        os.system(cmd)
    else:
        out_h = open(fasta2, "w")
        for record in FastaParser(fasta1):
            desc = str(record.description)
            desc = desc.lstrip(">")
            newTag = TagInfor[desc]
            seq = str(record.sequence)
            out_h.write(">%s\n%s\n" % (newTag, seq))
        out_h.close()




##################################################


################################# Merge statistics of multiple samples ########################
def merge_multiple_sample_stats(bamSums, out_file):
    for i in range(len(bamSums)):
        bamStat = bamSums[i]
        if i == 0:
            cmd  = "cat %s > %s" % (bamStat, out_file)
        else:
            cmd = "sed '1d' %s >> %s" % (bamStat, out_file)
        os.system(cmd)

###############################################################################################


###########################################################
def VCF_variant_to_bed(vcf_file, bed_file):
    """
    vcf_file:
    1   10440   .   C   A   70.57   .   AC=1;AF=0.25;AN=4;BaseQRankSum=2.19;ClippingRankSum=0;DP=75;ExcessHet=3.0103;FS=8.451;MLEAC=2;MLEAF=0.5;MQ=16.92;MQRankSum=-1.006;QD=10.08;ReadPosRankSum=-1.133;SOR=2.788  0/1:4,3:7:96:0|1:10440_C_A:96,0,145
    1   10447   .   C   A   67.57   .   AC=1;AF=0.25;AN=4;BaseQRankSum=1.98;ClippingRankSum=0;DP=28;ExcessHet=3.9794;FS=3.68;MLEAC=2;MLEAF=0.5;MQ=40.29;MQRankSum=-0.792;QD=9.65;ReadPosRankSum=-1.068;SOR=2.258    0/1:5,2:7:93:0|1:10440_C_A:93,0,170
    1   13838   .   C   T   30.05   .   AC=1;AF=0.167;AN=6;BaseQRankSum=1.38;ClippingRankSum=0;DP=11;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.167;MQ=22.77;MQRankSum=-1.383;QD=7.51;ReadPosRankSum=-0.674;SOR=0.693    0/1:2,2:4:60:0|1:13813_T_G:60,0,101

    bed_file:
    1_10440 1   10440   C   A
    1_10447 1   10447   C   A
    1_13838 1   13838   C   T
    """
    in_h = open(vcf_file, "r")
    out_h = open(bed_file, "w")
    for line in in_h:
        lines = line.strip().split("\t")
        if lines[0].startswith("#"):
            continue
        Chr, Pos = lines[:2]
        Ref, Alt = lines[3:5]
        ID = Chr + "_" + Pos
        out_h.write("%s\t%s\t%s\t%s\t%s\n" % (ID, Chr, Pos, Ref, Alt))
    in_h.close()
    out_h.close()
##############################################################



#################################### get sample from group information ######################################

def parser_group(group_infor, column):
    """
    args:
    group_infor (file):
    Sample  Group
    M625-0  M625
    M625-1  M625
    M625-2  M625
    
    column (string):
    "Group"

    
    returns:
    GroupSamples (dict):
    {"M625": ["M625-0", "M625-1", "M625-2"]}
    """
    GroupSamples = collections.defaultdict(list)

    column = str(column).strip()
    group_h = open(group_infor, "r")
    ### get the column of group information of samples
    headers = group_h.readline().strip().split("\t")
    targetIndex = None
    for i in range(len(headers)):
        title = headers[i]
        if title == column:
            targetIndex = i
            break
    if targetIndex == None:
        print("Please check whether the target column of group information is %s." % column)
        sys.exit(1)

    ### get the information of group and samples.
    for line in group_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        group = lines[targetIndex]
        GroupSamples[group].append(sample)
    group_h.close()

    return GroupSamples



def get_samples(group_infor, column, fileString, group):
    """
    args:
    fileString (string):
    "M625-0.gvcf,M625-1.gvcf,M625-2.gvcf"

    group (string):
    "M625"

    returns:
    Targets (list):
    file_name
    """
    GroupSamples = parser_group(group_infor, column)
    ## GroupSamples: {'M625': ['M625-0', 'M625-1', 'M625-2'], 'M628': ['M628-0', 'M628-1', 'M628-2']}

    files = fileString.split(",")
    files = [f.strip() for f in files]
    ### get base name of the file
    fileBases = [f.split("/")[-1].split(".")[0].split("_")[0] for f in files]
    ## fileBases: ['M625-0', 'M625-1', 'M625-2', 'M628-0', 'M628-1', 'M628-2']

    ### get the samples of the group name based on the file group_infor
    try:
        samples = GroupSamples[group]
    except KeyError:
        print("Please check whether the target group name %s is in the dict %s." % (group, GroupSamples))
        sys.exit(1)

    Targets = []
    for i in range(len(fileBases)):
        fileName = files[i]
        if fileBases[i] in samples:
            Targets.append(fileName)
    # print(sorted(Targets))
    return sorted(Targets)

#################################################################################################################



################################## parser the meta information

class MetaInfor:
    """
    Parser the meta information file which contains the group and sample relationship.
    """
    def __init__(self, group_file):
        """
        argv:
            group_file:
            Sample  Group   Type
            M625-0  M625    Proband
            M625-1  M625    Parent
            M625-2  M625    Parent
        """
        self.group_file = group_file
        self.GroupSamples = collections.defaultdict(list)
        self.SampleGroup = {}
        self.GroupParents = collections.defaultdict(list)
        self.GroupProband = {}

    def parse_meta_information(self):
        in_h = open(self.group_file, "r")
        header = in_h.readline()
        for line in in_h:
            lines = line.strip().split("\t")
            try:
                sample, group, relationship = lines[:3]
                relationship = relationship.lower()
            except IndexError:
                print("Please check and make sure that the group file %s should contain three colunms, such as 'Sample', 'Geoup' and 'Type'." % self.group_file)
                sys.exit(1)

            ### store the information in the dict
            self.GroupSamples[group].append(sample)
            self.SampleGroup[sample] = group

            if relationship == "proband":
                self.GroupProband[group] = sample
            elif relationship == "parent":
                self.GroupParents[group].append(sample)
            else:
                print("Please make sure that the string in third column of file %s is 'proband' or 'parents'." % self.group_file)
                sys.exit(1)
        in_h.close()


    def group_proband(self, Trio_name):
        """
        argv:
            'M625' (string): group anme

        return:
            'M625-0' (string): proband name
        """
        self.parse_meta_information()
        try:
            proband = self.GroupProband[Trio_name]
        except KeyError:
            sys.exit("Please ckeck whether the trio group name %s is in the meta information file %s." % (Trio_name, self.group_file))
        return proband

    def group_parents(self, Trio_name):
        """
        argv:
            'M625' (string): group name

        return:
            ['M625-1', 'M625-2'] (list): parents' names
        """
        self.parse_meta_information()
        try:
            parents = self.GroupParents[Trio_name]
        except KeyError:
            sys.exit("Please ckeck whether the trio group name %s is in the meta information file %s." % (Trio_name, self.group_file))
        return parents

    def sample_group(self, sample_name):
        """
        argv:
            'M625-0' (string): proband name

        return:
            'M625' (string): group name
        """
        self.parse_meta_information()
        try:
            group = self.SampleGroup[sample_name]
        except KeyError:
            sys.exit("Please ckeck whether the sample name %s is in the meta information file %s." % (sample_name, self.group_file))
        return group

    def group_sample(self, trio_name):
        """
        argv:
            'M625' (string): group name

        return:
            ['M625-0', 'M625-1', 'M625-2'] (list): all sample names
        """
        self.parse_meta_information()
        try:
            samples = self.GroupSamples[trio_name]
        except KeyError:
            sys.exit("Please ckeck whether the group name %s is in the meta information file %s." % (sample_name, self.group_file))
        return samples

#################################################################################



############################## tumor and normal pair information ##########################

class PairInfor:
    def __init__(self, meta_file):
        """
        meta_file:

        Sample  Pair    Type
        SRR60   SRR6    tumor
        SRR62   SRR6    normal
        """
        self.meta_file = meta_file
        self.PairSample = collections.defaultdict(dict)
        self.SamplePair = {}
        self.SampleType = {}
        self.SampleSex = {}

    def parse_pair_meta_infor(self):
        in_h = open(self.meta_file, "r")
        headers = in_h.readline().strip().split("\t")
        headerNames = [h.lower() for h in headers]
        headerSet = set(headerNames)
        if headerSet == {"sample", "pair", "type", "sex"}:
            sampleIndex = headerNames.index("sample")
            pairIndex = headerNames.index("pair")
            typeIndex = headerNames.index("type")
            sexIndex = headerNames.index("sex")
        else:
            print("Please check that the columns, such as 'sample', 'pair', 'type' (insenstive of case) should be in the header line.")
            sys.exit(1)

        for line in in_h:
            line = line.strip()
            if line != "":
                lines = line.split("\t")
                Sample = lines[sampleIndex]
                Pair = lines[pairIndex]
                Type = lines[typeIndex]
                Sex = lines[sexIndex]

                if Type in ["tumor", "normal"]:
                    self.SamplePair[Sample] = Pair
                    self.SampleType[Sample] = Type
                    self.SampleSex[Sample] = Sex
                    self.PairSample[Pair][Type] = Sample

                else:
                    print("Please make sure that the column Type should be 'tumor' or 'normal'.")
                    sys.exit(1)
        in_h.close()

    def pair_tumor(self, pair_name):
        self.parse_pair_meta_infor()
        try:
            tumorSample = self.PairSample[pair_name]["tumor"]
        except KeyError:
            sys.exit("Please check whetehr the pair name %s is in the metafile %s." % (pair_name, self.meta_file))
        return tumorSample

    def pair_normal(self, pair_name):
        self.parse_pair_meta_infor()
        try:
            normalSample = self.PairSample[pair_name]["normal"]
        except KeyError:
            sys.exit("Please check whetehr the pair name %s is in the metafile %s." % (pair_name, self.meta_file))
        return normalSample

    def sample_pair(self, sample_name):
        self.parse_pair_meta_infor()
        try:
            pairName = self.SamplePair[sample_name]
        except KeyError:
            sys.exit("Please check whetehr the sample name %s is in the metafile %s." % (sample_name, self.meta_file))
        return pairName

    def sample_type(self, sample_name):
        self.parse_pair_meta_infor()
        try:
            typeName = self.SampleType[sample_name]
        except KeyError:
            sys.exit("Please check whetehr the sample name %s is in the metafile %s." % (sample_name, self.meta_file))
        return typeName

    def sample_sex(self, sample_name):
        self.parse_pair_meta_infor()
        try:
            sexName = self.SampleSex[sample_name]
        except KeyError:
            sys.exit("Please check whetehr the sample name %s is in the metafile %s." % (sample_name, self.meta_file))
        return sexName
    
  

###########################################################################################



