#!/usr/bin/env python
import collections
import argparse
import sys


#usage: python ~/github/NanoHub/src/NanoHub/parseMetaInfor.py --meta meta_information.txt --group M625 --sample M625-0

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.04.26"


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



def main():
    parser = argparse.ArgumentParser(description="Get the target information based on the meta information.")
    parser.add_argument("-m", "--meta", help="The input meta information file containing sample and group information.")
    parser.add_argument("-g", "--group", help="The group name.")
    parser.add_argument("-s", "--sample", help="The sample name.")
    args = parser.parse_args()


    proband_name = MetaInfor(args.meta).group_proband(args.group)
    print(proband_name)
    parent_names = MetaInfor(args.meta).group_parents(args.group)
    print(parent_names)
    group_name = MetaInfor(args.meta).sample_group(args.sample)
    print(group_name)
    sample_names = MetaInfor(args.meta).group_sample(args.group)
    print(sample_names)


if __name__ == "__main__":
    main()

