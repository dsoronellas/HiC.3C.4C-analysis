#!/usr/local/bin/python

import pandas
import numpy
import pybedtools
import argparse
import os
from operator import attrgetter

# # Class that sort command-line menu options
class SortingHelpFormatter(argparse.HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle


parser = argparse.ArgumentParser(description='Generate Interaction Maps of Given Size (bp)')
parser.add_argument('-i', '--inputFile', dest='inputFile', metavar='FILE',
                    help='Tab-delimited non-binary input file [Required]',
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument('-c', '--chrom', dest='chrom', required=False,
                    default='chr6', help='Chromosome to calculate Maps [chr6]')
parser.add_argument('-r', '--res', dest='res', required=False,
                    default=2000, help='Resolution (bp) to make the matrix [2000]')
parser.add_argument('-t', '--threads', dest='threads', required=False,
                    help='Number of processors [1]', default=1)
parser.add_argument('-v', '--version',
                    action='version', version='%(prog)s (myprog version 0.1)')
args = parser.parse_args()


def generate_matrix_from_reference(chromosome, res):
    totalChromSize = pybedtools.chromsizes('hg19')[chromosome][1]
    names = []
    for bp in range(0, totalChromSize, res):
        if bp + res < totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(bp + res))
        elif bp + res >= totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(totalChromSize))
    ### Create an empty matrix that will be filled with counts
    dim = len(names)
    matrix = pandas.DataFrame(numpy.zeros((dim, dim)))
    matrix.columns = names
    matrix.index = names
    return (matrix, names)


def compare_and_count(bin1, bin2, line):
    s_bin1 = bin1.replace(':', ' ').replace('-', ' ').split()
    s_bin2 = bin2.replace(':', ' ').replace('-', ' ').split()
    id1, c1, s1, e1, sc1, st1, c2, s2, e2, sc2, st2 = line.strip().split('\t')

    interaction = 0
    if bin1 == bin2:
        # same chromosome and in-bin interaction
        if (c1 == c2 == s_bin1[0]) and all([ s1 >= s_bin1[1], e2 <= s_bin2[2] ]):
            interaction += 1
    elif bin1 != bin2:
        if all([ s_bin1[0] == s_bin2[0], (c1 == c2 == s_bin1[0]) ]):    # read is in target chromosome
            if all([all([ s1 >= s_bin1[1], s1 <= s_bin1[2] ]), ([ s2 >= s_bin2[1], s2 <= s_bin2[2] ]) ]):
                interaction += 1
        else



def file_block(fp, number_of_blocks, block):
    '''
    A generator that splits a file into blocks and iterates
    over the lines of one of the blocks.
    http://xor0110.wordpress.com/2013/04/13/how-to-read-a-chunk-of-lines-from-a-file-in-python/
    '''

    assert 0 <= block and block < number_of_blocks
    assert 0 < number_of_blocks

    fp.seek(0, 2)
    file_size = fp.tell()

    ini = file_size * block / number_of_blocks
    end = file_size * (1 + block) / number_of_blocks

    if ini <= 0:
        fp.seek(0)
    else:
        fp.seek(ini - 1)
        fp.readline()

    while fp.tell() < end:
        yield fp.readline()


if __name__ == '__main__':
    matrix, names = generate_matrix_from_reference(str(args.chrom), int(args.res))

    fp = open(args.inputFile)
    number_of_chunks = int(args.threads)
    for chunk_number in range(number_of_chunks):
        for line in file_block(fp, number_of_chunks, chunk_number):
            for bin1 in names:
                for bin2 in names:
                    matrix[bin1][bin2] = compare_and_count(bin1, bin2, line)


