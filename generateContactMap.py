#!/usr/local/bin/python

import collections
import re
import pandas
import pybedtools
import argparse
import os
from operator import attrgetter

## Class that sort command-line menu options
class SortingHelpFormatter(argparse.HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)

## Function to check if file name exist
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

parser = argparse.ArgumentParser(description='Generate Interaction Maps of Given Size (bp)')
parser.add_argument('-i', '--inputFile', default='chr6.txt', dest='inputFile', metavar='FILE',
                    help='Tab-delimited non-binary input file [Required]')
parser.add_argument('-o', '--outFile', default='chr6.matrix', dest='outFile', metavar='FILE',
                    help='Tab-delimited output file without extension [Required]')
parser.add_argument('-c', '--chrom', dest='chrom', required=False,
                    default='chr6', help='Chromosome to calculate Maps [chr6]')
parser.add_argument('-r', '--res', dest='res', required=False,
                    default=2000, help='Resolution (bp) to make the matrix [2000]')
parser.add_argument('-t', '--threads', dest='threads', required=False,
                    help='Number of processors [1]', default=1)
parser.add_argument('-v', '--version',
                    action='version', version='%(prog)s (myprog version 0.1)')
args = parser.parse_args()


def sort_keys_alphanumeric(x):
    number = int(re.sub('[^0-9]', '', x[0]) )
    length = len(x[0])
    return length, number

def generate_dictionary_of_bins(chromosome, res):
    totalChromSize = pybedtools.chromsizes('hg19')[chromosome][1]
    names = []
    for bp in range(0, totalChromSize, res):
        if bp + res < totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(bp + res))
        elif bp + res >= totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(totalChromSize))

    ## Create a dictionary of bins
    bin = int()
    binDictionary = dict()
    for items in range(0,len(names)):
        bin += 1
        binDictionary[str(chromosome) + ":" + str(bin)] = names[items]

    ## sort dictionary
#    binDictionary = collections.OrderedDict(sorted(binDictionary.items(), key=sort_keys_alphanumeric))
    return binDictionary


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
    ## generate matrix ---------------------------------------------------------
    binDictionary = generate_dictionary_of_bins(str(args.chrom), int(args.res))

    ## variable initialization
    matrix = collections.defaultdict(int)
    interaction1 = str()
    interaction2 = str()

    ## open input file # by chunks -----------------------------------------------
    fp = open(args.inputFile)
#    number_of_chunks = int(args.threads)
#    for chunk_number in range(number_of_chunks):
#        for line in file_block(fp, number_of_chunks, chunk_number):
    for line in fp:
        id, c1, s1, e1, sc1, st1, c2, s2, e2, sc2, st2, flags = line.rstrip().split('\t')
        firstFlag = False
        secondFlag = False
        for key, value in binDictionary.iteritems():
            binSize = str(value).replace(':', ' ').replace('-', ' ').split()[2]
            i1 = int(((int(s1)+int(e1))/2)/int(binSize))+1
            i2 = int(((int(s2)+int(e2))/2)/int(binSize))+1
            if i1 == 1 and firstFlag == False:
                interaction1 = value
                firstFlag = True
            if i2 == 1 and secondFlag == False:
                interaction2 = value
                secondFlag = True
            if all([firstFlag == True, secondFlag == True]):
                matrix[interaction1,interaction2] += 1
                break
#                if len(matrix) >= 20: break

    matrix2 = pandas.DataFrame(matrix.values(), index=pandas.MultiIndex.from_tuples(matrix.keys())).unstack(1)
    matrix2.fillna(0, inplace=True)
    args.outFile = str(args.outFile) + ".matrix.txt"
    matrix2.to_csv(args.outFile, sep='\t')
