#!/usr/local/bin/python

import pandas
import numpy
import pybedtools
import argparse
from argparse import HelpFormatter
from operator import attrgetter

class SortingHelpFormatter(HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)

parser = argparse.ArgumentParser(description='Generate Interaction Maps of Given Size (bp)')
parser.add_argument('-i', '--inputFile', dest='inputFile', metavar='FILE',
                    help='Tab-delimited non-binary input file [Required]')
parser.add_argument('-c', '--chrom', dest='chrom', required=False,
                    default='chr6', help='Chromosome to calculate Maps [chr6]')
parser.add_argument('-r', '--res', dest='res', required=False,
                    default=2000, help='Resolution (bp) to make the matrix [2000]')
parser.add_argument('-t', '--threads', dest='threads', required=False,
                    help='Number of processors [1]', default=1)
parser.add_argument('-v', '--version',
                    action='version', version='%(prog)s (myprog version 0.1)')
args = parser.parse_args()

def generate_matrix_from_reference(chromosome,res):
    totalChromSize = pybedtools.chromsizes('hg19')[chromosome][1]

    names = []
    for bp in range(0,totalChromSize,res):
        if bp+res < totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(bp+res))
        elif bp+res >= totalChromSize:
            names.append(str(chromosome) + ":" + str(bp) + "-" + str(totalChromSize))

    ### Create an empty matrix that will be filled with counts
    dim = len(names)
    matrix = pandas.DataFrame(numpy.zeros((dim,dim)))

    matrix.columns = names
    matrix.index = names
    return(matrix)

def file_block(fp, number_of_blocks, block):
    '''
    A generator that splits a file into blocks and iterates
    over the lines of one of the blocks.
    http://xor0110.wordpress.com/2013/04/13/how-to-read-a-chunk-of-lines-from-a-file-in-python/
    '''

    assert 0 <= block and block < number_of_blocks
    assert 0 < number_of_blocks

    fp.seek(0,2)
    file_size = fp.tell()

    ini = file_size * block / number_of_blocks
    end = file_size * (1 + block) / number_of_blocks

    if ini <= 0:
        fp.seek(0)
    else:
        fp.seek(ini-1)
        fp.readline()

    while fp.tell() < end:
        yield fp.readline()

if __name__ == '__main__':
    matrix = generate_matrix_from_reference(str(args.chrom), int(args.res))
    print matrix.iloc[3:3,:]

    fp = open(args.inputFile)
    number_of_chunks = int(args.threads)
    for chunk_number in range(number_of_chunks):
        print chunk_number, 100 * '='
        for line in file_block(fp, number_of_chunks, chunk_number):
            readName, chr1, start1, end1, strand1, score1, chr2, start2, end2, strand2, score2, flags = line.strip().split('\t')
            print flags


