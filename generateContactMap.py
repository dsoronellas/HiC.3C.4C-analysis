#!/usr/local/bin/python

import collections
import pybedtools
import argparse
import os
from time import gmtime, strftime
from operator import attrgetter

# # Class that sort command-line menu options
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

## Menu
parser = argparse.ArgumentParser(description='Generate Interaction Maps of Given Size (bp)')
parser.add_argument('-i', '--inputFile', default='chr6.txt', dest='inputFile', metavar='FILE',
                    help='Tab-delimited non-binary input file [Required]')
parser.add_argument('-o', '--outFile', default='chr6.matrix', dest='outFile', metavar='FILE',
                    help='Tab-delimited output file without extension [Required]')
parser.add_argument('-c', '--chrom', dest='chrom', required=False,
                    default='chr6', help='Chromosome to calculate Maps [chr6]')
parser.add_argument('-r', '--res', dest='res', required=False,
                    default=2000, help='Resolution (bp) to make the matrix [2000]')
parser.add_argument('-s', '--subset', dest='subset', required=False,
                    default='149900000-160000000',
                    help='Region to subset the matrix, needed for high-res matrices. [145000000-165000000]')
parser.add_argument('-t', '--threads', dest='threads', required=False,
                    help='Number of processors [1]', default=1)
parser.add_argument('-v', '--version',
                    action='version', version='%(prog)s (myprog version 1.0)')
args = parser.parse_args()


def file_block(fp, number_of_blocks, block):
    '''
   A generator that splits a file into blocks and iterates
   over the lines of one of the blocks.
   http://xor0110.wordpress.com/2013/04/13/how-to-read-a-chunk-of-lines-from-a-file-in-python/
   '''

    assert 0 <= block < number_of_blocks
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

    args.outFile = str(args.outFile) + ".matrix.txt"

    print '\nContact Map Generation Script v1.0'
    print '----------------------------------'
    print 'Input File = ' + str(args.inputFile)
    print 'Output File = ' + str(args.outFile)
    print 'Resolution = ' + str(args.res) + 'bp'
    print 'Interaction between = ' + str(args.chrom) + ' vs ' + str(args.chrom)
    print 'Region subset is = ' + str(args.chrom) + ':' + str(args.subset)
    print 'Threads = ' + str(args.threads) + '\n'

    ## variable initialization
    matrix = collections.defaultdict(int)
    interaction1 = str()
    interaction2 = str()

    ## open input file # by chunks -----------------------------------------------
    print 'Finding interactions ...\n'

    # Log every N lines.
    LOG_EVERY_N = 1000000
    factor = 1
    counter = 0

    fp = open(args.inputFile)
    # number_of_chunks = int(args.threads)
    # for chunk_number in range(number_of_chunks):
    # for line in file_block(fp, number_of_chunks, chunk_number):

    for line in fp:
        ## Count lines
        counter += 1

        ## Get fields into variables
        id, c1, s1, e1, sc1, st1, c2, s2, e2, sc2, st2, flags = line.rstrip().split('\t')

        ## Map each mate at bins ...
        i1 = int(((int(s1) + int(e1)) / 2) / int(args.res)) + 1
        i2 = int(((int(s2) + int(e2)) / 2) / int(args.res)) + 1

        ## Generate a dictionary of bins with interaction mapped ...
        matrix[str(args.chrom) + ":" + str(i1), str(args.chrom) + ":" + str(i2)] += 1

        ## Log for number of lines processed ...
        if (counter % LOG_EVERY_N) == 0:
            print '[INFO ' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' ] ' + str(
                LOG_EVERY_N * factor) + ' fragments processed ...'
            factor += 1
    else:
        lastLine = LOG_EVERY_N * (factor - 1)
        if counter != lastLine:
            ## Last Line processed log
            print '[INFO ' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' ] ' + str(counter) + ' fragments processed ...'

    fp.close()

    ## generate a list of bins
    totalBins = int(pybedtools.chromsizes('hg19')[args.chrom][1]) / int(args.res)
    bins = list()
    for splits in range(0, totalBins):
        bins.append(str(args.chrom) + ":" + str(splits + 1))

    ## generate a list of bins
    print '\nSaving results to a matrix-like tab-delimited file'
    start = int(args.subset.split('-')[0]) / int(args.res)
    stop = (int(args.subset.split('-')[1]) / int(args.res)) + 1

    ## Build matrix ------
    fp = open(args.outFile, 'w')

    for z in range(start, stop):  ## Header print
        if z == start:
            fp.write('%s\t%s\t' % ('bins', bins[z]))
        elif z < stop:
            fp.write('%s\t' % bins[z])
        elif z == stop:
            fp.write('%s\n' % bins[z])

    for i in range(start, stop):  ## Row iterator
        fp.write('%s' % bins[i])
        for j in range(start, stop):  ## Column iterator
            fp.write('\t%s' % matrix[bins[i], bins[j]])
        fp.write('\n')

    fp.close()

