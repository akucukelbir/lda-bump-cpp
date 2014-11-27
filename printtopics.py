#!/usr/local/bin/python

# printtopics.py: Prints the words that are most prominent in a set of
# topics.
#
# Copyright (C) 2014  Alp Kucukelbir (modified)
# Copyright (C) 2010  Matthew D. Hoffman
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, os, re, random, math, urllib2, time, cPickle
import numpy

#import onlineldavb

def main():
    """
    Displays topics fit by onlineldavb.py. The first column gives the
    (expected) most prominent words in the topics, the second column
    gives their (expected) relative prominence.
    """
    vocab      = str.split(file(sys.argv[1]).read())
    testlambda = numpy.loadtxt(sys.argv[2])

    for k in xrange(0, len(testlambda), 3):
        lambdak = list(testlambda[k, :])
        lambdak = lambdak / sum(lambdak)
        temp = zip(lambdak, range(0, len(lambdak)))
        temp = sorted(temp, key = lambda x: x[0], reverse=True)

        if k+1 < len(testlambda):
            lambdak2 = list(testlambda[k+1, :])
            lambdak2 = lambdak2 / sum(lambdak2)
            temp2 = zip(lambdak2, range(0, len(lambdak2)))
            temp2 = sorted(temp2, key = lambda x: x[0], reverse=True)

        if k+2 < len(testlambda):
            lambdak3 = list(testlambda[k+2, :])
            lambdak3 = lambdak3 / sum(lambdak3)
            temp3 = zip(lambdak3, range(0, len(lambdak3)))
            temp3 = sorted(temp3, key = lambda x: x[0], reverse=True)

        print '|topic %d|' % (k),
        if k+2 >= len(testlambda):
            if k+1 >= len(testlambda):
                print
            else:
                print '%35s' % (''),
                print '|topic %d|' % (k+1)
        else:
                print '%35s' % (''),
                print '|topic %d|' % (k+1),
                print '%35s' % (''),
                print '|topic %d|' % (k+2)


        for i in range(0, 10):
            print '%20s --- %.4f' % (vocab[temp[i][1]], temp[i][0]),
            if k+2 >= len(testlambda):
                if k+1 >= len(testlambda):
                    print
                else:
                    print '%10s' % (''),
                    print '%20s --- %.4f' % (vocab[temp2[i][1]], temp2[i][0])
            else:
                print '%10s' % (''),
                print '%20s --- %.4f' % (vocab[temp2[i][1]], temp2[i][0]),
                print '%10s' % (''),
                print '%20s --- %.4f' % (vocab[temp3[i][1]], temp3[i][0])

        print

if __name__ == '__main__':
    main()
