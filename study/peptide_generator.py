#! /usr/bin/python

import itertools
import sys
from datetime import datetime


def main():
    symbolsStr = 'ACDEFGHIKLMNPQRSTVWY'

    perms = itertools.product(symbolsStr, repeat=9)
    
    start = 5 * 10E5
    count = 10E5
    
    if (len(sys.argv) >= 3):
        start = int(sys.argv[1])
        count = int(sys.argv[2])

    slice_of_product = itertools.islice(perms, start, start+count)
    for p in slice_of_product :
        print '>'
        print ''.join(p)


if __name__ == "__main__":
    main()