#! /usr/bin/python



def local_permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = range(n)
    cycles = range(n, n-r, -1)
    idx = 0
    #print idx, tuple(pool[i] for i in indices[:r])

    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i

            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]

                idx += 1
                #print idx, tuple(pool[i] for i in indices[:r])

                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

def local_product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    #print len(list(pools))
    #print pools
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def local_islice(iterable, *args):
    # islice('ABCDEFG', 2) --> A B
    # islice('ABCDEFG', 2, 4) --> C D
    # islice('ABCDEFG', 2, None) --> C D E F G
    # islice('ABCDEFG', 0, None, 2) --> A C E G
    s = slice(*args)
    it = iter(xrange(s.start or 0, s.stop or sys.maxint, s.step or 1))
    nexti = next(it)
    for i, element in enumerate(iterable):
        if i == nexti:
            yield element
            nexti = next(it)



import itertools
from datetime import datetime

symbolsStr = 'ACDEFGHIKLMNPQRSTVWY'

#benchmark : product
print "itertools"
print datetime.now().time()
perms = itertools.product(symbolsStr, repeat=9)
#iStart = 4000000000
#iEnd = 4000000010
#slice_of_product = itertools.islice(perms, iStart, iEnd)
#slice_of_product = itertools.islice(perms, 8000000000, 8000000010)
slice_of_product = local_islice(perms, 200000, 200010)
#slice_of_product = local_islice(perms, 4000000000, 4000000010)

print list(slice_of_product)
#list(perms)
#print len(list(perms))
print datetime.now().time()

print "local"
print datetime.now().time()
perms = local_product(symbolsStr, repeat=4)
#list(perms)
#print len(list(perms))
print datetime.now().time()


#benchmark : permutation
"""
print "itertools"
print datetime.now().time()
p = itertools.permutations(symbolsStr, 5)
print len(list(p))
print datetime.now().time()

print "\nlocal"
print datetime.now().time()
p = local_permutations(symbolsStr, 5)
print len(list(p))
print datetime.now().time()
"""

#demo
"""
print "itertools : "
print list(itertools.permutations(symbolsStr, 3))

print "local : "
print list(local_permutations(symbolsStr, 3))
"""

