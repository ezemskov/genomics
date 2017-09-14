#! /usr/bin/python

import itertools

def genFragment(size):
    return list(itertools.permutations('ACDEFGHIKLMNPQRSTVWY', size))

p6 = genFragment(6)

for i in range(3000000):
  p = p6[i] + p6[i]
  permStr = ""
  for c in p[0:9]:
    permStr +=c
  print ">", i
  print permStr
  i += 1
