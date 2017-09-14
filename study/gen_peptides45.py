import itertools

def genFragment(size):
    return list(itertools.permutations('ACDEFGHIKLMNPQRSTVWY', size))

p4 = genFragment(4)
p5 = genFragment(5)

for i in range(100000):
  p = p4[i] + p5[i]
  permStr = ""
  for c in p:
    permStr +=c
  print ">", i
  print permStr
  i += 1
