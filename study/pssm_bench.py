import pytest
import time
import random
import math

"""
weight_matrix = []
peptide = []
peptide_len = 0;

score_max = 0;
score_min = 0;
score_range = 0;
"""

alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
alph_len = len(alphabet)
peptide = "CLMPCGRRQ"
peptide_len = len(peptide)

def random_row():
    return dict(zip(alphabet, [random.uniform(-2, 2) for x in range(alph_len)]))

weight_matrix = [random_row() for y in range(9)]

score_max = 0.8;
score_min = 0.8 * (1 - math.log(50000) / math.log(500));
score_range = score_max - score_min;


def setup():
    pass

def teardown():
    pass

def setup_module(module):
    pass
 
def teardown_module(module):
    pass
 
def setup_function(function):
    pass
 
def teardown_function(function):
    pass

def score_one_peptide():
    score = 0
    for ch_pos in range(len(peptide)):
        score += weight_matrix[ch_pos][peptide[ch_pos]]

    score = score / peptide_len;
    score = max(min(score, score_max), score_min)
    ic50 = 50000 **((score_max - score)/score_range);

    return ic50

def score_one_peptide_loop():

    with open("res_py.txt", "w") as f:
        for i in range(10000000):
            ic50 = score_one_peptide();          

        for i in range(1000000):
            f.write("{0} {1}\n".format(peptide, ic50))
    

def test_score_one_peptide(benchmark):
    result = benchmark(score_one_peptide_loop)
