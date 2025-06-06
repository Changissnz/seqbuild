from seqgen.rch_gen import * 
from intigers.idt_proc import * 
from morebs2.numerical_generator import prg__LCG 

import unittest

### lone file test 
"""
python -m tests.test_extraneous
"""
###
class ExtraneousMethods(unittest.TestCase):

    def test__extraneous__prg__integer_sets__mult__case1(self):
        # case 1 
        n = 10
        c = 2
        r = 0.2
        crange = (-5,6)
        mrange = (2,20)
        prg = prg__LCG(3,8,11,99)
        num_attempts_per_nc = 100

        pm = prg__integer_sets__mult(n,c,r,crange,\
        mrange,prg,num_attempts_per_nc)

        assert pm == [[-2, -20, -20], \
        [-3, -30, -57, -57, -30], \
        [-49, -37]]

        c = 0
        for pm_ in pm: c += len(pm_)
        assert c == 10 

        # case 2 
        n = 90
        c = 5
        r = 0.2
        crange = (-15,60)
        mrange = (-3,22)
        prg = prg__LCG(7,9,11,98)
        num_attempts_per_nc = 100

        pm2 = prg__integer_sets__mult(n,c,r,crange,\
        mrange,prg,num_attempts_per_nc)

        assert pm2 == [\
            [-8, -64, 0, 8, -64, -64, 0, 8, \
            -64, -64, 0, 8, -64, -64, 0], \
            [59, 1239, 767, 295, -177, 1239, \
            767, 295, -177, 1239, 767, 295, \
            -177, 1239, 767, 295], \
            [-1, -2, -18, -17, -10, -2, -18, \
            -17, -10, -2, -18, -17, -10], \
            [13, 91, -39, 208, 182, 91, -39, \
            208, 182, 91, -39, 208, 182], \
            [52, 572, 988, 156, 208, 572, 572, \
            988, 156, 208, 572, 572, 988, 156, \
            208, 572], \
            [-163, -138, -109, -142, -145, -172, \
            -121, -137, -100, -159, -102, -166, \
            -165, -173, -135, -82, -95]]

        c = 0
        for pm_ in pm2: c += len(pm_)
        assert c == 90 
        return

if __name__ == '__main__':
    unittest.main()