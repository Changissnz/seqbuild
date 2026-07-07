from seqgen.fit22_gen import * 
from morebs2.numerical_generator import * 

import unittest

### lone file test 
"""
py -m tests.test_fit22_gen 
"""
###
class Fit22GenMethods(unittest.TestCase):

    def test__Fit22Gen__next__case1(self):

        # subcase 1 

        prg = prg__LCG(1,3.2,4,101.1) 
        fg = Fit22Gen(prg,None,None,None) 

        L = [] 
        for _ in range(1000): 
            L.append(next(fg))
        assert len(set(L)) == 1000 

        # subcase 2 
        prg = prg__LCG(1,3,4,101) 
        fg2 = Fit22Gen(prg,None,None,None) 

        L2 = [] 
        for _ in range(1000): 
            L2.append(next(fg2))
        assert len(set(L2)) == 138, "got {}".format(len(set(L2)))

    def test__Fit22Gen__next__case2(self): 
        # subcase 1 
        prg = prg__n_ary_alternator(-50,50,-50)
        prg2 = prg__LCG(4.5,23,7.1,11.7)

        fg2 = Fit22Gen(deepcopy(prg),None,None,None) 

        L2 = [] 
        for _ in range(1000): 
            L2.append(next(fg2))
        assert len(set(L2)) == 83 

        # subcase 2 

        fg3 = Fit22Gen(prg2,prg,prg2,prg) 

        L3 = [] 
        for _ in range(1000): 
            L3.append(next(fg3))
        assert len(set(L3)) == 1000 

        fg4 = Fit22Gen(prg,prg2,prg,prg2) 

        # subcase 3  

        L4 = [] 
        for _ in range(1000): 
            L4.append(next(fg4))
        assert len(set(L4)) == 69 

if __name__ == '__main__':
    unittest.main()
