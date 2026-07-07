from seqgen.lps_gen import * 
from morebs2.numerical_generator import * 
import unittest

### lone file test 
"""
py -m tests.test_lps_gen 
"""
###
class LPSGenMethods(unittest.TestCase):

    def test__LPSGen__next__case1(self):

        prg = prg__LCG(1,3,2,101) 
        prg2 = None 
        prg3 = None 

        lg = LPSGen(prg,prg2,prg3)

        L = [] 
        for _ in range(1000): 
            L.append(next(lg))
        assert len(set(L)) == 41 

        prg = prg__LCG(1,3.2,4,101.1) 
        lg = LPSGen(prg,prg2,prg3)

        L2 = [] 
        for _ in range(1000): 
            L2.append(next(lg))
        assert len(set(L2)) == 1000 

        prg = prg__LCG(1,3.2,4.1,11.1) 
        lg = LPSGen(prg,prg2,prg3)

        L3 = [] 
        for _ in range(1000): 
            L3.append(next(lg))
        assert len(set(L3)) == 15  


if __name__ == '__main__':
    unittest.main()
