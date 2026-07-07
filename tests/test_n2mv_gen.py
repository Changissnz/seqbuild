from seqgen.n2mv_gen import *
#from morebs2.numerical_generator import prg__LCG
from morebs2.numerical_generator import prg__n_ary_alternator

import unittest

def N2MVGen__sample_QP(prg,prg2=None): 
    nm_range = (3,13)
    #prg = prg__LCG(45,-15,177,166)
    #prg = prg__LCG(155.66,83.455,-321.45,10000 + 41/43)
    #prg = prg__constant(0) 
    #prg = prg__n_ary_alternator(-100,100,-100)

    function_update_frequency_range = (75,150)
    return N2MVGen(nm_range,prg,prg2,function_update_frequency_range)

### lone file test 
"""
py -m tests.test_n2mv_gen 
"""
###
class N2MVGenMethods(unittest.TestCase):

    '''
    '''
    def test__N2MVGen__next__case1(self):

        # subcase 1: 1000 samples, 979 unique 
        prg = prg__n_ary_alternator(-100,100,-100)
        ng = N2MVGen__sample_QP(prg) 

        L = [] 
        for _ in range(1000): 
            L.append(next(ng))
        assert len(set(L)) == 979 

        # subcase 2: 2000 samples, 1866 unique 
        prg = prg__n_ary_alternator(-100,100,-100)
        ng = N2MVGen__sample_QP(prg) 

        L = [] 
        for _ in range(2000): 
            L.append(next(ng))
        assert len(set(L)) == 1866  

    """
    base PRNG is constant => output has 1 unique value. 
    """
    def test__N2MVGen__next__case2(self): 

        # subcase 1 
        prg = prg__constant(0) 
        ng = N2MVGen__sample_QP(prg) 
        L = [] 
        for _ in range(1000): 
            L.append(next(ng))
        assert len(set(L)) == 1 

        # subcase 2 
        prg = prg__constant(5) 
        ng = N2MVGen__sample_QP(prg) 
        L = [] 
        for _ in range(1000): 
            L.append(next(ng))
        assert len(set(L)) == 1 

    def test__N2MVGen__next__case3(self): 

        # subcase 1 
        prg = prg__LCG(155.66,83.455,-321.45,10000 + 41/43)
        ng = N2MVGen__sample_QP(prg) 

        L = [] 
        for _ in range(10000): 
            L.append(next(ng))
        assert len(set(L)) == 9979

        # subcase 2 
        prg = prg__LCG(155.66,83.455,-321.45,10000 + 41/43)
        prg2 = prg__n_ary_alternator(-100,100,-100)
        ng2 = N2MVGen__sample_QP(prg,prg2)  

        L2 = [] 
        for _ in range(10000): 
            L2.append(next(ng2))
        assert len(set(L2)) == 9979, "got {}".format(len(set(L2)))

        # check that the two cases produce entirely different sequences, 
        # due to the second PRNG instantiated with a non-null `prg2`. 
        L,L2 = np.array(L), np.array(L2) 

        D = np.round(L - L2,5) 
        Q = np.where(D != 0)[0] 
        assert len(Q) == 10000, "got {}".format(len(Q)) 

if __name__ == '__main__':
    unittest.main()