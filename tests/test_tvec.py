from intigers.tvec import * 
from morebs2.numerical_generator import prg__LCG
from intigers.mod_prng import * 
from mini_dm.minmax_freq import * 
import unittest

### lone file test 
"""
python -m tests.test_tvec
"""
###
class TrinaryVectorMethods(unittest.TestCase):

    def test__TrinaryVector__one_instance__v1(self):
        prg = prg__LCG(43,22,31,3000)
        tv = TrinaryVec.one_instance__v1(0,20,0.5,prg)

        prg2 = prg__iterable([700,700,200,500])
        tv2 = TrinaryVec.one_instance__v1(0,20,0.5,prg2)

        assert np.all(tv.l == [0,0,0,1,1,0,0,1,0,1,1,0,1,0,1,1,1,1,0,0]) 
        assert np.all(tv2.l == [-1,0,0,-1,0,0,-1,-1,0,-1,-1,0,0,-1,-1,-1,0,-1,0,0])
        return 

    def test__TrinaryVector__one_instance__v2(self):

        prg = prg__LCG(43,22,31,3000)
        fm = {-1:0.333333,0:0.3333333,1:0.3333333}
        tv = TrinaryVec.one_instance__v2(20,fm,prg)
        q1 = vec_to_frequency_map(np.array(tv.l,dtype=np.int32))
        assert q1 == {1:7,-1:7,0:6}

        prg2 = prg__iterable([700,700,200,500])
        fm = {-1:0.2,0:0.1,1:0.7}
        tv2 = TrinaryVec.one_instance__v2(20,fm,prg2)
        q2 = vec_to_frequency_map(np.array(tv2.l,dtype=np.int32))
        assert q2 == {-1:4,1:14,0:2}

if __name__ == '__main__':
    unittest.main()
