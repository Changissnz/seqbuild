from intigers.tvec import * 
from morebs2.numerical_generator import prg__LCG,prg__constant
from intigers.mod_prng import * 
from mini_dm.minmax_freq import * 
import unittest

### lone file test 
"""
py -m tests.test_tvec
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

    def test__generate_m_unique_trinary_vectors__case_1(self): 

        prg = prg__constant(1) 
        prg2 = prg__LCG(-57,12,312,-311) 

        ndim = 5  
        m = 15 

        T = generate_m_unique_trinary_vectors(ndim,m,prg,attempt_ratio=2.0)
        T2 = generate_m_unique_trinary_vectors(ndim,m,prg2,attempt_ratio=2.0)

        assert len(list(T)) == len(list(T2)) == m 

        T3 = generate_m_unique_trinary_vectors(ndim,m,prg,attempt_ratio=1.0) 
        T4 = generate_m_unique_trinary_vectors(ndim,m,prg2,attempt_ratio=1.0)

        assert len(list(T3)) == 11 
        assert len(list(T4)) == 15  

if __name__ == '__main__':
    unittest.main()
