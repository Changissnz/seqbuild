from seqgen.shadow_gen import * 
from morebs2.numerical_generator import * 
from intigers.extraneous import * 
import unittest

### lone file test 
"""
python3 -m tests.test_shadow_gen 
"""
###
class QualVecMethods(unittest.TestCase):

    def test__QualVec__apply_noise__case1(self):
        its = IntSeq([5,4,3])
        vec = np.array([2,6,171,5141,141,41,1,414,17352,13416])
        qual = "tvec"
        prg = prg__LCG(87,51,-311,40040)
        prg = wrap_ranged_modulo_over_generator(prg,(-40000,40000.1))  

        qv = QualVec(deepcopy(vec),qual,qual_op=sub)
        xx = qv.apply_noise(prg) 

        q = np.sum(np.abs(xx - vec))
        assert q == np.int64(45574) 

    def test__QualVec__apply_noise__case2(self):  
        vec = np.array([2,6,171,5141,141,41,1,414,17352,13416])
        qual = "fvec"
        prg = prg__LCG(87,51,-311,40040)

        qv2 = QualVec(deepcopy(vec),qual,qual_op=sub) 
        xx = qv2.apply_noise(prg) 

        xx_sol = np.array([2.0000e+00, 6.0000e+00, 1.7100e+02, 5.1410e+03,\
            1.4100e+02,4.1000e+01,np.nan, 4.1400e+02, 1.7352e+04, 1.3416e+04]) 
        assert np.isnan(xx[6]) 

        xx[6] = 4 
        xx_sol[6] = 4 
        assert np.all(xx == xx_sol) 

    def test__QualVec__apply_noise__case3(self):  
        vec = np.array([2,6,171,5141,141,41,1,414,17352,13416])
        qual = "optri"
        prg = prg__LCG(87,51,-311,40040)

        qv2 = QualVec(deepcopy(vec),qual,qual_op=sub,inverted_qual_op=add) 
        xx = qv2.apply_noise(prg) 

        err = np.sum(np.abs(np.array(xx) - vec))
        err_sol = np.float64(536055.0)
        assert err == err_sol

if __name__ == '__main__':
    unittest.main()
