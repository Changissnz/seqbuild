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

    def test__QualVec__apply_noise(self):
        its = IntSeq([5,4,3])
        vec = np.array([2,6,171,5141,141,41,1,414,17352,13416])
        qual = "tvec"
        prg = prg__LCG(87,51,-311,40040)
        prg = wrap_ranged_modulo_over_generator(prg,(-40000,40000.1))  

        qv = QualVec(deepcopy(vec),qual,qual_op=sub)
        xx = qv.apply_noise(prg) 

        q = np.sum(np.abs(xx - vec))
        assert q == np.int64(45574)

if __name__ == '__main__':
    unittest.main()
