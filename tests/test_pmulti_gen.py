from seqgen.pmulti_gen import * 
from morebs2.numerical_generator import prg__LCG,prg_to_prg__LCG_sequence 
from morebs2.aprng_gauge import coverage_of_sequence
import unittest

### lone file test 
"""
py -m tests.test_pmulti_gen 
"""
###
class PartitionedMultiGenMethodTests(unittest.TestCase):

    def test__PartitionedMultiGen__next__case1(self): 
        prg = prg__LCG(7334.2542,252.522,-6752.784,112675.66)
        prg_seq = prg_to_prg__LCG_sequence(prg,10,3.42)
        super_range = [-2000.5,10000.67] 
        segment_size_range = [20,112]

        pmg = PartitionedMultiGen(prg_seq,super_range,segment_size_range)

        L = np.round([next(pmg) for _ in range(9000)],5) 

        X = coverage_of_sequence(L,super_range,max_radius=0.5)
        assert X == 0.50921,"got {}".format(X)

        X2 = coverage_of_sequence(L,super_range,max_radius=3.)
        assert X2 == 0.98102

        assert len(set(L)) == 9000 

if __name__ == '__main__':
    unittest.main()