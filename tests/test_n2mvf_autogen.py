from intigers.n2mvf_autogen import *
from morebs2.numerical_generator import prg__LCG
import unittest

### lone file test 
"""
python -m tests.test_n2mvf_autogen
"""
###
class LCPVectorMap__TypeCShiftMethods(unittest.TestCase):

    '''
    test demonstrates sensitivity to changes in input, 
    as demonstrated by the two very different output 
    vectors. These changes could not have been possible 
    under linear (normal,index-to-index) vector multiplication. 
    '''
    def test__LCPVectorMap__TypeCShift__one_LCPVectorMap__case1(self):
        nm = (5,9)
        subvec_size_shifter = prg__iterable([3,2,2])
        prg1 = prg__LCG(43,100,31,511)
        prg2 = prg__LCG(49,93,5,499) 
        lmap = LCPVectorMap__TypeCShift.one_LCPVectorMap(nm,subvec_size_shifter,prg1,prg2)

        s0 = np.array([1,2,3,4,5])
        s1 = np.array([1,2,3,4,50])

        y0 = lmap.fit(s0)
        y1 = lmap.fit(s1)
        assert not np.any(np.round(y0 - y1,5) == 0)



if __name__ == '__main__':
    unittest.main()

