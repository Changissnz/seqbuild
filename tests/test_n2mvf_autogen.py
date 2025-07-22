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
    via comparison of the two "very" different output vectors, 
    despite the only difference in index-to-index comparison of 
    input vectors being index 4. These changes could not have been 
    possible under linear vector multiplication (normal,index-to-index). 
    '''
    def test__LCPVectorMap__TypeCShift__one_LCPVectorMap__case1(self):
        nm = (5,9)
        subvec_size_shifter = prg__iterable([3,2,2])
        prg1 = prg__LCG(43,100,31,511)
        prg2 = prg__LCG(49,93,5,499) 
        lmap = LCPVectorMap__TypeCShift.one_LCPVectorMap(nm,subvec_size_shifter,prg1,prg2)

        s0 = np.array([1,2,3,4,5])
        s1 = np.array([1,2,3,4,50])

        y0 = lmap.apply(s0)
        y1 = lmap.apply(s1)
        assert not np.any(np.round(y0 - y1,5) == 0)

class ModulatedN2MVectorMapMethods(unittest.TestCase):

    def test__ModulatedN2MVectorMap__apply__case1(self):
        vx = np.array([1,2,30,46,15,7,9,10])
        modmap = ModulatedN2MVectorMap(vx,add) 

        sol1 = np.array([ 2,  4, 31, 48, 16,  9, 10, 12])
        x1 = modmap.apply(np.array([1,2])) 
        assert np.all(sol1 == x1) 

        sol2 = np.array([ 2,  4, 34, 47, 17, 11, 10, 12])
        x2 = modmap.apply(np.array([1,2,4]))
        assert np.all(sol2 == x2) 

        sol3 = np.array([11, 22, 40, 66, 25, 27, 19, 30])
        x3 = modmap.apply(np.array([10,20]))
        assert np.all(sol3 == x3) 

if __name__ == '__main__':
    unittest.main()

