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

class N2MVectorFunctionGenMethods(unittest.TestCase): 

    """
    demonstrates deterministic quality of <N2MVectorFunction> 
    and difference in values between modes `replace` and 
    `accumulate`. 
    """
    def test__N2MVectorFunctionGen__one_N2MVF__case1(self): 
        nm = (5,12) 
        prg = prg__LCG(45,60,100,312)
        prg2 = prg__LCG(145,68,-100,411)

        nvfg = N2MVectorFunctionGen(nm,prg,prg2,mode="replace")
        nvfg_ = deepcopy(nvfg) 
        q = nvfg.one_N2MVF() 

        x = np.array([4,14,3,13,32])
        r1 =q.apply(x)
        sol1 = np.array([ 300.,596.,2316.,596.,4980.,\
            2316.,424.,596.,443.,808.,596.,596.])
        assert np.all(r1 == sol1)

        nvfg_.mode = "accumulate" 
        q2 = nvfg_.one_N2MVF() 
        r2 = q2.apply(x) 
        assert not np.any(r1 == r2) 

    """
    observes the degrees of the m-space indices. 
    """
    def test__N2MVectorFunctionGen__one_N2MVF__case2(self): 

        # case 1 
        prg = prg__LCG(23,3,4,1200) 
        prg2 = prg__LCG(15,-68,100,4110)

        nm = (2,10) 
        nvfg = N2MVectorFunctionGen(nm,prg,prg2,mode="replace")
        q = nvfg.one_N2MVF() 

        mindex_q = {0: 2, 1: 2, 2: 2, 3: 2, \
                4: 0, 5: 0, 6: 2, 7: 2, 8: 3, \
                9: 2}
        assert q.mindex == mindex_q 

        x0 = np.array([4,13])
        r0 =q.apply(x0)

        x1 = np.array([4,7]) 
        r1 =q.apply(x1)

        assert np.all(np.where(r0 == r1)[0] == [4,5])

        x2 = np.array([14,13]) 
        r2 =q.apply(x2)
        assert np.all(r0 == r2)

        # case 2 
        prg = prg__LCG(22,41,170,20010) # performs well! 
        nm = (4,9) 
        nvfg2 = N2MVectorFunctionGen(nm,prg,prg2,mode="replace")
        q2 = nvfg2.one_N2MVF() 

        mindex_q2 = {0: 6, 1: 12, 2: 18, 3: 12, 4: 12, 5: 6, 6: 6, 7: 13, 8: 6}
        assert mindex_q2 == q2.mindex 

if __name__ == '__main__':
    unittest.main()

