from intigers.poly_output_fitter_ import * 
import unittest

### lone file test 
"""
python3 -m tests.test_poly_output_fitter_
"""
###
class PolyOutputFitterMethods(unittest.TestCase):

    def test__PolyOutputFitterVar2__solve(self):
        n = 5 
        x1,x2 = -3,7
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 

        n = 5
        x1,x2 = 3,7
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 

        n = 4
        x1,x2 = -6,9
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 

        n = 3
        x1,x2 = 2,5
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 

        n = 7 
        x1,x2 = -17,5
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=2,prng=None)
        pofv.solve()
        assert pofv.apply(x1) == pofv.apply(x2) 

        n = 8 
        x1,x2 = -22,4
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=3,prng=None)
        pofv.solve()
        assert pofv.apply(x1) == pofv.apply(x2) 

if __name__ == '__main__':
    unittest.main()