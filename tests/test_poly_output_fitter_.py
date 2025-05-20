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

    def test__PolyOutputFitterVar2__resolve(self):

        n = 8 
        x1,x2 = 3,7
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=3,prng=None)
        pofv.solve()
        pofv.resolve(4,21)
        assert pofv.apply(x1) == pofv.apply(x2)
        assert pofv.poly[n-4] == 21

        n = 6
        x1,x2 = -3,5
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=2,prng=None)
        pofv.solve()
        pofv.resolve(3,30)
        assert pofv.apply(x1) == pofv.apply(x2)
        assert pofv.poly[n-3] == 30

        n = 5
        x1,x2 = -13,8
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=4,prng=None)
        pofv.solve()
        pofv.resolve(3,30)
        assert pofv.apply(x1) == pofv.apply(x2)
        assert pofv.poly[n-3] == 30

class UDLinSysSolverMethods(unittest.TestCase):

    def test__UDLinSysSolver__apply_case1(self):
        # case 1 
        M = np.array([\
            [5,7,1,1,1],\
            [10,14,2,3,2],\
            [15,20,2,3,2]])
        Y = np.array([10,20,38])

        ulss = UDLinSysSolver(M,Y)
        ulss.initial_eval()
        ulss.cancel() 
        ulss.postcancel_solve()

        fx = {2:3,4:4} 
        ulss.set_freevar_values(fx)

        q0 = ulss.apply(ulss.m[0])
        q1 = ulss.apply(ulss.m[1])
        q2 = ulss.apply(ulss.m[2])

        assert q0 == ulss.y[0] 
        assert q1 == ulss.y[1] 
        assert q2 == ulss.y[2] 

        # case 2 
        M = np.array([\
            [3,4,7,5,4],\
            [2,5,-4,3,-7],\
            [5,-12,6,-4,10]])
        Y = np.array([20,12,-36])

        ulss = UDLinSysSolver(M,Y)
        ulss.initial_eval()
        ulss.cancel() 
        ulss.postcancel_solve()

        fx = {0:5,1:30} 
        ulss.set_freevar_values(fx)

        q0 = ulss.apply(ulss.m[0])
        q1 = ulss.apply(ulss.m[1])
        q2 = ulss.apply(ulss.m[2])

        assert q0 == ulss.y[0] 
        assert q1 == ulss.y[1] 
        assert q2 == ulss.y[2] 

    def test__UDLinSysSolver__apply_case2(self):

        M = np.array([\
            [3,4,7,5,4,14],\
            [2,5,-4,3,-7,17],\
            [5,-12,6,-4,10,31],\
            [31,10,-4,6,-12,5]])
        Y = np.array([20,12,-36,120])


        ulss = UDLinSysSolver(M,Y)
        ulss.solve()

        fx = {1:5,3:30}
        ulss.set_freevar_values(fx)

        q0 = ulss.apply(ulss.m[0])
        q1 = ulss.apply(ulss.m[1])
        q2 = ulss.apply(ulss.m[2])
        q3 = ulss.apply(ulss.m[3])

        assert q0 == ulss.y[0] 
        assert q1 == ulss.y[1] 
        assert q2 == ulss.y[2] 
        assert q3 == ulss.y[3] 


if __name__ == '__main__':
    unittest.main()