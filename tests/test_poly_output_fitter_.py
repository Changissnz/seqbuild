from intigers.poly_output_fitter_ import * 
from morebs2.numerical_generator import prg__constant,prg__n_ary_alternator

import unittest

### lone file test 
"""
python -m tests.test_poly_output_fitter_
"""
###
class PolyOutputFitterMethods(unittest.TestCase):

    def test__DEFAULT_MAXPOW4BASE(self):
        x = 0
        assert DEFAULT_MAXPOW4BASE(x) == DEFAULT_POWER_RANGE[1] 

        x = 1
        assert DEFAULT_MAXPOW4BASE(x) == DEFAULT_POWER_RANGE[1] 

        x = 481
        assert DEFAULT_MAXPOW4BASE(x) == 3

        x = 1290
        assert DEFAULT_MAXPOW4BASE(x) == 3

        x = 200
        assert DEFAULT_MAXPOW4BASE(x) == 4 


    def test__PolyOutputFitterVar1__solve(self):
        n = 7 
        x1 = 100
        c = 10 ** 7 + 3 * 10 ** 6 + 2 * 10 **5 + \
            6 * 10 ** 3 + 7 * 10 ** 2 + 4 * 10 ** 1 + 8

        prg = prg__n_ary_alternator(s0=2,s1=59,start=-5)

        pofv1 = PolyOutputFitterVar1(n,x1,c,prg,default_sizemod=True)
        pofv1.solve()
        x1,c = pofv1.x1,pofv1.c
        assert pofv1.apply(x1) == c 

        cep = pofv1.to_CEPoly()
        scep = str(cep)
        assert scep == "-5983x^7  -36009x^6  -216165x^5  " +\
            "-1297162x^4  -7782973x^3  -46697843x^2  " +\
            "-280187064x^1 + 1681122397x^0" 
        assert cep.apply(pofv1.x1) == pofv1.apply(pofv1.x1)

    """
    demonstrates size limitations for np.int32 fitting 
    """
    def test__PolyOutputFitterVar1__solve__case2(self):
        l = 24 
        x = 43215678
        prng = prg__constant(0)

        pofv1 = PolyOutputFitterVar1(7,l//2,x,prng,default_sizemod=True)
        pofv1.solve()
        assert pofv1.apply(l//2) != pofv1.c
        assert pofv1.apply(pofv1.x1) == pofv1.c

        pofv2 = PolyOutputFitterVar1(7,l,x,prng,default_sizemod=False)
        pofv2.solve()
        try: 
            pofv2.apply(l) 
            assert False
        except: 
            assert True 

    def test__PolyOutputFitterVar2__solve(self):
        n = 5 
        x1,x2 = -3,7
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 1x^5  -9x^4 + 11x^3  -7x^2 + 4x^1'

        n = 5
        x1,x2 = 3,7
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 1x^5  -9x^4 + 14x^3  -4x^2 + 13x^1'

        n = 4
        x1,x2 = -6,9
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 1x^4  -7x^3 + 31x^2  -3x^1'

        n = 3
        x1,x2 = 2,5
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=1,prng=None)
        pofv.solve() 
        assert pofv.apply(x1) == pofv.apply(x2) 
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 1x^3  -7x^2 + 10x^1'

        n = 7 
        x1,x2 = -17,5
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=2,prng=None)
        pofv.solve()
        assert pofv.apply(x1) == pofv.apply(x2) 
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 2x^7  -36x^6  -1189x^5 + 33x^4 + 1080x^3  -15x^2  -337x^1'

        n = 8 
        x1,x2 = -22,4
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=3,prng=None)
        pofv.solve()
        assert pofv.apply(x1) == pofv.apply(x2)
        cep = pofv.to_CEPoly()
        assert str(cep) == ' 3x^8 + 66x^7  -2x^6  -46x^5 + 39x^4 + 1701x^3  -24x^2  -836x^1'
        
        
    def test__PolyOutputFitterVar2__solve__case2(self):

        n=3
        x1,x2 = 7000,32000 
        coeff = 3#992
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=coeff,prng=None,default_sizemod=False)
        pofv.solve() 
        x1,x2 = pofv.x1,pofv.x2 
        assert pofv.apply(x1) == pofv.apply(x2) 

        n=5
        x1,x2 = 59,12
        coeff = 31
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=coeff,prng=None,default_sizemod=False)
        pofv.solve() 
        x1,x2 = pofv.x1,pofv.x2 
        assert pofv.apply(x1) == pofv.apply(x2) 

    def test__PolyOutputFitterVar2__solve__case3(self):

        n=3
        x1,x2 = 7000,320#17,13 
        coeff = 3#992
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=coeff,prng=None,default_sizemod=False,\
            order_pair=False)
        pofv.solve() 
        cep = pofv.to_CEPoly()
        assert pofv.apply(x1) == 4216440480000
        assert pofv.apply(x2) == 100830489600
        assert not pofv.is_solved()

    """
    case that uses PRNG 
    """
    def test__PolyOutputFitterVar2__solve__case4(self):
        
        n=8
        x1,x2 = 7000,320 
        coeff = 3

        prng = prg__constant(123)
        pofv = PolyOutputFitterVar2(n,x1,x2,coeff=coeff,prng=prng,default_sizemod=True,\
            order_pair=True)
        pofv.solve() 
        cep = pofv.to_CEPoly()

        x1,x2 = pofv.x1,pofv.x2 
        assert pofv.apply(x1) == pofv.apply(x2) 

        pofv.resolve(3,6)
        x1,x2 = pofv.x1,pofv.x2 
        assert pofv.apply(x1) == pofv.apply(x2) 
        
        return 

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
        ulss.solve() 

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
        ulss.solve()

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

    '''
    case of stealthy inconsistent system 
    '''
    def test__UDLinSysSolver__apply_case3(self):
        M = np.array([\
            [1,7,1,7],\
            [7,1,7,1]])

        Y = np.array([10,20])

        ulss = UDLinSysSolver(M,Y)
        ulss.solve() 

        d = {1:10,3:-5}
        ulss.set_freevar_values(d)
        ulss.solve_missing_reps(ulss.varvec) 
        d = {2:23}
        ulss.set_freevar_values(d,2) 
        assert not ulss.constat


if __name__ == '__main__':
    unittest.main()