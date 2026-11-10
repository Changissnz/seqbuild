from intigers.mod_prng import * 
from intigers.prng_pw_op import * 
from morebs2.numerical_generator import * 
import unittest

### lone file test 
"""
py -m tests.test_mod_prng
"""
###
class ModPRNGClasses(unittest.TestCase):
 
    def test__NoiseAddedIterableContainer__next__case1(self):
        L = [414,54,76,91,231,451,-26,-96,64]
        v = 23 
        P = NoiseAddedIterableContainer(L,v,op_first=False) 

        X0 = [int(next(P)) for _ in range(len(L))] 
        X1 = [int(next(P)) for _ in range(len(L))] 
        X2 = [int(next(P)) for _ in range(len(L))] 

        assert X0 == L 
        assert X1 == list(array(L) + v)
        assert X2 == list(array(L) + 2 * v)

        # check for correct reset, w/ index 
        v2 = 45
        P.reset(reset_index=True,v=v2) 

        assert next(P) == L[0] + v2 
        assert next(P) == L[1] + v2 

        # check for correct reset, w/o index 
        v3 = 66
        P.reset(reset_index=False,v=v3)
        assert next(P) == L[2] + v2 + v3 

        # check for correct reset, w/ index 
        P.reset(reset_index=True,v=v3)
        assert next(P) == L[0] + v3 

    """
    same parameters and function calls as in case #1, except 
    `op_first` is set to True in the beginning. 
    """
    def test__NoiseAddedIterableContainer__next__case2(self):

        L = [414,54,76,91,231,451,-26,-96,64]
        v = 23 
        P = NoiseAddedIterableContainer(L,v,op_first=True) 

        X0 = [int(next(P)) for _ in range(len(L))] 
        X1 = [int(next(P)) for _ in range(len(L))] 
        X2 = [int(next(P)) for _ in range(len(L))] 

        assert X0 == list(array(L) + v) 
        assert X1 == list(array(L) + 2 * v)
        assert X2 == list(array(L) + 3 * v)

        next(P) 
        next(P) 

        # check for correct reset, w/ index 
        v2 = 45
        P.reset(reset_index=True,v=v2) 

        assert next(P) == L[0] + v2 
        assert next(P) == L[1] + v2 

        # check for correct reset, w/o index 
        v3 = 66
        P.reset(reset_index=False,v=v3)
        assert next(P) == L[2] + v2 + v3 

        # check for correct reset, w/ index 
        P.reset(reset_index=True,v=v3)
        assert next(P) == L[0] + v3 

    """
    uses <PRNGPairwiseOp> instead of standard pairwise op. 
    """
    def test__NoiseAddedIterableContainer__next__case3(self):

        L = [414,54,76,91,231,451,-26,-96,64]
        v = 23 

        q = prg__LCG(5625.24,-241.51,6252.532,45156.33)
        F = PRNGPairwiseOp([add,mul,sub],q) 

        P = NoiseAddedIterableContainer(L,v,f=F.apply,op_first=True) 

        X = [next(P) for _ in range(20)] 
        P.reset(reset_index=True) 
        X2 = [next(P) for _ in range(20)] 
        assert X == X2 
        return 

if __name__ == '__main__':
    unittest.main()