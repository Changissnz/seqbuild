from intigers.mod_prng import * 
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
        P = NoiseAddedIterableContainer(L,v) 

        X0 = [int(next(P)) for _ in range(len(L))] 
        X1 = [int(next(P)) for _ in range(len(L))] 
        X2 = [int(next(P)) for _ in range(len(L))] 

        assert X0 == L 
        assert X1 == list(array(L) + v)
        assert X2 == list(array(L) + 2 * v)

if __name__ == '__main__':
    unittest.main()