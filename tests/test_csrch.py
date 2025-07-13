from mini_dm.csrch import * 
import unittest

### lone file test 
"""
python -m tests.test_csrch
"""
###
class ClosestCombinationSearchMethods(unittest.TestCase):

    def test__ClosestCombinationSearch__find__case1(self):
            
        k = 3 
        seq = [4,100,21,34,100,6,7,17,70,31] 

        t = [100,6,7,17] 
        cfunc = lambda x: x == t 

        cs = ClosestCombinationSearch(k,seq,cfunc,num_attempts = 10 ** 5)
        r,s = cs.find() 
        assert r == t 
        assert cs.num_attempts_ == 99704 

        t3 = [10070316] 
        t3 = [70,31]
        t3 = [-7]
        cfunc = lambda x: set(x) == set(t3) 
        cs3 = ClosestCombinationSearch(k,seq,cfunc,num_attempts = 10 ** 5)
        stat = False
        r,s = cs3.find() 
        assert type(r) == type(None) and not s 

        t2 = [100,70,31,6] 
        cfunc = lambda x: set(x) == set(t2) 
        cs2 = ClosestCombinationSearch(k,seq,cfunc,num_attempts = 10 ** 5)
        r,s = cs2.find() 
        assert r == [100,6,70,31] 
        assert cs2.num_attempts_ == 99699

if __name__ == '__main__':
    unittest.main()
